#include <random>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream

using namespace std;
using namespace sdsl;


namespace gedmap_sample{

std::vector<std::mt19937> g1;

uint32_t rand_ui(){
	return g1[omp_get_thread_num()]();
}


/** @brief return a random unsigned int
 * rand() returns an positiv integer with MAX_VAL = 2147483647
 * This function return an unsigned integer with MAX_VAL 2*2147483647
 */
//uint32_t rand_ui(){
//	if(rand()%2)	return (uint32_t) rand() + (uint32_t) rand();
//	else			return (uint32_t) rand();
//}


string read_EDS_substring(uint32_t position, string & EDS, adjacency & adj, uint32_t length);

uint32_t get_rand_position_before_edge(adjacency & adj){
	uint64_t node_pos;
	while(true){
		uint32_t idx = rand_ui() % adj.backward_targets.size();
		node_pos = adj.backward_targets[idx];
		const auto targets = adj(node_pos);
		if( targets.size() > 1 ) break;		
	}
	return node_pos - 50;
}

tuple<string,string,uint32_t> add_errors(string read, uint32_t target_length);


string a_sample(string & EDS, pos_EDS_to_FA_type & p2FA, adjacency & adj){
	string read = "";
	std::stringstream sample_stream;
	uint32_t pos = -1;
	uint32_t length = (uint32_t)((double)LENGTH * 2 );
	
	while(read.size() < length){
		pos = rand_ui() % (EDS.size()-LENGTH);
		if(over_edges) pos = get_rand_position_before_edge(adj);
		if(!p2FA.empty() && !p2FA.ref_ind[pos]) continue; //when transform is used then dont start at variant positions
		read = read_EDS_substring(pos, EDS, adj, length);
		//when read_EDS_substring returns an error or a too short read we just try again
	}
	
	string 	sample;
	string 	CIGAR;
	uint32_t 	D;
	tie		(sample,CIGAR,D) = add_errors(read, LENGTH);
	
	bool rc = RC_RATE && !(rand_ui() % RC_RATE);
	if(rc) sample = gedmap_encode::rev_complement(sample);
	
	CIGAR = gedmap_encode::RL_encode(CIGAR);
			
	string chrom_name = "";
	string off;
	pos++;// 1-indexed
	
	uint32_t chrom_pos = 0;
	
	if(!p2FA.empty()) tie(chrom_name,chrom_pos,off) = p2FA(pos);
	else{chrom_name = '?'; chrom_pos = pos;}
	
	sample_stream << "@ref_" << chrom_name 
		<< "_pos_" << chrom_pos
		<< "_DIST_" << D
		<< "_CIGAR_"<< CIGAR 
		<< (rc?"_rc":"")
		<< '\n';
	
	sample_stream << sample << "\n+\n" << string(sample.size(),'?') << '\n';
	return sample_stream.str();
}

std::pair<std::string,std::string> mp_sample(string & EDS, pos_EDS_to_FA_type & p2FA, adjacency & adj) {
	string read;
	uint32_t pos = -1;
	uint32_t length = FRAGMENT_LENGTH * 2;
	
	while (read.size() < length){
		pos = rand_ui() % (EDS.size() - FRAGMENT_LENGTH);
		if(over_edges) pos = get_rand_position_before_edge(adj);
		if (!p2FA.empty() && !p2FA.ref_ind[pos]) continue; //when transform is used then dont start at variant positions
		read = read_EDS_substring(pos, EDS, adj,length);
		//when read_EDS_substring returns an error or a too short read we just try again
	}
	
	auto[sample, CIGAR, D] = add_errors(read, FRAGMENT_LENGTH);
	
	bool rc = RC_RATE && !(rand_ui() % RC_RATE);
	std::string sample_rev = gedmap_encode::rev_complement(sample);
	if(rc) std::swap(sample, sample_rev);
	
	//CIGAR = gedmap_encode::RL_encode(CIGAR);
			
	string chrom_name = "";
	string off;
	pos++; // 1-indexed
	
	uint32_t chrom_pos = 0;
	
	if (!p2FA.empty()) tie(chrom_name, chrom_pos, off) = p2FA(pos);
	else chrom_name = '?', chrom_pos = pos;
	
	const auto gen = [&](bool ot, const std::string& s, size_t len) {
		std::stringstream str;
		str << "@ref_" << chrom_name 
			<< "_pos_" << chrom_pos
			<< (rc?"_rc":"")
			<< '/' << char('1' + ot)
			<< '\n'
			<< s.substr(0, len)
			<< "\n+\n"
			<< std::string(std::min(len,s.size()), '?')
			<< '\n';
		return str.str();
	};
	return std::make_pair( gen(false, sample, LENGTH), gen(true, sample_rev, LENGTH) );
}

int main(int argc,  char** argv){	
	
	using namespace gedmap_io;
	using namespace std;
	
	dotline();
	print_row("GEDMAP SAMPLE");
	dotline();

	
	vector<sys_timer> tv = vector<sys_timer>(1);
	tv[0].start();
	
	flush_row( "LOADING...");
	string EDS; 
	adjacency adj;
	pos_EDS_to_FA_type p2FA;
	ofstream fastq;	
	string fname_out;
	
	std::ofstream fastq_mp;
	
	handle_input(argc,argv,EDS,adj,p2FA,fastq,fastq_mp,fname_out);

	if (fastq_mp.is_open()) flush_row("Running in paired-end mode");

	{
		omp_get_max_threads();
		for (int i = 0; i < omp_get_max_threads(); i++)
			g1.emplace_back(std::mt19937(i * COUNT * 123456LL));
	}
	
	
	flush_row("Generating: 0%");
	
	vector<string> reads(COUNT), reads2;
	uint32_t count = 0;
	if (fastq_mp.is_open())
	{
		reads2.resize(COUNT);
		#pragma omp parallel for 
		for(uint32_t i = 0; i < COUNT;i++){
			std::tie(reads[i], reads2[i]) = mp_sample(EDS, p2FA, adj);
			
			#pragma omp critical
			{
				count++;
				if ((count*100)%COUNT == 0) flush_row( "Generating: " + to_string((count*100)/COUNT) + "%");
			}
		}
		for (const auto& rd2 : reads2)
			fastq_mp << rd2;
		fastq_mp.close();
	}
	else
	{
		#pragma omp parallel for 
		for(uint32_t i = 0; i < COUNT;i++){
			reads[i] = a_sample(EDS, p2FA, adj);
			
			#pragma omp critical
			{
				count++;
				if((count*100)%COUNT == 0) flush_row( "Generating: " + to_string((count*100)/COUNT) + "%");
			}
		}
	}
	
	for(uint32_t i = 0; i < COUNT;i++) fastq << reads[i];
	
	flush_row( to_string(COUNT) + "/" + to_string(COUNT));
	fastq.close();	
	
	print_row("Generation completed");
	print_row("Reads stored in:",fname_out);
	dotline();
	print_row("Time overall:",tv[0].stop_and_get(), " s");
	dotline();
	
	return 0;
}



char random_base(){
	switch(rand_ui()%4){
		case 0: return 'A';
		case 1: return 'T';
		case 2: return 'G';
		default: return 'C';
	}
}


uint32_t next_special(uint32_t position, string & EDS){
	char c = EDS[position];
	while( c == 'A'|| c == 'T'|| c == 'G'|| c == 'C'){
		position++;
		c = EDS[position];
	}
	return position;
}



/** @brief return a random alternative of the alternative scope starting at position*/
string choose_alt(uint32_t & position, string & EDS){
	uint32_t end = position;
	
	deque<uint32_t> boundaries;	
	boundaries.push_back(position);
	
	while(EDS[end]!= ')'){
		end++;
		if(EDS[end] == '|'){
			boundaries.push_back(end);
		}
		if(end == EDS.size()) throw runtime_error("invalid EDS");
	}
	boundaries.push_back(end);
	
	uint32_t alt =  rand_ui() % (boundaries.size()-1); 
	
	position = end+1;
	
	string alt_s = EDS.substr(boundaries[alt]+1,boundaries[alt+1]-boundaries[alt]-1);	
	
			
	std::size_t found = alt_s.find('#');
	if (found!=std::string::npos){ throw runtime_error("invalid EDS??");}
	
	return alt_s;
}


/** @brief reads a substring of length SAMPLE_LENGTH   starting at position start_pos*/
string read_EDS_substring( uint32_t start_pos, string & EDS, adjacency & adj, uint32_t length){
	string EDS_buf = EDS.substr(start_pos, 3*length);
	
	uint32_t 	position = 0;
	string 	read = "";
	uint32_t 	maxN = 5; //DONT ALLOW MORE THAN 5 N IN READ

	while(read.size() < length && position < EDS_buf.size()){
		char c = EDS_buf[position];
		
		switch(c){
			case 'A':
			case 'T':
			case 'G':
			case 'C':{
				//uint32_t next_p = next_special(position, EDS);
				uint32_t next_p = next_special(position, EDS_buf);
				if(next_p - position > length - read.size())
					next_p = position + length - read.size();		
				read = read + EDS_buf.substr(position,next_p-position);
				position = next_p;
				break;
			}
			case 'N':
				if(maxN-- == 0) return "";
				read += random_base();
				position++;
				break;
			case '(':
				try{
					read += choose_alt(position,EDS_buf);
				}catch(std::runtime_error & e){
					//other problems
					return ""; //Whatever... try again
				}
				break;
			case ')':
				position++;
				break;
			case '|':
				//skip to ')'
				while(EDS[position] != ')')
					position ++;
				position++;
				break;
			case '#':
				try{
					if(!adj.initialised) return ""; //no ADJ-file given... dont go over #
					const auto targets = adj(position+start_pos);
					if(targets.size() == 0) return ""; //(END OF CHROMOSOM) ... try again 
					uint32_t random_target;
					random_target = targets[ rand_ui()%targets.size() ];
					
					int tmp = 0;
					while(over_edges && targets.size() > 1 && random_target==(position+start_pos) && tmp++ < 100)
						random_target = targets[ rand_ui()%(targets.size())];

					
					start_pos = random_target+1;
					position = 0;
					EDS_buf = EDS.substr(start_pos, 3*length);
					
					break;
				}catch(std::runtime_error & e){
					//is chromosom border or other problems
					return ""; //Whatever... try again
				}
				
			default:
				throw runtime_error("invalid EDS ???, position="+to_string(position)+ " c=" + c);
		}
	}
	return read;
}

tuple<string,string,uint32_t> add_errors(string read, uint32_t target_length){
	stringstream  sample;
	stringstream CIGAR;
	uint32_t D = 0;
	uint32_t written = 0;
	uint32_t i = 0;
	
	while (written < target_length){
		char c = read[i++];
		
		if( IN_RATE &&  !(rand_ui()%IN_RATE)){
		//INSERTION
			sample << random_base();
			D++;
			written++;
			CIGAR << 'I';
			i--;
			continue;
		}
		
		if( DEL_RATE && !(rand_ui()%DEL_RATE) ){
		//DELETION
			D++;
			CIGAR << 'D';
			continue;
		}
		
		if( ERROR_RATE && !(rand_ui()%ERROR_RATE)){
			//REPLACEMENT
			char c2 = c;
			while(c2 == c)
				c2 = random_base();			
			D++;
			written++;
			sample << c2;
			CIGAR << 'M';
			continue;
		}
		
		written++;
		sample << c;
		CIGAR << 'M';
	}
	
	return make_tuple(sample.str(),CIGAR.str(),D);
}

}

