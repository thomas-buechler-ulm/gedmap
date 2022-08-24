namespace gedmap_align_min{

vector<sys_timer> tv(3);

const uint32_t tv_ALL  = 0;
const uint32_t tv_LOAD = 1;
const uint32_t tv_MAP  = 2;


template<uint8_t t_width, class int_type>
void paralell_processor( gedmap_mini::minimizer_index & eoc, const std::string & EDS, const adjacency & ADJ, const pos_EDS_to_FA_type & p2FA, std::ifstream & fastq_s, std::ofstream& o_s);

template<uint8_t t_width, class int_type>
void paralell_processor_batch( gedmap_mini::minimizer_index & eoc, const std::string & EDS, const adjacency & ADJ, const pos_EDS_to_FA_type & p2FA, std::ifstream & fastq_s, std::ofstream& o_s);

template<uint8_t t_width, class int_type>
std::vector<fasta_read<t_width,int_type>> read_batch(std::ifstream & fastq_s);

template<uint8_t t_width, class int_type>
fasta_read<t_width,int_type> get_read_from_stream(std::ifstream & fastq_s);

int main(int argc,  char** argv){
	using namespace gedmap_io;
	using namespace std;
	using namespace sdsl;
	print_prog_headline("GEDMAP ALIGN");

	tv[tv_ALL].start();
	tv[tv_LOAD].start();

	gedmap_mini::minimizer_index mini;
	string EDS;
	adjacency ADJ = adjacency();
	pos_EDS_to_FA_type p2FA;

	ifstream fastq_in_stream;
	ofstream out_stream;

	handle_input(argc, argv, mini, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
	assert(ADJ.initialised);

	
	print_row("Search", argv[3]);
	print_row("in", argv[2]);
	print_row("using", argv[4]);
	dotline();

	tv[tv_LOAD].stop();
	tv[tv_MAP].start();

	const uint32_t int_width = bitsneeded(EDS.size());

	// DIFFERENTIATE BETWEEN INTEGER SIZES
	if		(int_width <=  8 ) paralell_processor< 8, uint8_t>(mini, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
	else if	(int_width <= 16 ) paralell_processor<16,uint16_t>(mini, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
	else if	(int_width <= 32 ) paralell_processor<32,uint32_t>(mini, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
	else if	(int_width <= 64 ) paralell_processor<64,uint64_t>(mini, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
	else throw runtime_error (" in align main: EDS too long");

	tv[tv_MAP].stop();

	print_row("Time for loading:",		tv[tv_LOAD]		.get(), " s");
	print_row("Time for mapping:",		tv[tv_MAP]		.get(), " s");
	print_row("Time for overall:", 		tv[tv_ALL].stop_and_get(), " s");
	dotline();
	return 0;
}

/** @brief 
 * read read from fastq file
 * map read in the index
 * write alignment to sam file
 */
template<uint8_t t_width, class int_type>
void paralell_processor( gedmap_mini::minimizer_index & mini, const std::string & EDS, const adjacency & ADJ,const pos_EDS_to_FA_type & p2FA, std::ifstream & fastq_s, std::ofstream& o_s){
	
	if(IN_ORDER){
		paralell_processor_batch<t_width,int_type>(mini,EDS,ADJ,p2FA,fastq_s,o_s);
		return;
	}
	
	string s;
	uint32_t read_count = 0;
	while( getline( fastq_s, s ) ) read_count++;
	read_count /= 4;
	fastq_s.clear();
	fastq_s.seekg(0);
	if(MAP_RC) read_count *= 2;
	

	gedmap_io::print_row("MAPPING:");
	uint32_t results	= 0;
	uint32_t mapped 	= 0;
	uint32_t count  	= 0;
	uint32_t chances_count	= 0;
	

	if(THREAD_COUNT) omp_set_num_threads(THREAD_COUNT);
	
	bool run = true;
	while(run){
		if(!MAP_RC){
			#pragma omp parallel for  schedule(dynamic,10)
			for(uint32_t r = 0; r < read_count; r++){
				
				fasta_read<t_width,int_type> read;
				bool ok = false;
				
				#pragma omp critical
				{
					if(fastq_s.good()){
						try{
							read = fasta_read<t_width,int_type>(fastq_s);
							ok = true;
						}catch(runtime_error & e){
							//NO READ LEFT
							run = false;
						}
					}
				}
				
				if(!ok) continue;
				
				uint32_t my_searches_chances = 0;
				
				for( uint32_t i = 0;
				 i < FRAGMENT_COUNT.size()
				 && (read.alignments.size() == 0 || get<1>(read.alignments[0]) >= DOUBT_DIST); i++){
					my_searches_chances++;
					read.get_fragments(FRAGMENT_COUNT[i], mini);
					read.get_positions(mini);
					read.find_hotspots(SPOT_SIZE, SPOT_HITS[i],CHECK_COLLI);
					read.start_aligner(EDS,ADJ,MAX_DIST[i],MAX_ALIGNS_C[i],MAX_ALIGNS_T[i],MAX_ALIGNS_O);
				}
				#pragma omp critical
				{
					uint32_t res;
					res = read.write_alignment(o_s,WRITE_FAILURE,p2FA);
					results += res;
					if(res) mapped++;
					count++;
					if(!(count%1000))
						gedmap_io::flush_row("Searched reads", to_string(count));
					chances_count += my_searches_chances;
				}
			}
		}else{
			#pragma omp parallel for  schedule(dynamic,10)
			for(uint32_t r = 0; r < read_count; r+=2){
				
				fasta_read<t_width,int_type> read;
				bool ok = false;
				
				#pragma omp critical
				{
					
					if(fastq_s.good()){
						try{
							read = fasta_read<t_width,int_type>(fastq_s);
							ok = true;
						}catch(runtime_error & e){
							//NO READ LEFT
							run = false;
						}
					}
				}
				
				if(!ok) continue;
				uint32_t my_searches_chances = 0;
				
				for( uint32_t i = 0;
				 i < FRAGMENT_COUNT.size()
				 && (read.alignments.size() == 0 || get<1>(read.alignments[0]) >= DOUBT_DIST) ; i++){
					my_searches_chances++;
					fasta_read<t_width,int_type> read_rev = fasta_read<t_width,int_type> (read.id +"_rev",gedmap_encode::rev_complement(read.sequence), "");
					read.get_fragments(FRAGMENT_COUNT[i], mini);
					read_rev.get_fragments(FRAGMENT_COUNT[i], mini); 
					read.get_positions(mini);
					read_rev.get_positions(mini);
					read.find_hotspots(SPOT_SIZE, SPOT_HITS[i],CHECK_COLLI);
					read_rev.find_hotspots(SPOT_SIZE, SPOT_HITS[i],CHECK_COLLI);
					read.start_aligner(EDS,ADJ,MAX_DIST[i],MAX_ALIGNS_C[i],MAX_ALIGNS_T[i],MAX_ALIGNS_O,read_rev);
				}
				#pragma omp critical
				{
					uint32_t res;
					res = read.write_alignment(o_s,WRITE_FAILURE,p2FA);
					results += res;
					if(res) mapped++;
					count++;
					if(!(count%1000))
						gedmap_io::flush_row("Searched reads", to_string(count));
					chances_count += my_searches_chances;
				}
			}
		}
	}
	fastq_s.close();
	o_s.close();

	gedmap_io::print_row("Searched reads", count);
	gedmap_io::print_row("Searches", chances_count);
	gedmap_io::print_row("Mapped reads", mapped);
	gedmap_io::print_row("Results", results);
	gedmap_io::dotline();
}


/** @brief 
 * read read from fastq file
 * map read in the index
 * write alignment to sam file
 */
template<uint8_t t_width, class int_type>
void paralell_processor_batch( gedmap_mini::minimizer_index & mini, const std::string & EDS, const adjacency & ADJ,const pos_EDS_to_FA_type & p2FA, std::ifstream & fastq_s, std::ofstream& o_s){
	fastq_s.clear();
	fastq_s.seekg(0);
	uint32_t batchsize = 10000;
	

	gedmap_io::print_row("MAPPING:");
	uint32_t results	= 0;
	uint32_t mapped 	= 0;
	uint32_t count  	= 0;
	uint32_t chances_count	= 0;
	

	if(THREAD_COUNT) omp_set_num_threads(THREAD_COUNT);
	
	bool run = true;
	while(run){
		vector<fasta_read<t_width,int_type>> reads;
		reads.reserve(batchsize);
		
		for(uint32_t i = 0; i < batchsize; i++){
			fasta_read<t_width,int_type> read;
			if(fastq_s.good()){
				try{
					reads.push_back(fasta_read<t_width,int_type>(fastq_s));
				}catch(runtime_error & e){
					//NO READ LEFT
					run = false;
					break;
				}
			}
		}
		
		uint32_t my_searches_chances = 0;	
		
		if(!MAP_RC){
			#pragma omp parallel for  schedule(dynamic,10)
			for(uint32_t r = 0; r < reads.size(); r++){				
				for( uint32_t i = 0;
				 i < FRAGMENT_COUNT.size()
				 && (reads[r].alignments.size() == 0 || get<1>(reads[r].alignments[0]) >= DOUBT_DIST); i++){
					my_searches_chances++;
					reads[r].get_fragments(FRAGMENT_COUNT[i], mini);
					reads[r].get_positions(mini);
					reads[r].find_hotspots(SPOT_SIZE, SPOT_HITS[i],CHECK_COLLI);
					reads[r].start_aligner(EDS,ADJ,MAX_DIST[i],MAX_ALIGNS_C[i],MAX_ALIGNS_T[i],MAX_ALIGNS_O);
				}
			}			
		}else{
			#pragma omp parallel for  schedule(dynamic,10)
			for(uint32_t r = 0; r < reads.size(); r++){			
				for( uint32_t i = 0;
				 i < FRAGMENT_COUNT.size()
				 && (reads[r].alignments.size() == 0 || get<1>(reads[r].alignments[0]) >= DOUBT_DIST) ; i++){
					my_searches_chances++;
					fasta_read<t_width,int_type> read_rev = fasta_read<t_width,int_type> (reads[r].id +"_rev",gedmap_encode::rev_complement(reads[r].sequence), "");
					reads[r].get_fragments(FRAGMENT_COUNT[i], mini);
					read_rev.get_fragments(FRAGMENT_COUNT[i], mini); 
					reads[r].get_positions(mini);
					read_rev.get_positions(mini);
					reads[r].find_hotspots(SPOT_SIZE, SPOT_HITS[i],CHECK_COLLI);
					read_rev.find_hotspots(SPOT_SIZE, SPOT_HITS[i],CHECK_COLLI);
					reads[r].start_aligner(EDS,ADJ,MAX_DIST[i],MAX_ALIGNS_C[i],MAX_ALIGNS_T[i],MAX_ALIGNS_O,read_rev);
				}
			}
		}
		
		for(uint32_t r = 0; r < reads.size(); r++){
			uint32_t res;
			res = reads[r].write_alignment(o_s,WRITE_FAILURE,p2FA);
			results += res;
			if(res) mapped++;
			count++;
			if(!(count%1000))
				gedmap_io::flush_row("Searched reads", to_string(count));
			chances_count += my_searches_chances;
		}
	}
	fastq_s.close();
	o_s.close();

	gedmap_io::print_row("Searched reads", count);
	gedmap_io::print_row("Searches", chances_count);
	gedmap_io::print_row("Mapped reads", mapped);
	gedmap_io::print_row("Results", results);
	gedmap_io::dotline();
}

}

