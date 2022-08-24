/**
 * @brief parses FA and VCF file to a EDS representation
 * 
 */

#include <utility> //pair
#include <set>
#include "sdsl/select_support.hpp"
#include <map>
#include <algorithm> //sort

namespace gedmap_parse{
using namespace std;
using namespace sdsl;

void parse2EDS(ifstream& ifs_FA, ifstream& ifs_VCF, ofstream& ofs_EDS, ofstream& log, uint32_t limit);

void include_node_boundaries(string & fname_eds, string & fname_adj);


vector<sys_timer> tv(3);

const uint32_t tv_ALL  	= 0;
const uint32_t tv_PARSE  = 1;
const uint32_t tv_UPDATE  = 2;


int_vector_buffer<1>	ref_ind_buf;
vector<string> 		chrom_names;
vector<uint64_t> 		chrom_starts;


//edges: CHROM -> PAIR = FROM TO 
typedef pair <uint64_t,uint64_t > edge_type;
vector< vector < edge_type > > edges;
vector< set<uint64_t> > positions;
int_vector_buffer<1>	alt_ind_buf;

bool LOG_FILE_USED;

int main(int argc, char** argv){
	using namespace gedmap_io;
	print_prog_headline("GEDMAP PARSE");
	tv[tv_ALL].start();
	
	ifstream 	ifs_FA;
	ifstream 	ifs_VCF;
	ofstream 	ofs_EDS;
	string 	fname_eds;
	string 	fname_adj; //ADJ	
	uint32_t 	limit;
	
	handle_input(argc, argv, ifs_FA, ifs_VCF,fname_eds,limit);
	
	
	
	
	
	//TODO LOGFILE
	string log_filename = string(argv[4]) + FEX_LOG; //TODO
	ofstream 	log(log_filename,ios_base::out|ios_base::trunc);
	
	//OUTPUT INFO
	print_row("Parse FA: "	, argv[2]);
	print_row("and  VCF: "	, argv[3]);
	print_row("to   EDS: "	, fname_eds);
	dotline();
	
	
	// FOR EDS2FA
	string 	fname_eds2fa = fname_eds + "." + FEX_2FA;
	string	fname_ref_ind_buf = TMP_DIR + TMP_PARSE_REF_IND;
	string	fname_alt_ind_buf = TMP_DIR + "/alt_ind_buf";
	ref_ind_buf = int_vector_buffer<1> (fname_ref_ind_buf, ios::out|ios_base::trunc);
	alt_ind_buf = int_vector_buffer<1> (fname_alt_ind_buf, ios::out|ios_base::trunc);
	chrom_names = vector<string>(0);
	chrom_starts = vector<uint64_t>(0);
	
	// DO PARSING
	tv[tv_PARSE].start();
	print_row("Read VCF and include short variations");
	ofs_EDS.open(fname_eds,ios_base::out|ios_base::trunc);	
	parse2EDS(ifs_FA, ifs_VCF,ofs_EDS,log,limit);	
	ref_ind_buf.close();
	alt_ind_buf.close();
	dotline();
	tv[tv_PARSE].stop();
	
	
	
	//INCLUDE POSITIONS
	if(INCLUDE_CNV){
		tv[tv_UPDATE].start();
		print_row("Include structural variations and genarate adjacency file");
		fname_adj = fname_eds + "." + FEX_ADJ;
		include_node_boundaries(fname_eds, fname_adj);
		tv[tv_UPDATE].stop();
	}
	
	
	//LOAD BUILD EDS2FA
	print_row("Generate pos2FA file");
	pos_EDS_to_FA_type pos2FA(fname_ref_ind_buf, chrom_names, chrom_starts);	
	store_to_file(pos2FA, fname_eds2fa);
	file_remove(fname_ref_ind_buf);
	file_remove(fname_alt_ind_buf);
	dotline();	
	
	
	tv[tv_ALL].stop();
	//WRITE INFO
	print_row("EDS stored to: ", fname_eds);
	print_row("ADJ stored to: ", fname_adj);
	print_row("POS2FA stored to: ", fname_eds2fa);	
	if(0 != (uint32_t) log.tellp()) print_row("Please notice logfile:",log_filename);
	else	file_remove(log_filename);
	dotline();	
	print_row("Time for parsing short vars:",	tv[tv_PARSE]	.get(), " s");
	if(INCLUDE_CNV)
	print_row("Time for including sv:",		tv[tv_UPDATE]	.get(), " s");
	print_row("Time in total:",			tv[tv_ALL]		.get(), " s");
	dotline();	
	return 0;
}

struct ref_ind_copy{
	sdsl::bit_vector* ref_ind; //input bitvector
	sdsl::int_vector_buffer<1>* ref_ind_out;//output bitvector stream
	uint64_t	next_pos;
	
	ref_ind_copy(sdsl::bit_vector*  ref_ind, sdsl::int_vector_buffer<1>* ref_ind_out):ref_ind(ref_ind),ref_ind_out(ref_ind_out){next_pos = 0;};
	
	void copy(uint64_t count){
		uint64_t last_pos = next_pos + count;
		while(next_pos < last_pos)
			*ref_ind_out << (*ref_ind)[next_pos++];
	};
	
	void copy_all(){
		while(next_pos < ref_ind->size())
			*ref_ind_out << (*ref_ind)[next_pos++];
	};
};


struct stream_to_stream_copy{
	ifstream* in;
	ofstream* out;
	char* buffer;
	uint32_t buffer_size;
	
	stream_to_stream_copy(ifstream* in, ofstream* out, uint32_t buffer_size):in(in),out(out),buffer_size(buffer_size){
		buffer = new char [buffer_size];
	};
	
	/**
	* @brief writes specified amount of data from one stream to another
	* @param in input stream
	* @param out output stream
	* @param buffer buffer for data
	* @param buffer_size
	* @param count number of chars to copy
	*/
	void copy(uint64_t count){
		while(count > buffer_size){
			in ->read (buffer,buffer_size);
			out->write(buffer,buffer_size);
			count -= buffer_size;
		}
		in->read (buffer,count);
		out->write(buffer,count);
	};

	/**
	* @brief writes all data left in one stream to another
	* @param in input stream
	* @param out output stream
	* @param buffer buffer for data
	* @param buffer_size
	*/
	void copy_all(){
		uint64_t pos =  in->tellg();
		in->seekg (0, in->end);
		uint64_t pos_end =  in->tellg();
		in->seekg(pos);
		copy(pos_end - pos);
	};
};




/**
 * @brief writes all data left in one stream to another
 * 
 * adds markers at the positions in the position array
 * generates new adjascency list from edges and the new positions
 * updates chrom_starts
 * updates ref_ind
 */
void include_node_boundaries(string & fname_eds, string & fname_adj){
	ifstream EDS_in;
	ofstream EDS_with_sv_out;
	string fname_eds_sv = fname_eds + "_with_sv";
	EDS_in.open(fname_eds ,ios_base::in);
	EDS_with_sv_out.open(fname_eds_sv,ios_base::out|ios_base::trunc);
	
	
	
	//Calculate positions of markers in finished EDS
	vector<uint64_t> 		new_chrom_starts(chrom_starts.size());
	map<uint64_t, uint64_t> pos_map;
	sdsl::bit_vector 		ref_ind;
	
	//UPDATE POSITIONS
	{
		sdsl::bit_vector 		alt_ind;
		load_from_file(ref_ind,	TMP_DIR + TMP_PARSE_REF_IND);
		load_from_file(alt_ind,	TMP_DIR + "/alt_ind_buf");	
		sdsl::select_support_mcl<1,1>  ref_ind_ss;
		sdsl::rank_support_v5<> ref_ind_rs;	
		sdsl::util::init_support(ref_ind_ss,&ref_ind);
		sdsl::util::init_support(ref_ind_rs,&ref_ind);
		
		uint64_t offset = 0; 
		for(uint64_t i = 0; i < chrom_names.size(); i++){
			new_chrom_starts[i] = chrom_starts[i] + offset;
			for(auto p = positions[i].begin(); p != positions[i].end(); p++){
				//find position p of chrom i in EDS:
				uint64_t eds_pos = ref_ind_ss(   *p +  ref_ind_rs(chrom_starts[i]) );
				while(alt_ind[eds_pos]) 
					//shift, because markers shoundnt be in alterative scope
					eds_pos++;
				//Add ofset because of previouse inserted markers
				eds_pos += offset++;
				pos_map[chrom_starts[i] + *p] = eds_pos;
			}
			
		}
	}
	
	//UPDATE EDGES
	{
		
		uint32_t edge_count = 0;	
		for(uint32_t i = 0; i < chrom_names.size(); i++){
			edge_count += edges[i].size();
		}
		
		vector<edge_type> new_edges(edge_count);
		uint32_t j = 0;
		for(uint32_t i = 0; i < chrom_names.size(); i++)
			for(auto e : edges[i])
				new_edges[j++] = make_pair(pos_map[chrom_starts[i] + e.first],pos_map[chrom_starts[i] + e.second]) ;
			
		adjacency A(new_edges);
		store_to_file(A,fname_adj);
	}
	
	
	//Write new EDS2FA AND EDS
	ref_ind_buf = int_vector_buffer<1> (TMP_DIR + TMP_PARSE_REF_IND, ios::out|ios_base::trunc);
	ref_ind_copy rc(&ref_ind,&ref_ind_buf);
	stream_to_stream_copy ssc(&EDS_in, &EDS_with_sv_out, 100*KB);
	
	uint64_t last_pos = -1;//in first block no marker was inserted. Compensates the -1 below.
	for (std::map<uint64_t,uint64_t>::iterator it=pos_map.begin(); it!=pos_map.end(); ++it){
		uint64_t count = it->second - last_pos -1; //-1 because of inserted marker
		last_pos = it->second;
		rc.copy(count);
		ssc.copy(count);
		EDS_with_sv_out.put('#');
		ref_ind_buf << 0;
	}
	ssc.copy_all();
	rc.copy_all();
	delete[] ssc.buffer;
	
	//close streams and overwrite old file
	EDS_in.close();
	EDS_with_sv_out.close();
	gedmap_io::file_move(fname_eds_sv,fname_eds);
	
	//update chrom start
	chrom_starts.swap(new_chrom_starts);
	
	//ref ind is already overwritten
	ref_ind_buf.close();
	
	
	
	
	return;
}

/*
void write_snp(vcfl vc, ofstream& out){
	out << "[";	
	ref_ind_buf << 0;
	out << vc.ref;
	ref_ind_buf << 1;
	for(ui i = 0; i < vc.alt.size(); i++)
		if(vc.alt[i] == 'A'||vc.alt[i] == 'C'||vc.alt[i] == 'T'||vc.alt[i] == 'G'){
			out << vc.alt[i];
			ref_ind_buf << 0;
		}
	out << "]";
	ref_ind_buf << 0;
}*/

void write_indel(vcfl vc, ofstream& out){	
// 	if(vc.ref.size() + vc.alt.size() == 3){//SNI
// 		out << vc.ref[0] << "?" << (vc.ref.size()>1?vc.ref[1]:vc.alt[1]);
// 		ref_ind_buf << 1<< 0 << (vc.ref.size()>1?1:0);
// 	}
// 	else{
		out << vc.ref[0] << "(|";
		ref_ind_buf << 1 << 0 << 0;
		alt_ind_buf << 0 << 1 << 1;
		out << (vc.ref.size()>1?vc.ref.substr(1):vc.alt.substr(1)) << ")";
		if(vc.ref.size()>1)
			for(unsigned int i = 0; i < vc.ref.length() - 1 ; i++){
				ref_ind_buf << 1;
				alt_ind_buf << 1;				
			}
		else 
			for(unsigned int i = 0; i < vc.alt.length() - 1 ; i++){
				ref_ind_buf << 0;
				alt_ind_buf << 1;							
			}
		ref_ind_buf << 0;
		alt_ind_buf << 1;
// 	}
}

void write_alt(vcfl vc, ofstream& out){
	deque<string> alts = vc.genaltQ();
	out << "(" << vc.ref;
	ref_ind_buf << 0;
	alt_ind_buf << 1;
	for(unsigned int i = 0; i < vc.ref.length(); i++){
		ref_ind_buf << 1;		
		alt_ind_buf << 1;
	}
	
	for(auto it = alts.begin(); it!=alts.end(); ++it){
		out << '|' << *it;
		ref_ind_buf << 0;		
		alt_ind_buf << 1;
		for(unsigned int i = 0; i < (*it).length(); i++){
			ref_ind_buf << 0;			
			alt_ind_buf << 1;
		}
	}
	out << ")";	
	ref_ind_buf << 0;
	alt_ind_buf << 1;
}


void copy_till_pos(uint32_t pos, uint32_t& current_fa_pos, ifstream& ifs_FA, ofstream& ofs_EDS, bool skip = false){
	char c;  
	while(current_fa_pos < pos){
		if(!ifs_FA.get(c)) 
			throw invalid_argument( "FA STREAM ENDED BEVOR POSITION " + to_string(pos) +" OF VARIATION WAS REACHED");
		
		if(c == '\n')
			continue;
		
		if(!(c == 'A' || c == 'T' || c == 'G'|| c == 'C'|| c == 'N')){
			gedmap_io::print_error("UNEXPECTED SYMBOL: ASCII=(" + to_string((uint32_t) c) + ")='" + c + "' IN FA" + " at pos " + to_string(ifs_FA.tellg()) + " replaced by 'N'");
			c = 'N';
		} 
		
		if(!skip){
			ofs_EDS << c;
			ref_ind_buf << 1;			
			alt_ind_buf << 0;
		}
		current_fa_pos++;
	}
}


void copy_till_end(ifstream& ifs_FA, ofstream& ofs_EDS, uint32_t& current_fa_pos){
	char c;
	while(ifs_FA.get(c)){
		
		if(c == '>'){
			//read gene name
			string next_chrom;
			getline(ifs_FA, next_chrom);
			next_chrom.substr(0,next_chrom.find('|')); //cut additional infos
			
			ofs_EDS << EDS_NODE_BOUNDARY;
			ref_ind_buf << 0;					
			alt_ind_buf << 1;
			chrom_names.push_back(next_chrom);
			chrom_starts.push_back(ofs_EDS.tellp());
			edges.push_back( vector<edge_type>(0) );
			positions.push_back( set<uint64_t>() );
			
			gedmap_io::print_row("SEQUENCE: '" + next_chrom + "' (at pos " + to_string( ofs_EDS.tellp()  ) +" of FA)");  
			continue;
		} 
		
		if (c == '\n')
			continue;
		
		if ( !(c == 'A' || c == 'T' || c == 'G'|| c == 'C'|| c == 'N') ){
			gedmap_io::print_error("UNEXPECTED SYMBOL: ASCII=(" + to_string((uint32_t) c) + ")='" + c + "' IN FA" + " at pos " + to_string(ifs_FA.tellg()) + "replaced by 'N'");
			c = 'N';
		} 
		ofs_EDS << c;
		ref_ind_buf << 1;		
		alt_ind_buf << 0;
		current_fa_pos++;
	}
}


string copy_till_next_chrom(ifstream& ifs_FA, ofstream& ofs_EDS, uint32_t& current_fa_pos,  string & target){
	char c;
	while(true){
		
		if(!ifs_FA.get(c)) 
			throw invalid_argument( "FA STREAM ENDED BEVOR SEQUENCE '" + target +"' WAS REACHED");
		
		if(c == '>'){
			//read gene name
			string next_chrom;
			getline(ifs_FA, next_chrom);
			return next_chrom.substr(0,next_chrom.find('|')); //cut additional infos
		} 
		
		if (c == '\n')	
			continue;
		
		if ( !(c == 'A' || c == 'T' || c == 'G'|| c == 'C'|| c == 'N') ){
			gedmap_io::print_error("UNEXPECTED SYMBOL: ASCII=(" + to_string((uint32_t) c) + ")='" + c + "' IN FA" + " at pos " + to_string(ifs_FA.tellg()) + "replaced by 'N'");
			c = 'N';
		} 
		ofs_EDS << c;
		ref_ind_buf << 1;
		alt_ind_buf << 0;
		current_fa_pos++;
	}
}

void parse2EDS(ifstream& ifs_FA, ifstream& ifs_VCF, ofstream& ofs_EDS, ofstream& log, uint32_t limit){
	string line;
	//FA
	getline(ifs_FA,line);//chromline
	
	string current_fa_chrom = line.substr(1,line.find('|'));
	uint32_t current_fa_pos = 0;
	
	//VCF
	vcfl last_v = vcfl();
	conflict_group cogru = conflict_group();
	
	uint32_t read_lines = 0;
	uint32_t applied_alt = 0;
	uint32_t applied_cnv = 0;
	
	chrom_names.push_back(current_fa_chrom);
	chrom_starts.push_back(ofs_EDS.tellp());
	edges.push_back( vector<edge_type>() );
	positions.push_back( set<uint64_t>() );
	
	gedmap_io::print_row("SEQUENCE: '" + current_fa_chrom + "' (at pos " + to_string( ofs_EDS.tellp()  ) +" of FA)"); 
	
	while(getline(ifs_VCF, line)){
		
		if(++read_lines % 100000 == 0) gedmap_io::flush_row("Read vcf lines", to_string(read_lines));
		
		
		if(line.size() == 0 || line[0] == '#') continue;			//skip comments
		
		vcfl v = vcfl(line);		
		if(v.equal(last_v))	//same variation as in line before
			continue;
		

				
		
		if(v.isALT() && ( limit==0 || v.ref.length() + v.alt.length() <= limit ) ){	
// applied_alt++; RESOLVE CONFLICTS
// if(cogru.conflict(v)){ //when conflict then added to cogru
// cogru.add(v);
// //log << linecount << CONFLICT_TXT <<endl;                            
// }else{ 
			
			if(!cogru.conflict(v)){ //when no conflict, write cogru
				//default case: write variation of cogru to file
				if(!cogru.empty()){
					vcfl vc = cogru.solve();
					copy_till_pos(vc.pos-1,current_fa_pos,ifs_FA,ofs_EDS);
					copy_till_pos(current_fa_pos+vc.ref.size(),current_fa_pos,ifs_FA,ofs_EDS,true);
// 					if(vc.isSNP()) 		write_snp(vc,ofs_EDS);						
// 					else 
					if (vc.isINDEL())	write_indel(vc,ofs_EDS);
					else	write_alt(vc,ofs_EDS);
				}

				
				while(v.chrom != current_fa_chrom){
					string next_chrom = copy_till_next_chrom(ifs_FA,ofs_EDS, current_fa_pos, v.chrom);
					current_fa_chrom =  next_chrom;
					current_fa_pos	= 0;
					ofs_EDS << EDS_NODE_BOUNDARY;
					ref_ind_buf << 0;					
					alt_ind_buf << 1;
					chrom_names.push_back(current_fa_chrom);
					chrom_starts.push_back(ofs_EDS.tellp());
					edges.push_back( vector<edge_type>() );
					positions.push_back( set<uint64_t>() );
					gedmap_io::print_row("SEQUENCE: '" + current_fa_chrom + "' (at pos " + to_string( ofs_EDS.tellp()  ) +" of FA)");  
				}
				cogru.clear();
				applied_alt++;
				cogru.add(v);
			}
		}
		else if (v.isCNV() && INCLUDE_CNV){
			while(v.chrom != current_fa_chrom){
				if(!cogru.empty()){
					vcfl vc = cogru.solve();
					copy_till_pos(vc.pos-1,current_fa_pos,ifs_FA,ofs_EDS);
					copy_till_pos(current_fa_pos+vc.ref.size(),current_fa_pos,ifs_FA,ofs_EDS,true);
// 					if(vc.isSNP()) 		write_snp(vc,ofs_EDS);						
// 					else 
					if (vc.isINDEL())	write_indel(vc,ofs_EDS);
					else	write_alt(vc,ofs_EDS);
					cogru.clear();
				}
				string next_chrom = copy_till_next_chrom(ifs_FA,ofs_EDS, current_fa_pos, v.chrom);
				current_fa_chrom =  next_chrom;
				current_fa_pos	= 0;
				ofs_EDS << EDS_NODE_BOUNDARY;
				ref_ind_buf << 0;					
				alt_ind_buf << 1;
				chrom_names.push_back(current_fa_chrom);
				chrom_starts.push_back(ofs_EDS.tellp());
				edges.push_back( vector<edge_type>() );
				positions.push_back( set<uint64_t>() );
				gedmap_io::print_row("SEQUENCE: '" + current_fa_chrom + "' (at pos " + to_string( ofs_EDS.tellp()  ) +" of FA)");
			}
			applied_cnv++;
			uint32_t type	= v.getCN();
			uint64_t pos	= v.pos;
			uint64_t end	= v.sv_end();
			
			edges.back().push_back( make_pair(pos,pos));
			if(type & vcfl::COPY_0) edges.back().push_back( make_pair(pos,end) );
			if(type & vcfl::COPY_2) edges.back().push_back( make_pair(end,pos) );
			edges.back().push_back( make_pair(end,end));
			positions.back().insert(pos);
			positions.back().insert(end);
			
		}else {
			//IGNORE
			log << "UNSUPPORTED VCF ENTRY -> IGNORE line " << read_lines << ": " << v << endl;
		}
	}
	
	if(!cogru.empty()){
		vcfl vc = cogru.solve();
		copy_till_pos(vc.pos-1,current_fa_pos,ifs_FA,ofs_EDS);
		copy_till_pos(current_fa_pos+vc.ref.size(),current_fa_pos,ifs_FA,ofs_EDS,true);
// 		if(vc.isSNP()) 		write_snp(vc,ofs_EDS);						
// 		else 
		if (vc.isINDEL())	write_indel(vc,ofs_EDS);
 		else	write_alt(vc,ofs_EDS);
	}
	
	try{copy_till_end(ifs_FA,ofs_EDS,current_fa_pos);}
	catch(invalid_argument &ia){}
	
	ofs_EDS << endl;
	ofs_EDS.close();
	
	
	ifs_FA.close();
	ifs_VCF.close();
	gedmap_io::print_row("Number of included short variations:",applied_alt);
	gedmap_io::print_row("Number of CNVs to include:",applied_cnv);
	
}

}
