#pragma once
//#define COUT(var) std::cout << #var << "=" << var <<std::endl
//#define TIME double(clock() - time) / CLOCKS_PER_SEC
//#define TAB '\t'

template<uint8_t t_width>
sdsl::int_vector_buffer<t_width>& operator <<(sdsl::int_vector_buffer<t_width>& ivb, const uint32_t x){ 
	ivb.push_back(x);
	return ivb;
}

/** @brief prints vector for debugging*/
template <uint8_t t_width>
std::ostream& operator <<(std::ostream& os, const sdsl::int_vector<t_width> vec){
	os << '<';
	for (uint32_t i = 0; i < vec.size(); i++ )
		os << ((uint32_t) vec[i]) << " ";
	os << '>';
	return os;
}

/** @brief prints vector for debugging*/
//template <class T>
//std::ostream& operator <<(std::ostream& os, const std::vector<T> & vec){
//	os << '<';
//	for (auto it = vec.begin() ; it != vec.end(); ++it)
//		os << (it==vec.begin()?"":",") << *it;
//	os << '>';
//	return os;
//}

template<class T>
string to_string(const std::vector<T> & vec){
	std::string out = "";
	for (auto it = vec.begin() ; it != vec.end(); ++it)
		out += (it==vec.begin()?"":",") + to_string(*it);
	return out;
}

std::vector<uint32_t> parse_uint32_vector(std::string s){
	std::vector<uint32_t> v(0);

	int begin 	= 0;
	int length	= 0;
	for(unsigned int i = 0; i < s.size(); i++)
		if(s[i] != ',')
			length++;
		else{    //if "," add number to v
			v.push_back( stoi(s.substr(begin,length)));
			begin 	= i+1;
			length 	= 0;
		}
	v.push_back(stoi(s.substr(begin,length)));
	return v;
}

namespace gedmap_io{
	using namespace std;
//line width for output
const int IO_LINE_WIDTH = 55;


void dotline(){
	cout << '\r' << setfill('.') << setw(IO_LINE_WIDTH) << '.' << setfill(' ')  << endl;
}
void print_error(string message){
	cerr << "\rERROR: "  << message << endl;
}

void print_row(string left, double rigth, string unit = ""){
	cout << '\r' << left
	<< setw(IO_LINE_WIDTH-left.length()-unit.length()) << fixed << setprecision(2) << rigth << unit << endl;
}

void print_row(string left, unsigned int rigth){
	cout << '\r' << left	<< setw(IO_LINE_WIDTH-left.length())	<< rigth << endl;
}

void print_row(string left, unsigned long int rigth){
	cout << '\r' << left	<< setw(IO_LINE_WIDTH-left.length())	<< rigth << endl;
}

void print_row(string left, string rigth=""){
	cout << '\r' << left 	<< setw(IO_LINE_WIDTH-left.length()) 	<< rigth << endl;
}

void flush_row(string left, string rigth=""){
	cout << '\r' << left 	<< setw(IO_LINE_WIDTH-left.length())	<< rigth << flush;
}

/** @brief flushes a progress bar*/
void flush_progress(double percent){
	int count = (IO_LINE_WIDTH-2)*percent;
	cout << "\r[" << string(count,'#') << string(IO_LINE_WIDTH-count-2,' ') << "]" << flush;
}




// ostream& operator <<(ostream& os, const vector<string> vec){
// 	os << '<';
// 	for (auto it = vec.begin() ; it != vec.end(); ++it)
// 		os << (it==vec.begin()?"":",") << *it;
// 	os << '>';
// 	return os;
// };

/**
 *	@brief checks if file exist
 * 	return false if it exists
 * 	returns true if not and prints Error message
 */
bool file_access(string filename){
	if(access((filename).c_str(), F_OK)){
		print_error("File " + filename + " not found");
		return true;
	}
	return false;
};

/**
 *	@brief removes file
 * 	return false if remove works
 * 	returns true if not and prints Error message
 */
bool file_remove(string filename){
	if(access((filename).c_str(), F_OK)){
		print_error("File " + filename + " not found");
		return true;
	}
	if(remove(filename.c_str())){
		print_error("File " + filename + " could not be removed");
		return true;
	}
	return false;
};

/**
 *	@brief removes file
 * 	return false if remove works
 * 	returns true if not
 */
bool file_remove_silent(string filename){
	if(access((filename).c_str(), F_OK)) return true;
	if(remove(filename.c_str())) return true;
	return false;
};

/**
 *	@brief copies file
 * 	return false if no error
 * 	returns true if not
 */
bool file_copy(string filename_from, string filename_to){
	std::error_code ec;
	std::filesystem::copy(filename_from, filename_to,ec);
	if(ec.value()){
		print_error("Copy from " + filename_from + " to " +filename_to + " (filesystem::copy error_code="+std::to_string(ec.value())+")");
		return true;
	}
	return false;
};

/**
 *	@brief moves file
 * 	return false if no error
 * 	returns true if not
 */
bool file_move(string filename_from, string filename_to){
	std::error_code ec;
	std::filesystem::rename(filename_from, filename_to,ec);
	if(ec.value()){
		print_error("Move from " + filename_from + " to " +filename_to + " (filesystem::copy error_code="+std::to_string(ec.value())+")");
		return true;
	}
	return false;
};




/**
 *	@brief loads sdsl object
 * 	return false if everithing is ok
 * 	returns true if not and prints Error message
 */
template <class T>
bool sdsl_load(string filename, T & obj){
	bool error = false;

	if(file_access(filename)){
		error = true;
	}else	if( !sdsl::load_from_file<T>(obj,filename)){
		print_error(" Object could not be loaded from " + filename);
		error = true;
	}
	return error;
};

/**
 *	@brief loads string from file
 * 	return false if everithing is ok
 * 	returns true if not and prints Error message
 */
bool load_string_from_file(string  filename, string & buffer){
	if(file_access(filename))
		return true;

	//load jis from file
	ifstream stream;
	stream.open(filename);
	stream >> buffer;
	stream.close();
	return false;
};

bool ends_with(std::string const & s, std::string const & ending){
	if (ending.size() > s.size()) return false;
	return(s.compare(s.size()-ending.size(),ending.size(),ending)==0);
}


void print_prog_headline(string prog){
	dotline();
	print_row(prog);
	dotline();
}

//################### ERROR MESSAGES ##########################
const string missing_arguments	= "Missing arguments";
const string not_loaded	= "Input could not be loaded";
string unknown_argument (string arg){
	return "Unknown argument: " + arg;
}
string missing_value (string arg){
	return "Value for argument '" + arg + "' missing";
}

//################### HELP MESSAGES ##########################
void print_short_help(string bin_path){
	dotline();
	print_row("GEDMAP (Graph-ED-string-MAPper) is a program collection around read mapping to a graph of ED stings.");
	print_row("version 1.0");
	dotline();
	cout
	 << "The containd programs are: \n - gedmap parse\n - gedmap index\n - gedmap align\n - gedmap sample\n";
	dotline();
	 cout << "call '"  << bin_path << " program_name -h' to print the manual of the program."
	 << "\ncall '"  << bin_path << " -h' to print the manuel of all programs.\n";
}

string expected_arguments(vector<string> &args, vector<string> & parameters){
	stringstream ss;
	ss  << args.size() << " Arguments expected:"  << endl ;
	for(uint32_t i = 0; i < args.size(); i++)
		ss << "[" << (i+1) << "] " << args[i] << endl ;
	ss << '\r' << setfill('.') << setw(IO_LINE_WIDTH) << '.' << setfill(' ')  << endl;
	if(parameters.size() != 0)
		ss << "Optional parameters: " << endl;
	for(uint32_t i = 0; i < parameters.size(); i++)
		ss << parameters[i] << endl;
	ss << endl;
	return ss.str();
}

void print_call_for_help(string bin_path){
	cout << "call '" <<  bin_path << " -h' for help message" << endl;
}

//######################## INPUT HANDLER #############################





}



namespace gedmap_parse{
using namespace gedmap_io;
void print_help(){
	vector<string> args { "filename of FA", "filename of VCF" , "output filename"};
	vector<string> params {
		"-nosv, do not include structural variants, i.e. copy number variation and other large variation (Then only and EDS is generated)",
		"-tmp tmp_dir , to set tmp direcoty (DEFAULT tmp_dir="+ TMP_DIR_DEFAULT +")",
		"-lim l , variants with ref or alt larger than l are handled as structural variants (DEFAULT l="+ to_string(PLAIN_ALT_LIMIT_DEFAULT) +")"};

		cout << "'gedmap parse' parses a FA file and a VCF file to an EDS graph." << endl;
		dotline();
		cout << expected_arguments(args, params);
}


/**
 * \param [2] filename of FA
 * \param [3] filename of VCF
 * \param [4] filename of output EDS
 */
void handle_input(int argc, char**& argv, string& fname_fa, string& fname_vcf, string& fname_eds){
	try{

		if(argc < 3) throw invalid_argument( "in gedmap parse: " + missing_arguments);
		if(string(argv[2]) == "-h" || string(argv[2]) == "--help"){
			print_help();
			exit(0);
		}

		if(argc < 5) throw invalid_argument( "in gedmap parse: " + missing_arguments);

		fname_fa  = argv[2];
		fname_vcf = argv[3];
		fname_eds = argv[4];

		if(!ends_with(fname_eds, FEX_EDS))
			fname_eds += "." + FEX_EDS;
		bool input_error = false;
		input_error |= file_access(fname_fa);
		input_error |= file_access(fname_vcf);

		if(input_error) throw invalid_argument("in gedmap parse: " + not_loaded);


		for(int i = 5; i < argc;){
			string param = argv[i++];
			if("-nosv" == param){
				gedmap_parse::INCLUDE_SV = false;
				continue;
			}

			if(i >= argc)	throw invalid_argument("in gedmap parse: " + missing_value(param));
			string value = argv[i++];
			if("-tmp" == param)
				TMP_DIR = value;
			else if("-lim" == param)
				PLAIN_ALT_LIMIT = stoi(value);
			else
				throw invalid_argument("in gedmap parse: " + unknown_argument(param));
		}
	}
	catch(exception & e){
		print_error(e.what());
		print_call_for_help(string(argv[0]) + " " + argv[1]);
		exit(1);
	}
}
}

namespace gedmap_index_min{
using namespace gedmap_io;
void print_help(){
	vector<string> args { "filename of GEDS"};
	vector<string> params {
		"-a fname,  file name of the adijacency file",
		"-2fa fname,  file of the 2fa file",
		"-k k,  kmer/minimizer size k (DEFAULT k="+ to_string(KMER_SIZE_DEFAULT) +")",
		"-w w,  window size w (DEFAULT w="+ to_string(WINDOW_SIZE_DEFAULT) +")",
		"-o fname,  file name of minimizer index (DEFAULT fname= geds_fname.min)",
		"-n n,  maximum number of N in minimizer (DEFAULT n="+ to_string(MAX_N_IN_SEED_DEFAULT) +")",
		"-tc x , maximum number of threads used (DEFAULT uses as many as avaiable)",
		"-trim x , removes sets greater than x from index (DEFAULT x="+ to_string(TRIM_DEFAULT) + " (x=0 for no trim))"};
	cout << "'gedmap index' calculates a minimizer index of a given EDS graph." << endl;
	dotline();
	cout << expected_arguments(args, params);
}

/**
 * \param [2] filename of EDS
 */
void handle_input(int argc, char**& argv, string & eds, string &  fname_min, adjacency & adj, pos_EDS_to_FA_type & p2FA){
	try{
		if(argc < 3) throw invalid_argument( "in gedmap index: " + missing_arguments);
		if(string(argv[2]) == "-h" || string(argv[2]) == "--help"){
			print_help();
			exit(0);
		}

		if(load_string_from_file(argv[2],eds)) throw invalid_argument("in gedmap index: geds" + not_loaded);
		fname_min = argv[2];

		for(int i = 3; i < argc;){
			string param = argv[i++];
			if(i >= argc)	throw invalid_argument("in gedmap index: " + missing_value(param));
			string value = argv[i++];
			if("-a" == param){
				if(sdsl_load(value, adj)) throw invalid_argument("in gedmap index: (-a)" + not_loaded);
			}
			else if("-2fa" == param){
				if(sdsl_load(value, p2FA)) throw invalid_argument("in gedmap index: 2fa" + not_loaded);}
			else if("-k" == param)
				KMER_SIZE = stoi(value);
			else if("-w" == param)
				WINDOW_SIZE = stoi(value);
			else if("-o" == param)
				fname_min = value;
			else if("-n" == param)
				MAX_N_IN_SEED = stoi(value);
			else if("-tc" == param)
				MAX_THREAD_COUNT = stoi(value);
			else if("-trim" == param)
				TRIM = stoi(value);
			else
				throw invalid_argument("in gedmap index: " + unknown_argument(param));
		}

		if(!ends_with(fname_min, FEX_MIN))
			fname_min += "." + FEX_MIN;

		if(KMER_SIZE > 32) invalid_argument("in gedmap index: k has a max value of 32)");

	}catch(std::exception& e){
		print_error(e.what());
		print_call_for_help(string(argv[0]) + " " + argv[1]);
		exit(1);
	}
}
}

namespace gedmap_align_min{
using namespace gedmap_io;
void print_help(){
	using namespace gedmap_align;
	vector<string> args { "filename of GEDS", "filename of FASTQ" , "filename of MINI"};
	vector<string> params {
	"-o            , fname, output will be stored in file fname (DEFAULT = [2]."+FEX_SAM+")",
	"-mp           , filename of FASTQ containing the mates (optional, presence indicates paired-end mode)",
	"-fragment-mean, mean length of fragment in paired-end mode (ignored if -mp is not present, DEFAULT=" + std::to_string(PE_FRAGMENT_LENGTH) + ")",
	"-fallback     , fallback for paired-end mapping (ignored if -mp is not present)",
	"-fmat         , maximum number of alignments tried for fallback (ignored if -mp and -fallback are not present, DEFAULT=" + std::to_string(MAX_ALIGNS_T_FALLBACK_DEFAULT) + ")",
	"-mam x        , max number of alignments used for pairing (DEFAULT x="	+to_string(MAX_ALIGNS_M_DEFAULT)+")",
	"-2fa          , .2fa-file , if given this is used to transform GEDS-positions to FA positions",
	"-a fname      , file name of the adijacency file",
	"-rc           , reversed complement of pattern will be searched, too",
	"-oa           , only aligned reads will be reported",
	//"-nc           , don't check for collinearity, when rankin seeds",
	"-io           , output reads in the same order as in the input (may be a bit slower and with higher memory)",
	"-mc x         , minimizer count, x minimizers will be looked up per read(DEFAULT x="+to_string(FRAGMENT_COUNT_DEFAULT)+")",
	"-ws x         , window size of hotspot (DEFAULT x="+to_string(SPOT_SIZE_DEFAULT)+")",
	"-wh x         , minimum minimizer score (DEFAULT x="+to_string(SPOT_HITS_DEFAULT)+")",
	"-dd x         , doubt distance, when best alignment has a distance x or greater, go to next round (DEFAULT x="+to_string(DOUBT_DIST_DEFAULT)+")",
	"-mac x        , max number of alignments completely calculated (DEFAULT x="	+to_string(MAX_ALIGNS_C_DEFAULT)+")",
	"-mat x        , max number of alignments tried to calculate (DEFAULT x="	+to_string(MAX_ALIGNS_T_DEFAULT)+")",
	"-mao x        , max number of alignments in output (DEFAULT x="		+to_string(MAX_ALIGNS_O_DEFAULT)+")",
	"-d x          , max distance in alignment (DEFAULT x="			+to_string(MAX_DIST_DEFAULT)+")",
	"-tmp tmp_dir  , to change DEFAULT tmp direcoty from /tmp to tmp_dir",
	"-tc x         , maximum number of threads used per index copy (DEFAULT uses as many as avaiable)"
	};

	cout << "'gedmap align' algings reads to the given GEDS and MINIMIZER INDEX" << endl;
	dotline();
	cout << expected_arguments(args, params);
	dotline();
	cout << "The program can run multiple rounds with different parameters per read." << endl;
	cout << "Therefore, the parametes mc,d,mac,mat can be given as a comma seperated list of equal length." << endl;
	cout << "If there was no alignment with distance smaller then -dd was calculated in the current round, the program performs the next round." << endl;
	cout << "F.e.: -d 3,10 means the alignment algorithm allows only 3 errors in the first round and 10 errors in the second round."  << endl;	
	dotline();
}

/**
 * \param [2] filename of EDS
 * \param [3] filename of FASTQ
 * \param [4] filename of EOC
 */
void handle_input(int argc,  char**& argv, gedmap_mini::minimizer_index & eoc, string & EDS, adjacency & ADJ, pos_EDS_to_FA_type & p2FA, ifstream & fastq_s, ifstream& fastq_s2, ofstream & o_s){
	using namespace std;

	string 	fname_EDS;
	string 	fname_fastq;
	string 	fname_eoc;
	string 	fname_sam = "";

	try{
		if(argc < 3) throw invalid_argument( "in gedmap align_min: " + missing_arguments);
		if(string(argv[2]) == "-h" || string(argv[2]) == "--help"){
			print_help();
			exit(0);
		}

		if(argc < 5) throw invalid_argument( "in gedmap align_min: " + missing_arguments);

		fname_EDS		= argv[2];
		fname_fastq		= argv[3];
		fname_eoc		= argv[4];
		bool input_error = false;

		for(int i = 5; i < argc;) {
			string param = argv[i++];
			if("-rc" == param){
				MAP_RC = true;
				continue;
			}else if("-oa" == param) {
				WRITE_FAILURE = false;
				continue;
			}else if("-nc" == param) {
				CHECK_COLLI = false;
				continue;
			}else if("-io" == param) {
				IN_ORDER = true;
				continue;
			}else if("-fallback" == param) {
				FALLBACK = true;
				continue;
			}

			if(i >= argc) throw invalid_argument("in gedmap index: " + missing_value(param));

			string value = argv[i++];

			if ("-mp" == param) { // paired_end
				fastq_s2.open(value);
				if (!fastq_s2.is_open()) {
					throw invalid_argument("in gedmap align: could not open mate pair files -mp " + param);
				}
			}
			else if ("-fragment-mean" == param)
				PE_FRAGMENT_LENGTH = stoi(value);
			else if ("-tmp" == param)
				TMP_DIR = value;
			else if ("-o" == param)
				fname_sam = value;
			else if("-mc" == param)
				FRAGMENT_COUNT = parse_uint32_vector(value);
			else if("-ws" == param)
				SPOT_SIZE = stoi(value);
			else if("-wh" == param)
				SPOT_HITS = parse_uint32_vector(value);
			else if("-dd" == param)
				DOUBT_DIST = stoi(value);
			else if("-mac" == param)
				MAX_ALIGNS_C = parse_uint32_vector(value);
			else if("-mam" == param)
				MAX_ALIGNS_M = parse_uint32_vector(value);
			else if("-mat" == param)
				MAX_ALIGNS_T = parse_uint32_vector(value);
			else if("-fmat" == param)
				MAX_ALIGNS_T_FALLBACK = stoi(value);
			else if("-mao" == param)
				MAX_ALIGNS_O = stoi(value);
			else if("-d" == param)
				MAX_DIST = parse_uint32_vector(value);
			else if("-2fa" == param)
				input_error |= sdsl_load(value, p2FA);
			else if("-tc" == param)
				THREAD_COUNT = stoi(value);
			else if("-a" == param)
				input_error |= sdsl_load(value, ADJ);
			else
				throw invalid_argument("in gedmap align: " + unknown_argument(param));
		}

		input_error |= load_string_from_file(fname_EDS, EDS);
		input_error |= sdsl_load(fname_eoc, eoc);
		input_error |= file_access(fname_eoc);

		if(input_error) throw invalid_argument("in gedmap align: " + not_loaded);

		if(SPOT_HITS.size() != FRAGMENT_COUNT.size()
			|| SPOT_HITS.size() != MAX_DIST.size()
			|| SPOT_HITS.size() != MAX_ALIGNS_C.size()
			|| SPOT_HITS.size() != MAX_ALIGNS_T.size())
			throw invalid_argument("conficting params: -mc -d -wh -mat -mac (none equal list size)");

		if(fname_sam.size() == 0)
			fname_sam = fname_fastq + "." +  FEX_SAM;
		fastq_s.open(fname_fastq);
		o_s.open(fname_sam);

	}catch(std::exception& e){
		print_error(e.what());
		print_call_for_help(string(argv[0]) + " " + argv[1]);
		exit(1);
	}
}
}

namespace gedmap_sample{
using namespace gedmap_io;

void print_help(){
	using namespace gedmap_sample;

	vector<string> args { "filename of GEDS"};
	vector<string> params {
		"-c x , sample x reads (count, DEFAULT x=" + to_string(COUNT_DEFAULT) +")",
		"-l x , reads length = x (length, DEFAULT x=" + to_string(LENGTH_DEFAULT) +")",
		"-e x , probability of a base to be false = 1/x (error rate, DEFAULT x=" + to_string(ERROR_RATE_DEFAULT) +")",
		"-e_d x , probability of a base indel = 1/x, (insertion rate, DEFAULT x=" + to_string(DEL_RATE_DEFAULT) +")",
		"-e_i x , probability of a base indel = 1/x, (deletion rate, DEFAULT x=" + to_string(IN_RATE_DEFAULT) +")",
		"-rc x , probability of a read to be a reverse complement = 1/x, (rev complement rate, DEFAULT x=" + to_string(RC_RATE_DEFAULT) +")",
		"-s x , seed for rng=x, (seed, DEFAULT x=0)",
		"-o fname , output will be written to file fname (DEFAULT  fname=[1]."+ FEX_SAMPLE +")",
		"-fragment-mean x , enables paired-end mode, x is length of fragment from which the reads-pairs come from (second file will be {argument to \"-o\"}.2)",
		"-2fa "+FEX_2FA+"-file , if given this is used to transform GEDS-positions to FA positions",
		"-a fname,  file name of the adijacency file",
		"-ae fname,  file name of the adijacency file, only sample reads that go over edges in the graph"
	};

	cout << "'gedmap sample' samples reads from the given GEDS." << endl;
	dotline();
	cout << expected_arguments(args, params);
	print_row("for all probabilities above:");
	print_row(" x : prob (1/x)");
	print_row(" 0 : 0");
	print_row(" 1 : 1");
	print_row(" 2 : .5");
	print_row(" 3 : .33");
	print_row(" etc.");
	dotline();
}


/**
 * \param [2] filename of EDS
 */
void handle_input(int argc,  char**& argv, std::string & EDS, adjacency & adj, pos_EDS_to_FA_type & p2FA, ofstream & o_s, ofstream& o_s_mp, string & fname_out){
	using namespace std;
	using namespace gedmap_sample;

	string 	fname_EDS;

	try{
		if(argc < 3) throw invalid_argument( "in gedmap sample: " + missing_arguments);
		if(string(argv[2]) == "-h" || string(argv[2]) == "--help"){
			print_help();
			exit(0);
		}

		string fname_EDS = argv[2];
		fname_out = fname_EDS + "." +  FEX_SAMPLE;
		bool input_error = false;

		for(int i = 3; i < argc;){
			std::string param = argv[i++];
			if(i >= argc) throw invalid_argument("in gedmap sample: " + missing_value(param));
			std::string value = argv[i++];

			if ("-fragment-mean" == param)
				FRAGMENT_LENGTH = stoi(value);
			else if ("-c" == param)
				COUNT = stoi(value);
			else if ("-l" == param)
				LENGTH = stoi(value);
			else if ("-e" == param)
				ERROR_RATE = stod(value);
			else if ("-e_d" == param)
				DEL_RATE = stoi(value);
			else if ("-e_i" == param)
				IN_RATE = stoi(value);
			else if ("-rc" == param)
				RC_RATE = stoi(value);
			else if ("-s" == param)
				SEED = stoi(value);
			else if ("-o" == param)
				fname_out = value;
			else if("-2fa" == param)
				input_error |= sdsl_load(value, p2FA);
			else if("-a" == param)
				input_error |= sdsl_load(value, adj);
			else if("-ae" == param){
				input_error |= sdsl_load(value, adj);
				over_edges = true;
			}
			else
				throw std::runtime_error("Unkown argument: " +  param);

		}
		input_error |= load_string_from_file(fname_EDS, EDS);
		o_s.open(fname_out);
		if (FRAGMENT_LENGTH < std::numeric_limits<uint32_t>::max())
			o_s_mp.open(fname_out + ".2");
		if(input_error) throw invalid_argument("in gedmap align: " + not_loaded);

	}catch(std::exception& e){
		print_error(e.what());
		print_call_for_help(string(argv[0]) + " " + argv[1]);
		exit(1);
	}

}
}

namespace gedmap_io{
	void print_long_help(string bin_path){
		print_short_help(bin_path);
		print_row("");
		print_row("");
		gedmap_parse::print_help();
		print_row("");
		print_row("");
		gedmap_index_min::print_help();
		print_row("");
		print_row("");
		gedmap_align_min::print_help();
		print_row("");
		print_row("");
		gedmap_sample::print_help();
	}
}
