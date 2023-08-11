#ifndef EDSM_DEFAULT_VAL_DEF
#define EDSM_DEFAULT_VAL_DEF

//##############  FILE EXTENSIONS ####################
const std::string FEX_EDS		= "geds";
const std::string FEX_2FA 		= "2fa";
const std::string FEX_2GFA 		= "2gfa";
const std::string FEX_ADJ 		= "adj";
const std::string FEX_EOC 		= "eoc";
const std::string FEX_MIN 		= "min";
const std::string FEX_SAM 		= "sam";
const std::string FEX_LOG 		= "log";
const std::string FEX_SAMPLE 		= "sample.fastq";


namespace gedmap_io{ const std::string VERSION="1.1"; }

//##############  TEMPORY FILES ####################
const std::string TMP_DIR_DEFAULT	= "/tmp";
std::string TMP_DIR 			= TMP_DIR_DEFAULT;


//##############  PARAMETERS ####################
namespace gedmap_parse{
	const bool INCLUDE_SV_DEFAULT	= true;
	bool INCLUDE_SV = INCLUDE_SV_DEFAULT;
	const std::string TMP_PARSE_REF_IND	= "/ref_ind.ivb.sdsl";
	const uint32_t PLAIN_ALT_LIMIT_DEFAULT = 50;
	uint32_t PLAIN_ALT_LIMIT = PLAIN_ALT_LIMIT_DEFAULT;
}

namespace gedmap_parse_gfa{
	bool print = false;
	uint32_t BUB_max_length = 50;
	uint32_t BUB_max_path_c = 20;

}

namespace gedmap_index_min{
	const uint32_t WINDOW_SIZE_DEFAULT		= 5;
	const uint32_t KMER_SIZE_DEFAULT		= 20;
	const uint32_t MAX_N_IN_SEED_DEFAULT 	= 2;
	const uint32_t TRIM_DEFAULT		 	= 1000;

	uint32_t WINDOW_SIZE				= WINDOW_SIZE_DEFAULT;
	uint32_t KMER_SIZE				= KMER_SIZE_DEFAULT;
	uint32_t MAX_N_IN_SEED 				= MAX_N_IN_SEED_DEFAULT;
	uint32_t TRIM	 				= TRIM_DEFAULT;
	uint32_t MAX_THREAD_COUNT			= omp_get_max_threads();

}
//ALIGN
//SETTINGS with Default values
namespace gedmap_align_min{
	const std::vector<uint32_t> FRAGMENT_COUNT_DEFAULT{80};
	const uint32_t SPOT_SIZE_DEFAULT		= 500;
	const std::vector<uint32_t> SPOT_HITS_DEFAULT{1};
	const uint32_t DOUBT_DIST_DEFAULT		= 7;
	const std::vector<uint32_t> MAX_DIST_DEFAULT{30};
	const std::vector<uint32_t> MAX_ALIGNS_C_DEFAULT = {5};
	const std::vector<uint32_t> MAX_ALIGNS_M_DEFAULT = {1,10};
	const std::vector<uint32_t> MAX_ALIGNS_T_DEFAULT= {10};
	const uint32_t MAX_ALIGNS_T_FALLBACK_DEFAULT = 100;
	const uint32_t MAX_ALIGNS_O_DEFAULT		= 1;
	uint32_t BATCH_SIZE_IN_ORDER			= 10000;
	uint32_t MAX_PATHS_IN_ALIGNMENT		= 1024;
	const uint32_t THREAD_COUNT_DEFAULT		= omp_get_max_threads();
	std::vector<uint32_t> FRAGMENT_COUNT	= FRAGMENT_COUNT_DEFAULT;
	//uint32_t BATCH_SIZE				= BATCH_SIZE_DEFAULT;
	uint32_t SPOT_SIZE				= SPOT_SIZE_DEFAULT;
	std::vector<uint32_t> SPOT_HITS		= SPOT_HITS_DEFAULT;
	uint32_t DOUBT_DIST				= DOUBT_DIST_DEFAULT;
	std::vector<uint32_t>  MAX_DIST		= MAX_DIST_DEFAULT;
	std::vector<uint32_t> MAX_ALIGNS_C		= MAX_ALIGNS_C_DEFAULT;
	std::vector<uint32_t> MAX_ALIGNS_M		= MAX_ALIGNS_M_DEFAULT; // number of alignments used for pairing
	std::vector<uint32_t> MAX_ALIGNS_T		= MAX_ALIGNS_T_DEFAULT;
	uint32_t MAX_ALIGNS_T_FALLBACK		= MAX_ALIGNS_T_FALLBACK_DEFAULT;
	uint32_t MAX_ALIGNS_O				= MAX_ALIGNS_O_DEFAULT;
	uint32_t THREAD_COUNT				= THREAD_COUNT_DEFAULT;
	bool IN_ORDER					= false;
	bool WRITE_FAILURE				= true;
	bool MAP_RC						= false;
	bool CHECK_COLLI					= true;
	uint32_t PE_FRAGMENT_LENGTH = 700; // default fragment length for paired end mapping
	bool FALLBACK = false;
}

//SAMPLE
namespace gedmap_sample{
	const uint32_t COUNT_DEFAULT 		= 100;
	const uint32_t LENGTH_DEFAULT 	= 100;
	const uint32_t ERROR_RATE_DEFAULT	= 100;
	const uint32_t DEL_RATE_DEFAULT 	= 1000;
	const uint32_t IN_RATE_DEFAULT 	= 1000;
	const uint32_t RC_RATE_DEFAULT 	= 2;
	const uint32_t SEED_DEFAULT 		= 0;
	
	bool 	over_edges				= false;

	uint32_t FRAGMENT_LENGTH = std::numeric_limits<uint32_t>::max(); // -fragment-mean
	uint32_t COUNT 		= COUNT_DEFAULT; //-c
	uint32_t LENGTH		= LENGTH_DEFAULT; //-l
	uint32_t ERROR_RATE	= ERROR_RATE_DEFAULT; //-e
	uint32_t DEL_RATE 	= DEL_RATE_DEFAULT; //-e_d
	uint32_t IN_RATE 		= IN_RATE_DEFAULT; //-e_i
	uint32_t RC_RATE 		= RC_RATE_DEFAULT;
	uint32_t SEED 		= SEED_DEFAULT;
}


//LEVINSTIN
namespace edsm_levinstein{	
	const uint32_t STRING_BUFFER_SIZE = 500;
	const bool BACKWARD 	= true;
	const bool FORWARD 	= false;
}

//ALIGN
namespace gedmap_align{	
	const uint32_t DEFAULT_POSPAIR_CAPACITY = 800000;
}

// seperates chromosoms in eds
const char EDS_NODE_BOUNDARY	= '#';
// seperates positions in EOC
const int EOC_SEPERATOR = -1;

const uint32_t KB = 1024;
const uint32_t MB = 1024*1024;
const uint32_t GB = 1024*1024*1024;


/**
 * @brief caclulates how many bits are needed to represent x
 * @param x a integer number  
 * @return  floor(log_2(x))+1 (1 if x = 0)
 */
uint32_t bitsneeded (uint64_t x) {
	uint32_t ret = 1; 
	while( 1 < x){
		x  >>= 1;
		ret  += 1;    
	}
	return ret;
};






#endif

