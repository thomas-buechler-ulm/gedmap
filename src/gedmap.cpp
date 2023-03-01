//STD
#include <string>
#include <vector>
#include <sstream>  
#include <deque>
#include <tuple>
#include <fstream>
#include <bits/stdc++.h> //explizit vector print help
#include <chrono>
#include <filesystem> 	//COPY FILE
#include <algorithm>	// std::sort
#include <omp.h>

//SDSL
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v5.hpp>

//GEDMAP HEADERS
#include "default_values.hpp"
#include "lib/include.hpp"

//EDSM PROGRAMS
#include "programs/include.hpp"

using namespace std;
using namespace sdsl;

int main(int argc,  char** argv){
	try{
		if(argc == 1) {gedmap_io::print_short_help(argv[0]); exit(1);}
		string prog =  argv[1];
		if(prog == "-h" || prog == "--help"){
			gedmap_io::print_long_help(argv[0]);
			exit(0);
		}

		if(prog == "parse")	return gedmap_parse	::main(argc, argv);
 		if(prog == "index")	return gedmap_index_min	::main(argc, argv);
		if(prog == "align")	return gedmap_align_min	::main(argc, argv);
		if(prog == "sample")	return gedmap_sample	::main(argc, argv);
		
		throw std::invalid_argument("in gedmap: " + gedmap_io::unknown_argument(prog));
	}catch(std::exception& e){		
		gedmap_io::print_error(e.what());
		gedmap_io::print_call_for_help(argv[0]);
		exit(1);
	}
	return 0;
}
