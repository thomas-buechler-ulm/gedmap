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

using namespace sdsl;


size_t med(vector<uint32_t> & v, int count){
	size_t med = 0;
	while(count > 0) {
		if(med >= v.size()) return -1;
		count -= v[med++];
	}
	if(med) med--;
	return med;
}

void analyse_kmi(char **argv, gedmap_mini::minimizer_index & mini){
	cout << argv[1] << '\t'
	     << mini.k << '\t'
	     << mini.w << '\t'
	     << mini.pos_sing.size() << '\t'
	     << mini.pos_mult.size() << '\t'
	     << mini.mult.size() << '\t'
	     << (mini.mult.size()/(double)mini.contained.size())*100 << '\t'
	     << ((mini.pos_sing.size()+mini.pos_sing.size())/ (double) (mini.mult.size())) << '\t'
	     << (size_in_bytes(mini.pos_sing)) << '\t'
	     << (size_in_bytes(mini.pos_mult)) << '\t'
	     << (size_in_bytes(mini.start)) << '\t'
	     << (size_in_bytes(mini.contained)+size_in_bytes(mini.contained_rs)) << '\t'
	     << (size_in_bytes(mini.mult)+size_in_bytes(mini.mult_rs)) << '\t'
	     << (size_in_bytes(mini.pos_sing) + size_in_bytes(mini.pos_mult) + size_in_bytes(mini.start)+size_in_bytes(mini.contained) + size_in_bytes(mini.mult)+size_in_bytes(mini.contained_rs)+size_in_bytes(mini.mult_rs));
}


vector<string> string_split(string & s, char delim){
	vector<string> out(0);
	size_t l = 0;
	size_t r;
	while(l < s.size()){
		r = s.find(delim,l);
		if(r == string::npos){
			out.push_back(s.substr(l));
			return out;
		}else{
			out.push_back(s.substr(l,r-l));
			l = r+1;
		}
	}
	throw runtime_error ("Strange things in my string_split implementation\n");
}

pair<uint32_t,uint32_t> count_hints(gedmap_mini::minimizer_index & mini, fasta_read<uint64_t>  & read  ){
	uint32_t correct = 0;
	uint32_t all = 0;
	vector<string> fields = string_split(read.id, '_');
	uint64_t starting_pos = stol(fields[3]);
	//read.get_fragments(read.sequence.size(), mini);
	auto pos_pairs = read_processor::get_positions(read_processor::get_fragments(read, 100, mini),mini);
	for(auto pp : pos_pairs){
		all++;
		uint64_t est_start = pp.eds_pos - pp.read_pos + 1;
		if( (starting_pos <= est_start) && (est_start - starting_pos) < 100)
			correct++;
	}
	return make_pair(correct,all);
}

void analyse_kmi_samples(char** argv, gedmap_mini::minimizer_index & mini){
	vector<uint32_t> c_hints;
	vector<uint32_t> hints;

	uint32_t read_count = 0;
	uint32_t c_sum = 0;
	uint64_t h_sum = 0;

	ifstream fq_in;
	fq_in.open(argv[2],ios_base::in);

	while(fq_in.good()){
		try{
			fasta_read<uint64_t>  read(fq_in);
			uint32_t ch,h;
			tie(ch,h) = count_hints(mini, read);
			c_sum += ch;
			h_sum += h;
			c_hints.push_back(ch);
			hints.push_back(h);
			read_count++;
		}catch(runtime_error & e){
		}
	}
	std::sort (c_hints.begin(), c_hints.end());
	std::sort (hints.begin(), hints.end());
	uint32_t no_correct = 0;
	uint32_t no_hint = 0;
	while(!c_hints[no_correct] && no_correct < c_hints.size()) no_correct++;
	while(!hints[no_hint] && no_hint < hints.size() ) no_hint++;

	cout << string(argv[2]) << '\t'
	<< c_hints[c_hints.size()/2] << '\t'
	<< hints[hints.size()/2] << '\t'
	<< (c_sum/read_count) << '\t'
	<< (h_sum/read_count) << '\t'
	<< no_correct << '\t'
	<< no_hint;
}

int main(int argc, char** argv){

	if(argc < 2){
		cout << "Argument expected" << endl << "[1] minimizer index" << endl << "([2] samples)" << endl; 
		return 0;
	}

	gedmap_mini::minimizer_index mini;
	
	if( argc == 2){
		cout << "in_file\tk\tw\tcount_positions_sing\tcount_positions_mult\tcount_kmers\tkmer_density\tavg_position_per_kmer\tsize_pos_sing\tsize_pos_mult\tsize_start\tsize_contained\tsize_mult\tsize_sum\n";
		if(!sdsl::load_from_file<gedmap_mini::minimizer_index>(mini,string(argv[1]))){
			cout << "could not load " << string(argv[1]) << '\t';
			return 0;
		}
		analyse_kmi(argv,mini);
		cout << endl;
	}
	else if( argc == 3){
		cout << "in_file\tk\tw\tcount_positions_sing\tcount_positions_mult\tcount_kmers\tkmer_density\tavg_position_per_kmer\tsize_pos_sing\tsize_pos_mult\tsize_start\tsize_contained\tsize_mult\tsize_sum\t";
		cout << "fq_fname\tmed_correct_hints\tmed_hints\tavg_correct_hints\tavg_hints\tcount_no_correct_hints\tcount_no_hint\n";
		if(!sdsl::load_from_file<gedmap_mini::minimizer_index>(mini,string(argv[1]))){
			cout << "could not load " << string(argv[1]) << '\t';
			return 0;
		}
		analyse_kmi(argv,mini);
		cout << '\t';
		analyse_kmi_samples(argv,mini);
		cout << endl;
	}
	else
		cout << "too many Arguments" << endl;
	return 0;
}
