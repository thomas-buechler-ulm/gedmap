#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>      // std::setw

using namespace std;


uint32_t abs(uint32_t x, uint32_t y){
	return (x < y) ? (y-x) : (x-y);
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

int main(int argc, char *argv[]){

	if(argc != 2)	cout << "1 argument (sam-file) expected" << endl; 

	ifstream sam_in;
	sam_in.open(argv[1],ios_base::in);

	//FOR FLAG
	uint32_t UNMAPPED 	= 4;
	uint32_t REV_COMP	= 16;
	uint32_t NOT_PRIM	= 256;
	uint32_t SUPPLEMENTARY	= 2048;

	uint32_t CORRECT_POSITION_TRESHOLD 	= 600;

	//COUNTERS
	uint32_t count_all 			= 0;
	uint32_t count_mapped 			= 0;
	uint32_t count_correct_pos		= 0;
	uint32_t count_correct_pos_better_eq	= 0;
	uint32_t count_correct_pos_worse	= 0;
	
	uint32_t count_other_pos_better 	= 0;
	uint32_t count_other_pos_equal 		= 0;
	uint32_t count_other_pos_worse 		= 0;

	uint32_t count_unmap			= 0;
	
	uint32_t NM_commulated			= 0;
	uint32_t NM_count			= 0;

	string entry;

	while(getline(sam_in, entry)){
		if(entry[0] == '@' || entry[0] == '[') continue;
		count_all++;

		vector<string> fields  = string_split(entry, '\t');

		uint32_t FLAG = stoi(fields[1]);

		if(FLAG & (NOT_PRIM + SUPPLEMENTARY)) continue; //ONLY LOOK AT PRIMARIES

		if(FLAG & UNMAPPED){ //UMAPPED
			count_unmap++;
			continue;
		}
		count_mapped++;
		
		vector<string> query_fields = string_split(fields[0],'_');
		
		uint64_t mapping_pos = stol( fields[3] );
		uint64_t mapping_pos_mate = fields[7]=="*"?0:stol( fields[7] );
		uint32_t mapping_dist = 9999;
		bool     mapping_dist_known = false;		
		//find dist
		for(uint32_t i = 0; i < fields.size(); i++)
			if(fields[i].compare(0,5,"NM:i:") == 0){
				mapping_dist = stoi (fields[i].substr(5));
				mapping_dist_known = true;
				
				NM_commulated += mapping_dist;
				NM_count++;
				break;
			}		

		bool correct_pos;
		if(query_fields.size() > 3){
			uint64_t sample_pos = stol( query_fields [3]);
			correct_pos = (abs(sample_pos,mapping_pos) <= CORRECT_POSITION_TRESHOLD) && (fields[2] == query_fields[1]); //Compare positions && Chromosome names
			if(!correct_pos && mapping_pos_mate) correct_pos = (abs(mapping_pos_mate,sample_pos) <= CORRECT_POSITION_TRESHOLD) && (fields[2] == query_fields[1]); //Compare positions && Chromosome names
			if(correct_pos) count_correct_pos++;
		}
		
		if(query_fields.size() > 5){
			uint32_t sample_dist =  stoi( query_fields [5]);
			
			if(correct_pos){
				if(mapping_dist_known && mapping_dist <= sample_dist)
					count_correct_pos_better_eq++;
				else
					count_correct_pos_worse++;
			}else{
				if(mapping_dist_known && mapping_dist < sample_dist)
					count_other_pos_better++;
				else if(mapping_dist_known && mapping_dist == sample_dist)
					count_other_pos_equal++;
				else
					count_other_pos_worse++;			
			}
		}
	}
	
	double avg_NM = 99;
	if(NM_count) avg_NM = (double) NM_commulated / NM_count;
	
	
	string tab = "\t";
	cout << "file"
	<< tab << "count_all"
	<< tab << "count_mapped" 
	<< tab << "count_unmap"
	<< tab << "count_correct_pos" 
	<< tab << "count_correct_pos_better_eq" 
	<< tab << "count_correct_pos_worse" 
	<< tab << "count_other_pos" 
	<< tab << "count_other_pos_better" 
	<< tab << "count_other_pos_equal"
	<< tab << "count_other_pos_worse" 
	<< tab << "avg_NM" 
	<< endl;
	
	
	cout << string(argv[1]) 
	<< tab << count_all 
	<< tab << count_mapped 
	<< tab << count_unmap 
	<< tab << count_correct_pos
	<< tab << count_correct_pos_better_eq 
	<< tab << count_correct_pos_worse 
	<< tab << count_other_pos_better +  count_other_pos_equal + count_other_pos_worse
	<< tab << count_other_pos_better 
	<< tab << count_other_pos_equal 
	<< tab << count_other_pos_worse 
	<< tab << avg_NM 
	<< endl;
}
