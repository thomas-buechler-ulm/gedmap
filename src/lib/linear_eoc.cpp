#pragma once
/**
 * @brief representation of the kmer index (exact occurences) that uses two 1d-vectors
 * 
 * Let k_id be a number that represents a kmer, than
 * 
 * k is the lenght of k
 * and the exact occurences of the kmer are:
 * occurences[ table[k_id], table[k_id+1]-1]
 */

struct linear_eoc_type{
	uint32_t k;
	sdsl::int_vector<1> indicator;//true if positions are stored for kmer k (table[k_id] < table[k_id+1]-1)
	sdsl::int_vector<0> occurences;
	sdsl::int_vector<0> table;
	
	linear_eoc_type(){};
	
	linear_eoc_type(uint32_t k,sdsl::int_vector<1> indicator, sdsl::int_vector<0> occurences, sdsl::int_vector<0> table):k(k),indicator(indicator),occurences(occurences),table(table){}; //TODO MAYME NOT NEEDED
	
	linear_eoc_type(uint32_t k, sdsl::int_vector<0> occurences, sdsl::int_vector<0> table):k(k),occurences(occurences),table(table){
		set_indicator();
	};
	
	
	
	
	/**
	 * @brief generates the two vectors from an 2D vector
	 * 
	 */
	template<class int_type>
	linear_eoc_type( vector<vector<int_type>> & two_dim_eoc , uint32_t pos_int_width,  uint32_t k):k(k){
		
		//count total positions to store
		uint64_t occ_size = 0;		
		for(uint64_t kmer_id = 0; kmer_id < two_dim_eoc.size(); kmer_id++)
			occ_size += two_dim_eoc[kmer_id].size();
		
		//allocate space for the two vectors
		occurences = sdsl::int_vector<0>(occ_size,0,pos_int_width);
		uint32_t occpointer_int_width = bitsneeded(occ_size);
		if(occpointer_int_width % 8 != 0) occpointer_int_width = occpointer_int_width - (occpointer_int_width % 8) + 8;
		table = sdsl::int_vector<0>(two_dim_eoc.size()+1,0,occpointer_int_width);
		
		//fill vectors
		uint64_t occ_pointer = 0;
		table[0] = 0;
		for(uint64_t kmer_id = 0; kmer_id < two_dim_eoc.size(); kmer_id++){
			for(uint64_t pos = 0; pos < two_dim_eoc[kmer_id].size(); pos++) 
				occurences[occ_pointer++] = two_dim_eoc[kmer_id][pos];
			table[kmer_id+1] = occ_pointer;
		}
		
		set_indicator();		
	}
	
	typedef sdsl::int_vector<>::size_type size_type;
	
	//serialization  
	/** \brief serealize to file*/  
	size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const{  
		sdsl::structure_tree_node *child 			= sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
		size_type written_bytes 	= 0;
		
		sdsl::int_vector<32> k_vec = sdsl::int_vector<32>(1,k,32);
		
		written_bytes += k_vec		.serialize(out, child, "k");
		written_bytes += indicator	.serialize(out, child, "ind");
		written_bytes += occurences	.serialize(out, child, "occ");
		written_bytes += table		.serialize(out, child, "table");

		sdsl::structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}
	
	/** \brief load from file*/
	void 	load(std::istream& in){
		sdsl::int_vector<32> k_vec;
		k_vec		.load(in);
		indicator	.load(in);
		occurences	.load(in);
		table		.load(in);
		k = k_vec[0];
		set_indicator();
	} 
	
	void set_indicator(){
		indicator = sdsl::int_vector<1>(table.size()-1,0,1);
		for(uint64_t kmer_id = 0; kmer_id < indicator.size(); kmer_id++)
			indicator[kmer_id] = (table[kmer_id] < table[kmer_id+1]);
	}
	
	
	/** @brief removes all mall postiones for kmers with more than max_positions occurences*/
	vector<uint64_t> trim(uint32_t max_positions){
		
		uint64_t old_occsize = occurences.size();
		uint64_t old_emptycount = 0;
		uint64_t new_emptycount = 0;
		uint64_t greates_set_kmer = 0;
		uint64_t greates_set_val = 0;
		
		
		uint64_t occ_pointer = 0;
		
		for(uint64_t kmer_id = 0; kmer_id < table.size()-1; kmer_id++){
			
			uint32_t positions_count = table[kmer_id+1] - table[kmer_id];
			uint32_t old_positions_begin = table[kmer_id];
			
			{//STATS
				if(positions_count > greates_set_val){
					greates_set_val = positions_count;
					greates_set_kmer = kmer_id;
				}
			}
			
			
			table[kmer_id] = occ_pointer;
			
			if( positions_count > max_positions || positions_count == 0){
				//IGNORE THIS KMER
				if(positions_count == 0) old_emptycount++;
				new_emptycount++;
			}else{
				//LET THIS KMER STAY
				//COPY VALS
				for(uint32_t i = 0; i < positions_count; i++)
					occurences[ occ_pointer++] = occurences[old_positions_begin + i];
			}
		}
		occurences.resize(occ_pointer);
		table[ table.size()-1 ] = occ_pointer; //set right baundary correct
		set_indicator();
		
		vector<uint64_t> stats = {old_occsize, occ_pointer, old_emptycount, new_emptycount, greates_set_kmer, greates_set_val};
		return stats;
		
	}
};


/** \brief for debugging */
std::ostream& operator<< (std::ostream& os, const linear_eoc_type& eoc){
	os << "--linear_eoc_type:" << endl;
	os << "k = " << eoc.k << endl;
	
	for(size_t kmer = 0; kmer < eoc.indicator.size(); kmer++){
		KMER<uint64_t> km(kmer,eoc.k);
		cout << km.get_kmer() << ": " ;
		for(size_t i = eoc.table[kmer]; i < eoc.table[kmer+1]; i++)
			cout << eoc.occurences[i] << " ";
		cout << endl;
	}
	os << "linear_eoc_type--" << endl;
	return os;	
}
