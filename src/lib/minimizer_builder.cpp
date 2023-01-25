#pragma once
namespace gedmap_mini{
using namespace std;
using namespace sdsl;

// Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences: code http://bit.ly/invihgi
uint64_t hash(uint64_t key) {
// 	return key;
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}

uint64_t hash_inverse(uint64_t key) {
// 	return key;
	uint64_t tmp;

	// Invert key = key + (key << 31)
	tmp = key-(key<<31);
	key = key-(tmp<<31);

	// Invert key = key ^ (key >> 28)
	tmp = key^key>>28;
	key = key^tmp>>28;

	// Invert key *= 21
	key *= 14933078535860113213u;

	// Invert key = key ^ (key >> 14)
	tmp = key^key>>14;
	tmp = key^tmp>>14;
	tmp = key^tmp>>14;
	key = key^tmp>>14;

	// Invert key *= 265
	key *= 15244667743933553977u;

	// Invert key = key ^ (key >> 24)
	tmp = key^key>>24;
	key = key^tmp>>24;

	// Invert key = (~key) + (key << 21)
	tmp = ~key;
	tmp = ~(key-(tmp<<21));
	tmp = ~(key-(tmp<<21));
	key = ~(key-(tmp<<21));

	return key;
}
	
/**
 * @brief representation of the kmer minimizer index 1d-vectors
 * 
 * minimizers are of lenght k
 * the window has size w
 * 
 * if inidicator[kmer], then positions for kmer are stored
 * 
 * let kmer be the x'th kmer for wich a position is stored, then: 
 *  1) ind_rs(kmer) = x
 *  2) positions[table[x]] contains the first postions for this kmer
 *  3) positions[table[x+1]-1] contains the last postions for this kmer
 */
struct minimizer_index{
	uint32_t k;
	uint32_t w;
	sdsl::int_vector<0> positions;
	sdsl::int_vector<0> table;
	sdsl::bit_vector    indicator;
	sdsl::rank_support_v<> ind_rs;
	
	
	template<typename kmer_int_type>
	//<begin,size>
	std::pair<size_t, size_t> boundaries(kmer_int_type kmer){
		if(! indicator[kmer] ) return make_pair(0,0);
		
		uint64_t x = ind_rs(kmer);
		return make_pair(table[x] , table[x+1] - table[x]);
	}
	
	template<typename kmer_int_type>
	sdsl::int_vector<0> operator()(kmer_int_type kmer) const{
		if(! indicator[kmer] ) 
			return sdsl::int_vector<0>(0,0,positions.width());
		
		uint64_t x = ind_rs(kmer);
		uint64_t first = table[x];
		uint64_t last = table[x+1];
		sdsl::int_vector<0> out(last-first,0,positions.width());
		uint64_t i = 0;
		while(first < last) out[i++] = positions[first++];
		return out;
	}
	
	
	typedef sdsl::int_vector<>::size_type size_type;
	//serialization  
	/** \brief serealize to file*/  
	size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const{  
		sdsl::structure_tree_node *child 			= sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
		size_type written_bytes 	= 0;
		
		sdsl::int_vector<32> constants = sdsl::int_vector<32>(2,k,32);
		constants[1] = w;
		
		written_bytes += constants	.serialize(out, child, "k_w");
		written_bytes += positions	.serialize(out, child, "positions");
		written_bytes += table		.serialize(out, child, "table");
		written_bytes += indicator	.serialize(out, child, "ind");
		written_bytes += ind_rs		.serialize(out, child, "ind_rs");

		sdsl::structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}
	
	/** \brief load from file*/
	void 	load(std::istream& in){
		sdsl::int_vector<32> constants;
		constants	.load(in);
		k = constants[0];
		w = constants[1];
		positions	.load(in);
		table		.load(in);
		indicator	.load(in);
		ind_rs	.load(in,&indicator);
	} 
	
	vector<uint64_t> trim(uint32_t max_positions){	
		vector<uint64_t> stats(6,0);
		
		stats[0] = positions.size();
		stats[2] = table.size()-1;
		
		
		uint64_t position_pointer = 0;
		uint64_t table_pointer_new = 0;		
		uint64_t table_pointer_old = -1;
		for(uint64_t kmer = 0; kmer < indicator.size(); kmer++){
			if(!indicator[kmer]) continue;
			
			table_pointer_old++;
			uint64_t kmer_begin = table[table_pointer_old];
			uint64_t kmer_size  = table[table_pointer_old+1] - kmer_begin;
			
			//STATS 
			if(kmer_size > stats[5]){
				stats[4] = kmer;
				stats[5] = kmer_size;
			}
			
			if( kmer_size > max_positions ){
				//ignore kmer
				indicator[kmer] = 0;
			}else{
				//copy positions[kmer_begin, kmer_begin+kmer_size[ to positions[position_pointer,position_pointer+kmer_size[
				table[table_pointer_new++] = position_pointer;
				for(uint64_t i = 0; i < kmer_size; i++)
					positions[position_pointer++] = positions[kmer_begin++];
			}			
		}
		table[table_pointer_new++] = position_pointer; // last entry contains size of positions
		positions.resize(position_pointer);
		table.resize(table_pointer_new);
		util::init_support(ind_rs,&indicator); 
		
		stats[1] = positions.size();
		stats[3] = table.size()-1;
		
		return stats;
	}	
};

struct minimizer_builder{
	
	uint64_t offset;	
	//min[i] stores the minimum kmers starting at position i + string_offset
	sdsl::int_vector<64> min;
	sdsl::bit_vector min_ind;
	
	minimizer_builder(){};
	
	minimizer_builder(size_t size, uint64_t offset):offset(offset){
		min = sdsl::int_vector<64>(size,-1,64);
		min_ind = sdsl::int_vector<1>(size,0,1);
	}
	
	//window size w
	void minimize(uint32_t w){
		min_ind = sdsl::int_vector<1>(min.size(),0,1);		
		for(size_t i = 0; i + w -1 < min.size(); i++){
			//check window starting at position i
			uint64_t m = min[i];
			uint64_t pos = i;
			for(size_t j = 1; j < w; j++){
				if(min[i+j] < m){
					m = min[i+j];
					pos = i+j;
				}
			}
			if(min[pos] != (uint64_t) -1)
				min_ind[pos] = 1;
		}
	}
	/* Instead of sorting -> count and bucket -> random access slows down
	minimizer_index index2(uint32_t k, uint32_t w, uint32_t t){
		if(min.size() != min_ind.size()) throw std::runtime_error("minimizer_index index: vectors do not have same size");
		if(t  > 65534 || t == 0) throw std::runtime_error("minimizer_index index: trim-value has to be in range [1,65534]");
		
		minimizer_index index = minimizer_index();
		index.k = k;
		index.w = w;
		index.indicator = sdsl::int_vector<1>(KMER<uint64_t>::number_of_kmers(k),0,1);
		
		sdsl::int_vector<32> count_per_kmer = sdsl::int_vector<32>(index.indicator.size(),0,32);
		
		cout << "count kmers/ calc indicator" << endl;
		
		
		//count kmers, trim and fill in indicator
		uint64_t count_of_kmers = 0;
		uint64_t count_of_positions = 0;
		uint64_t count_of_kmers_befor_trim = 0;
		uint64_t count_of_positions_befor_trim = 0;
		for(uint64_t i = 0; i < min.size(); i++){
			if( min_ind[i]){ 
				uint64_t kmer = hash_inverse(min[i]); //kmer is a minimizer
				count_of_positions_befor_trim++;
				if(count_per_kmer[kmer] == 0){ //first seen -> mark in index
					index.indicator[kmer] = 1;
					count_of_kmers++;
					count_per_kmer[kmer]++;
					count_of_positions++;
					count_of_kmers_befor_trim++;
				} 
				else if( count_per_kmer[kmer] == t){ //trim
					count_per_kmer[kmer] = t+1; //not use count_per_kmer[kmer] anymore
					index.indicator[kmer] = 0; //trimed, therefor not present in index
					count_of_kmers--;
					count_of_positions -= t;
				}
				else if( count_per_kmer[kmer] < t ){ //add kmer
					count_per_kmer[kmer]++;
					count_of_positions++;
				}
				
			}
		}		
		util::init_support(index.ind_rs,&index.indicator); 
		//transform count_per_kmer to be comulative / tmp table
		cout << "tmp table" << endl;
		uint32_t count_smaller_i = 0;
		for(uint32_t i = 0; i < count_per_kmer.size(); i++ ){
			if(count_per_kmer[i] > t) count_per_kmer[i] = 0;
			uint32_t count_i = count_per_kmer[i];
			count_per_kmer[i] =  count_smaller_i;
			count_smaller_i += count_i;
		}
		cout << "count_smaller_i " << count_smaller_i << endl;
		cout << "count_per_kmer[count_per_kmer.size()-1] " << count_per_kmer[count_per_kmer.size()-1] << endl;
		
		
// 		cout << "calc table" << endl;
// 		
// 		sdsl::int_vector<0>(count_of_kmers+1,0,occpointer_int_width);
// 		
// 		
		
// 		uint32_t occpointer_int_width = bitsneeded(count_of_positions);
// 		if(occpointer_int_width % 8 != 0) occpointer_int_width = occpointer_int_width - (occpointer_int_width % 8) + 8;
// 		index.table = sdsl::int_vector<0>(count_of_kmers+1,0,occpointer_int_width);
// 		
// 		//build table
// 		uint64_t next_kmer_in_table = 1;
// 		for(uint64_t kmer = 0; kmer < count_per_kmer.size(); kmer++){
// 			if(index.indicator[kmer]){ //not trimmed and not empty
// 				index.table[ next_kmer_in_table ] = index.table[ next_kmer_in_table - 1 ] + count_per_kmer[kmer];
// 				next_kmer_in_table++;
// 			}
// 		}
		
		
		cout 	<< "count_of_positions " << count_of_positions << endl
			<< "count_of_kmers " << count_of_kmers  << endl
			<< "count_of_positions_befor_trim " << count_of_positions_befor_trim  << endl
			<< "count_of_kmers_befor_trim"<< count_of_kmers_befor_trim << endl;
		
		cout << "fill in positions why?" << endl;
		
		
		//build positions
// 		count_per_kmer.resize(0);		
		uint32_t occ_width = bitsneeded(min.size());
		if(occ_width % 8) occ_width = occ_width + 8 - (occ_width % 8);
		index.positions = sdsl::int_vector<0>(count_of_positions,0,occ_width);
		
		
		for(uint64_t i = 0; i < min.size(); i++){
			if( min_ind[i]){ 
				uint64_t kmer = hash_inverse(min[i]); //kmer is a minimizer
				if(kmer >= index.indicator.size() || kmer >= count_per_kmer.size()){
					cout << "kmer " << kmer << " index.indicator.size()  " << index.indicator.size() << " count_per_kmer.size() " << count_per_kmer.size()<< endl;
				}
				else
				if(index.indicator[kmer]){ //not trimmed
					if(count_per_kmer[kmer] >= index.positions.size() ){
						cout << "kmer " << kmer << " count_per_kmer[kmer] " << count_per_kmer[kmer]  << " index.positions.size()  " << index.positions.size() <<endl;
					}else{
					index.positions[ count_per_kmer[kmer] ] = i;
					count_per_kmer[kmer] = count_per_kmer[kmer]+1;}
				}
			}
		}
		
		exit(0);
		
		cout << "adjust table" << endl;
		
		
		for(uint32_t i = index.table.size(); i > 0; i--)
			index.table[i] = index.table[i-1];
		index.table[0] = 0;
		
		return index;
		
	}*/
	
	minimizer_index index(uint32_t k, uint32_t w, uint32_t t){
		string fn_key_val_buf = TMP_DIR+ "/key_val_buf";
		string fn_tmp_positions = TMP_DIR+ "/tmp_positions";
		string fn_tmp_table = TMP_DIR+ "/tmp_table";
		
		
		if(min.size() != min_ind.size()) throw std::runtime_error("minimizer_index index: vectors do not have same size");
		
		uint32_t occ_width = bitsneeded(min.size());
		cout << "Generate key-value pairs" << endl;
// 		std::vector<std::pair<uint64_t,uint64_t>> kmer_pos;
		
		sdsl::int_vector_buffer<64> key_val_buf_out(fn_key_val_buf, ios::out|ios_base::trunc, 100*MB, 64 );
		uint32_t count_pairs = 0;
		for(uint64_t i = 0; i < min.size(); i++){
			if( min_ind[i]){
				key_val_buf_out.push_back(hash_inverse(min[i]));
				key_val_buf_out.push_back(i);
				count_pairs++;
			}
		}
		key_val_buf_out.close();
		
		
// 		std::vector<uint64_t>().swap(min); //DELETE MIN / MIN_IND
		min= sdsl::int_vector<64>(0,0,64);
		min_ind = sdsl::int_vector<1>(0,0,1);
		
		sdsl::int_vector_buffer<64> key_val_buf_in(fn_key_val_buf, ios::in, 100*MB, 64 );
		std::vector<std::pair<uint64_t,uint64_t>> kmer_pos(count_pairs);
		for(uint32_t i = 0; i < count_pairs; i++){
			kmer_pos[i] = std::make_pair( key_val_buf_in[2*i] ,key_val_buf_in[2*i+1] );
		}
		key_val_buf_in.close();
// 		std::sort(kmer_pos.begin(), kmer_pos.end());
		cout << "Sort key-value pairs"<< endl;
		mergeSortParallel(kmer_pos); //sort parallel
		
		
		cout << "Write to index"<< endl;
		
		minimizer_index index = minimizer_index();
		
		index.k = k;
		index.w = w;
		
		
		
		
		if(occ_width % 8) occ_width = occ_width + 8 - (occ_width % 8); //Round up to next multiple of 8
		
// 		index.positions = sdsl::int_vector<0>(kmer_pos.size(),0,occ_width);
		sdsl::int_vector_buffer<0> positions_buf(fn_tmp_positions, ios::out|ios_base::trunc, 100*MB, occ_width );
		
		uint32_t occpointer_int_width = bitsneeded(kmer_pos.size());
		if(occpointer_int_width % 8 != 0) occpointer_int_width = occpointer_int_width - (occpointer_int_width % 8) + 8;
		
// 		index.table = sdsl::int_vector<0>(KMER<uint64_t>::number_of_kmers(k)+1,0,occpointer_int_width);
		sdsl::int_vector_buffer<0> table_buf(fn_tmp_table, ios::out|ios_base::trunc, 100*MB, occpointer_int_width );
		
		index.indicator = sdsl::int_vector<1>(KMER<uint64_t>::number_of_kmers(k),0,1);
		
		
		uint64_t occ_pointer = 0;
// 		uint64_t tab_pointer = 0;
		
		
		uint64_t i = 0;
		
		uint32_t trimmed_kmers = 0;
		uint32_t trimmed_positions = 0;
			
		while(i < kmer_pos.size()){
			uint64_t kmer = kmer_pos[i].first;
			
			uint32_t count = 0;
			while( (i+count) < kmer_pos.size() && kmer == kmer_pos[(i+count)].first) count++;
			
			if(t && count > t){ // trim
				trimmed_kmers++;
				trimmed_positions += count;
				i += count;
			}else{
				index.indicator[kmer] = 1;
// 				index.table[tab_pointer++] = occ_pointer;
				table_buf.push_back(occ_pointer);
				while(count-- > 0) positions_buf[occ_pointer++] =  kmer_pos[i++].second;
			}
				

				
// 			while(i < kmer_pos.size() && kmer == kmer_pos[i].first){
// 				index.positions[occ_pointer++] =  kmer_pos[i++].second;
// 				positions_buf[occ_pointer++] =  kmer_pos[i++].second;
// 			}
		}
		
		cout << "(trimmed kmers " << trimmed_kmers << ')'<< endl;
		cout << "(trimmed positions " << trimmed_positions << ')' << endl;
		
// 		index.table[tab_pointer++] = occ_pointer;
		table_buf.push_back(occ_pointer);
// 		index.table.resize(tab_pointer);
		table_buf.close();
		positions_buf.close();
		
		std::vector<std::pair<uint64_t,uint64_t>>().swap(kmer_pos);
		
		load_from_file(index.positions,fn_tmp_positions);
		load_from_file(index.table,fn_tmp_table);
		
		remove(fn_tmp_positions.c_str());
		remove(fn_tmp_table.c_str());
		remove(fn_key_val_buf.c_str());
		
		util::init_support(index.ind_rs,&index.indicator); 
		return index;
	}
};


/** \brief for debugging */
std::ostream& operator<< (std::ostream& os, const minimizer_index& m){
	os << "--minimizer index:" << endl;
	os << "k = " << m.k << endl;
	os << "w = " << m.w << endl;
	
	uint64_t x = -1;
	for(size_t kmer = 0; kmer < m.indicator.size(); kmer++){
		KMER<uint64_t> km(kmer,m.k);
		if(m.indicator[kmer]) {
			x++;
			cout << km.get_kmer() << ":" << m(kmer) << endl;
			
		}
	}
	os << "minimizer index--" << endl;
	return os;	
}



}
