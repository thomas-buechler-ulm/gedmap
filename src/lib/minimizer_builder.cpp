#pragma once
#include <sdsl/sd_vector.hpp>
#include <sdsl/bit_vectors.hpp>
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
	
	sd_vector<> contained; //EF CODING
	bit_vector mult; 
	
	int_vector<0> start; //EF CODING??
	
	rank_support_sd<> contained_rs;
	rank_support_v<> mult_rs;
	
	int_vector<0> pos_sing;
	int_vector<0> pos_mult;
	
	
	template<typename kmer_int_type>
	bool is_in_index(kmer_int_type kmer) const{
		return (contained[kmer]);
	}
	
	template<typename kmer_int_type>
	//<begin,size>
	std::pair<size_t, size_t> get_position(kmer_int_type kmer) const{
		if(!contained[kmer]) return make_pair(0,0);
		
		uint64_t x = contained_rs(kmer);
		uint64_t y = mult_rs(x);
		
		if(!mult[x]) 
			return make_pair(pos_sing[x-y],1);
		else 
			return make_pair(start[y],start[y+1]-start[y]);
	}
	
	template<typename kmer_int_type>
	sdsl::int_vector<0> operator()(kmer_int_type kmer) const{
		size_t begin, count, i = 0;
		tie(begin,count) = get_position(kmer);
		
		if(count == 0)	return sdsl::int_vector<0>(0,0,1);
		if(count == 1)	return sdsl::int_vector<0>(1,begin,pos_sing.width());
	
		sdsl::int_vector<0> out(count,0,pos_mult.width());	
		while(count --> 0) out[i++] = pos_mult[begin++];
		return out;
	}
	
	
	typedef sdsl::int_vector<>::size_type size_type;
	//serialization  
	/** \brief serealize to file*/  
	size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const{  
		sdsl::structure_tree_node *child 			= sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
		size_type written_bytes 	= 0;
		
		sdsl::int_vector<32> constants = sdsl::int_vector<32>(2,k,32);
		constants[0] = k;
		constants[1] = w;
				
		written_bytes += constants	.serialize(out, child, "k_w");
		written_bytes += contained	.serialize(out, child, "contained");
		written_bytes += mult		.serialize(out, child, "mult");
		written_bytes += start		.serialize(out, child, "start");
		written_bytes += contained_rs	.serialize(out, child, "contained_rs");
		written_bytes += mult_rs	.serialize(out, child, "mult_rs");
		written_bytes += pos_sing	.serialize(out, child, "pos_sing");
		written_bytes += pos_mult	.serialize(out, child, "pos_mult");
		
		sdsl::structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}
	
	/** \brief load from file*/
	void 	load(std::istream& in){
		sdsl::int_vector<32> constants;
		constants	.load(in);
		k		= constants[0];
		w		= constants[1];
		contained	.load(in);
		mult		.load(in);
		start		.load(in);
		contained_rs	.load(in, &contained);
		mult_rs		.load(in, &mult);
		pos_sing	.load(in);
		pos_mult	.load(in);
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
	
	minimizer_index index(uint32_t k, uint32_t w, uint32_t t){
		string fn_key_val_buf = TMP_DIR+ "/key_val_buf";
		string fn_tmp_pos_sing = TMP_DIR+ "/tmp_positions_single";
		string fn_tmp_pos_mult = TMP_DIR+ "/tmp_positions_mult";
		string fn_tmp_start = TMP_DIR+ "/tmp_start";
		
		
		if(min.size() != min_ind.size()) 
			throw std::runtime_error("minimizer_index index: vectors do not have same size");
		
		uint32_t pos_width = bitsneeded(min.size());
		cout << "Generate key-value pairs" << endl;
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
		uint32_t start_width = bitsneeded(count_pairs);
				
		
		
// 		//DELETE MIN / MIN_IND
		min= sdsl::int_vector<64>(0,0,64);
		min_ind = sdsl::int_vector<1>(0,0,1);
		
		cout << "Read key-value pairs" << endl;
		sdsl::int_vector_buffer<64> key_val_buf_in(fn_key_val_buf, ios::in, 100*MB, 64 );
		std::vector<std::pair<uint64_t,uint64_t>> kmer_pos(count_pairs);
		for(uint32_t i = 0; i < count_pairs; i++)
			kmer_pos[i] = std::make_pair( key_val_buf_in[2*i] ,key_val_buf_in[2*i+1] );
		key_val_buf_in.close();
		
		cout << "Sort key-value pairs"<< endl;
		mergeSortParallel(kmer_pos); //sort parallel
		
		
		cout << "Analyse key-value pairs"<< endl;
		
		uint64_t n = KMER<uint64_t>::number_of_kmers(k);
		uint64_t kmers = 0, mult_kmer = 0, i = 0;
		
		while(i < kmer_pos.size()){
			uint64_t kmer = kmer_pos[i].first;
			
			uint32_t count = 0;
			while( (i+count) < kmer_pos.size() && kmer == kmer_pos[(i+count)].first) count++;
			
			if(t && count > t){ 
				// trim -> dont count			
			}else{
				kmers++;
				if(count > 1) mult_kmer++;
			}
			i += count;
		}		
		
		
		cout << kmers << "/" << n << " " << mult_kmer  << "/" << kmers << endl;
		
		cout << "Generate index"<< endl;
		
		minimizer_index index = minimizer_index();
		
		index.k = k;
		index.w = w;
		if(pos_width % 8) 	pos_width = pos_width + 8 - (pos_width % 8); //Round up to next multiple of 8
		if(start_width % 8)	start_width = start_width - (start_width % 8) + 8;
		
		sdsl::int_vector_buffer<0> pos_sing(fn_tmp_pos_sing, ios::out|ios_base::trunc, 100*MB, pos_width );
		sdsl::int_vector_buffer<0> pos_mult(fn_tmp_pos_mult, ios::out|ios_base::trunc, 100*MB, pos_width );
		sdsl::int_vector_buffer<0> start_buf(fn_tmp_start,   ios::out|ios_base::trunc, 100*MB, start_width );
		
		sd_vector_builder contained_b(n,kmers);
		index.mult = int_vector<1>(kmers,0,1);
		
		i = 0;
		uint64_t pos_mul_pointer	= 0;
		uint64_t contained_counter	= 0;
		uint64_t trimmed_kmers		= 0;
		uint64_t trimmed_positions	= 0;
		
		while(i < kmer_pos.size()){
			uint64_t kmer = kmer_pos[i].first;
			
			uint32_t count = 0;
			while( (i+count) < kmer_pos.size() && kmer == kmer_pos[(i+count)].first) count++;
			
			if(t && count > t){ // trim
				trimmed_kmers++;
				trimmed_positions += count;
				i += count;
			}else{
				if(count == 1){
					pos_sing.push_back(kmer_pos[i++].second);
				}else{
					index.mult[contained_counter] = 1;
					start_buf.push_back(pos_mul_pointer);
					while(count-- > 0) pos_mult[pos_mul_pointer++] =  kmer_pos[i++].second;
				}
				contained_b.set(kmer);
				contained_counter++;
			}
		}
		
		cout << "(trimmed kmers " << trimmed_kmers << ')'<< endl;
		cout << "(trimmed positions " << trimmed_positions << ')' << endl;
		
		start_buf.push_back(pos_mul_pointer);
		start_buf.close();
		pos_sing.close();
		pos_mult.close();
		
		std::vector<std::pair<uint64_t,uint64_t>>().swap(kmer_pos);
		
		load_from_file(index.start,fn_tmp_start);
		load_from_file(index.pos_sing,fn_tmp_pos_sing);
		load_from_file(index.pos_mult,fn_tmp_pos_mult);
		
		remove(fn_tmp_start.c_str());
		remove(fn_tmp_pos_sing.c_str());
		remove(fn_tmp_pos_mult.c_str());
		remove(fn_key_val_buf.c_str());
		
		index.contained	= sd_vector<>(contained_b);
		
		util::init_support(index.contained_rs,&index.contained);
		util::init_support(index.mult_rs,&index.mult);
		
		return index;
	}
};


/** \brief for debugging */
std::ostream& operator<< (std::ostream& os, const minimizer_index& m){
	os << "--minimizer index:" << endl;
	os << "k = " << m.k << endl;
	os << "w = " << m.w << endl;
	
	uint64_t x = -1;
	for(size_t kmer = 0; kmer < m.contained.size(); kmer++){
		KMER<uint64_t> km(kmer,m.k);
		if(m.contained[kmer]) {
			x++;
			cout << km.get_kmer() << ":" << m(kmer) << endl;
			
		}
	}
	os << "minimizer index--" << endl;
	return os;	
}



}
