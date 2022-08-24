
#include <string_view>

namespace gedmap_index_min{

using namespace std;
using namespace sdsl;
using namespace gedmap_mini;

vector<sys_timer> tv(7);

const uint32_t tv_ALL  	= 0;
const uint32_t tv_LOAD 	= 1;
	//INDEX NODES CALC MINIMIZER BUILD INDEX
const uint32_t tv_INDEX_NODES  	= 2;
const uint32_t tv_COMBINE_NODES	= 3;
const uint32_t tv_BUILD_INDEX	= 4;
const uint32_t tv_STORING  	= 5;
const uint32_t tv_TRIM  	= 6;

minimizer_index calculate_index( string & eds, const adjacency & adj, const pos_EDS_to_FA_type & p2FA);


int main(int argc,  char** argv){
	using namespace gedmap_io;
	print_prog_headline("GEDMAP INDEX");

	tv[tv_ALL].start();

	//LOAD
	string eds;
	string min_fname;
	adjacency adj;
	pos_EDS_to_FA_type p2FA;

	tv[tv_LOAD].start();
	handle_input(argc, argv,  eds, min_fname, adj, p2FA);
	tv[tv_LOAD].stop();

	print_row("GEDS:", argv[2]);
	if(adj.initialised) print_row("including the graph file");
	print_row("K:", KMER_SIZE);
	print_row("W:", WINDOW_SIZE);
	print_row("INDEX:", min_fname);
	dotline();

	// CALCULATE INDEX
	minimizer_index index = calculate_index(eds, adj, p2FA);
	
	tv[tv_STORING].start();
	print_row("Store index");
	file_remove_silent(min_fname);
	if(!store_to_file<minimizer_index>(index, min_fname))	print_error("Storing of index failed!");
	else	print_row("Index written. Number of positions in index ", index.positions.size());
	tv[tv_STORING].stop();

	//FINISH
	dotline();
	print_row("Time for loading:",		tv[tv_LOAD]			.get(), " s");
	print_row("Time for indexing nodes:",	tv[tv_INDEX_NODES]	.get(), " s");
	print_row("Time for combining nodes",	tv[tv_COMBINE_NODES]	.get(), " s");
	print_row("Time for building index",	tv[tv_BUILD_INDEX]	.get(), " s");
	if(TRIM)print_row("Time for trimming index",	tv[tv_TRIM]		.get(), " s");
	print_row("Time for storing:",		tv[tv_STORING]		.get(), " s");
	print_row("Time in total:", 	 		tv[tv_ALL].stop_and_get(), " s");
	dotline();
	return 0;
}


void index_out(const string_view & eds, minimizer_builder & mb, uint64_t offset, std::map<uint64_t,vector<string>> & node_K);
void include_kmers_over_edges_for_node( const string_view & node, minimizer_builder & mb, vector<string> K);
vector<string> edge_set(uint64_t i, const adjacency & adj, std::map<uint64_t,vector<string>> & node_K);

template<class E>
void remove_duplettes(vector<E> & A);

minimizer_index calculate_index( string & eds, const adjacency & adj, const pos_EDS_to_FA_type & p2FA){
	
// 	if(k >= 16){ gedmap_io::print_error("k is to big?"); exit(0);}
	
	vector<uint64_t> boundaries; //boundaries are chromosome bounderies or node boundaries
	{
		set<uint64_t> boundaries_set; 
		std::set<uint64_t>::iterator hint = boundaries_set.begin();
		//include chromosome boundaries
		for(uint32_t i = 1; i < p2FA.chrom_starts.size(); i++ )
			hint = boundaries_set.insert( hint, p2FA.chrom_starts[i]-1 );
		
		hint = boundaries_set.begin();
		//include node boundaries
		if(adj.initialised)
			for(auto f_e : adj.forward_edges)
				hint = boundaries_set.insert( hint, f_e.first );
		
// 		hint = boundaries_set.begin()
// 		for(std::map<uint64_t,uint64_t>::iterator b_e = adj.backward_edges.begin(); b_e != adj.backward_edges.end(); b_e++){
// 			hint = boundaries_set.insert( hint, b_e->first );
// 		}
		boundaries = vector<uint64_t>(boundaries_set.size());
		uint32_t i = 0;
		for(set<uint64_t> ::iterator it = boundaries_set.begin(); it != boundaries_set.end(); it++) boundaries[i++] = *it;
	}
	
	uint32_t node_count = boundaries.size()+1;
	
	vector<minimizer_builder > node_indexes(node_count);
	vector<pair<uint64_t,uint64_t>> node_bounds(node_count);
	std::map<uint64_t,vector<string>> node_K_map;
	
	tv[tv_INDEX_NODES].start();
	gedmap_io::print_row("Calculate index per node");
	#pragma omp parallel for 
	for(uint64_t i = 0; i < node_count ; i++){
		uint64_t start = i==0?0:boundaries[i-1];
		uint64_t end = i==(boundaries.size())?eds.size():boundaries[i];
		uint64_t size = end-start;

		string_view node(eds.data() + start,size);
		
		node_indexes[i] 	= minimizer_builder(size,start);
		node_bounds[i]	= make_pair(start,end);
		
		index_out(node, node_indexes[i], start, node_K_map);
	}
	
	if(adj.initialised) {
		gedmap_io::print_row("Include edges");
		#pragma omp parallel for
		for(uint64_t i = 0; i < node_count; i++){
			
			uint64_t start = node_bounds[i].first;
			uint64_t end = node_bounds[i].second;
			uint64_t size = end-start;
			string_view node(eds.data() + start,size);
			
			vector<string> K = edge_set( end , adj, node_K_map);
			
			if(!K.empty())
				include_kmers_over_edges_for_node(node, node_indexes[i], K);
		}
	}
	
	vector<pair<uint64_t,uint64_t>>().swap(node_bounds);
	
	gedmap_io::print_row("Calculate minimizers");	
	#pragma omp parallel for
	for(uint64_t i = 0; i < node_count; i++)
		node_indexes[i].minimize(WINDOW_SIZE);
	
	
	tv[tv_INDEX_NODES].stop();
	tv[tv_COMBINE_NODES].start();
	
	
	gedmap_io::print_row("Combine node minimizers");
	
	
	minimizer_builder mb = minimizer_builder(); 
	{ /* //Try to save space, did not work well
	sdsl::int_vector_buffer<64> tmp_min_buf("tmp_min_buf", ios::out|ios_base::trunc, 100*MB, 64 );
	sdsl::int_vector_buffer<1> tmp_min_ind_buf("tmp_min_ind_buf", ios::out|ios_base::trunc, 100*MB, 1 );
	for(uint64_t i = 0; i < node_count ; i++){
		for(uint64_t j = 0; j < node_indexes[i].min.size(); j++ ){
			tmp_min_buf[node_indexes[i].offset + j] = node_indexes[i].min[j];
			tmp_min_ind_buf[node_indexes[i].offset + j] = node_indexes[i].min_ind[j];
		}
	}
	tmp_min_buf.close();
	tmp_min_ind_buf.close();
	vector<minimizer_builder >().swap(node_indexes);
	load_from_file(mb.min,"tmp_min_buf");
	load_from_file(mb.min_ind,"tmp_min_ind_buf");
	*/
	}
	{
	mb = minimizer_builder (eds.size(),0); //TODO THIS COPYING IS MEMORY PEAK. 2x 30GB INDEX 
	#pragma omp parallel for
	for(uint64_t i = 0; i < node_count ; i++){
		for(uint64_t j = 0; j < node_indexes[i].min.size(); j++ ){
			mb.min[node_indexes[i].offset + j] = node_indexes[i].min[j];
			mb.min_ind[node_indexes[i].offset + j] = node_indexes[i].min_ind[j];
		}
	}	
	vector<minimizer_builder >().swap(node_indexes);
	}
	
	tv[tv_COMBINE_NODES].stop();
	tv[tv_BUILD_INDEX].start();
	
	gedmap_io::print_row("Build minimizer index");
	// DIFFERENTIATE BETWEEN INTEGER SIZES
	
	
	string().swap(eds);
// 	minimizer_index mini = mb.index2(KMER_SIZE,WINDOW_SIZE,TRIM);
	minimizer_index mini = mb.index(KMER_SIZE,WINDOW_SIZE,TRIM);
	
	
	tv[tv_BUILD_INDEX].stop();
	
	return mini;
}


//###############DECLARATION OF HELP METHODS##################
void add_letter_at_front(uint32_t k, vector<string> & K, char c);

void fill_in_kmers_for_position(uint64_t pos, uint32_t k, vector<string> & K, minimizer_builder & mb);

vector<string> merge( vector<string> & K1, vector<string> & K2);

void add_to_node_map(uint64_t pos, vector<string> & K, std::map<uint64_t,vector<string>> & node_map);
//#############################################################

/**
 * @brief calculates kmer index
 */
void index_out(const string_view & eds, minimizer_builder & mb, uint64_t offset, std::map<uint64_t,vector<string>> & node_K){
	uint64_t i = eds.size();
	char c;
	//strings in K are allways sorted!
	vector<string> K;		//kmers
	vector<string> K_i;	//kmers before entering alt scope
	vector<string> K_o; 	//kmers when hitting boarder to next alt
	K.push_back("");

	while( i > 0 ){
		c = eds[--i];
		switch(c){
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'N':{
				add_letter_at_front(KMER_SIZE,K,c);
				fill_in_kmers_for_position(i,KMER_SIZE,K, mb);
				break;
			}
			case ')':
				K_i = vector<string>(K);
				break;
			case '|':
				K_o = merge(K_o, K);
				K = vector<string>(K_i);
				break;
			case '(':
				K = merge(K_o, K);
				K_o.clear();
				K_i.clear();
				break;
			case EDS_NODE_BOUNDARY: //BARRIER TO NEXT CHROM
				#pragma omp critical
				{
					add_to_node_map(i+offset,K,node_K);
				}
				K.clear();
				K_o.clear();//should already be empty, but to be sure
				K_i.clear();//should already be empty, but to be sure
				K.push_back("");
			default:
				;//IGNORE LETTER
		}
	}
}

/**
 * @brief lookup set K for position i
 * 
 * first, lookup all edges starting from i in adj
 * then , iterate over edges and union all sets K stored in node_K
 */
vector<string> edge_set(uint64_t i, const adjacency & adj, std::map<uint64_t,vector<string>> & node_K){
	vector<string> K; 	//kmers on position i / output
	
	// lookup all edges starting from i in adj
	std::map<uint64_t,uint64_t>::const_iterator  node = adj.forward_edges.find(i);
	uint64_t  first_edge, last_edge;
	first_edge = node->second;
	node++;
	if(node != adj.forward_edges.cend()) last_edge = node->second-1;
	else last_edge = adj.forward_targets.size()-1;

	//iterate over edges and union all sets stored in node_K
	for(uint64_t e = first_edge; e <= last_edge; e++){ 
		uint64_t t = adj.forward_targets[e];
		std::map<uint64_t,vector<string>>::iterator target_set = node_K.find(t); //find target nodes
		if(target_set == node_K.end()) throw runtime_error("index_eoc, include_kmers_over_edges: target_set not found - maybe eds and adj file do not match");
		K = merge(K,target_set->second);
	}
	
	return K;
}


/**
 * @brief adds kmers that are going over edges
 * iterates over nodes in adj
 * for each node all conected nodes will be loaded
 * possible kmers vor each connected node will be added to possible kmers of this node
 * makes k steps and adds kmers to eoc
 */
void include_kmers_over_edges_for_node(const string_view & node, minimizer_builder & mb, vector<string> K){
	uint64_t i = node.size();
	vector<string> K_i;	//kmers before entering alt scope
	vector<string> K_o; 	//kmers when hitting boarder to next alt

	//start at pos i with k steps to go.
	uint32_t k = KMER_SIZE-1;
	uint32_t k_i = 0;
	uint32_t k_o = 0;
	while( i > 0 && (k > 0 || k_i > 0 || k_o > 0) ){
		char c = node[--i];
		switch(c){
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'N':{
				add_letter_at_front(KMER_SIZE,K,c);
				if(k > 0){
					fill_in_kmers_for_position(i,KMER_SIZE,K, mb);
					k--;
				}
				break;
			}
			case ')':
				K_i = vector<string>(K);
				k_i = k;
				k = 0;
				break;
			case '|':
				K_o = merge(K_o, K);
				if (k_o < k) k_o = k; // k_o = max(k_o,k)
				K = vector<string>(K_i);
				break;
			case '(':
				K = merge(K_o, K);
				if (k < k_o) k = k_o; // k = max(k_o,k)
				k_o = 0;
				k_i = 0;
				K_o.clear();
				K_i.clear();
				break;
			case EDS_NODE_BOUNDARY: //BARRIER TO NEXT CHROM
				break; //do not go over two edges
			default:
				;//IGNORE LETTER
		}
	}

}

//################################# HELP METHODS ##################

/** @brief returns the number of N is the given string T*/
uint32_t count_N(string & T){
	uint32_t c = 0;
	for(uint32_t i = 0; i < T.size(); i++)
		if(T[i] == 'N') c++;
	return c;
}


/**
 * @brief returns a vector of all possible kmers represented by the given string T
 * T can include N and therefor represent several kmers
 *
 * kmer_type:  uint32_t for k <= 16; uint64_t for k <= 32
 */
template<typename kmer_type>
vector<kmer_type> possible_kmers(string & T){
	uint32_t cN = count_N(T);
	if(cN > MAX_N_IN_SEED) return vector<kmer_type>(0); //DO NOT USE THIS T, TO MANY N

	// K stores all kmers
	vector<kmer_type> K(1);
	K[0]=0;

	for(uint32_t k = 0; k < T.size(); k++){
		if(T[k] != 'N'){ //add T[k] to the back of all kmers
			for(uint32_t i = 0; i < K.size(); i++)
				K[i] = KMER<kmer_type>::add_back(K[i],T[k]);
		}
		else{ //add each nucleotide to the back of all kmers
			vector<kmer_type> K_tmp = vector<kmer_type>(K.size()*4);
			uint32_t j = 0;
			for(uint32_t i = 0; i < K.size(); i++)
				for(char c : "ACGT")
					if(c != 0) K_tmp[j++] = KMER<kmer_type>::add_back(K[i],c);
			K.swap(K_tmp);
		}
	}
	return K;
}

/** @brief removes duplettes in a sorted vector A*/
template<class E>
void remove_duplettes(vector<E> & A){
	if(A.size() == 0) return;
	uint32_t idx = 1;
	bool duplette_occured = false;
	for(uint32_t i = 1; i < A.size(); i++){
		if(A[i] == A[i-1])
			duplette_occured = true;
		else if (duplette_occured)
			A[idx++] = A[i];
		else
			idx++;
	}
	A.resize(idx);
}

/** @brief transforms  vector<string>  to vector<kmer_type>
 * ignores strings with to many N
 * replaces N by A,C,G and T
 * gets the kmer-Identifier
 * @return a vector with kmers (sorted and unique)
 */
template<typename kmer_type>
vector<kmer_type> transform_string_vec_to_kmer_type_vec(vector<string> & K, uint32_t k){
	//COUNT RESULTING KMERS ( sum(4^numer_of_N) )
	uint32_t count = 0;
	for(string kmer : K){
		uint32_t cN = count_N(kmer);
		if(cN <= MAX_N_IN_SEED && kmer.size() == k) {
			uint32_t c = 1;
			c  <<= (2*cN);
			count += c;
		}
	}
	//ENTER NEW KMERS TO new_kmers
	vector<kmer_type> new_kmers(count);
	uint32_t it = 0;
	for(string kmer : K){
		if(kmer.size() < k)
			continue;
		vector<kmer_type> tmp = possible_kmers<kmer_type>(kmer);
		for(kmer_type k : tmp)
			new_kmers[it++] = k;
	}
	//AVOID  DUPLETTES
	std::sort(new_kmers.begin(), new_kmers.end());
	remove_duplettes(new_kmers);
	return new_kmers;
}

/** @brief add kmers an beginning of a node to the map*/
void add_to_node_map(uint64_t pos, vector<string> & K, std::map<uint64_t,vector<string>> & node_map){
	if( node_map.find(pos) != node_map.end() )
		throw runtime_error("INDEX_EOC > add_to_node_map > node already has values");
	node_map[pos] = K;
}

/**
 * @brief adds minimum kmer starting at position pos to minimizer_builder
 */
void fill_in_kmers_for_position(uint64_t pos, uint32_t k, vector<string> & K, minimizer_builder & mb){
	vector<uint64_t> K_new = transform_string_vec_to_kmer_type_vec<uint64_t>(K, k);
	
	uint64_t min = -1;
	for( auto kmer : K_new){
		uint64_t h = gedmap_mini::hash(kmer);
		if(h < min) min = h;
	}
	
	if(min != (uint64_t) -1) if(min < mb.min[pos]) mb.min[pos] = min;
}

/**
 * @brief merges two sorted string vectors to one sorted string vector (without duplettes)
 */
vector<string> merge( vector<string> & K1, vector<string> & K2){
	vector<string> K(K1.size() + K2.size());
	vector<string>::iterator k1 = K1.begin();
	vector<string>::iterator k2 = K2.begin();
	uint32_t k = 0;

	while(k1 != K1.end() && k2 != K2.end()){
		if(*k1 < *k2){
			K[k++] = *k1;
			k1++;
		}else if (*k1 == *k2){
			K[k++] = *k1;
			k1++;
			k2++;
		} else{
			K[k++] = *k2;
			k2++;
		}
	}

	while(k1 != K1.end()){
		K[k++] = *k1;
		k1++;
	}

	while(k2 != K2.end()){
		K[k++] = *k2;
		k2++;
	}
	K.resize(k);
	return K;
}

/**
 * @brief adds c to the front of each kmer in K
 * ensures the the size of a kmer is at maximum k (by dropping the letter at the back)
 * ensures K is a sorted vector with unique elements (when K was alerady sorted)
 */
void add_letter_at_front(uint32_t k, vector<string> & K, char c){
	for(vector<string>::iterator kmer = K.begin(); kmer != K.end(); kmer++){
		if(kmer->size() < k) 	*kmer = c + *kmer;
		else 				*kmer = c + kmer->substr(0,k-1);
	}
	remove_duplettes(K);
}


}
