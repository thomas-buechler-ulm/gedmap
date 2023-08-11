
#include <string_view>

namespace gedmap_index_kmi{

using namespace std;
using namespace sdsl;

vector<sys_timer> tv(6);

const uint32_t tv_ALL  	= 0;
const uint32_t tv_LOAD 	= 1;
const uint32_t tv_INDEX_OUT  	= 2;
const uint32_t tv_INDEX_IN  	= 3;
const uint32_t tv_STORING  	= 4;
const uint32_t tv_TRIM  	= 5;

template<class int_type, uint32_t int_width>
linear_eoc_type calculate_index(const uint32_t k, const string & eds, const adjacency & adj, const pos_EDS_to_FA_type & p2FA);


int main(int argc,  char** argv){
	using namespace gedmap_io;
	print_prog_headline("GEDMAP INDEX");

	tv[tv_ALL].start();

	//LOAD
	string eds;
	uint32_t k;
	string eoc_fname;
	adjacency adj;
	pos_EDS_to_FA_type p2FA;

	tv[tv_LOAD].start();
	handle_input(argc, argv,  eds, k, eoc_fname, adj, p2FA);
	tv[tv_LOAD].stop();

	print_row("EDS:", argv[2]);
	if(adj.initialised) print_row("including the graph file");
	print_row("K:", k);
	print_row("EOC:", eoc_fname);
	dotline();

	// CALCULATE INDEX
	// DIFFERENTIATE BETWEEN INTEGER SIZES
	const uint32_t int_width = bitsneeded(eds.size());
	linear_eoc_type index;
	if		(int_width <= 8  ) index = calculate_index<uint8_t ,  8>(k, eds, adj, p2FA);
	else if	(int_width <= 16 ) index = calculate_index<uint16_t, 16>(k, eds, adj, p2FA);
	else if	(int_width <= 24 ) index = calculate_index<uint32_t, 24>(k, eds, adj, p2FA);
	else if	(int_width <= 32 ) index = calculate_index<uint32_t, 32>(k, eds, adj, p2FA);
	else if	(int_width <= 48 ) index = calculate_index<uint64_t, 48>(k, eds, adj, p2FA);
	else if	(int_width <= 64 ) index = calculate_index<uint64_t, 64>(k, eds, adj, p2FA);
	else throw runtime_error (" in index main: EDS too long");

	if(TRIM){
		tv[tv_TRIM].start();
		vector<uint64_t> trim_stats =  index.trim(TRIM);
		gedmap_io::dotline();
		gedmap_io::print_row("TRIM REPORT");
		gedmap_io::print_row("position count before TRIM", trim_stats[0]);
		gedmap_io::print_row("position count after  TRIM", trim_stats[1]);
		gedmap_io::print_row("emptyset count before TRIM", trim_stats[2]);
		gedmap_io::print_row("emptyset count after  TRIM", trim_stats[3]);
		KMER<uint64_t> most = KMER(trim_stats[4], k);
		gedmap_io::print_row("most frequent kmer", most.get_kmer());
		gedmap_io::print_row("frequency", trim_stats[5]);
		tv[tv_TRIM].stop();
		gedmap_io::dotline();
	}
	
	tv[tv_STORING].start();
	print_row("Store index");
	file_remove_silent(eoc_fname);
	if(!store_to_file<linear_eoc_type>(index, eoc_fname))	print_error("Storing of index failed!");
	else	print_row("Index written. Number of positions in index ", index.occurences.size());
	tv[tv_STORING].stop();

	//FINISH
	dotline();
	print_row("Time for loading:",		tv[tv_LOAD]		.get(), " s");
	print_row("Time for indexing:",		tv[tv_INDEX_OUT]	.get(), " s");
	print_row("Time for rearrange index",	tv[tv_INDEX_IN]	.get(), " s");
	if(TRIM)print_row("Time for trimming index",	tv[tv_TRIM]		.get(), " s");
	print_row("Time for storing:",		tv[tv_STORING]	.get(), " s");
	print_row("Time in total:", 	 		tv[tv_ALL].stop_and_get(), " s");
	dotline();
	return 0;
}

template<typename kmer_type, class int_type, class off_stream_type>
void index_out(const uint32_t k, const string_view & eds,  off_stream_type & off_stream, vector<uint32_t>& count, uint64_t offset, std::map<uint64_t,vector<string>> & node_K);

template<typename kmer_type, class int_type, class off_stream_type>
void include_kmers_over_edges( const uint32_t k_c, const string & eds, off_stream_type & off_stream, vector<uint32_t>& count, const adjacency & adj, std::map<uint64_t,vector<string>> & node_K);

template<class E>
void remove_duplettes(vector<E> & A);

template<class int_type, uint32_t int_width>
linear_eoc_type calculate_index(const uint32_t k, const string & eds, const adjacency & adj, const pos_EDS_to_FA_type & p2FA){
	
	if(k >= 16){ gedmap_io::print_error("k is to big?"); exit(0);}
	
	gedmap_io::print_row("Calculate index");
	linear_eoc_type eoc;
	uint64_t number_of_kmers = KMER<uint64_t>::number_of_kmers(k);

	tv[tv_INDEX_OUT].start();

	uint32_t chrom_count = p2FA.chrom_starts.size();
	uint32_t calc_count  = chrom_count + (adj.initialised?1:0);

	//array that counts the kmers for each chromosome
	vector<vector<uint32_t>> counts(calc_count, vector<uint32_t>(number_of_kmers,0));
	//map that stores the value of K, for each node boundary of the EDS
	std::map<uint64_t,vector<string>> node_K_map;


	gedmap_io::print_row("Calculate exact occurences");
	vector<vector<uint64_t>> key_vals(calc_count,vector<uint64_t>(0));
	uint64_t chrom_ready = 0;
	omp_set_num_threads(chrom_count < MAX_THREAD_COUNT? chrom_count : MAX_THREAD_COUNT); //max of chrom_count and MAX_THREAD_COUNT

	#pragma omp parallel for
	for(uint64_t i = 0; i < chrom_count ; i++){
		uint64_t	chrom_begin = p2FA.chrom_starts[i];
		uint64_t	chrom_length = ( (i < chrom_count-1) ? p2FA.chrom_starts[i+1] : eds.size() ) - p2FA.chrom_starts[i];
 		string_view chrom_i(eds.data() + chrom_begin,chrom_length);

		index_out<uint32_t,int_type>(k, chrom_i,  key_vals[i], counts[i] , p2FA.chrom_starts[i], node_K_map);

		key_vals[i].shrink_to_fit();
		#pragma omp critical
		{
			chrom_ready++;
			gedmap_io::flush_row("Chromosomes ready", to_string(chrom_ready) + "/" + to_string(chrom_count) );
		}
	}
	
	gedmap_io::print_row("Calculate exact occurences, that span over edges");
	if(adj.initialised) include_kmers_over_edges<uint32_t,int_type>(k, eds, key_vals[chrom_count], counts[chrom_count], adj, node_K_map);
	key_vals[chrom_count].shrink_to_fit();
	
	tv[tv_INDEX_OUT].stop();
	tv[tv_INDEX_IN].start();

	gedmap_io::print_row("Count number of occurences per kmer");

	//array that contains all occurences
	sdsl::int_vector<0> occurences;
	uint64_t occurences_size = 0;
	//array that stores the boundaries of each kmer in the occurences-array //see linear_eoc_type
	sdsl::int_vector<0> table;

	omp_set_num_threads(MAX_THREAD_COUNT);
	#pragma omp parallel for
	for(uint64_t j = 0; j < number_of_kmers; j++){
		for(uint64_t i = 1; i < counts.size(); i++) counts[i][j] += counts[i-1][j];
		#pragma omp atomic
		occurences_size += counts[ counts.size()-1 ][j];
	}

	{	// initialise table
		uint32_t occpointer_int_width = bitsneeded(occurences_size);
		if(occpointer_int_width % 8 != 0) 
			occpointer_int_width = occpointer_int_width - (occpointer_int_width % 8) + 8;
		table		= sdsl::int_vector<0>(number_of_kmers+1,-1,occpointer_int_width);
		table[0]	= 0;
		for(uint64_t i = 1; i < table.size(); i++) 
			table[i] = table[i-1] + counts[ counts.size()-1 ] [i-1];
	}

	gedmap_io::print_row("Calculate array positions of occurence values");
	chrom_ready = 0;
	omp_set_num_threads(calc_count < MAX_THREAD_COUNT? calc_count : MAX_THREAD_COUNT); //max of calc_count and MAX_THREAD_COUNT
 	#pragma omp parallel for
	for(unsigned int i = 0; i < calc_count; i++){
		for(uint64_t j = 0; j < key_vals[i].size();j+=2){
			uint64_t kmer = key_vals[i][j];
			counts[i][kmer]--;
			key_vals[i][j] = table[kmer] + counts[i][kmer];
		}
		#pragma omp critical
		{
			chrom_ready++;
			gedmap_io::flush_row("Chromosomes ready", to_string(chrom_ready) + "/" + to_string(calc_count) );
		}
	}
	
	gedmap_io::print_row("Build index");
	occurences = int_vector<0>(occurences_size,-1,int_width);
	uint32_t tries = 0;
	while(true){
		tries++;
		uint64_t cc = 0;
		chrom_ready = 0;
		
		if(tries <= 3){
			#pragma omp parallel for
			for(unsigned int i = 0; i < calc_count; i++){
				uint64_t ci = 0;
				for(uint64_t j = 0; j < key_vals[i].size();){
					uint64_t target_pos = key_vals[i][j++];
					uint64_t pos  = key_vals[i][j++];
					occurences[target_pos] = pos; 
					ci++;
				}
				#pragma omp critical
				{
					cc += ci;
					chrom_ready++;
					gedmap_io::flush_row("Chromosomes ready", to_string(chrom_ready) + "/" + to_string(calc_count) );
				}
			}
		
		}else{
			gedmap_io::print_row("Try restoring using only one thread");
			for(unsigned int i = 0; i < calc_count; i++){
				uint64_t ci = 0;
				for(uint64_t j = 0; j < key_vals[i].size();){
					uint64_t target_pos = key_vals[i][j++];
					uint64_t pos  = key_vals[i][j++];
					occurences[target_pos] = pos; 
					ci++;
				}
				cc += ci;
				chrom_ready++;
				gedmap_io::flush_row("Chromosomes ready", to_string(chrom_ready) + "/" + to_string(calc_count) );
			}
		}
		
		gedmap_io::flush_row("Check index");
		
		if(cc != occurences.size()) throw runtime_error("Build Index: Less postions written, than expected");
		
		int_vector<0> tmp = int_vector<0>(1,-1,int_width);
		uint64_t empty = tmp[0];
		uint64_t empty_cells = 0;
		for(uint64_t i = 1; i < occurences_size; i++)
			if((uint64_t) occurences[i] == empty) empty_cells++;
			

		
		if(empty_cells)
			if(tries <= 4)
				gedmap_io::flush_row(to_string(empty_cells) + " entries missing. Start restoring:");
			else
				throw runtime_error("Not able to restore all entries of the index. May try using less threads (parameter -tc).");
		else {
			gedmap_io::print_row("Check index, OK");
			break;
		}
	}

	if(calc_count != chrom_count){
		gedmap_io::print_row("Merge occurences with occurences, that span over edges");
		omp_set_num_threads(MAX_THREAD_COUNT);
		#pragma omp parallel for schedule(dynamic,10)
		for(uint64_t i = 0; i < number_of_kmers; i++){
			uint64_t start	= table[i];
			uint64_t mid	= counts[ counts.size()-1 ] [i] + table[i];
			uint64_t end	= table[i+1];
			if(mid != end) std::inplace_merge(occurences.begin() + start, occurences.begin() + mid , occurences.begin() + end );
		}
	}

	eoc = linear_eoc_type(k,occurences,table);

	tv[tv_INDEX_IN].stop();
	return eoc;
}


//###############DECLARATION OF HELP METHODS##################
void add_letter_at_front(uint32_t k, vector<string> & K, char c);

template<typename kmer_type,class int_type, class off_stream_type>
void fill_in_kmers_for_position(uint64_t pos, uint32_t k, vector<string> & K, off_stream_type& off_stream, vector<uint32_t>& counts);

vector<string> merge( vector<string> & K1, vector<string> & K2);

void add_to_node_map(uint64_t pos, vector<string> & K, std::map<uint64_t,vector<string>> & node_map);
//#############################################################

/**
 * @brief calculates kmer index
 * kmer_type:  uint32_t for k <= 16; uint64_t for k <= 32
 */
template<typename kmer_type, class int_type, class off_stream_type>
void index_out(const uint32_t k, const string_view & eds, off_stream_type &off_stream, vector<uint32_t>& counts, uint64_t offset, std::map<uint64_t,vector<string>> & node_K){
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
				add_letter_at_front(k,K,c);
				fill_in_kmers_for_position<kmer_type,int_type>(i+offset,k,K, off_stream, counts);
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
 * @brief adds kmers that are going over edges
 * iterates over nodes in adj
 * for each node all conected nodes will be loaded
 * possible kmers vor each connected node will be added to possible kmers of this node
 * makes k steps and adds kmers to eoc
 */
template<typename kmer_type, class int_type, class off_stream_type>
void include_kmers_over_edges( const uint32_t k_c, const string & eds, off_stream_type& off_stream, vector<uint32_t>& counts, const adjacency & adj, std::map<uint64_t,vector<string>> & node_K ){
	if(!adj.initialised) return;
	for(std::map<uint64_t,uint64_t>::const_iterator node = adj.forward_edges.cend(); node!= adj.forward_edges.cbegin(); ){
		node--;
		uint64_t i, first_edge, last_edge;

		i = node->first;
		vector<string> K; 	//kmers on position i
		vector<string> K_i;	//kmers before entering alt scope
		vector<string> K_o; 	//kmers when hitting boarder to next alt

		//get kmers on #
		first_edge = node->second;
		std::map<uint64_t,uint64_t>::const_iterator next = node;
		next++;
		if(next != adj.forward_edges.cend()) last_edge = next->second-1;
		else last_edge = adj.forward_targets.size()-1;

		for(uint64_t j = first_edge; j <= last_edge; j++){ //for all predecessors
			uint64_t pred = adj.forward_targets[j];
			std::map<uint64_t,vector<string>>::iterator node_pred = node_K.find(pred); //find predecessors
			if(node_pred == node_K.end()) throw runtime_error("index_eoc, include_kmers_over_edges: predecessor not found - maybe eds and adj file dont match");
			K = merge(K,node_pred->second);
		}

		//start at pos i with k steps to go.
		uint32_t k = k_c-1;
		uint32_t k_i = 0;
		uint32_t k_o = 0;
		while( i > 0 && (k > 0 || k_i > 0 || k_o > 0) ){
			char c = eds[--i];
			switch(c){
				case 'A':
				case 'C':
				case 'G':
				case 'T':
				case 'N':{
					add_letter_at_front(k_c,K,c);
					if(k > 0){
						fill_in_kmers_for_position<kmer_type,int_type>(i,k_c,K, off_stream, counts);
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
 * @brief adds all kmers starting at position pos to index eoc
 *
 * kmer_type:  uint32_t for k <= 16; uint64_t for k <= 32
 */
template<typename kmer_type,class int_type, class off_stream_type>
void fill_in_kmers_for_position(uint64_t pos, uint32_t k, vector<string> & K, off_stream_type & off_stream, vector<uint32_t>& counts){
	vector<kmer_type> K_new = transform_string_vec_to_kmer_type_vec<kmer_type>(K, k);
	for(kmer_type kmer : K_new){
		off_stream.push_back(kmer);
		off_stream.push_back(pos);
		counts[kmer]++;
	}
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
