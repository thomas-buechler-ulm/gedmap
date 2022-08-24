#include "kmer.cpp"
#include "align_dp.cpp"

#include <algorithm>    // std::sort
#include <utility>      // std::pair, std::make_pair
#include <tuple>        // std::tuple, std::get, std::tie, std::ignore

#include "rsort.cpp"

using namespace std;
using namespace sdsl;


template <uint8_t t_width, typename int_type>
struct fasta_read{
	typedef tuple<int_type,int_type,int_type> int_type_triple;
	typedef tuple<string,uint32_t,int_type,bool> alignment_type;

	string	id;			//name of seq
	string	sequence;		//
	string	qual;

	//[EOC]QUERY KMERS
	//first  = id of kmer
	//second = start of this kmer in pattern
	vector<pair<int_type,int_type>> kmer_pairs;
	bool filled_1;

	//QUERY POSITIONS
	//first  = position in EDS
	//second = position in read
	vector<pair<int_type,int_type>> pos_pairs;
	bool filled_3;

	//QUERY HOTSPOTS
	vector< int_type_triple > hotspots;
	//0 positions in EDS that are in a hot spot
	//1 coresponding positions in pattern
	//2 number of seeds at this hotspot
	bool filled_4;

	//ALIGNMENTS
	vector<alignment_type> alignments;
	bool filled_5;

	fasta_read(){
		filled_1 = filled_3 = filled_4 = filled_5 = false;
	}
	
	fasta_read(string id ,string sequence, string qual):id(id),sequence(sequence), qual(qual){
		filled_1 = filled_3  = filled_4 = filled_5 = false;
	}

	fasta_read(ifstream & fastq_s){
		using namespace std;
		getline(fastq_s, id);
		if(!fastq_s.good()) throw runtime_error ("NO READ IN STREAM");
		getline(fastq_s, sequence);
		getline(fastq_s, qual);  //plus
		getline(fastq_s, qual); //qual
		if(!fastq_s.good()){
			gedmap_io::print_error("READ " + id + " incomplete");
			throw runtime_error ("NO READ IN STREAM");
		}
	}

	bool operator==(const fasta_read<t_width,int_type> & fr){
		return id == fr.id;
	}

	void delete_kmer_vectors(){
		vector<pair<int_type,int_type>>(0).swap(kmer_pairs);
		filled_1 = false;
	}

	void delete_pos_vectors(){
		vector<pair<int_type,int_type>>(0).swap(pos_pairs);
		filled_3 = false;
	}

	void delete_hints(){
		vector< tuple<int_type,int_type,int_type> >(0).swap(hotspots);
		filled_4 = false;
	}

	void delete_alignment(){
		vector<alignment_type>(0).swap(alignments);
		filled_5 = false;
	};

	/**
	 * @brief gets all seeds
	 *
	void get_seed_kmers(const uint32_t shift, uint32_t k){
		assert(k <= sequence.length());

		uint32_t seed_count = (sequence.length()-k) / shift+1;

		kmer_pairs	= vector<pair<int_type,int_type>>(seed_count);

		uint32_t i = 0;
		for(	int_type  seed_begin =  0;
			seed_begin 	+  k <= sequence.size();
			seed_begin 	+= shift)
		{
			kmer_pairs[i++] = make_pair<int_type,int_type>((int_type) KMER<int_type>::to_id(sequence.substr(seed_begin, k)),(int_type) seed_begin);
		}
		filled_1 = true;
	};*/
	
	
	void get_fragments(const uint32_t fragment_max_count, gedmap_mini::minimizer_index & mini){
		assert(mini.w + mini.k <= sequence.length());

		kmer_pairs	= vector<pair<int_type,int_type>>();
		kmer_pairs.reserve(fragment_max_count);

		vector<uint64_t> window(mini.w,-1);
		
		uint32_t kmer_start 	= 0;
		uint32_t next_c 		= mini.k;
		KMER<uint64_t> kmer(sequence.substr(0, mini.k));
		
		//initial window
		for(uint32_t i = 0; i < mini.w; i++){
			window[i] = gedmap_mini::hash(kmer.content);
			//go to next kmer
			kmer.rm_front();
			kmer.add_back_f(sequence[next_c++]);
			kmer_start++;
		}
		
		uint32_t wsp = 0; //window start pointer
		uint32_t ksp = 0; //kmer start pointer
		uint32_t last_kmer_pos = -1;
		
		while(true){
			//calc min
			uint32_t min_pos 	= 0;
			for(uint32_t i = 1; i < mini.w; i++){
				if(window[i] < window[min_pos]) 
					min_pos = i;
			}
			
				
			// min = window[min_pos]
			uint32_t kmer_pos = (mini.w + min_pos - wsp) % mini.w + ksp;
			int_type kmer_id = gedmap_mini::hash_inverse(window[min_pos]);
			
			//add to kmer pairs
			if(last_kmer_pos != kmer_pos && mini.indicator[kmer_id]){
				last_kmer_pos = kmer_pos;
				kmer_pairs.push_back(make_pair<int_type,int_type>((int_type) kmer_id ,(int_type) kmer_pos));
			}
			
			if(kmer_pairs.size() >= fragment_max_count || next_c >= sequence.size())
				break;
			
			//add next kmer to window
			window[wsp] = gedmap_mini::hash(kmer.content);
			wsp = (wsp+1) % mini.w;
			ksp++;	
			//go to next kmer
			kmer.rm_front();
			kmer.add_back_f(sequence[next_c++]);
			kmer_start++;
		}
		filled_1 = true;
	}
	
	/**
	 * @brief generate fragments from sequence
	 */
	void get_fragments(const uint32_t fragment_count, const uint32_t k){
		assert(k <= sequence.length());
		
		/*if( k*fragment_count >=  sequence.length()){
			kmer_pairs	= vector<pair<int_type,int_type>>(fragment_count);
			uint32_t fragment_begin = 0;
			for(uint32_t i = 0; i < fragment_count; i++){
				kmer_pairs[i] = make_pair<int_type,int_type>((int_type) KMER<int_type>::to_id(sequence.substr(fragment_begin, k)),(int_type) fragment_begin);
				fragment_begin += k;
			}
			return;
		}*/
		
		uint32_t this_fragment_count = fragment_count;
		
		if(  (sequence.length() - k)/ (double) fragment_count   < 1){
			gedmap_io::print_error("Not possible to generate that many fragments of the sequence. Set fragment_count value to (read.length - k)");
			this_fragment_count = (sequence.length() - k);
		}
		
		kmer_pairs	= vector<pair<int_type,int_type>>(this_fragment_count);
		
		for(uint32_t i = 0; i < this_fragment_count; i++){			
			uint32_t fragment_begin = i* ((sequence.length() - k)/ (double) this_fragment_count);
			kmer_pairs[i] = make_pair<int_type,int_type>((int_type) KMER<int_type>::to_id(sequence.substr(fragment_begin, k)),(int_type) fragment_begin);
		}
		
		filled_1 = true;
	};
	
	/**
	 * @brief generate fragments from sequence
	 */
	void get_fragments(const uint32_t fragment_count, const uint32_t k, sdsl::bit_vector & indicator){
		assert(k <= sequence.length());
		
		
		uint32_t this_fragment_count = 0;
		kmer_pairs	= vector<pair<int_type,int_type>>(fragment_count);
		
		for(uint32_t i = 0; i < k && this_fragment_count < fragment_count; i++){
			for(uint32_t fragment_begin = i; fragment_begin < (sequence.length() - k) && this_fragment_count < fragment_count; fragment_begin += k){
				int_type kmer_id = KMER<int_type>::to_id(sequence.substr(fragment_begin, k));
				if( indicator[kmer_id] ){
					kmer_pairs[this_fragment_count++] = make_pair<int_type,int_type>( (int_type) kmer_id , (int_type) fragment_begin);
				}
			}
		}
		
		if( this_fragment_count < fragment_count){
			kmer_pairs.resize(this_fragment_count);
// 			gedmap_io::print_error("Only generated " + to_string(this_fragment_count) + " fragment for read " + id);
		}
		
		filled_1 = true;
	};

	/**
	 * @brief queries the fragment positions
	 * @param mini the minimizer index
	 */
	void get_positions(gedmap_mini::minimizer_index & mini){
		pos_pairs.reserve(2*kmer_pairs.size());
		vector<uint32_t> boundaries(kmer_pairs.size()+1);
		uint32_t b_i = 0;
		boundaries[b_i++] = 0;
		for(auto it = kmer_pairs.begin(); it != kmer_pairs.end(); it++){
			int_type k = get<0>(*it);
			size_t begin, count;
			tie(begin,count) = mini.boundaries(k);
			if(count > 0){
				uint64_t old_size = pos_pairs.size();
				pos_pairs.resize(old_size + count);

				//COPY VALUES FROM KMI
				for( uint32_t j  = 0; j < count; j++)
					pos_pairs[old_size + j] = make_pair(mini.positions[begin + j], get<1>(*it));
				assert(std::is_sorted(pos_pairs.begin() + old_size, pos_pairs.begin() + old_size + count));

				boundaries[b_i++] = pos_pairs.size();
			}
		}
// 		std::sort(pos_pairs.begin(),pos_pairs.end());
		boundaries.resize(b_i);
		multi_merge(pos_pairs.begin(),boundaries);
		assert(std::is_sorted(pos_pairs.begin(), pos_pairs.end()));
		delete_kmer_vectors();
		filled_3 = true;
	};
	
	
	/**
	 * @brief queries the fragment positions
	 * @param eoc the EOC-index
	 */
	void get_positions(linear_eoc_type & eoc){
		vector<uint32_t> boundaries(kmer_pairs.size()+1);
		uint32_t b_i = 0;
		boundaries[b_i++] = 0;
		for(auto it = kmer_pairs.begin(); it != kmer_pairs.end(); it++){
			int_type k = get<0>(*it);
			uint32_t begin = eoc.table[k];
			uint32_t count = eoc.table[k+1] - begin;
			if(count > 0){
				uint64_t old_size = pos_pairs.size();
				pos_pairs.resize(old_size + count);

				//COPY VALUES FROM KMI
				for( uint32_t j  = 0; j < count; j++) pos_pairs[old_size + j] = make_pair(eoc.occurences[begin + j], get<1>(*it)) ;

				boundaries[b_i++] = pos_pairs.size();
			}
		}
		boundaries.resize(b_i);
		multi_merge(pos_pairs.begin(),boundaries);
		delete_kmer_vectors();
		filled_3 = true;
	};
	/**
	 * @brief queries the seed positions in the batch
	 * @param eoc 		the EOC vector, eoc[1] is seperator, that seperates kmers
	 * @param batch	vector containing the patterns
	 */
	//template<class eoc_type>
	//static void get_positions_in_batch_eoc(vector<fasta_read<t_width,int_type>> & batch, eoc_type & eoc);


	/**
	* @brief performs hotspot finding
	* 
	* @param window_size
	* @param window_hits 
	* @param check_col check for colliniarity
	*/
	void find_hotspots(const uint32_t window_size, const uint32_t window_hits, bool check_col);


	/**
	* @brief aligns hotspots to EDS
	* 
	* @param EDS	pangenome
	* @param D	max Dist 
	* @param max_a_c maximum number of completed alignments
	* @param max_a_t maximum number of started alignments
	* @param max_a_c maximum number of alignments in output
	* 
	*/
	void start_aligner(const string & EDS, const adjacency & adj, uint32_t D, uint32_t max_a_c, uint32_t max_a_t, uint32_t max_a_o){
		fasta_read<t_width,int_type> DUMMY_FR;
		start_aligner(EDS,adj,D,max_a_c,max_a_t,max_a_o,DUMMY_FR);
	};

	/**
	* @brief aligns hotspots to EDS with considering a reverse complement
	* 
	* @param EDS	pangenome
	* @param D	max Dist 
	* @param max_a_c maximum number of completed alignments
	* @param max_a_t maximum number of started alignments
	* @param max_a_c maximum number of alignments in output
	* @param rev_compl	reverse complement of read
	* 
	*/
	void start_aligner(const string & EDS, const adjacency & adj, uint32_t D, uint32_t max_a_c, uint32_t max_a_t, uint32_t max_a_o, fasta_read<t_width,int_type> & rev_compl);	

	/**
	* @brief writes alignment in file 
	* @param ofs	ofstream
	* @param write_failure	if true a line with name and stars in written if read could not be aligned
	* @param transform		transforms position in EDS to position in FA
	* @return number of alignments in output
	*/
	uint32_t write_alignment(ofstream & ofs, bool write_failure, const pos_EDS_to_FA_type & transform){
		using namespace std;
		assert(filled_5);
		
		uint32_t res = alignments.size();
		
		
		string 	QNAME	= id.substr(1);
		uint32_t 	FLAG		= 0;
		string	RNAME	= "*";
		uint32_t 	POS		= 0;
		uint32_t	MAPQ 	= 255; //TODO
		string 	CIGAR	= "*";
		string 	RNEXT 	= "*";
		uint32_t 	PNEXT 	= 0;
		uint32_t 	TLEN		= 0;
		string	SEQ		= sequence;
		string 	QUAL		= qual;
		
		//FOR FLAG
		uint32_t UNMAPPED = 4;
		uint32_t REV_COMP	= 16;
		uint32_t NOT_PRIM	= 256;
		
		for(uint32_t i = 0; i < alignments.size(); i++){
			RNAME 	= "*";
			POS 		= get<2>(alignments[i])+1;
			CIGAR 	= get<0>(alignments[i]);
			SEQ		= sequence;
			QUAL		= qual;
			
			FLAG = 0;
			if( i > 0 )	
				FLAG += NOT_PRIM;
			
			if(get<3>(alignments[i])){
			 	FLAG += REV_COMP;
				SEQ = gedmap_encode::rev_complement(SEQ);
				QUAL = gedmap_encode::rev_complement(QUAL,false);
			}
			uint32_t off = 0;
			if(!transform.empty()) tie(RNAME,POS,off) = transform(POS);
			
				
			ofs 	<< QNAME << '\t' << FLAG << '\t' << RNAME << '\t' << POS << '\t' << MAPQ << '\t' << CIGAR << '\t' << RNEXT << '\t' << PNEXT << '\t' << TLEN << '\t' << SEQ << '\t' << QUAL << '\t'
				<< "NM:i:" << ((uint32_t)get<1>(alignments[i])) <<  (off?(" XO:i:" + to_string(off)):"") << endl;
		}
		
		if(!alignments.size() && write_failure)
			ofs 	<< QNAME << '\t' << UNMAPPED << '\t' << RNAME << '\t' << POS << '\t' << MAPQ << '\t' << CIGAR << '\t' << RNEXT << '\t' << PNEXT << '\t' << TLEN << '\t' << SEQ << '\t' << QUAL <<  endl;
		delete_alignment();
		return res;
	};
};

/**
template <uint8_t t_width, class int_type>
template<class eoc_type>
void fasta_read<t_width,int_type>::get_positions_in_batch_eoc(vector<fasta_read<t_width,int_type>> & batch, eoc_type & eoc){
	sys_timer batch_time;
	batch_time.start();
	//calc numer of intervals
	uint32_t kmer_count = 0;
	for(auto r : batch){
		assert(r.filled_1);
		kmer_count += r.kmer_pairs.size();
	}

	vector<int_type_triple> all_kmers = vector<int_type_triple>(kmer_count);

	uint32_t k = 0;
	for(uint32_t i = 0; i < batch.size(); i++){
		for(uint32_t j = 0; j < batch[i].kmer_pairs.size(); j++){
			all_kmers[k++] =(  make_tuple<>(batch[i].kmer_pairs[j].first, batch[i].kmer_pairs[j].second , i )   );
		}
		batch[i].delete_kmer_vectors();
	}

	for(size_t i = 0; i < batch.size(); i++){
		batch[i].pos_pairs.reserve(edsm_align::DEFAULT_POSPAIR_CAPACITY);
	}
	
	sort(all_kmers.begin(),
		all_kmers.end(), 
		[] ( int_type_triple const& t1, int_type_triple const& t2) -> bool{ return get<0>(t1) < get<0>(t2); } //SORT BY HITS ID
	);

	
	int_type sep = eoc[1];
	size_t eoc_p = 1;
	int_type eoc_id = 0;
	int_type i = 0;
	while(true){
		
		eoc_p++;		
		int_type p = eoc[eoc_p];
		while(p == sep && eoc_p < eoc.size()){
			eoc_id++;
			eoc_p++;
			p = eoc[eoc_p];
		}
		if(eoc_p == eoc.size()) break;			
			
		
		for(; get<0>(all_kmers[i]) < eoc_id && i < all_kmers.size() ; i++);
		
		for(size_t j = i; get<0>(all_kmers[j]) == eoc_id; j++)
			batch[ get<2>(all_kmers[j]) ].pos_pairs.push_back( make_pair(p,get<1>(all_kmers[j]) ));
	}
}*/


template <uint8_t t_width, class int_type>
void fasta_read<t_width,int_type>::find_hotspots(const uint32_t window_size, const uint32_t window_hits, bool check_col){
 	assert(filled_3);//positions must be availablee in vector
	 
	// if (window_hits <= 1) check_col = false;

	//boundaries of search window
	uint32_t l = 0;
	uint32_t r = 0;
	uint32_t np = 0; // next position than can be included in hotspot
	
	// remove duplicates from pos_pairs (TODO: why is this necessary?)
	pos_pairs.erase( std::unique(pos_pairs.begin(), pos_pairs.end()), pos_pairs.end() );

#ifndef NDEBUG
	assert(std::is_sorted(pos_pairs.begin(), pos_pairs.end(),
			[] (auto const& t1, auto const& t2) { return get<0>(t1) < get<0>(t2); }));
	for(uint32_t i = 1; i < pos_pairs.size(); i++)
		assert( get<0>(pos_pairs[i]) != get<0>(pos_pairs[i-1]) || get<1>(pos_pairs[i]) != get<1>(pos_pairs[i-1])  );
#endif

	hotspots.reserve(pos_pairs.size());

	std::vector<uint32_t> lis;
	while( l < pos_pairs.size() && r < pos_pairs.size()) {
		//shift r as far to the right as possible, i.e. [pos_i[l],pos_i[r]] < window_size
		while(r + 1 < pos_pairs.size() && pos_pairs[r+1].first < pos_pairs[l].first + window_size) 
			r++;

		//if interval is big enough
		if(r-l+1 >= window_hits) {
			if(check_col){
#ifndef NDEBUG
				// cl_v[i] := number of smaller positions in pos_s[l..l+i]
				std::vector<uint32_t> cl_v(r-l+1,0);
				//if (l >= np && window_hits <= 1)
				//	hotspots.emplace_back(pos_pairs[l].first, pos_pairs[l].second, 1);
				for(uint32_t i = 1; i < cl_v.size(); i++){
					for(uint32_t j = 0; j < i; j++)
						if(pos_pairs[l+j].second < pos_pairs[l+i].second && cl_v[i] <= cl_v[j] )
							cl_v[i] = cl_v[j] + 1;
					//if(i + l >= np && cl_v[i] + 1 >= window_hits)
					//	hotspots.push_back( make_tuple<>(pos_pairs[l+i].first, pos_pairs[l+i].second, cl_v[i]+1) );
				}
#endif
				uint32_t s = 0;
				if (lis.size() < r-l+1) lis.resize(r-l+1);
				for (uint32_t i = l; i <= r; i++) {
					const auto it = std::lower_bound(lis.begin(), lis.begin() + s, pos_pairs[i].second);
					const uint32_t p = std::distance(lis.begin(), it);
					assert(p == cl_v[i - l]);
					*it = pos_pairs[i].second;
					if (p == s) s++;
					if (i >= np && p + 1 >= window_hits)
						hotspots.emplace_back( pos_pairs[i].first, pos_pairs[i].second, p + 1 );
				}
			}else {
				if (np < l) np = l;
				for(uint32_t i = np; i <= r; i++)
					hotspots.push_back( make_tuple<>(pos_pairs[i].first, pos_pairs[i].second, r-l+1) );
			}
		}

		//shift r by one
		r++;
		np = r;
		if(r >= pos_pairs.size())
			break;

		//shift l so that [pos_i[l],pos_i[r]] < window_size holds
		while(pos_pairs[l].first + window_size  < pos_pairs[r].first) 
			l++;	
	}
	delete_pos_vectors();
	
#ifndef NDEBUG 
	assert(std::is_sorted(hotspots.begin(), hotspots.end(),
			[] (auto const& t1, auto const& t2) { return get<0>(t1) < get<0>(t2); }));
	for(uint32_t i = 1; i < hotspots.size(); i++)
		assert( get<0>(hotspots[i]) != get<0>(hotspots[i-1]) || get<1>(hotspots[i]) != get<1>(hotspots[i-1])  );
#endif
	//INVERSE SORT BY HITS
//	rsort<true>(hotspots.data(), hotspots.data() + hotspots.size(), [] (const auto& v) { return get<2>(v); });
	sort(hotspots.begin(), hotspots.end(), 
		[] ( int_type_triple const& t1, int_type_triple const& t2) -> bool{ return get<2>(t1) > get<2>(t2); } );
	
	filled_4 = true;
};

template <uint8_t t_width, class int_type>
void fasta_read<t_width,int_type>::start_aligner(const string & EDS, const adjacency & adj,  uint32_t D, uint32_t max_a_c, uint32_t max_a_t, uint32_t max_a_o, fasta_read<t_width,int_type> & r_c) {
	assert(D <= std::numeric_limits<uint8_t>::max());
	assert(filled_4);
	
	uint32_t aligns 	= 0;
	uint32_t tries	= 0;
	uint32_t i 	 	= 0;
	uint32_t i_r_c	= 0;
	
	std::vector<std::tuple<uint32_t, int_type, int_type, bool>> tmp_alignments;
	//TODO AUTOMATIC DECREASE OF D
	if(!alignments.empty() && get<1>(alignments[0]) < D ) D = get<1>(alignments[0]);

	//	map storing pair of position and error count
// 	std::map<int_type,uint32_t> pos_d; //TODO

	while(aligns < max_a_c	// only calculate max_a_c alignments
		&& tries < max_a_t	// only start aligner max_a_t times
		&& (i < hotspots.size() || i_r_c < r_c.hotspots.size() ) // another hint left
		)
	{		
		bool align_r_c = (i == hotspots.size()) // no further read in h_pos_i
					|| (i_r_c < r_c.hotspots.size()  && get<2>(hotspots[i]) < get<2>(r_c.hotspots[i_r_c])); // or next rev komp hint has more window hits 
		{
			tries++;
			const auto [eds_pos, read_pos] = align_r_c // position in eds and read, respectively
					? std::make_pair(std::get<0>(r_c.hotspots[i_r_c]), std::get<1>(r_c.hotspots[i_r_c]))
					: std::make_pair(std::get<0>(    hotspots[i    ]), std::get<1>(    hotspots[i    ]));
			const uint32_t min_dist = align_r_c
				? align<false, uint8_t, uint32_t>(EDS, adj, eds_pos, r_c.sequence, read_pos, D)
				: align<false, uint8_t, uint32_t>(EDS, adj, eds_pos,     sequence, read_pos, D);

			if constexpr (false) {
				const auto [cig,dist,offset] = align_r_c
					? edsm_levinstein::levinstein<int_type>(EDS, adj, r_c.sequence, D, eds_pos, read_pos)
					: edsm_levinstein::levinstein<int_type>(EDS, adj ,sequence, D, eds_pos, read_pos);

				if ((offset == EDS.size()) != (min_dist > D) || (min_dist <= D && min_dist != dist)) {
					#pragma omp critical
					{
						const string_view read = (align_r_c ? r_c.sequence : sequence);
						cerr << "D = " << D << " eds_pos = " << eds_pos << " read_pos = " << read_pos << endl;
						cerr << "old says " << (offset == EDS.size() ? "impossible"s : to_string(dist))
							<< ", new says " << (min_dist > D ? "impossible"s : to_string(min_dist)) << endl;
						cerr << (align_r_c ? "reverse" : "not reverse") << endl;
						cerr << "READ = " << read << endl;
						cerr << "old cig = " << cig << endl;

						cerr << "FORWARD:\n ";
						cerr << "\tread: " << read.substr(read_pos) << endl;
						cerr << "\teds : " << EDS.substr(eds_pos, 2 * read.size()) << endl;

						cerr << "BACKWARD:\n ";
						cerr << "\tread: " << read.substr(0, read_pos) << " | " << read.substr(read_pos, 10) << endl;
						cerr << "\teds : " << EDS.substr(eds_pos - 2 * read.size(), 2 * read.size()) << " | " << EDS.substr(eds_pos, 10)  << endl;

						exit(1);
					}
				}
				if (min_dist <= D) {
					const auto[cigar, pos] = align_r_c
							? align<true, uint8_t, uint32_t>(EDS, adj, eds_pos, r_c.sequence, read_pos, min_dist)
							: align<true, uint8_t, uint32_t>(EDS, adj, eds_pos,     sequence, read_pos, min_dist);
					const uint32_t dist2 = align_r_c
						? align<false, uint8_t, uint32_t>(EDS, adj, pos, r_c.sequence, 0, D)
						: align<false, uint8_t, uint32_t>(EDS, adj, pos,     sequence, 0, D);
					if (dist2 > min_dist) {
						#pragma omp critical
						{
							cerr << "eds_pos = " << (uint64_t)eds_pos << " read_pos = " << (uint64_t)read_pos << " pos = " << (uint64_t)pos << " min_dist = " << min_dist << " dist2 = " << dist2 << endl;

							exit(1);
						}
					}
				}
			}
			
			if (min_dist <= D) {
				// ali found
				aligns++;
				tmp_alignments.emplace_back(min_dist, eds_pos, read_pos, align_r_c);

					
				//TODO AUTOMATIC DECREASE OF D
				if( min_dist < D) D = min_dist;
			}
		}
		
		if(align_r_c)
			i_r_c++;
		else
			i++;
	}
	//rsort(alignments.data(), alignments.data() + alignments.size(),  [] (const auto& v) { return get<1>(v); });
	//SORT BY D
 	sort(tmp_alignments.begin(), tmp_alignments.end(),
 		[] ( auto const& t1, auto const& t2) { return get<0>(t1) < get<0>(t2); } );
	
	const auto report_count = min(aligns, max_a_o);
	tmp_alignments.resize(report_count);

	alignments.resize(report_count);

	for (size_t i = 0; i < report_count; i++) {
		const auto& [min_dist, eds_pos, read_pos, align_r_c] = tmp_alignments[i];
		const auto[cigar, pos] = align_r_c
				? align<true, uint8_t, uint32_t>(EDS, adj, eds_pos, r_c.sequence, read_pos, min_dist)
				: align<true, uint8_t, uint32_t>(EDS, adj, eds_pos,     sequence, read_pos, min_dist);
		alignments[i] = std::make_tuple(cigar, min_dist, pos, align_r_c);
	}
	tmp_alignments.clear();

	delete_hints();

	filled_5 = true;
};
