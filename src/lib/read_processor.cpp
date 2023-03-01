#include "kmer.cpp"
#include "align_dp.cpp"

#include <algorithm>    // std::sort
#include <utility>      // std::pair, std::make_pair
#include <tuple>        // std::tuple, std::get, std::tie, std::ignore

#include "rsort.cpp"

using namespace std;
using namespace sdsl;

//template<typename int_type>
//using kmer_pair = named_type< std::pair<int_type,int_type>, struct KmerPairTag >;
template<typename int_type>
struct kmer_pair {
	int_type id = 0; // id of kmer
	int_type start = 0; // start of this kmer in pattern
	constexpr kmer_pair(int_type id, int_type start) : id(id), start(start) {}
	constexpr kmer_pair() = default;
};
template<typename int_type>
struct query_position {
	int_type eds_pos = 0; // position in EDS
	int_type read_pos = 0; // position in read
	constexpr query_position(int_type eds_pos, int_type read_pos)
		: eds_pos(eds_pos), read_pos(read_pos) {}
	constexpr query_position() = default;
	auto operator<=>(const query_position<int_type>& rhs) const = default;
};
template<typename int_type >
struct hotspot {
	int_type eds_pos = 0; // position in eds
	int_type read_pos = 0; // position in read
	int_type quality = 0; // higher is better. usually the length of LIS ending here
	constexpr hotspot(int_type eds_pos, int_type read_pos, int_type quality)
		: eds_pos(eds_pos), read_pos(read_pos), quality(quality) {}
	constexpr hotspot() = default;
};
// result of running the alignment algorithm on a read. Contains necessary data
// to (quickly) reconstruct the alignment incl. cigar string
template<typename int_type>
struct temp_alignment {
	std::pair<uint32_t,uint32_t> dist; // (dist_left, dist_right)
	int_type eds_pos; // eds pos
	int_type read_pos; // read_pos
	uint32_t source; // source
	bool reverse_compl; // reverse_compl
	constexpr uint32_t get_dist() const { return dist.first + dist.second; }
	constexpr temp_alignment() = default;
	constexpr temp_alignment(std::pair<uint32_t,uint32_t> dist, int_type eds_pos, int_type read_pos, bool reverse_compl, uint32_t source)
		: dist(dist), eds_pos(eds_pos), read_pos(read_pos), source(source), reverse_compl(reverse_compl) {}
};
template<typename int_type>
struct alignment {
	std::string cigar;
	uint32_t dist;
	int_type pos;
	bool align_r_c;
	constexpr alignment(
		std::string&& cigar,
		uint32_t dist,
		int_type pos,
		bool align_r_c)
		: cigar(std::move(cigar))
		, dist(dist)
		, pos(pos)
		, align_r_c(align_r_c)
	{ }
	constexpr alignment() = default;
};

template <typename int_type>
struct fasta_read;

namespace read_processor {
	using dist_type = uint8_t;

	template< typename int_type >
	vector<kmer_pair<int_type>>
	get_fragments(
		const fasta_read<int_type>& read,
		const uint32_t fragment_max_count,
		const gedmap_mini::minimizer_index & mini
	) {
		assert(mini.w + mini.k <= read.sequence.length());

		vector<kmer_pair<int_type>> kmer_pairs;
		kmer_pairs.reserve(fragment_max_count);

		uint32_t next_c 		= mini.k;
		KMER<uint64_t> kmer(std::string_view(read.sequence).substr(0, mini.k));
		
		uint32_t ksp = 0; //kmer start pointer
		uint32_t last_kmer_pos = -1;
		// Assume we have two elements v_i > v_j in the current window with i < j
		// Then we can discard v_i, since v_j will be in the window at least as
		// long as v_i, and v_j is better than v_i
		// -> We maintain a deque that is sorted in increasing order according
		// to both the index and the value (i.e. the minimum is always at the
		// front). When adding a new element to the back, we first remove all
		// elements that are worse that it. They are of course all at the back
		// of the deque, so we have amortised constant time complexity per
		// add-operation.
		std::deque<std::pair<uint64_t, uint32_t>> cur;
		const auto add = [&] (uint64_t val, uint32_t i) {
			while (!cur.empty() and cur.front().second + mini.w <= i) cur.pop_front();
			while (!cur.empty() and cur.back().first > val) cur.pop_back();
			cur.emplace_back(val, i);
		};
		
		//initial window
		for(uint32_t i = 0; i < mini.w; i++){
			add(gedmap_mini::hash(kmer.content), i);
			//go to next kmer
			kmer.rm_front();
			kmer.add_back_f(read.sequence[next_c++]);
		}
		while(true){
			//calc min
			
			const auto[kmer_id_, kmer_pos] = cur.front();
			const int_type kmer_id = gedmap_mini::hash_inverse(kmer_id_);
			
			//add to kmer pairs
			if(last_kmer_pos != kmer_pos && mini.is_in_index(kmer_id)){
				last_kmer_pos = kmer_pos;
				kmer_pairs.emplace_back((int_type) kmer_id, (int_type) kmer_pos);
			}
			
			if (kmer_pairs.size() >= fragment_max_count or next_c >= read.sequence.size())
				break;
			
			//add next kmer to window
			add(gedmap_mini::hash(kmer.content), ksp + mini.w);
			ksp++;
			//go to next kmer
			kmer.rm_front();
			kmer.add_back_f(read.sequence[next_c++]);
		}
		return kmer_pairs;
	}
	
	template< typename int_type >
	vector<kmer_pair<int_type>>
	get_fragments(
		const fasta_read<int_type>& read,
		const uint32_t fragment_count,
		const uint32_t k
	) {
		assert(k <= read.sequence.length());
		
		/*if( k*fragment_count >=  read.sequence.length()){
			kmer_pairs	= vector<pair<int_type,int_type>>(fragment_count);
			uint32_t fragment_begin = 0;
			for(uint32_t i = 0; i < fragment_count; i++){
				kmer_pairs[i] = make_pair<int_type,int_type>((int_type) KMER<int_type>::to_id(read.sequence.substr(fragment_begin, k)),(int_type) fragment_begin);
				fragment_begin += k;
			}
			return;
		}*/
		
		uint32_t this_fragment_count = fragment_count;
		
		if(  (read.sequence.length() - k) < fragment_count) [[unlikely]] {
			gedmap_io::print_error("Not possible to generate that many fragments of the sequence. Set fragment_count value to (read.length - k)");
			this_fragment_count = (read.sequence.length() - k);
		}
		
		vector<kmer_pair<int_type>> kmer_pairs(this_fragment_count);
		
		for(uint32_t i = 0; i < this_fragment_count; i++){			
			uint32_t fragment_begin = i * (read.sequence.length() - k) / this_fragment_count;
			kmer_pairs[i] = kmer_pair<int_type>
				( (int_type) KMER<int_type>::to_id(std::string_view(read.sequence).substr(fragment_begin, k))
				, (int_type) fragment_begin
				);
		}
		
		return kmer_pairs;
	}
	
	template< typename int_type >
	vector<kmer_pair<int_type>>
	get_fragments(
		const fasta_read<int_type>& read,
		const uint32_t fragment_count,
		const uint32_t k,
		const sdsl::bit_vector & indicator
	) {
		assert(k <= read.sequence.length());
		
		uint32_t this_fragment_count = 0;
		vector<kmer_pair<int_type>> kmer_pairs(fragment_count);
		
		for(uint32_t i = 0; i < k && this_fragment_count < fragment_count; i++){
			for(uint32_t fragment_begin = i; fragment_begin < (read.sequence.length() - k) && this_fragment_count < fragment_count; fragment_begin += k){
				int_type kmer_id = KMER<int_type>::to_id(std::string_view(read.sequence).substr(fragment_begin, k));
				if( indicator[kmer_id] ){
					kmer_pairs[this_fragment_count++] = kmer_pair<int_type>( (int_type) kmer_id , (int_type) fragment_begin);
				}
			}
		}
		
		if( this_fragment_count < fragment_count){
			kmer_pairs.resize(this_fragment_count);
			// 			gedmap_io::print_error("Only generated " + to_string(this_fragment_count) + " fragment for read " + id);
		}
		
		return kmer_pairs;
	}
	
	/**
	 * @brief queries the fragment positions
	 * @param mini the minimizer index
	 */
	template< typename int_type >
	vector<query_position<int_type>>
	get_positions(
		const vector<kmer_pair<int_type>>& kmer_pairs,
		const gedmap_mini::minimizer_index & mini
	) {
		vector<query_position<int_type>> pos_pairs;
		pos_pairs.reserve(2*kmer_pairs.size());
		vector<uint32_t> boundaries(kmer_pairs.size()+1);
		size_t b_i = 0;
		boundaries[b_i++] = 0;
		for (const auto&[id, start] : kmer_pairs) {
			auto[begin,count] = mini.get_position(id);
			if(count > 0){
				const auto old_size = pos_pairs.size();
				pos_pairs.resize(old_size + count);

				if (count == 1)
					pos_pairs[old_size] = query_position<int_type>(begin, start);
				else
					for (size_t j  = 0; j < count; j++)
						pos_pairs[old_size + j] = query_position<int_type>(mini.pos_mult[begin + j], start);
				assert(std::is_sorted(pos_pairs.begin() + old_size, pos_pairs.begin() + old_size + count));

				boundaries[b_i++] = pos_pairs.size();
			}
		}
		boundaries.resize(b_i);
		multi_merge(pos_pairs.begin(), boundaries);
		assert(std::is_sorted(pos_pairs.begin(), pos_pairs.end()));
		return pos_pairs;
	}

	/**
	* @brief performs hotspot finding
	* 
	* @param window_size
	* @param window_hits 
	* @param check_col check for colliniarity
	*/
	template < typename int_type >
	vector<hotspot<int_type>>
	find_hotspots(
		const vector<query_position<int_type>>& pos_pairs,
		const uint32_t window_size,
		const uint32_t window_hits,
		bool check_col
	) {
		// if (window_hits <= 1) check_col = false;
	
		//boundaries of search window
		uint32_t l = 0;
		uint32_t r = 0;
		uint32_t np = 0; // next position than can be included in hotspot
	
	#ifndef NDEBUG
		// assume pos_pairs is sorted (after first component)
		assert(std::is_sorted(pos_pairs.begin(), pos_pairs.end(),
				[] (auto const& t1, auto const& t2) { return t1.eds_pos < t2.eds_pos; }));
		// assume there are no duplicates
		for(uint32_t i = 1; i < pos_pairs.size(); i++)
			assert( pos_pairs[i] != pos_pairs[i-1] );
	#endif
	
		vector<hotspot<int_type>> hotspots;
		hotspots.reserve(pos_pairs.size());
	
		std::vector<uint32_t> lis;
		while( l < pos_pairs.size() && r < pos_pairs.size()) {
			//shift r as far to the right as possible, i.e. [pos_i[l],pos_i[r]] < window_size
			while(r + 1 < pos_pairs.size() && pos_pairs[r+1].eds_pos < pos_pairs[l].eds_pos + window_size) 
				r++;
	
			//if interval is big enough
			if(r-l+1 >= window_hits) {
				if(check_col){
	#define FAST_LIS
	#if !defined(NDEBUG) || !defined(FAST_LIS)
					// cl_v[i] := number of smaller positions in pos_s[l..l+i]
					std::vector<uint32_t> cl_v(r-l+1,0);
		#ifndef FAST_LIS
					if (l >= np && window_hits <= 1)
						hotspots.emplace_back(pos_pairs[l].eds_pos, pos_pairs[l].read_pos, 1);
		#endif
					for(uint32_t i = 1; i < cl_v.size(); i++){
						for(uint32_t j = 0; j < i; j++)
							if(pos_pairs[l+j].read_pos < pos_pairs[l+i].read_pos && cl_v[i] <= cl_v[j] )
								cl_v[i] = cl_v[j] + 1;
		#ifndef FAST_LIS
						if(i + l >= np && cl_v[i] + 1 >= window_hits)
							hotspots.emplace_back( pos_pairs[l+i].eds_pos, pos_pairs[l+i].read_pos, cl_v[i]+1 );
		#endif
					}
	#endif
	#ifdef FAST_LIS
					uint32_t s = 0;
					if (lis.size() < r-l+1) lis.resize(r-l+1);
					for (uint32_t i = l; i <= r; i++) {
						const auto it = std::lower_bound(lis.begin(), lis.begin() + s, pos_pairs[i].read_pos);
						const uint32_t p = std::distance(lis.begin(), it);
						assert(p == cl_v[i - l]);
						*it = pos_pairs[i].read_pos;
						if (p == s) s++;
						if (i >= np && p + 1 >= window_hits)
							hotspots.emplace_back( pos_pairs[i].eds_pos, pos_pairs[i].read_pos, p + 1 );
					}
	#endif
				}else {
					if (np < l) np = l;
					for(uint32_t i = np; i <= r; i++)
						hotspots.emplace_back( pos_pairs[i].eds_pos, pos_pairs[i].read_pos, r-l+1 );
				}
			}
	
			//shift r by one
			r++;
			np = r;
			if(r >= pos_pairs.size())
				break;
	
			//shift l so that [pos_i[l],pos_i[r]] < window_size holds
			while(pos_pairs[l].eds_pos + window_size  < pos_pairs[r].eds_pos) 
				l++;	
		}
		
	#ifndef NDEBUG 
		assert(std::is_sorted(hotspots.begin(), hotspots.end(),
				[] (auto const& t1, auto const& t2) { return t1.eds_pos < t2.eds_pos; }));
		for(uint32_t i = 1; i < hotspots.size(); i++)
			assert( hotspots[i].eds_pos != hotspots[i-1].eds_pos || hotspots[i].read_pos != hotspots[i-1].read_pos  );
	#endif
		//INVERSE SORT BY HITS
		//	rsort<true>(hotspots.data(), hotspots.data() + hotspots.size(), [] (const auto& v) { return get<2>(v); });
		sort(hotspots.begin(), hotspots.end(), 
			[] ( const auto& t1, const auto& t2) { return t1.quality > t2.quality; } );
		return hotspots;
	}
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
	template < typename int_type >
	std::vector<temp_alignment<int_type>>
	align(
		const std::vector<hotspot<int_type>>& hotspots,
		const fasta_read<int_type>& read,
		const std::vector<hotspot<int_type>>& hotspots_r_c,
		const fasta_read<int_type>& read_r_c,
		const std::string & EDS,
		const adjacency & adj,
		uint32_t D,
		const uint32_t max_a_c, // max alignments calculated
		const uint32_t max_a_t // max tries
	) {
		assert(D <= std::numeric_limits<dist_type>::max());
		
		uint32_t aligns 	= 0;
		uint32_t tries	= 0;
		uint32_t i 	 	= 0;
		uint32_t i_r_c	= 0;
		
		std::vector<temp_alignment<int_type>> tmp_alignments;

		//if(!alignments.empty() && get<1>(alignments[0]) < D ) D = get<1>(alignments[0]);
	
		while (aligns < max_a_c	// only calculate max_a_c alignments
			&& tries < max_a_t	// only start aligner max_a_t times
			&& (i < hotspots.size() || i_r_c < hotspots_r_c.size() ) // another hint left
		) {		
			bool align_r_c = (i == hotspots.size()) // no further read in h_pos_i
						|| (i_r_c < hotspots_r_c.size() && hotspots[i].quality < hotspots_r_c[i_r_c].quality); // or next rev komp hint has more window hits 
			{
				tries++;
				const auto [eds_pos, read_pos] = align_r_c // position in eds and read, respectively
						? std::make_pair(hotspots_r_c[i_r_c].eds_pos, hotspots_r_c[i_r_c].read_pos)
						: std::make_pair(hotspots    [i    ].eds_pos, hotspots    [i    ].read_pos);
				const auto dist = align_r_c
					? align<false, dist_type, uint32_t>(EDS, adj, eds_pos, read_r_c.sequence, read_r_c.qual, read_pos, D)
					: align<false, dist_type, uint32_t>(EDS, adj, eds_pos,     read.sequence,     read.qual, read_pos, D);
	
				if (dist.first + dist.second <= D) {
					// ali found
					aligns++;
					tmp_alignments.emplace_back(dist, eds_pos, read_pos, align_r_c, align_r_c ? i_r_c : i);
	
					D = dist.first + dist.second;
				}
			}
			
			if(align_r_c)
				i_r_c++;
			else
				i++;
		}
		//rsort(alignments.data(), alignments.data() + alignments.size(),  [] (const auto& v) { return get<1>(v); });
		//SORT BY D

		// TODO: we actually only need the best report_count
	 	sort(tmp_alignments.begin(), tmp_alignments.end(), [](const auto & t1, const auto & t2) { return t1.get_dist() < t2.get_dist(); } );
		
		return tmp_alignments;
	}
	// finalizes the first max_out alignments
	template<typename int_type>
	std::vector<alignment<int_type>>
	finalize_alignments(
		const fasta_read<int_type>& read,
		const fasta_read<int_type>& read_r_c,
		const std::string& EDS,
		const adjacency & adj,
		const std::vector<temp_alignment<int_type>>& tmp,
		const size_t max_out
	) {
		const size_t num_out = std::min(max_out, tmp.size());
		std::vector<alignment<int_type>> res(num_out);
		for (size_t i = 0; i < num_out; i++) {
			auto[cigar, pos] = tmp[i].reverse_compl
					? align<true, dist_type, uint32_t>(EDS, adj, tmp[i].eds_pos, read_r_c.sequence, read_r_c.qual, tmp[i].read_pos, tmp[i].dist)
					: align<true, dist_type, uint32_t>(EDS, adj, tmp[i].eds_pos,     read.sequence,     read.qual, tmp[i].read_pos, tmp[i].dist);
			res[i] = alignment<int_type>(std::move(cigar), tmp[i].get_dist(), pos, tmp[i].reverse_compl);
		}
		return res;
	}

	template<typename int_type>
	std::vector<alignment<int_type>>
	start_aligner(
		const std::vector<hotspot<int_type>>& hotspots,
		const fasta_read<int_type>& read,
		const std::vector<hotspot<int_type>>& hotspots_r_c,
		const fasta_read<int_type>& read_r_c,
		const std::string & EDS,
		const adjacency & adj,
		uint32_t D,
		const uint32_t max_a_c, // max alignments calculated
		const uint32_t max_a_t, // max tries
		const uint32_t max_a_o
	) {
		return finalize_alignments(
			read, read_r_c,
			EDS,
			adj,
			align(hotspots, read, hotspots_r_c, read_r_c, EDS, adj, D, max_a_c, max_a_t),
			max_a_o);
	}
	template<typename int_type>
	std::vector<alignment<int_type>>
	start_aligner(
		const std::vector<hotspot<int_type>>& hotspots,
		const fasta_read<int_type>& read,
		const std::string & EDS,
		const adjacency & adj,
		uint32_t D,
		const uint32_t max_a_c, // max alignments calculated
		const uint32_t max_a_t, // max tries
		const uint32_t max_a_o
	) {
		fasta_read<int_type> dummy;
		std::vector<hotspot<int_type>> hotspots_dummy;
		return finalize_alignments(
			read, dummy,
			EDS,
			adj,
			align(hotspots, read, hotspots_dummy, dummy, EDS, adj, D, max_a_c, max_a_t),
			max_a_o);
	}


	// just chains get_fragments, get_positions and find_hotspots
	template < typename int_type >
	vector<hotspot<int_type>>
	find_hotspots(
		const fasta_read<int_type>& read,
		const gedmap_mini::minimizer_index & mini,
		uint32_t fragment_max_count, // get_fragments
		uint32_t window_size, // find_hotspots
		uint32_t window_hits, // find_hotspots
		bool check_col         // find_hotspots
	) {
		auto fragments = get_fragments(read, fragment_max_count, mini);
		auto pos_pairs = get_positions(std::move(fragments), mini);
		return find_hotspots(std::move(pos_pairs), window_size, window_hits, check_col);
	}
	
	template< typename int_type >
	void
	append_alignments(
		std::vector< alignment<int_type> >& alignments,
		std::vector< alignment<int_type> >&& tmp
	) {
		const auto b = tmp.size();
		tmp.insert(
			tmp.end(),
			std::make_move_iterator(alignments.begin()),
			std::make_move_iterator(alignments.end()));
		std::inplace_merge(
			tmp.begin(),
			tmp.begin() + b,
			tmp.end(),
			[](const auto& lhs, const auto& rhs) {
				return lhs.dist < rhs.dist;
			});
		alignments = std::move(tmp);
	}

	using SAM_FLAG = uint32_t;
	namespace SAM_FLAGS {
		constexpr SAM_FLAG PROPER_PAIR = 2;
		constexpr SAM_FLAG UNMAPPED = 4;
		constexpr SAM_FLAG REV_COMP = 16;
		constexpr SAM_FLAG NOT_PRIMARY = 256;
	} // namespace SAM_FLAGS
	
	template< typename int_type >
	uint32_t
	write_alignment(
		const std::vector< alignment<int_type> >& alignments,
		const fasta_read<int_type>& read,
		const fasta_read<int_type>& read_r_c,
		ostream & ofs,
		bool write_failure,
		const pos_EDS_to_FA_type & transform
	) {
		string 	QNAME	= read.id.substr(1);
		string	RNAME	= "*";
		uint32_t 	POS		= 0;
		uint32_t	MAPQ 	= 255; //TODO
		string 	CIGAR	= "*";
		string 	RNEXT 	= "*";
		uint32_t 	PNEXT 	= 0;
		uint32_t 	TLEN		= 0;
		const string& SEQ		= read.sequence;
		const string& R_SEQ = read_r_c.sequence;
		const string& QUAL		= read.qual;
		const string& R_QUAL	= read_r_c.qual;
		
		if(alignments.empty()) {
			if (write_failure)
				ofs 	<< QNAME << '\t' << SAM_FLAGS::UNMAPPED << '\t' << RNAME << '\t' << POS << '\t' << MAPQ << '\t' << CIGAR << '\t' << RNEXT << '\t' << PNEXT << '\t' << TLEN << '\t' << SEQ << '\t' << QUAL << '\n';
			return 0;
		}
		
		for (size_t i = 0; i < alignments.size(); i++){
			RNAME 	= "*";
			POS 		= alignments[i].pos+1;
			const auto& CIGAR 	= alignments[i].cigar;
			
			SAM_FLAG FLAG = 0;
			if ( i > 0 )	
				FLAG |= SAM_FLAGS::NOT_PRIMARY;
			
			if(alignments[i].align_r_c)
				FLAG |= SAM_FLAGS::REV_COMP;
			uint32_t off = 0;
			if(!transform.empty()) tie(RNAME,POS,off) = transform(POS);
			
				
			ofs 	<< QNAME << '\t' << FLAG << '\t' << RNAME << '\t' << POS << '\t' << MAPQ << '\t' << CIGAR << '\t' << RNEXT << '\t' << PNEXT << '\t' << TLEN << '\t' 
				<< (alignments[i].align_r_c ? R_SEQ : SEQ) << '\t'
				<< (alignments[i].align_r_c ? R_QUAL : QUAL) << '\t'
				<< "NM:i:" << ((uint32_t)alignments[i].dist) <<  (off?(" XO:i:" + to_string(off)):"") << '\n';
		}
		
		return alignments.size();
	}

	/**
	* @brief writes alignment
	* @param ofs	ostream
	* @param write_failure	if true a line with name and stars in written if read could not be aligned
	* @param transform		transforms position in EDS to position in FA
	* @return number of alignments in output
	*/
	template< typename int_type >
	uint32_t
	write_alignment(
		const std::vector< alignment<int_type> >& alignments,
		const fasta_read<int_type>& read,
		ostream & ofs,
		bool write_failure,
		const pos_EDS_to_FA_type & transform
	) {
		return write_alignment(
			alignments,
			read,
			read.get_rev_compl(),
			ofs,
			write_failure,
			transform);
	}

	template<typename int_type>
	uint32_t
	write_aligned_mates(
		const std::vector< alignment<int_type> >& alignments_l,
		const std::vector< alignment<int_type> >& alignments_r,
		const fasta_read<int_type>& read_1,
		const fasta_read<int_type>& read_2,
		ostream & ofs,
		const pos_EDS_to_FA_type & transform
	) {
		assert(alignments_l.size() == alignments_r.size());
		// TODO: write failed mates?
		for (size_t i = 0; i < alignments_l.size(); i++) {
			assert(alignments_l[i].align_r_c == alignments_r[i].align_r_c);
			const bool align_r_c = alignments_l[i].align_r_c;

			std::array< const alignment<int_type>*, 2 > alignments{{ &alignments_l[i], &alignments_r[i] }};
			std::array< const fasta_read<int_type>*, 2 > reads{{ &read_1, &read_2 }};
			if (align_r_c) {
				std::swap(alignments[0], alignments[1]);
				std::swap(reads[0], reads[1]);
				for (size_t i = 0; i < 2; i++)
					reads[i] = &reads[i]->get_rev_compl();
			}
			std::string_view QNAME = std::string_view(reads[0]->id);
			QNAME.remove_prefix(1); // remove "@"
			QNAME.remove_suffix(2); // remove "/1" (or "/2")

			// write l
			SAM_FLAG flag = SAM_FLAGS::PROPER_PAIR; // we only produce proper pairs
			if (align_r_c) flag |= SAM_FLAGS::REV_COMP;
			if (i > 1) flag |= SAM_FLAGS::NOT_PRIMARY;


			std::array< uint32_t, 2 > POS, off{{0, 0}};
			std::array< std::string, 2 > RNAME;
			for (size_t i = 0; i < 2; i++)
			{
				POS[i] = alignments[i]->pos + 1;
				if (!transform.empty())
					std::tie(RNAME[i], POS[i], off[i]) = transform(POS[i]);
				else
					RNAME[i] = "*";
			}
			
			for (size_t i = 0; i < 2; i++)
			{
				const uint32_t MAPQ = 255; //TODO
				constexpr std::string_view RNEXT = "*";
				constexpr uint32_t TLEN = 0;

				ofs << QNAME << '\t' << flag << '\t' << RNAME[i] << '\t' << POS[i] << '\t' << MAPQ << '\t' << alignments[i]->cigar << '\t' << RNEXT << '\t' << POS[1-i] /* PNEXT */ << '\t' << TLEN << '\t' 
					<< reads[i]->sequence << '\t'
					<< reads[i]->qual << '\t'
					<< "NM:i:" << (alignments[i]->dist) <<  (off[i]?(" XO:i:" + to_string(off[i])):"") << '\n';
			}
		}
		return alignments_l.size();
	}
} // namespace read_processor

template <typename int_type>
struct fasta_read{
	// cigar, min_dist, pos, align_r_c

	string	id;			//name of seq
	string	sequence;		// TODO: this should probably not be a string, at least not for seeding
	string	qual;

	mutable fasta_read<int_type>* m_rev_compl = nullptr;

	fasta_read() { }

	const fasta_read<int_type>& get_rev_compl() const {
		if (!m_rev_compl) {
			m_rev_compl = new fasta_read<int_type>(
				std::string(id),// + "_rev",
				gedmap_encode::rev_complement(sequence),
				gedmap_encode::rev_complement<false>(qual));
		}
		return *m_rev_compl;
	}
	
	fasta_read(string&& id, string&& sequence, string&& qual)
		: id(std::move(id))
		, sequence(std::move(sequence))
		, qual(std::move(qual))
	{ }

	fasta_read(istream & fastq_s){
		getline(fastq_s, id);
		if (!fastq_s.good()) [[unlikely]] throw runtime_error ("NO READ IN STREAM");
		getline(fastq_s, sequence);
		getline(fastq_s, qual);  //plus
		getline(fastq_s, qual); //qual
		if (!fastq_s.good()) [[unlikely]] {
			gedmap_io::print_error("READ " + id + " incomplete");
			throw runtime_error ("NO READ IN STREAM");
		}
	}
	~fasta_read() {
		delete m_rev_compl;
	}
};
