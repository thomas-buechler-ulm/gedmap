#pragma once

#include "paired_end/edsm_frontier.hpp"

namespace gedmap_align_min {

using pe_event_type = std::tuple<int64_t, uint32_t, std::optional<int64_t>>;

template< typename int_type, typename F, typename GET_TMP_ALIGNMENTS >
inline
std::vector< temp_alignment<int_type> >
try_mate(
	const ED_Graph<uint32_t>& msq,
	std::vector< pe_event_type >& events,
	std::vector< pe_event_type >& events_rev, // TODO: maybe dont copy these
	std::vector<uint32_t> (&right_source)[2],
	const std::vector< temp_alignment<int_type> >& tmp_alignments_l,
	size_t ali_begin, size_t ali_end,
	F mk_hotspot,
	GET_TMP_ALIGNMENTS get_tmp_alignments,
	double& opt
) {
	const size_t events_pre_size = events.size(), events_rev_pre_size = events_rev.size();
	for ( ; ali_begin != ali_end; ali_begin++) {
		const auto& ali = tmp_alignments_l[ali_begin];
		//tmp_alignments_l[i].source = std::numeric_limits<uint32_t>::max();
		// TODO: find the rightmost matched symbol
		(ali.reverse_compl ? events_rev : events)
			.emplace_back(
				ali.eds_pos - ali.read_pos,
				ali_begin,
				-(int64_t)ali.get_dist());
	}

	// rank candidates according to distance to already aligned "tmp_alignments_l"
	std::vector<std::tuple<uint32_t,bool,double>> opt_mate // (mate index in hotspots_* / pos_pairs_*, r_c, dist)
		( tmp_alignments_l.size()
		, std::make_tuple(std::numeric_limits<uint32_t>::max(), false, std::numeric_limits<double>::infinity()));
	paired_end<int64_t, int64_t>([&opt_mate](auto a, auto b, auto val) {
			if (val < std::get<2>(opt_mate[b]))
				opt_mate[b] = std::make_tuple(a, false, val);
		},
		events,
		EDG_view<uint32_t>{&msq},
		msq,
		PE_FRAGMENT_LENGTH);
	paired_end<int64_t, int64_t>([&opt_mate](auto a, auto b, auto val) {
			if (val < std::get<2>(opt_mate[b]))
				opt_mate[b] = std::make_tuple(a, true, val);
		},
		events_rev,
		EDG_view<uint32_t>{&msq},
		msq,
		PE_FRAGMENT_LENGTH);

	events.resize(events_pre_size), events_rev.resize(events_rev_pre_size);
	
	opt = std::numeric_limits<double>::infinity();

	for (size_t i = 0; i < 2; i++) {
		right_source[i].clear();
		right_source[i].reserve(opt_mate.size());
	}
	return [&] {
		std::vector<hotspot<int_type>> final_hotspots, final_hotspots_r_c;
		for (size_t i_l = 0; i_l < opt_mate.size(); i_l++) {
			const auto&[i_r, r_c, val] = opt_mate[i_l];
			if (val == std::numeric_limits<double>::infinity()) 
				continue;
			opt = std::min(opt, val);
			right_source[r_c].emplace_back(i_l);

			(r_c ? final_hotspots_r_c : final_hotspots)
				.emplace_back(mk_hotspot(r_c, i_r, val));
		}
		std::sort(final_hotspots.begin(), final_hotspots.end(),
				[&](const auto& lhs, const auto& rhs) {
					return lhs.quality > rhs.quality;
				});
		std::sort(final_hotspots_r_c.begin(), final_hotspots_r_c.end(),
				[&](const auto& lhs, const auto& rhs) {
					return lhs.quality > rhs.quality;
				});
		return get_tmp_alignments(std::move(final_hotspots), std::move(final_hotspots_r_c));
	}();
}

// NOTE: we assume FR-reads, i.e. one read aligns to the forward strand and the
// other aligns to the reverse strand at a later position so that they point
// towards each other. We of course also consider the case that the fragment is
// the reverse-complement
template< class int_type >
void
map_pairs(
	gedmap_mini::minimizer_index & mini,
	const ED_Graph<uint32_t>& msq, // minimum shift query
	const std::string & EDS,
	const adjacency & adj,
	const pos_EDS_to_FA_type & p2FA,
	std::istream& fastq_l,
	std::istream& fastq_r,
	std::ofstream& o_s)
{
	// set large buffers for fastq_{l,r}. Mostly useful for many threads
	std::vector<char> outbuffer(1<<23);
	o_s.rdbuf()->pubsetbuf(outbuffer.data(), outbuffer.size());

	gedmap_io::print_row("Fragment length ", PE_FRAGMENT_LENGTH);
	if(THREAD_COUNT) omp_set_num_threads(THREAD_COUNT);

	const auto num_reads = [&] {
		std::string tmp;
		size_t num_lines = 0;
		while (std::getline(fastq_l, tmp))
			num_lines++;
		fastq_l.clear();
		fastq_l.seekg(0, std::ios::beg);
		return num_lines / 4;
	} ();

	std::cerr << "max_aligns_t: " << MAX_ALIGNS_T[0]
		<< " max_aligns_c: " << MAX_ALIGNS_C[0]
		<< " max_aligns_o: " << MAX_ALIGNS_O
		<< " max_aligns_m: " << MAX_ALIGNS_M
		<< std::endl;

	size_t num_failed = 0;
	double avg = 0; size_t avg_cnt = 0;
	size_t num_aligns = 0, num_mate_candidates = 0;
	size_t max_aligns = 0;
	size_t sum_sol = 0;
	size_t num_tried_twice = 0;

	const auto get_fragments = [&](const fasta_read<int_type>& read) {
		return read_processor::get_fragments<int_type>(read, FRAGMENT_COUNT[0], mini);
	};
	const auto get_positions = [&](const vector<kmer_pair<int_type>>& kmer_pairs) {
		return read_processor::get_positions<int_type>(kmer_pairs, mini);
	};
	const auto get_tmp_alignments = [&]
		( const fasta_read<int_type>& read_l, vector<hotspot<int_type>>&& hotspots_l
		, const fasta_read<int_type>& read_c, vector<hotspot<int_type>>&& hotspots_c
		)
	{
		return read_processor::align
			( std::move(hotspots_l), read_l
			, std::move(hotspots_c), read_c
			, EDS, adj, MAX_DIST[0]
			, MAX_ALIGNS_C[0], MAX_ALIGNS_T[0]
			);
	};

	size_t count = 0;
	#pragma omp parallel for schedule(dynamic,10)
	for (size_t read_i = 0; read_i < num_reads; read_i++) {
		fasta_read<int_type> read_l, read_r;
		#pragma omp critical
		{
			read_l = fasta_read<int_type>(fastq_l);
			read_r = fasta_read<int_type>(fastq_r);
		}
		const auto& read_l_r_c = read_l.get_rev_compl();
		const auto& read_r_r_c = read_r.get_rev_compl();


		// TODO: use different FRAGMENT_COUNT, SPOT_HIT, ...
		
		auto hotspots_l = read_processor::find_hotspots(read_l, mini,
				FRAGMENT_COUNT[0],
				SPOT_SIZE, SPOT_HITS[0], CHECK_COLLI);
		auto hotspots_r = read_processor::find_hotspots(read_r, mini,
				FRAGMENT_COUNT[0],
				SPOT_SIZE, SPOT_HITS[0], CHECK_COLLI);
		auto pos_pairs_r_r_c = get_positions(get_fragments(read_r_r_c));
		auto pos_pairs_l_r_c = get_positions(get_fragments(read_l_r_c));
		std::vector< hotspot<int_type> > hotspots_r_r_c, hotspots_l_r_c;
		const bool pair_hotspots = false;

		// TODO: hotspots_l.size() < hotspots_r.size()) {
		// align left part of fragment to reference ("left" according to the reference)
		// NOTE: align_r_c = true means that the entire fragment is the reverse-complement
		auto tmp_alignments_l = get_tmp_alignments
			( read_l, std::move(hotspots_l)
			, read_r, std::move(hotspots_r));

		// build events for mate-pairing
		std::vector< std::tuple<int64_t, uint32_t, std::optional<int64_t>> > events, events_rev;
		if (pair_hotspots)
		{
			hotspots_r_r_c = read_processor::find_hotspots(std::move(pos_pairs_r_r_c), SPOT_SIZE, SPOT_HITS[0], CHECK_COLLI);
			for (uint32_t i = 0; i < hotspots_r_r_c.size(); i++)
				events.emplace_back(hotspots_r_r_c[i].eds_pos - hotspots_r_r_c[i].read_pos, i, std::nullopt);

			hotspots_l_r_c = read_processor::find_hotspots(std::move(pos_pairs_l_r_c), SPOT_SIZE, SPOT_HITS[0], CHECK_COLLI);
			for (uint32_t i = 0; i < hotspots_l_r_c.size(); i++)
				events_rev.emplace_back(hotspots_l_r_c[i].eds_pos - hotspots_l_r_c[i].read_pos, i, std::nullopt);
		}
		else
		{
			for (uint32_t i = 0; i < pos_pairs_r_r_c.size(); i++)
				events.emplace_back(pos_pairs_r_r_c[i].eds_pos - pos_pairs_r_r_c[i].read_pos, i, std::nullopt);
			for (uint32_t i = 0; i < pos_pairs_l_r_c.size(); i++)
				events_rev.emplace_back(pos_pairs_l_r_c[i].eds_pos - pos_pairs_l_r_c[i].read_pos, i, std::nullopt);
		}
		
		double opt = std::numeric_limits<double>::infinity();

		bool tried_twice = false;

		std::vector<uint32_t> right_source[2]; // not-rev, rev
		std::vector< temp_alignment<int_type> > tmp_alignments_r;
		for (size_t i = 0, num_try = 0; i < tmp_alignments_l.size() && num_try < MAX_ALIGNS_M.size(); num_try++) {
			tried_twice = num_try;
			size_t j = std::min(tmp_alignments_l.size(), i + MAX_ALIGNS_M[num_try]);
			tmp_alignments_r = try_mate<int_type>(
				msq,
				events, events_rev,
				right_source,
				tmp_alignments_l, i, j,
				[&] (bool r_c, uint32_t i_r, double val) {
					if (pair_hotspots)
						return (r_c ? hotspots_l_r_c : hotspots_r_r_c)[i_r];
					else {
						const auto& [eds_pos, read_pos] = (r_c ? pos_pairs_l_r_c : pos_pairs_r_r_c)[i_r];
						return hotspot<int_type>(eds_pos, read_pos, (int_type) std::max(0.0, PE_FRAGMENT_LENGTH - (double)val));
					}
				},
				[&] (std::vector<hotspot<int_type>>&& hotspots_l, std::vector<hotspot<int_type>>&& hotspots_c) {
					return get_tmp_alignments(read_r_r_c, std::move(hotspots_l), read_l_r_c, std::move(hotspots_c));
				},
				opt);
			if (!tmp_alignments_r.empty()) [[likely]] break;
			i = j;
		}

		// mate
		auto[alignments_l, alignments_r] = [&] {
			std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> mates; // score, tmp_l, tmp_r
			for (size_t i = 0; i < tmp_alignments_r.size(); i++) {
				const auto src = right_source[ tmp_alignments_r[i].reverse_compl ][ tmp_alignments_r[i].source ];
				assert(tmp_alignments_l[src].reverse_compl == tmp_alignments_r[i].reverse_compl);

				mates.emplace_back( tmp_alignments_l[src].get_dist() + tmp_alignments_r[i].get_dist(),
						src, i);
			}
			std::sort(mates.begin(), mates.end(), [](const auto& lhs, const auto& rhs) { return std::get<0>(lhs) < std::get<0>(rhs); });

			if (mates.size() > MAX_ALIGNS_O)
				mates.resize(MAX_ALIGNS_O);
			std::vector<temp_alignment<int_type>> tmp_l(mates.size()), tmp_r(mates.size());
			for (size_t i = 0; i < mates.size(); i++) {
				tmp_l[i] = tmp_alignments_l[std::get<1>(mates[i])];
				tmp_r[i] = tmp_alignments_r[std::get<2>(mates[i])];
			}
			return std::make_pair(
				read_processor::finalize_alignments(read_l, read_r, EDS, adj, tmp_l, MAX_ALIGNS_O),
				read_processor::finalize_alignments(read_r_r_c, read_l_r_c, EDS, adj, tmp_r, MAX_ALIGNS_O));
		}();

		#pragma omp critical
		{
			const auto cnt = read_processor::write_aligned_mates(
				alignments_l, alignments_r,
				read_l, read_r,
				o_s,
				p2FA);

			num_tried_twice += tried_twice;

			max_aligns = std::max(max_aligns, alignments_l.size());
			num_aligns += alignments_l.size();
			if (pair_hotspots)
				num_mate_candidates += hotspots_l_r_c.size() + hotspots_r_r_c.size();
			else
				num_mate_candidates += pos_pairs_l_r_c.size() + pos_pairs_r_r_c.size();
			if (cnt > 0) {
				avg += opt, avg_cnt++;
				sum_sol += cnt;
			} else
				num_failed++;
			count++;
			if (count % 1000 == 0)
				gedmap_io::flush_row("Searched read pairs", to_string(count));
		}
		
	}
	std::cerr << std::endl << "avg dist: " << avg / avg_cnt << " (cnt: " << avg_cnt << ")" << std::endl;
	std::cerr << "avg sols (for pairs with solutions): " << sum_sol / (double)avg_cnt << std::endl;
	std::cerr << num_failed << " could not pair" << std::endl;
	std::cerr << num_aligns/(double)count << " " << num_mate_candidates/(double)count << std::endl;
	std::cerr << "max_aligns: " << max_aligns << std::endl;
	std::cerr << "num_tried_twice: " << num_tried_twice << std::endl;
}

} // namespace gedmap_align_min
