#pragma once

#include "fast_frontier.hpp"
#include "ed_graph.hpp"
//#include "timer.hpp"

#include <iostream> // TODO: rm
#include <algorithm>
#include <cassert>

// only search in forward direction
// query on nullopt
// f is called with (query payload, payload of min, value)
template<typename T, typename R, typename Payload, typename F, typename EDG, typename SKQ>
void paired_end(F f, std::vector<std::tuple<T, Payload, std::optional<R>>> ts, EDG g, const SKQ& skq, const T& dist) {
	const R max_err = dist * 5 + 1;

	//measure("sorting", [&] {
	std::sort(ts.begin(), ts.end(), [] (const auto& lhs, const auto& rhs) {
		return std::get<0>(lhs) < std::get<0>(rhs) or (std::get<0>(lhs) == std::get<0>(rhs) and std::get<2>(lhs) and !std::get<2>(rhs));
	});
	//});

	using frontier_type = fast_frontier<T, double, Payload>;

	struct before_open_s {
		T open, close, min_dist;
		frontier_type frontier;
	};

	frontier_type fr(dist);
	std::optional<before_open_s> before_open = std::nullopt;
	frontier_type accumulate(dist);

	const auto process_event = [&] (const T& t, const bool add_ev) {
		assert(g);
		fr.commit(g.next_position());
		
		if (!fr.down.empty() and g.next_event_type() == VAR_OPEN) {
			if (const auto[min_dist, close, all_eq, num_skip] = g.get_open_payload(); all_eq and close <= t) {
				// skip this variant site
				const auto open = g.next_position();
				
				fr.shift( close + 2 - open - min_dist );
				g.skip_events(num_skip + 1); // +1 since we only had .peek_next() above
				return;
			}
		}
		
		// TODO: used very rarely
		if (before_open
				and before_open->close < t // next query/add is after variant
				and before_open->frontier.max_t() < before_open->open + before_open->min_dist) { // frontier from before variant only contributes up

			assert(g.next_event_type() != VAR_OPEN);

			fr.shift( before_open->close + 1 - g.next_position() );
			accumulate.merge( std::move(fr) );

			//g.unsafe_seek(before_open->close + 1);
			while (g.next_event_type() != VAR_CLOSE) g.skip_events(1);
			g.skip_events(1); // skip close
			assert(g.next_event_type() == VAR_OPEN);
			
			before_open->frontier.commit( before_open->open + before_open->min_dist );
			before_open->frontier.shift( before_open->close + 2 - before_open->open - before_open->min_dist );
			accumulate.merge( std::move(before_open->frontier) );

			fr = std::move(accumulate);
			before_open = std::nullopt;
			
			return;
		}
		if (!before_open and add_ev and fr.down.empty()) {
			//	add: clear everything and restart
			if (auto q = fr.template query<double>(t); !q or *q > max_err) {
				if (g.seek(t))
					return;
			}
		}

		if (g.next_event_type() == VAR_OPEN and fr.down.empty()) {
			// TODO: maybe use the slow way for short distances
			const auto p = g.next_position();
			if (g.seek(t)) {
				fr.shift( skq.get_skip(p, g.next_position() - 1) );
				return;
			}
		}

		{
			const auto t = g.next_position();
			switch (g.next_event_type()) {
				case VAR_OPEN:
					{
						const auto[min_dist, close, all_eq, _] = g.get_open_payload();
						g.skip_events(1);
						before_open = before_open_s { .open = t, .close = close, .min_dist = min_dist, .frontier = fr };
						fr.shift(1); // ignore '('
						accumulate = frontier_type(dist);
					}
					break;
				case VAR_CLOSE:
					g.skip_events(1);
					fr.shift(1); // for ')'
					assert(before_open);

					// no extra shift needed, since terminal symb. of variant (')' or '|') is already accounted for in accumulate
					accumulate.merge( std::move(fr) );
					fr = std::move(accumulate);
					before_open = std::nullopt;

					break;
				case VAR_SEP:
					g.skip_events(1);
					assert(before_open);

					fr.shift( before_open->close + 1 - t );
					accumulate.merge( std::move(fr) );

					fr = before_open->frontier;
					fr.shift( t - before_open->open + 1 );

					break;
			}
			fr.commit(t);
		}
	};

	// TODO: don't add unless there is a chance it will be used
	// -> look forward to next query, if it has distance >= dist we
	// can simply note it and consider it upon query
	size_t ts_next_query_i = 0;
	T ts_next_query_open = 0;
	std::pair<R, Payload> ts_next_query_add_up;

	const auto advance_ts_next_query_i = [&] {
		while (ts_next_query_i < ts.size() && std::get<2>(ts[ts_next_query_i]))
			ts_next_query_i++;
		if (ts_next_query_i < ts.size())
		{
			auto tmp_g = g;
			tmp_g.seek(std::get<0>(ts[ts_next_query_i]));
			ts_next_query_open = std::min( (T)tmp_g.next_position() - 1, std::get<0>(ts[ts_next_query_i]) );
		}
		ts_next_query_add_up.first = std::numeric_limits<R>::max();
	};
	advance_ts_next_query_i();
	for (size_t ts_i = 0; ts_i < ts.size(); ts_i++) {
		auto[t, payload, cost] = ts[ts_i];

		bool add_up = false;

		//if (ts_i == ts_next_query_i and ts_next_query_add_up.first < std::numeric_limits<R>::max())
		//{ // we must add something at ts_next_query_open
		//	// TODO

		//	t = ts_next_query_open;
		//	cost = { ts_next_query_add_up.first };
		//	payload = ts_next_query_add_up.second;

		//	add_up = true;

		//	ts_i--;
		//	ts_next_query_add_up.first = std::numeric_limits<R>::max();
		//}
		//else if (cost and ts_next_query_i <= ts_i) {
		//	ts_next_query_i = ts_i + 1;
		//	advance_ts_next_query_i();
		//	if (ts_next_query_i >= ts.size()) break; // no next query -> stop
		//}

		//while (cost and before_open and g and g.next_position() <= t)
		//	process_event(t, !!cost);
		//if (cost and not before_open and (not g or g.next_position() > t))
		//{
		//	const auto d = ts_next_query_open - t - skq.get_skip(t, ts_next_query_open);
		//	if (d >= dist) {
		//		if (const R val = (d - dist) + *cost; val < ts_next_query_add_up.first)
		//		{
		//			ts_next_query_add_up = std::make_pair(val, payload);
		//		}
		//		continue;
		//	}
		//}

		if (g and g.next_event_type() == VAR_OPEN) {
			// TODO: maybe use faster rank-select data structure?
			assert(!before_open);

			auto tar = fr.max_t();
			bool ok = false; // fr.down will be empty before t

			// find min dist until t
			T total_shift = 0;
			for (auto tmp_g = g; tmp_g and tmp_g.next_position() <= t; ) {
				assert(tmp_g.next_event_type() == VAR_OPEN);
				const auto open = tmp_g.next_position();
				if (tar + total_shift <= open) {
					ok = true;
					g = tmp_g;
					fr.shift(total_shift);
					break;
				}
				const auto[min_dist, close, all_eq, num_skip] = tmp_g.get_open_payload();
				if (close > t) break;
				total_shift += close + 2 - open - min_dist;

				tmp_g.skip_events(num_skip + 1); // +1 since we only had .peek_next() above
			}
			
			if (ok)
				if (auto q = fr.template query<double>(t); !cost and (not q or *q > max_err)) {
					continue;
				}
		}

		while (g and g.next_position() <= t)
			process_event(t, !!cost);

		if (!cost)
		{ // perform query
			//std::cerr << "QUERY " << t << " : " << fr.template query<R>(t) << std::endl;
			if (const auto res = fr.template query<>(t); res and res->second <= max_err) {
				f(payload, res->first, res->second);
			}
		}
		else
		{ // add to frontier
			fr.commit(t); // for performance
			if (add_up) {
				//std::cerr << "add_up @" << t << " : " << *cost << std::endl;
				fr.add_up(t, *cost, payload);
			} else
				fr.add(t, *cost, payload);
		}
	}
}
