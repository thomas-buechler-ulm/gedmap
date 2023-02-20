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
	const R max_err = dist * 3 + 1;

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

	// TODO: don't add unless there is a chance it will be used
	for (const auto&[t, payload, cost] : ts) {
		if (g and g.peek_next().second.index() == VAR_OPEN) {
			assert(!before_open);
			// find min dist until t
			auto tar = fr.max_t();
			T total_shift = 0;
			bool ok = false; // fr.down will be empty before t
			for (auto tmp_g = g; tmp_g and tmp_g.peek_next().first <= t; ) {
				const auto[open, ev] = tmp_g.peek_next();
				assert(ev.index() == VAR_OPEN);
				if (tar + total_shift <= open) {
					ok = true;
					g = tmp_g;
					fr.shift(total_shift);
					break;
				}
				const auto&[min_dist, close, all_eq, num_skip] = std::get<VAR_OPEN>(ev);
				if (close > t) break;
				total_shift += close + 2 - open - min_dist;

				tmp_g.skip(num_skip + 1); // +1 since we only had .peek_next() above
			}
			if (ok)
				if (auto q = fr.template query<double>(t); !cost and (!q or *q > max_err)) {
					continue;
				}
		}

		while (g and g.peek_next().first <= t) {
			fr.commit(g.peek_next().first);

			if (!fr.down.empty() and g.peek_next().second.index() == VAR_OPEN) {
				const auto[open, ev] = g.peek_next();
				if (const auto&[min_dist, close, all_eq, num_skip] = std::get<VAR_OPEN>(ev); all_eq and close <= t) {
					// skip this variant site
					
					fr.shift( close + 2 - open - min_dist );
					g.skip(num_skip + 1); // +1 since we only had .peek_next() above
					//g.pop_next();
					//while (g.pop_next().second.index() != VAR_CLOSE) { }
					continue;
				}
			}
			
			// TODO: used very rarely
			if (before_open
					and before_open->close < t // next query/add is after variant
					and before_open->frontier.max_t() < before_open->open + before_open->min_dist) { // frontier from before variant only contributes up

				assert(g.peek_next().second.index() != VAR_OPEN);

				fr.shift( before_open->close + 1 - g.peek_next().first );
				accumulate.merge( std::move(fr) );

				//g.unsafe_seek(before_open->close + 1);
				while (g.pop_next().second.index() != VAR_CLOSE) { }
				
				before_open->frontier.commit( before_open->open + before_open->min_dist );
				before_open->frontier.shift( before_open->close + 2 - before_open->open - before_open->min_dist );
				accumulate.merge( std::move(before_open->frontier) );

				fr = std::move(accumulate);
				before_open = std::nullopt;
				
				continue;
			}
			if (!before_open and cost and fr.down.empty()) {
				//	add: clear everything and restart
				if (auto q = fr.template query<double>(t); !q or *q > max_err) {
					if (g.seek(t))
						continue;
				}
			}

			if (g.peek_next().second.index() == VAR_OPEN and fr.down.empty()) {
				// TODO: maybe use the slow way for short distances
				const auto p = g.peek_next().first;
				if (g.seek(t)) {
					fr.shift( skq.get_skip(p, g.peek_next().first - 1) );
					continue;
				}
			}

			{
				const auto[t, ev] = g.pop_next();
				switch (ev.index()) {
					case VAR_OPEN:
						{
							const auto&[min_dist, close, all_eq, _] = std::get<VAR_OPEN>(ev);
							before_open = before_open_s { .open = t, .close = close, .min_dist = min_dist, .frontier = fr };
							fr.shift(1); // ignore '('
							accumulate = frontier_type(dist);
						}
						break;
					case VAR_CLOSE:
						fr.shift(1); // for ')'
						assert(before_open);

						// no extra shift needed, since terminal symb. of variant (')' or '|') is already accounted for in accumulate
						accumulate.merge( std::move(fr) );
						fr = std::move(accumulate);
						before_open = std::nullopt;

						break;
					case VAR_SEP:
						assert(before_open);

						fr.shift( before_open->close + 1 - t );
						accumulate.merge( std::move(fr) );

						fr = before_open->frontier;
						fr.shift( t - before_open->open + 1 );

						break;
				}
				fr.commit(t);
			}
		}

		if (!cost)
		{ // perform query
			if (const auto res = fr.template query<>(t); res and res->second <= max_err) {
				f(payload, res->first, res->second);
			}
		}
		else
		{ // add to frontier
			fr.commit(t); // for performance
			fr.add(t, *cost, payload);
		}
	}
}
