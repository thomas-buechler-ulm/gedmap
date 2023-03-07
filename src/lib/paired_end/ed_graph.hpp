#pragma once

#include <string>
#include <vector>
#include <cassert>
#include <variant>
#include <tuple>
#include <algorithm>

#include <sdsl/bit_vectors.hpp>

enum EDG_EVENT {
	VAR_OPEN = 0,
	VAR_SEP = 1,
	VAR_CLOSE = 2
};

template<typename T>
class EDG_rev;

template<typename T, bool reverse>
struct EDG_view;

template<typename T>
class ED_Graph {
public:
	static_assert(not std::is_signed_v<T>);
	using skip_length_t = uint16_t;
	struct open_payload {
		T min_dist, close; // NOTE: min_dist is number of chars + 1
		bool all_eq;
		skip_length_t num_skip = 0;
	};
	using event_payload = std::variant
							< open_payload		// OPEN
							, std::monostate	// SEP
							, open_payload		// CLOSE, same as OPEN but with 'close' pointing to OPEN
							>;
	using event_type = std::pair
						< T // position
						, event_payload>;
private:
	// TODO: use bit rank/select vector 2FA instead of dist
	std::string edg;
	std::vector<event_type> events;
	sdsl::bit_vector ref_ind; // 1 iff in shortest variant (or outside variant scope)
	sdsl::rank_support_v5<0> rs_ref_ind;

	void build(std::string_view edg) {
		events.clear();
		bool all_eq = true; // are all distances in variant equal?
		size_t open_i = 0;
		std::pair<T,T> min_dist; // (dist, start pos)
		ref_ind.resize(edg.size());
		events.emplace_back(
			static_cast<T>(0),
			event_payload(std::in_place_index<VAR_CLOSE>, T{0}, T{0}, false));
		for (size_t i = 1; i <= edg.size(); i++)
			switch (edg[i-1]) {
				case '(':
				{
					for (size_t j = events.back().first + 1; j < i; j++)
						ref_ind[j-1] = true;
					open_i = events.size();
					min_dist = std::make_pair(std::numeric_limits<T>::max(), (size_t)0);
					all_eq = true;
					events.emplace_back(
						static_cast<T>(i),
						event_payload(std::in_place_index<VAR_OPEN>, T{0}, T{0}, false));
				}
				break;
				case ')':
				{
					assert(!events.empty());
					const auto dist = static_cast<T>(i - events.back().first);
					if (min_dist.first != std::numeric_limits<T>::max() && min_dist.first != dist)
						all_eq = false;
					min_dist = std::min(min_dist, std::make_pair(dist, events.back().first + 1));
	
					assert(min_dist.first != std::numeric_limits<T>::max());
	
					// mark in ref_ind
					for (T j = min_dist.second; j + 1 < min_dist.second + min_dist.first; j++)
					{
						ref_ind[j] = true;
					}
	
					const size_t skip = events.size() - open_i;
					assert(skip == static_cast<skip_length_t>(skip));
					assert(min_dist.first > 0);
					std::get<VAR_OPEN>(events[open_i].second) = open_payload
						{ min_dist.first
						, static_cast<T>(i)
						, all_eq
						, static_cast<skip_length_t>(skip)
						};
	
					events.emplace_back(
						static_cast<T>(i),
						event_payload(std::in_place_index<VAR_CLOSE>,
							min_dist.first,
							events[open_i].first,
							all_eq,
							static_cast<skip_length_t>(skip)));
				}
				break;
				case '|':
				{
					assert(!events.empty());
					const auto dist = static_cast<T>(i - events.back().first);
					if (min_dist.first != std::numeric_limits<T>::max() && min_dist.first != dist)
						all_eq = false;
					min_dist = std::min(min_dist, std::make_pair(dist, events.back().first + 1));

					events.emplace_back(
						static_cast<T>(i),
						event_payload(std::in_place_index<VAR_SEP>));
				}
				break;
			}
		for (size_t j = events.back().first + 1; j <= edg.size(); j++)
			ref_ind[j-1] = true;
		assert(events.empty() or events.back().second.index() == VAR_CLOSE);

		events.emplace_back(
			static_cast<T>(edg.size() + 1),
			event_payload(std::in_place_index<VAR_OPEN>, T{0}, T{0}, false));
	}
public:
	ED_Graph<T>& operator=(ED_Graph<T>&& rhs) = delete; // because sdsl::rank_support doesnt have move-assignment
	ED_Graph(const ED_Graph<T>& rhs) = delete;
	ED_Graph(ED_Graph<T>&& rhs) = default;
	ED_Graph() = default;
	explicit ED_Graph(std::string_view edg) {
		build(edg);
		rs_ref_ind = sdsl::rank_support_v5<0>(&ref_ind);
	}
	~ED_Graph() = default;
	friend class EDG_view<T, true>;
	friend class EDG_view<T, false>;
	friend class EDG_rev<T>;
	// hi is expected to be at ')' (or before)
	T get_skip(T lo, T hi) const { // [inc,inc]
		return rs_ref_ind(hi+1) - rs_ref_ind(lo);
	}
};

template<typename T>
class EDG_rev {
	const ED_Graph<T>* g;
	T max_pos;
public:
	EDG_rev(const ED_Graph<T>* g)
		: g(g)
		, max_pos(g->events.back().first - 2)
		{}
	inline T get_skip(T lo, T hi) const {
		return g->get_skip(max_pos - hi, max_pos - lo);
		std::tie(lo,hi) = std::make_pair(max_pos - hi, max_pos - lo);
	}
};

template<typename T, bool reverse = false>
class EDG_view {
public:
	using event_type = ED_Graph<T>::event_type;
	using open_payload = ED_Graph<T>::open_payload;
private:
	const event_type* m_events;
	size_t m_num_events;
	size_t i;
	T max_pos;
	inline T translate_pos(const T& t) const {
		if constexpr (reverse)
			return max_pos - t;
		else
			return t - 1;
	}
	inline T untranslate_pos(const T& t) const {
		if constexpr (reverse)
			return max_pos - t;
		else
			return t + 1;
	}
public:
	EDG_view(const ED_Graph<T>* g)
		: m_events(g->events.data())
		, m_num_events(g->events.size())
		, i(reverse ? m_num_events - 2 : 1)
		, max_pos(m_events[m_num_events-1].first - 1)
	{}

	size_t num_events() const { return m_num_events - 2; }

	operator bool() const
	{
		if constexpr (reverse)
			return i > 0;
		else
			return i + 1 < m_num_events;
	}

	// may land in an alternative scope, use with care. If you do not want to land in
	// an alternative scope, use seek
	void unsafe_seek(const T& pos) {
		if constexpr (reverse)
		{
			auto it = std::lower_bound(std::make_reverse_iterator(m_events + i), std::make_reverse_iterator(m_events), pos, [&](const auto& v, auto) {
				return v.first > untranslate_pos(pos);
			});
			i = std::distance(it, std::make_reverse_iterator(m_events));
		}
		else
		{
			auto it = std::lower_bound(m_events + i, m_events + m_num_events, pos, [&](const auto& v, auto) {
				return v.first < untranslate_pos(pos);
			});
			i = std::distance(m_events, it);
		}
	}
	inline void skip_events(ptrdiff_t n) {
		if constexpr (reverse) {
			assert((ptrdiff_t)i >= n);
			i -= n;
		} else {
			assert(i + n < m_num_events);
			i += n;
		}
	}
	bool seek(const T& pos) {
		const size_t l_i = i;
		unsafe_seek(pos);
		// we must not "land" in alternative scope
		if constexpr (reverse)
		{
			if (i > 0)
				while (m_events[i].second.index() != VAR_CLOSE)
					i++;
		}
		else
		{
			if (i != m_num_events)
				// we have very few alternatives in a scope, so a binary search
				// is probably not worth it
				while (m_events[i].second.index() != VAR_OPEN)
					i--;
		}
		return i != l_i;
	}

	T next_position() const {
		assert(i < m_num_events);
		return translate_pos(m_events[i].first);
	}
	EDG_EVENT next_event_type() const {
		assert(i < m_num_events);
		if constexpr (reverse)
			return static_cast<EDG_EVENT>(VAR_CLOSE - m_events[i].second.index());
		else
			return static_cast<EDG_EVENT>(m_events[i].second.index());
	}
	open_payload get_open_payload() const {
		assert(next_event_type() == VAR_OPEN);
		auto res = std::get<reverse ? VAR_CLOSE : VAR_OPEN>(m_events[i].second);
		res.close = translate_pos(res.close);
		return res;
	}
};
