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
	VAR_CLOSE,
	VAR_SEP
};

template<typename T>
struct EDG_view;

template<typename T>
class ED_Graph {
public:
	using skip_length_t = uint16_t;
	struct open_payload {
		T min_dist, close; // NOTE: min_dist is number of chars + 1
		bool all_eq;
		skip_length_t num_skip = 0;
	};
	using event_payload = std::variant
							< open_payload		// OPEN
							, std::monostate	// CLOSE
							, std::monostate	// SEP
							>;
	using event_type = std::pair
						< T // position
						, event_payload>;
private:
	std::vector<event_type> events;

	sdsl::bit_vector shortest_path; // 1 iff in shortest variant (or outside variant scope)
	sdsl::rank_support_v5<0> rs_ref_ind;

	sdsl::bit_vector event_pos;
	sdsl::rank_support_v5<1> rs_event_pos;

	void build(std::string_view edg) {
		events.clear();
		bool all_eq = true; // are all distances in variant equal?
		size_t open_i = 0;
		std::pair<T,T> min_dist; // (dist, start pos)
		shortest_path.resize(edg.size() + 1);
		event_pos.resize(edg.size() + 1);
		for (size_t i = 0; i < edg.size(); i++)
			if (edg[i] == '(') {
				for (size_t j = events.empty() ? 0 : events.back().first + 1; j < i; j++)
					shortest_path[j] = true;
				open_i = events.size();
				min_dist = std::make_pair(std::numeric_limits<T>::max(), (size_t)0);
				all_eq = true;
				events.emplace_back(
					static_cast<T>(i),
					event_payload(std::in_place_index<VAR_OPEN>, T{0}, T{0}, false));
			} else if (edg[i] == ')') {
				assert(!events.empty());
				const auto dist = static_cast<T>(i - events.back().first);
				if (min_dist.first != std::numeric_limits<T>::max() && min_dist.first != dist)
					all_eq = false;
				min_dist = std::min(min_dist, std::make_pair(dist, events.back().first + 1));

				// mark in shortest_path
				for (T i = min_dist.second; i + 1 < min_dist.second + min_dist.first; i++)
					shortest_path[i] = true;

				const size_t skip = events.size() - open_i;
				assert(skip == static_cast<skip_length_t>(skip));
				assert(min_dist.first > 0);
				std::get<VAR_OPEN>(events[open_i].second) = open_payload{ min_dist.first, static_cast<T>(i), all_eq, static_cast<skip_length_t>(skip) };

				events.emplace_back(
					static_cast<T>(i),
					event_payload(std::in_place_index<VAR_CLOSE>));
			} else if (edg[i] == '|') {
				assert(!events.empty());
				const auto dist = static_cast<T>(i - events.back().first);
				if (min_dist.first != std::numeric_limits<T>::max() && min_dist.first != dist)
					all_eq = false;
				min_dist = std::min(min_dist, std::make_pair(dist, events.back().first + 1));

				events.emplace_back(
					static_cast<T>(i),
					event_payload(std::in_place_index<VAR_SEP>));
			}
		for (size_t j = events.empty() ? 0 : events.back().first + 1; j < edg.size(); j++)
			shortest_path[j] = true;
		assert(events.empty() or events.back().second.index() == VAR_CLOSE);

		events.emplace_back(
			static_cast<T>(edg.size()),
			event_payload(std::in_place_index<VAR_OPEN>, T{0}, T{0}, 0));
		events.emplace_back(events.back());

		for (const auto& e : events)
			event_pos[ e.first ] = true;
		rs_event_pos = sdsl::rank_support_v5<1>(&event_pos);
	}
public:
	ED_Graph<T>& operator=(ED_Graph<T>&& rhs) = delete; // because sdsl::rank_support doesnt have move-assignment
	ED_Graph(const ED_Graph<T>& rhs) = delete;
	ED_Graph(ED_Graph<T>&& rhs) = default;
	ED_Graph() = default;
	explicit ED_Graph(std::string_view edg) {
		//std::cout << "edg.size() = " << edg.size() << std::endl;
		build(edg);
		rs_ref_ind = sdsl::rank_support_v5<0>(&shortest_path);
	}
	~ED_Graph() = default;
	friend class EDG_view<T>;
	// hi is expected to be at ')' (or before)
	T get_skip(T lo, T hi) const { // (excl,inc]
		return rs_ref_ind(hi+1) - rs_ref_ind(lo);
	}
};

template<typename T>
class EDG_view {
public:
	using event_type = ED_Graph<T>::event_type;
private:
	const event_type* m_events;
	const sdsl::rank_support_v5<1>* m_rs_event_pos;
	size_t m_num_events;
	T m_cur_start;
	size_t i;
	size_t pos;
public:
	EDG_view(const ED_Graph<T>* g)
		: m_events(g->events.data())
		, m_rs_event_pos(&(g->rs_event_pos))
		, m_num_events(g->events.size() - 1)
		, m_cur_start(0)
		, i(0)
		, pos(0)
	{}

	operator bool() const { return i + 1 < m_num_events; }

	// skip until after next VAR_CLOSE
	void seek_end() {
		while (i < m_num_events && m_events[i].second.index() != VAR_CLOSE)
			i++;
	}

	// may land in an alternative scope, use with care. If you do not want to land in
	// an alternative scope, use seek
	void unsafe_seek(const T& pos) {
		i = (*m_rs_event_pos)(pos);
	}
	inline void skip(ptrdiff_t n) {
		i += n;
		assert(i < m_num_events);
	}
	bool seek(const T& pos) {
		const size_t l_i = i;
		unsafe_seek(pos);
		// we must not "land" in alternative scope
		if (i != m_num_events)
			// we have very few alternatives in a scope, so a binary search
			// is probably not worth it
			while (m_events[i].second.index() != VAR_OPEN)
				i--;
		return i != l_i;
	}

	size_t remaining_events() const noexcept {
		return m_num_events - i;
	}
	std::string_view view_until_next(std::string_view edg) const {
		assert(i <= m_num_events);
		return edg.substr(m_cur_start, m_events[i].first - m_cur_start);
	}
	const event_type& peek_next() const {
		assert(i < m_num_events);
		return m_events[i];
	}
	const event_type& pop_next() {
		const auto& res = peek_next();
		i++;
		m_cur_start = res.first + 1;
		return res;
	}
	void undo_pop() {
		assert(i > 0);
		i--;
	}
};
