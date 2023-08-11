#pragma once

#include <optional>
#include <utility>
#include <queue>
#include <tuple>
#include <cassert>
#include <algorithm>
#include <ostream>

template<typename T, typename F = std::less<T>>
std::optional<T>
max_opt(
	const std::optional<T>& lhs,
	const std::optional<T>& rhs,
	const F& cmp = F{}
) {
	if (!lhs) return rhs;
	if (!rhs or cmp(*rhs, *lhs))
		return lhs;
	else
		return rhs;
}

template<typename T, typename R, typename Payload>
class fast_frontier {
	static_assert(std::is_signed_v<T>);
	T dist;
	T m_shift = 0;
	struct P {
		T min_x;
		R min_y;
		Payload payload;
		bool eq(const P& rhs) const { return min_x == rhs.min_x and min_y == rhs.min_y; }
		bool operator<=>(const P&) const = default;
		friend std::ostream& operator<<(std::ostream& out, const P& p) {
			return out << "P { " << p.min_x << ", " << p.min_y << ", " << p.payload << '}';
		}
	};
	std::deque<P> q;
	std::deque<P> down; // sorted by .min_x (t), down.front().min_x is smallest
	inline void commit_down(const T& t) {
		while (!down.empty() and down.front().min_x <= t)
			down.pop_front();
	}
	inline std::optional<P> get_down(const T& t) {
		commit_down(t);
		if (down.empty())
			return std::nullopt;
		else
			return { down.front() };
	}
	std::optional<P> up;
	inline void insert_up(const T& t, P&& tmp) {
		//std::cerr << " insert_up: " << t << " : ";
		//if (!up) std::cerr << "nullopt ";
		//else std::cerr << eval(*up, t) << ' ';
		//std::cerr << "   " << eval(tmp, t) << std::endl;
		if (!up or eval(*up, t) > eval(tmp, t))
			up = std::move(tmp);
	}
	inline void commit_up(const T& t) {
		while (!q.empty() and q.front().min_x <= t)
		{
			insert_up(t, std::move(q.front()));
			q.pop_front();
		}
	}
	inline std::optional<P> get_up(const T& t) {
		commit_up(t);
		return up;
	}
public:
	bool down_is_empty() const { return down.empty(); }
	T max_t() const {
		if (down.empty())
			return T{0};
		else
			return down.back().min_x + m_shift;
	}
	static R eval(const P& p, const T& t) {
		if (p.min_x <= t)
			return (t - p.min_x) + p.min_y;
		else
			return (p.min_x - t) + p.min_y;
	}

	size_t down_queue_size() const { return down.size(); }
	fast_frontier(const T& dist) : dist(dist) {}
	fast_frontier(const fast_frontier&) = default;
	~fast_frontier() = default;
	inline void shift(const T& t) {
		m_shift += t;
	}
	inline void commit(T t) {
		t -= m_shift;
		commit_down(t); commit_up(t);
	}
	using query_result_types = std::tuple<Payload, R, std::tuple<Payload, T, R>>;
	template< typename Res = std::tuple<Payload, T, R> >
	std::optional<Res> query(T t) {
		//std::cerr << " query(" << t << ")  m_shift: " << m_shift << std::endl;
		t -= m_shift;

		auto res = max_opt(get_down(t), get_up(t), [&](const auto& lhs, const auto& rhs) {
			return eval(lhs, t) > eval(rhs, t);
		});
		if (res) {
			if constexpr (std::is_same_v<Res, Payload>)
				return res->payload;
			else if constexpr (std::is_same_v<Res, R>)
				return eval(*res, t);
			else if constexpr (std::is_same_v<Res, std::tuple<Payload, T, R>>) {
				auto val = eval(*res, t);
				return std::make_tuple(res->payload, (t + dist) - res->min_x, std::move(val));
			} else if constexpr (std::is_same_v<Res, P>) {
				res->min_x -= dist - m_shift;
				return res;
			} else
				[]<bool flag = false>() {
					static_assert(flag, "illegal return type for query");
				}();
		} else {
			return std::nullopt;
		}
	}
	void add_up(T t, R r, Payload payload) {
		//std::cerr << " m_shift = " << m_shift << std::endl;
		commit(t);
		t -= m_shift;
		insert_up(t, P { t, r, payload });
	}
	void add(T t, R r, Payload payload) {
		t -= m_shift;

		const P p { t + dist, r, payload };
		const auto eval_down = [](const P& lhs, const P& rhs) {
			assert(lhs.min_x >= rhs.min_x);
			if (!(lhs.min_x >= rhs.min_x)) __builtin_unreachable();
			const auto t = rhs.min_x;//std::min(lhs.min_x, rhs.min_x);
			//return eval(lhs, t) <= eval(rhs, t);
			return eval(lhs, t) <= rhs.min_y;
		};
		while (!down.empty() and down.back().min_x <= t)
			down.pop_back();
		while (!down.empty() and eval_down(p, down.back()))
			down.pop_back();
		if (down.empty() or down.back().min_x < p.min_x or p.min_y < down.back().min_y) {
			q.emplace_back(p);
			down.emplace_back(p);
		}
		assert(std::is_sorted(down.begin(), down.end(), [](const auto& lhs, const auto& rhs) {
			return std::tie(lhs.min_x, lhs.min_y) < std::tie(rhs.min_x, rhs.min_y);
		}));
	}
	void merge(fast_frontier<T, R, Payload>and rhs) {
		assert(dist == rhs.dist);

		if (rhs.down.empty()) {
			if (rhs.up) {
				auto p = std::move(*rhs.up);
				p.min_x += rhs.m_shift - m_shift;
				bool ok = !up;
				if (!ok) {
					const auto t = std::max(p.min_x, up->min_x);
					ok = eval(*up, t) > eval(p, t);
				}
				if (ok)
					up = std::move(p);
			}

			rhs.up = std::nullopt;
			return;
		}

		const size_t mid = down.size();
		std::vector<P> a(mid + rhs.down.size());
		for (size_t i = 0; i < mid; i++)
			a[i] = P
				{ down[i].min_x + m_shift - dist
				, down[i].min_y
				, down[i].payload
				};
		down.clear();

		a.reserve(mid + rhs.down.size());

		for (size_t i = 0; i < rhs.down.size(); i++)
			a[i + mid] = P
				{ rhs.down[i].min_x - dist + rhs.m_shift
				, rhs.down[i].min_y
				, rhs.down[i].payload
				};
		rhs.down.clear();

		if (rhs.up) {
			auto[t,r,payload] = std::move(*rhs.up);
			P p { t + rhs.m_shift - m_shift, r, payload };
			bool ok = !up;
			if (!ok) {
				const auto t = std::max(p.min_x, up->min_x);
				ok = eval(*up, t) > eval(p, t);
			}
			if (ok)
				up = std::move(p);

			rhs.up = std::nullopt;
		}

		q.clear();

		//std::inplace_merge(a.begin(), a.begin() + mid, a.end());
		//a.erase(std::unique(a.begin(), a.end()), a.end());

		//for (const auto&[t,r] : a) {
		//	add(t, r);
		//}

		{ // set union
			size_t l = 0, r = mid;
			while (l < mid and r < a.size()) {
				if (a[l].min_x < a[r].min_x) {
					add(a[l].min_x, a[l].min_y, a[l].payload);
					while (r < a.size() and a[l].eq(a[r])) r++;
					l++;
				} else {
					add(a[r].min_x, a[r].min_y, a[r].payload);
					while (l < mid and a[l].eq(a[r])) l++;
					r++;
				}
			}
			for ( ; l < mid; l++) add(a[l].min_x, a[l].min_y, a[l].payload);
			for ( ; r < a.size(); r++) add(a[r].min_x, a[r].min_y, a[r].payload);
		}
	}
};

template<typename P, typename T>
auto eval_start(P p, const T& t, const T& dist) {
	p.min_x += dist;

	if (p.min_x <= t)
		return (t - p.min_x) + p.min_y;
	else
		return (p.min_x - t) + p.min_y;
}
