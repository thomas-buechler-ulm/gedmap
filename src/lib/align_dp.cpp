#pragma once

#include <string>
#include <vector>
#ifndef NDEBUG
#include <unordered_set>
#endif
#include <cassert>
#include <deque>
#include <tuple>
#include <numeric>
#include <limits>
#include <cstring>
#include <cassert>
#include <math.h>

using Read = std::string_view;
using ReadQual = Read;

namespace edsm_align {

constexpr bool try_no_gap_first = true;

constexpr uint8_t gap_continue_cost = 1;
constexpr uint8_t gap_start_cost = 6;
constexpr uint32_t max_mismatch_cost = 4;
constexpr uint32_t min_mismatch_cost = 0;
constexpr bool variable_gap_costs = false;

constexpr bool matches(char c, char e) {
	return e == 'N' || e == c || c == 'N';
}

constexpr uint8_t min_qual = 0x21, max_qual = 0x7e;
// compute the probability that the base is wrong using the
// Phread quality score
// c is the quality from the FASTQ file: 0x21 (lowest quality) to 0x7e (highest quality)
constexpr float wrong_base_prob(const char& c) {
	assert(min_qual <= c && c <= max_qual);
	const uint8_t d = c - min_qual;
	// d = -10 log_10 p
	// -d/10 = log_10 p
	// p = 10^{-d/10}
	//const auto res = phred_lookup[d];
	const auto res = exp(log(10.f) * d / -10.f);
	assert(res >= 0 && res <= 1);
	return res;
}

template<typename D>
std::array<D, max_qual + 1 - min_qual> gap_start_cost_lookup;
template<typename D>
inline constexpr D read_gap_start_cost(const char& b) {
	if (variable_gap_costs)
		return gap_start_cost_lookup<D>[b - min_qual];
	else
		return 6;
}
template<typename D>
constexpr D read_gap_continue_cost(const char&) {
	return 1;
}

template<typename D>
std::array<D, max_qual + 1 - min_qual> mismatch_cost_lookup;

template<typename D>
void m_init_lookuptables() {
	// mismatch cost
	for (uint8_t i = min_qual; i <= max_qual; i++)
		mismatch_cost_lookup<D>[i - min_qual] = static_cast<D>(
				min_mismatch_cost + (max_mismatch_cost - min_mismatch_cost) * (1.f - wrong_base_prob(i))
			);
	// gap start cost
	for (uint8_t i = min_qual; i <= max_qual; i++)
		gap_start_cost_lookup<D>[i - min_qual] = static_cast<D>(
				1 + 5 * (1.f - wrong_base_prob(i))
			);
}
void init_lookuptables() {
	m_init_lookuptables<uint8_t>();
	m_init_lookuptables<uint16_t>();
	m_init_lookuptables<uint32_t>();
	m_init_lookuptables<uint64_t>();
}

template<typename D, bool match>
constexpr D mismatch_cost(const char& b) {
	//return match ? 0 : 4;
	//return match ? 0 : static_cast<D>(4 * (1.f - wrong_base_prob(b)));
	return match ? 0 : mismatch_cost_lookup<D>[b - min_qual];
}
template<typename D>
constexpr D mismatch_cost(bool match, const char& b) {
	return match ? mismatch_cost<D,true>(b) : mismatch_cost<D,false>(b);
}


class RLE {
private:
	std::deque<std::pair<char, size_t>> q;
public:
	RLE() {}
	RLE& append(char c) {
		if (!q.empty() && c == q.back().first)
			q.back().second++;
		else
			q.emplace_back(c, 1);
		return *this;
	}
	RLE& prepend(char c) {
		if (!q.empty() && c == q.front().first)
			q.front().second++;
		else
			q.emplace_front(c, 1);
		return *this;
	}
	static std::string concat(RLE&& lhs, RLE&& rhs) {
		if (!rhs.q.empty() && !lhs.q.empty() && lhs.q.back().first == rhs.q.front().first)
		{
			lhs.q.back().second += rhs.q.front().second;
			rhs.q.pop_front();
		}
		std::stringstream res;
		for (const auto& q : {std::move(lhs.q), std::move(rhs.q)})
			for (const auto&[c,num] : q)
			{
				if (num > 1) res << num;
				res << c;
			}
		return res.str();
	}
};

using col_id = uint32_t;
using row_id = uint32_t;

static constexpr size_t pad(size_t num_bytes) {
	return ((num_bytes - 1u) | 7u) + 1u;
}
// Don't look at this, it contains evil pointer magic
template<typename D, typename I, bool traceback>
struct Table_ {
	struct Entry {
		D dist;
		char op;
		col_id col;
		I row;
	};
	using EType = std::conditional_t<traceback, Entry, D>;
	std::vector<uint8_t> m_data;
	static constexpr size_t s_size_size() { return pad(2 * sizeof(I)); }
	static constexpr size_t s_data_size(const I& size) { return pad(size * sizeof(EType)); }
	static constexpr size_t get_byte_size(const I& size) { return s_size_size() + s_data_size(size); }
	inline I get_lo(const col_id& col) const { assert(columns.count(col)); return reinterpret_cast<const I*>(m_data.data() + col)[0]; }
	inline I get_hi(const col_id& col) const { assert(columns.count(col)); return reinterpret_cast<const I*>(m_data.data() + col)[1]; }
protected:
#ifndef NDEBUG
	std::unordered_set<col_id> columns;
#endif
	// NOTE: this does not set any data except lo and hi
	inline col_id add_col_(I lo, I hi) {
		//assert (lo <= hi);
		const col_id last_added = m_data.size();
		m_data.resize(last_added + get_byte_size(lo <= hi ? hi + 1 - lo : 0));
		assert(columns.emplace(last_added).second);
		I* bounds = reinterpret_cast<I*>(m_data.data() + last_added);
		bounds[0] = lo, bounds[1] = hi;
		return last_added;
	}

	inline const EType* get_data(const col_id& col) const {
		assert(columns.count(col));
		return reinterpret_cast<const EType*>(m_data.data() + col + s_size_size());
	}
	inline EType* get_data(const col_id& col) { return const_cast<EType*>(std::as_const(*this).get_data(col)); }

	inline const D& get_(const col_id& col, const I& row) const {
		assert(row >= 0 && row <= get_hi(col) - get_lo(col));
		if constexpr (traceback)
			return get_data(col)[row].dist;
		else
			return get_data(col)[row];
	}
	inline D& get_(const col_id& col, I row) { return const_cast<D&>(std::as_const(*this).get_(col, row)); }
public:
	inline void extend(const col_id& col, const I& num, const D& val) {
		if (num == 0) return;
		row_id old_sz;
		{
			auto sz = reinterpret_cast<I*>(m_data.data() + col);
			old_sz = sz[1] + 1 - sz[0];

			assert(col + s_size_size() + s_data_size(old_sz) == m_data.size());
			sz[1] += num;
			m_data.resize(col + s_size_size() + s_data_size(old_sz + num));
		}
		EType* it = reinterpret_cast<EType*>(m_data.data() + col + s_size_size()) + old_sz;
		for (I i = 0; i < num; i++) {
			if constexpr (traceback)
				it[i].dist = val;
			else
				it[i] = val;
		}
		assert(col + s_size_size() + s_data_size(get_hi(col) + 1 - get_lo(col)) == m_data.size());
	}

	inline std::tuple<I,I> get_dim(const col_id& col) const {
		assert(columns.count(col));
		const auto p = reinterpret_cast<const I*>(m_data.data() + col);
		return std::tie(p[0], p[1]);
	}

	inline col_id add_col(I lo, I hi, D def) {
		const auto c_ix = add_col_(lo, hi);
		for (I i = lo; i <= hi; i++) {
			assert(get_lo(c_ix) == lo);
			assert(get_hi(c_ix) == hi);
			get_(c_ix, i - lo) = def;
		}
		return c_ix;
	}
	template<typename It>
	Table_(I lo, I hi, It qual) {
		m_data.reserve(get_byte_size(1000)); // TODO: make input-dependant?
		const auto c_ix = add_col_(lo, hi);
		get_(c_ix, 0) = 0;
		if (lo == hi) return;
		D cost = read_gap_start_cost<D>(qual[1]);
		if constexpr (traceback)
		{
			get_data(0)[1] = Entry { .dist = cost, .op = 'I', .col = 0, .row = 0 };
			for (I i = 2; i <= hi - lo; i++)
				get_data(0)[i] = Entry { .dist = cost += read_gap_continue_cost<D>(qual[i]), .op = 'I', .col = 0, .row = i - 1 };
		}
		else
		{
			get_(c_ix, 1) = cost;
			for (I i = 2; i <= hi - lo; i++)
				get_(c_ix, i) = cost += read_gap_continue_cost<D>(qual[i]);
		}
	}

	inline const D& get(const col_id& col, const I& row) const {
		return get_(col, row - get_lo(col));
	}

	inline const D& update(const col_id& col, const I& i, const D& v, const col_id&, const I&, const char&) {
		assert(i >= get_lo(col) && i <= get_hi(col));
		D& res = get_(col, i - get_lo(col));
		res = std::min(res, v);
		return res;
	}
	inline const D& update_noop(const col_id& col, const I& i, const col_id& ref_col, const I& ref_i) {
		D& res = get_(col, i - get_lo(col));
		const D& ref = get_(ref_col, ref_i - get_lo(ref_col));
		res = std::min(res, ref);
		return res;
	}

	inline void prune(const col_id& col, const D& res, I& lo, I& hi) {
		assert(get_lo(col) <= lo && hi <= get_hi(col));
		while (lo <= hi && get(col, lo) >= res) lo++;
		while (hi >= lo && get(col, hi) >= res) { assert(hi > 0); hi--; }
	}
};

template<typename D, typename I, bool traceback>
struct Table;
template<typename D, typename I>
struct Table<D, I, false> : public Table_<D, I, false> {
	template<typename It>
	Table(I lo, I hi, It read) : Table_<D, I, false>(lo, hi, read) {}
};
template<typename D, typename I>
struct Table<D, I, true> : public Table_<D, I, true> {
public:
	using Base = Table_<D, I, true>;
	using Base::get_hi;
	using Base::get_lo;
	using Base::get_dim;
	using Base::m_data;
	using Base::get_;
	using typename Base::EType;
	using Base::s_size_size;
	using Base::s_data_size;
	template<typename It>
	Table(I lo, I hi, It qual) : Base(lo, hi, qual) { }

	inline const EType& get_pred(const col_id& col, const I& row) const {
		return reinterpret_cast<const EType*>(m_data.data() + col + s_size_size())[row - get_lo(col)];
	}
private:
	inline EType& get_pred(const col_id& col, const I& row) {
		return const_cast<EType&>(std::as_const(*this).get_pred(col, row));
	}
public:

	inline const D& update(const col_id& col, const I& i, const D& v, const col_id& pre_col, const I& pre_i, const char& op) {
		assert(i >= get_lo(col) && i <= get_hi(col));
		auto& res = get_pred(col, i);
		if (v < res.dist) {
			res.dist = v;
			res.col = pre_col;
			res.row = pre_i;
			res.op = op;
		}
		return res.dist;
	}
	inline const D& update_noop(const col_id& col, const I& i, const col_id& ref_col, const I& ref_i) {
		auto& res = get_pred(col, i);
		const auto& ref = get_pred(ref_col, ref_i);
		if (ref.dist < res.dist) {
			res = ref;
		}
		return res.dist;
	}
};

template<
	typename D, // index type for error count
	typename PC, // iterator type for eds
	typename RP, // iterator type for read
	typename QP, // iterator type for quality
	typename JI> // jump index
void align_no_indel(const JI& jump_index,
		PC j, const PC j_end, const RP& read_b, const QP& qual,
		const size_t& read_size,
		std::vector<std::pair<uint32_t, D>> read_pos, // (pos, cost)
		D& res, PC& res_pos
) {
	assert(!read_pos.empty());
	std::vector<std::pair<uint32_t, D>> var_start, var_acc;

	for ( ; res > 0 && j != j_end; ++j) {
		const auto skip_to_var_end = [&] {
			// go to next ')'
			for (PC nx; (nx = std::next(j)) != j_end && *nx != ')'; )
				j = nx;
		};
		const auto skip_to_next_var = [&] {
			// go to next '|' or ')'
			for (PC nx; (nx = std::next(j)) != j_end && *nx != '|' && *nx != ')'; )
				j = nx;
		};

		const auto prune = [&res](std::vector<std::pair<uint32_t, D>>& p) {
			erase_if(p, [&res] (const auto& x) { return x.second >= res; });
		};

		const auto ed_c = *j;
		switch (ed_c) {
			case '#':
			{
				for (auto it : jump_index(j)) {
					align_no_indel<D, PC, RP, QP, JI>(
						jump_index,
						it, j_end, read_b, qual,
						read_size,
						read_pos,
						res, res_pos);
					if (res == 0) return; // optimal solution found
					prune(read_pos);
					if (read_pos.empty()) return;
				}
				return;
			}
			case '(':
			{
				assert(var_start.empty());
				assert(!read_pos.empty());
				var_start = read_pos;
				break;
			}
			case ')':
			{
				if (var_start.empty()) {
					// did not open variant
					//  -> just continue
					continue;
				}

				// terminate variant ending here
				var_acc.insert(var_acc.end(),
					std::make_move_iterator(read_pos.begin()),
					std::make_move_iterator(read_pos.end()));
				read_pos.clear();
				if (var_acc.empty()) {
					// no useful variant
					return;
				}

				// merge
				std::sort(var_acc.begin(), var_acc.end());
				for (size_t i = 0; i < var_acc.size(); ) {
					auto[p, opt] = var_acc[i];
					i++;
					while (i < var_acc.size() && var_acc[i].first == p) {
						//opt = std::min(opt, var_acc[i].second);
						i++;
					}
					if (opt < res) // prune
						read_pos.emplace_back(p, opt);
				}

				if (read_pos.empty()) return;

				var_acc.clear();
				var_start.clear();

				break;
			}
			case '|':
			{
				if (var_start.empty()) {
					// never opened an alternative -> skip to end
					skip_to_var_end();
					if (j != j_end) {
						j++;
						assert (j == j_end || *j == ')');
					} else return;
					break;
				}

				// merge
				var_acc.insert(var_acc.end(),
					std::make_move_iterator(read_pos.begin()),
					std::make_move_iterator(read_pos.end()));
				
				prune(var_start);
				if (var_start.empty()) return;
				read_pos = var_start;
				break;
			}
			default:
			{
				// ed_c is not a syntax symbol
				assert(ed_c == 'N' || ed_c=='A' || ed_c=='C' || ed_c=='G' || ed_c=='T');

				{
					std::vector<std::pair<uint32_t, D>> new_pos;
					for (auto[p,d] : read_pos) {
						assert(d < res);
						if (!matches(read_b[p], ed_c) && (d += mismatch_cost<D,false>(qual[p])) >= res)
							continue;
						if (++p == read_size) {
							assert(res > d);
							res = d;
							prune(new_pos);
							res_pos = j;
							if (res == 0) return; // optimal solution found
						} else
							new_pos.emplace_back(p, d);
					}
					read_pos = std::move(new_pos);
				}
				
				if (read_pos.empty()) {
					if (var_start.empty()) {
						// we are not in an alternatve
						return;
					} else {
						// try next alt.
						skip_to_next_var();
					}
				}
			}
		}
	}
}

struct nil {}; // because we cant have a variable of type void

// traceback = 0: no traceback, returns type pair<D, PC> // (dist, res_pos)
// traceback = 1: traceback, returns nothing
template<bool traceback = false,
	typename D, // index type for error count
	typename PC, // iterator type for eds
	typename RP, // iterator type for read
	typename QP, // iterator type for quality
	typename JI, // jump index
	typename F, // callable for traceback
	typename R = std::conditional_t< // return type
		not traceback,
		std::pair<D,PC>,
		nil>,
	typename Tab = Table<D, row_id, traceback>>
R align(const JI& jump_index, const PC& j, const PC& j_end, const RP& read_b, const RP& read_e, const QP& qual, D max_err, F f) {
	static_assert(std::is_unsigned_v<D>);

	// NOTE: we ignore j -- read_b[0] here and add the possible difference later
	
	const size_t read_size = std::distance(read_b, read_e);
	assert(read_size > 0);
	if (read_size == 1) {
		if constexpr (traceback)
			return nil{};
		else
			return R(0, j);
	}

	D res = max_err + 1;
	PC res_pos = j;

	if (try_no_gap_first and (not traceback || max_err == 0)) { // try alignment with no gaps
		align_no_indel<D, PC, RP, QP, JI>(
			jump_index,
			std::next(j), j_end, read_b, qual,
			read_size,
			{{(uint32_t) 1, (D) 0}},
			res, res_pos);
		if (res <= max_err) {
			if (res == 0) { // TODO: also works for res == mismatch_cost_ (mismatch_cost_ < gap_start_cost)
				// done
				if constexpr (not traceback)
					return R(res, res_pos);
				else {
					for (size_t i = 1; i < read_size; i++)
						f('=');
					return nil{};
				}
			} else {
				max_err = res;
				res++; // restore res = max_err + 1
			}
		} else if (max_err == 0) {
			// impossible
			if constexpr (not traceback) {
				return R(res, res_pos);
			} else {
				assert(false);
			}
		}
	}

	row_id lo = 0, hi;
	{ // determine hi
		//hi = std::min(read_size - 1, (size_t) max_err);
		if (max_err < read_gap_start_cost<D>(qual[1]))
			hi = 0; // we cant have a gap at the start
		else {
			if constexpr (variable_gap_costs) {
				hi = 1;
				D cost = read_gap_start_cost<D>(qual[1]);
				for (size_t i = 2; hi+1 < read_size && (cost + read_gap_continue_cost<D>(qual[i])) <= max_err; i++)
					cost += read_gap_continue_cost<D>(qual[i]), hi++;
			} else {
				hi = 1 + (max_err - read_gap_start_cost<D>(*qual)) / read_gap_continue_cost<D>(qual[0]); // 1+ is for gap start
				hi = std::min((size_t) hi, read_size - 1);
			}
		}
	}
	Tab columns(lo, hi, qual); // initial column
	
	col_id res_col;
	if (hi + 1 == read_size) { // initial col. already contains valid solution
		res = columns.get(0, hi);
		res_col = 0, res_pos = j;
		hi--;
	}

	const auto last_gap_col = columns.add_col(0, 0, max_err+1); // (empty) column for gaps

	align<traceback, D, PC, RP, QP, R, Tab, JI>(
		jump_index,
		std::next(j), j_end, read_b, qual,
		read_size,
		columns,
		{{lo, 0}}, {{hi, 0}},
		res, res_col, res_pos,
		{{0, last_gap_col}}); // last_col (points to col. before cur.)
	
	if constexpr (not traceback)
		return R(res, res_pos);
	else {
		assert(res <= max_err);

		if constexpr (traceback) assert(res == max_err);

		// generate CIGAR string

		{
			col_id res_row = read_size - 1;
			for ( ; res_col > 0 || res_row > 0; ) {
				const auto& e = std::as_const(columns).get_pred(res_col, res_row);
				f(e.op);
				res_col = e.col, res_row = e.row;
			}
		}

		return nil{};
	}
}

template<typename It>
struct eds_reverse_iterator : public It {
	using value_type = typename It::value_type;
	using difference_type = typename It::difference_type;
	constexpr eds_reverse_iterator() = default;
	constexpr explicit eds_reverse_iterator(const It& it) : It(it) {}
	constexpr eds_reverse_iterator(const eds_reverse_iterator<It>& rhs) : It(static_cast<const It&>(rhs)) {}
	constexpr eds_reverse_iterator(eds_reverse_iterator<It>&& rhs) : It(static_cast<It&&>(std::move(rhs))) {}
	constexpr eds_reverse_iterator<It>& operator=(const eds_reverse_iterator<It>& rhs) { It::operator=(rhs); return *this; }
	constexpr eds_reverse_iterator<It>& operator=(eds_reverse_iterator<It>&& rhs) { It::operator=(std::move(rhs)); return *this; }
	~eds_reverse_iterator() = default;

	constexpr value_type operator*() const {
		const auto res = It::operator*();
		switch (res) {
			case '(': return ')';
			case ')': return '(';
			default: return res;
		}
	}
	constexpr eds_reverse_iterator<It>& operator++() { It::operator++(); return *this; } // prefix
	constexpr eds_reverse_iterator<It>  operator++(int x) { return eds_reverse_iterator<It>(It::operator++(x)); } // postfix
	constexpr eds_reverse_iterator<It>& operator+=(difference_type n) { It::operator+=(n); return *this; }
	constexpr eds_reverse_iterator<It>  operator+(difference_type n) const { return eds_reverse_iterator<It>(It::operator+(n)); }
	constexpr eds_reverse_iterator<It>& operator--() { It::operator--(); return *this; } // prefix
	constexpr eds_reverse_iterator<It>  operator--(int x) { return eds_reverse_iterator<It>(It::operator--(x)); } // postfix
	constexpr eds_reverse_iterator<It>& operator-=(difference_type n) { It::operator+=(n); return *this; }
	constexpr eds_reverse_iterator<It>  operator-(difference_type n) const { return eds_reverse_iterator<It>(It::operator-(n)); }
};
template<typename It>
bool operator==(const eds_reverse_iterator<It>& lhs, const eds_reverse_iterator<It>& rhs) { return static_cast<const It&>(lhs) == static_cast<const It&>(rhs); }
template<typename It>
bool operator!=(const eds_reverse_iterator<It>& lhs, const eds_reverse_iterator<It>& rhs) { return static_cast<const It&>(lhs) != static_cast<const It&>(rhs); }

template<typename D, typename RP, typename QP, typename Tab>
inline void align_dp(
	Tab& columns,
	const std::array<col_id, 2>& last_cols,
	const std::array<col_id, 2>& new_cols,
	const std::array<row_id, 2>& lo, std::array<row_id, 2>& hi, // hi[0] may change when read_gap_start_cost is not constant
	const RP& read_b, const QP& qual, const char& ed_c,
	const D& res, const size_t& read_size
) {
	const auto[l_lo,l_hi] = columns.get_dim(last_cols[0]);

	std::vector<D> read_del(hi[0] + 1 - lo[0]);
	row_id max_hi = hi[0];

	for (row_id i = lo[0]; i <= hi[0]; i++) {
		if (l_lo <= i && i <= l_hi) // start gap, delete ed_c
			columns.update(new_cols[1], i, columns.get(last_cols[0], i) + gap_start_cost, last_cols[0], i, 'D');
		if (lo[1] <= i && i <= hi[1]) // end gap in eds
			columns.update_noop(new_cols[0], i, new_cols[1], i);


		if (l_lo <= i - 1 && i - 1 <= l_hi) {
			const bool match = matches(read_b[i], ed_c);
			columns.update(new_cols[0], i,
				columns.get(last_cols[0], i - 1) + mismatch_cost<D>(match, qual[i]),
				last_cols[0], i-1, match ? '=' : 'X');
		}
		if (i > lo[0]) {
			const auto v = read_del[i - lo[0] - 1] = columns.get(new_cols[0], i - 1) + read_gap_start_cost<D>(qual[i]); // start gap in read
			if (v < res)
				max_hi = std::max(max_hi, i + (res - v - 1));
		}
	}
	if (const auto i = hi[0] + 1; i < read_size) {
		const auto v = read_del[i - lo[0] - 1] = columns.get(new_cols[0], i - 1) + read_gap_start_cost<D>(qual[i]);
		if (v < res)
			max_hi = std::max(max_hi, i + (res - v - 1));
	}
	max_hi = std::min((size_t) max_hi, read_size - 1);

	columns.extend(new_cols[0], max_hi - hi[0], res);
	const auto read_delete_col = columns.add_col(lo[0] + 1, max_hi, res);

	for (row_id i = lo[0] + 1; i <= std::min(hi[0]+1, max_hi); i++) {
		columns.update(read_delete_col, i, read_del[i - (lo[0] + 1)], new_cols[0], i - 1, 'I'); // start gap in read
		if (i > lo[0] + 1)
			columns.update(read_delete_col, i, columns.get(read_delete_col, i - 1) + read_gap_continue_cost<D>(qual[i]), read_delete_col, i - 1, 'I'); // continue gap in read
		columns.update_noop(new_cols[0], i, read_delete_col, i);
	}
	if (max_hi > hi[0]) {
		for (row_id i = std::max(hi[0], lo[0]+1) + 1; i <= max_hi; i++) {
			columns.update(read_delete_col, i, columns.get(read_delete_col, i - 1) + read_gap_continue_cost<D>(qual[i]), read_delete_col, i - 1, 'I');
			columns.update_noop(new_cols[0], i, read_delete_col, i);
		}

		hi[0] = max_hi;
	}
}

template<bool traceback,
	typename D, // index type for error count
	typename PC, // iterator type for eds
	typename RP, // iterator type for read
	typename QP, // iterator type for quality
	typename R, // return type
	typename Tab = Table<D, row_id, traceback>,
	typename JI> // jump index
void align(const JI& jump_index,
		PC j, const PC& j_end, const RP& read_b, const QP& qual,
		const size_t& read_size,
		Tab& columns,
		std::array<row_id, 2> lo, std::array<row_id, 2> hi,
		D& res, col_id& res_col, PC& res_pos,
		std::array<col_id, 2> last_cols
) {
	constexpr size_t INVALID = std::numeric_limits<col_id>::max();
	std::array<row_id, 2> var_lo{{0,0}}, var_hi = var_lo; // (no_gap, gap), initialisation is only to make the warnings go away
	std::array<col_id, 2> var_start = {{INVALID, INVALID}}; // (no_gap, gap)
	std::vector<std::array<col_id, 2>> variants;

	for ( ; res > 0 && j != j_end; ++j) {
		const auto skip_to_var_end = [&] {
			// go to next ')'
			for (PC nx; (nx = std::next(j)) != j_end && *nx != ')'; )
				j = nx;
		};
		const auto skip_to_next_var = [&] {
			// go to next '|' or ')'
			for (PC nx; (nx = std::next(j)) != j_end && *nx != '|' && *nx != ')'; )
				j = nx;
		};
		const auto commit_variant = [&] {
			variants.emplace_back(last_cols);
			assert(lo[0] >= columns.get_lo(last_cols[0]));
			assert(hi[0] <= columns.get_hi(last_cols[0]));
			for (size_t k = 0; k < 2; k++)
				if (lo[k] <= hi[k]) {
					var_lo[k] = std::min(var_lo[k], lo[k]);
					var_hi[k] = std::max(var_hi[k], hi[k]);
				}
		};

		const auto ed_c = *j;

		switch (ed_c) {
			case '#':
			{
				for (auto it : jump_index(j)) {
					align<traceback, D, PC, RP, QP, R, Tab, JI>(
						jump_index,
						it, j_end, read_b, qual,
						read_size,
						columns,
						lo, hi,
						res, res_col, res_pos,
						last_cols);
					if (res == 0) return; // optimal solution found
					for (size_t k = 0; k < 2; k++) // res might have changed
						columns.prune(last_cols[k], res, lo[k], hi[k]);
					if (lo[0] > hi[0]) break;
				}
				return;
			}
			case '(':
			{
				assert(var_start[0] == INVALID);
				assert(variants.empty());
				var_lo.fill(read_size);
				var_hi.fill(0);
				var_start = last_cols;
				break;
			}
			case ')':
			{
				if (var_start[0] == INVALID) {
					// did not open variant
					//  -> just continue
					continue;
				}

				// terminate variant ending here
				if (lo[0] <= hi[0]) {
					commit_variant();
				} else if (variants.empty()) {
					// no useful variant
					return;
				}

				// merge
				lo = var_lo, hi = var_hi;
				const auto c_ix = columns.add_col(lo[0], hi[0], res);
				const auto c_gap_ix = columns.add_col(lo[1], hi[1], res);
				last_cols = {{c_ix, c_gap_ix}};
				for (const auto& var : variants) {
					for (size_t k = 0; k < 2; k++)
						for (row_id i = std::max(lo[k], columns.get_lo(var[k]));
								i <= std::min(hi[k], columns.get_hi(var[k]));
								i++)
							columns.update_noop(last_cols[k], i, var[k], i);
				}
				for (size_t k = 0; k < 2; k++)
					columns.prune(last_cols[k], res, lo[k], hi[k]);
				if (lo[0] > hi[0]) return;

				variants.clear();
				var_start[0] = INVALID;
				break;
			}
			case '|':
			{
				if (var_start[0] == INVALID) {
					// never opened an alternative -> skip to end
					skip_to_var_end();
					if (j != j_end) {
						j++;
						assert (j == j_end || *j == ')');
					} else return;
					break;
				}
				if (lo[0] <= hi[0]) {
					commit_variant();
				}

				// start new alt
				last_cols = var_start;
				for (size_t k = 0; k < 2; k++) {
					std::tie(lo[k], hi[k]) = columns.get_dim(last_cols[k]);
					columns.prune(last_cols[k], res, lo[k], hi[k]);
				}
				if (lo[0] > hi[0]) return; // no other alternative can beat optimum
				break;
			}
			default:
			{
				assert(lo[0] <= hi[0]);
				// ed_c is not a syntax symbol
				assert(ed_c == 'N' || ed_c=='A' || ed_c=='C' || ed_c=='G' || ed_c=='T');

				if (hi[0] + 1 < read_size) hi[0]++;

				const auto[l_lo_gap, l_hi_gap] = std::make_pair(lo[1], hi[1]);

				lo[1] = lo[0], hi[1] = hi[0];
				if (hi[1] + 1 < read_size) hi[1]++;

				std::array<col_id, 2> new_cols = [&] {
					auto eds_gap_col = columns.add_col(lo[1], hi[1], res);
					auto no_gap = columns.add_col(lo[0], hi[0], res);
					return std::array<col_id, 2>{{ no_gap, eds_gap_col }};
				}();

				for (row_id i = std::max(lo[1], l_lo_gap); i <= std::min(hi[1], l_hi_gap); i++) // continue gap in eds
					columns.update(new_cols[1], i, columns.get(last_cols[1], i) + gap_continue_cost, last_cols[1], i, 'D');

				align_dp<D,RP,QP,Tab>(columns, last_cols, new_cols, lo, hi, read_b, qual, ed_c, res, read_size);
				
				// prune
				columns.prune(new_cols[0], res, lo[0], hi[0]);
				columns.prune(new_cols[1], res, lo[1], hi[1]);
				
				if (lo[0] <= read_size-1 && read_size-1 <= hi[0]) {
					// we have a solution

					// this holds because of the pruning above
					assert(columns.get(new_cols[0], read_size - 1) < res);
					res = columns.get(new_cols[0], read_size - 1);
					res_pos = j;
					if constexpr (traceback)
						res_col = new_cols[0];
					if (res == 0) return; // optimal solution found
				}

				last_cols = new_cols;

				if (lo[0] > hi[0]) {
					if (var_start[0] == INVALID) {
						// we are not in an alternatve
						return;
					} else {
						// try next alt.
						skip_to_next_var();
					}
				}
			}
		}
	}
}

} // namespace edsm_align

// when traceback==true, it is assumed that the optimal alignment has exactly
// max_err errors
// returns ((dist in fwd dir, dist ind backw dir), start in eds) when traceback == false, and
// cigar otherwise
// NOTE: dist in backwards direction includes read[i], j
template<bool traceback = false,
	typename D, // error type
	typename PC, // index for eds
	typename ADJ, // adjacency
	typename R = std::conditional_t<traceback, std::string, std::pair<std::pair<D,D>, size_t>>,
	typename MD = std::conditional_t<traceback, std::pair<D,D>, D>>
R align_dp(std::string_view eds, const ADJ& adj, PC j, const Read& read, const ReadQual& qual, size_t i, MD max_err) {
	using namespace edsm_align;

	if (read.size() == 0)
		return R{};

	assert(eds[j] == 'N' || eds[j] == 'A' || eds[j] == 'G' || eds[j] == 'T' || eds[j] == 'C');
	assert(qual.size() == read.size());

	const eds_reverse_iterator backwards_start(eds.crbegin() + eds.size() - 1 - j);
	
	// check whether j and read[i] match
	char match_symb = '=';
	if (!edsm_align::matches(eds[j], read[i])) {
		match_symb = 'X';
		const D cost = edsm_align::mismatch_cost<D,false>(qual[i]);
		if constexpr (traceback) {
			assert(cost <= max_err.second);
			max_err.second -= cost;
		} else {
			if (cost > max_err) return std::make_pair(std::make_pair(0, max_err + 1), size_t{0});
			max_err -= cost;
		}
	}
	RLE cigar_backward{};

	auto err_b = align<traceback, D>(
		[&](auto it) {
			assert(*it == '#');
			const auto i = std::distance(it, eds_reverse_iterator(eds.crend())) - 1;
			assert(eds[i] == '#');
			std::vector<eds_reverse_iterator<std::string_view::const_reverse_iterator>> res;
			for (auto p : adj(i, ADJ::BACKWARD)) {
				auto p_it = eds_reverse_iterator(eds.crbegin() + eds.size() - p);
				assert(*(p_it-1) == '#');
				res.emplace_back(p_it);
			}
			return res;
		},
		backwards_start,
		eds_reverse_iterator(eds.crend()),
		read.crbegin() + read.size() - 1u - i, read.crend(),
		qual.crbegin() + read.size() - 1u - i,
		[&] { if constexpr (traceback) return max_err.second; else return max_err; }(),
		[&](char c) { cigar_backward.append(c); });

	if constexpr (traceback)
		cigar_backward.append(match_symb);

	if constexpr (!traceback)
		if (err_b.first > max_err)
			return std::make_pair(std::make_pair(0, err_b.first), size_t{0});

	RLE cigar_forward{};
	const auto err_f = align<traceback, D>(
		[&](auto it) {
			assert(*it == '#');
			const auto i = std::distance(eds.cbegin(), it);
			assert(eds[i] == '#');
			std::vector<std::string_view::const_iterator> res;
			for (auto p : adj(i, ADJ::FORWARD)) {
				auto p_it = eds.cbegin() + p + 1;
				assert(*(p_it-1) == '#');
				res.emplace_back(p_it);
			}
			return res;
		},
		eds.cbegin() + j,
		eds.cend(),
		read.cbegin() + i, read.cend(),
		qual.cbegin() + i,
		[&] { if constexpr (traceback) return max_err.first; else return max_err - err_b.first; }(),
		[&](char c) { cigar_forward.prepend(c); });
	if constexpr (not traceback) {
		const size_t res_pos = err_f.first + err_b.first > max_err
			? 0
			: j - std::distance(backwards_start, err_b.second);
		return std::make_pair(std::make_pair(err_f.first, err_b.first), res_pos);
	} else {
		return RLE::concat(std::move(cigar_backward), std::move(cigar_forward));
	}
}
