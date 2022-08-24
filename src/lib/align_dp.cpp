#pragma once

#include "adjacency.cpp"

#include <string>
#include <vector>
#include <cassert>
#include <tuple>
#include <numeric>
#include <iostream>
#include <limits>
#include <cstring>

using Read = std::string_view;

namespace edsm_align {

constexpr bool matches(char c, char e) {
	return e == 'N' || e == c || c == 'N';
}

class RLE {
private:
	std::string res;
	std::pair<char, size_t> cur;
	void commit() {
		if (cur.second > 1)
			res += std::to_string(cur.second) + cur.first;
		else if (cur.second == 1)
			res += cur.first;
		cur.second = 0;
	}
public:
	RLE& operator+=(char c) {
		if (cur.first == c)
			cur.second++;
		else {
			commit();
			cur = std::make_pair(c, 1);
		}
		return *this;
	}
	RLE& append(char c) {
		return *this += c;
	}
	template<typename It>
	RLE& append(It b, const It e) {
		while (b != e)
			append(*(b++));
		return *this;
	}
	std::string&& finalize() {
		commit();
		return std::move(res);
	}
};

static constexpr size_t pad(size_t num_bytes) {
	return ((num_bytes - 1u) | 7u) + 1u;
}
// Don't look at this, it contains evil pointer magic
template<typename D, typename I, bool traceback>
struct Table_ {
	std::vector<uint8_t> m_data;
	static constexpr size_t get_bytes(I size) {
		return pad(2*sizeof(I)) // lo, hi
			+ pad(size * sizeof(D)) // actual DP data
			+ (traceback ? pad(size * sizeof(std::pair<uint32_t,I>)) : 0u); // possible data for backtracking
	}
	inline D& get_(uint32_t col, I row) {
		return reinterpret_cast<D*>(m_data.data() + col + pad(2*sizeof(I)))[row];
	}
	inline const D& get_(uint32_t col, I row) const {
		return reinterpret_cast<const D*>(m_data.data() + col + pad(2*sizeof(I)))[row];
	}
	uint32_t last_added;
	inline uint32_t add_col_(I lo, I hi) {
		last_added = m_data.size();
		m_data.resize(last_added + get_bytes(hi + 1 - lo));
		return last_added;
	}
	inline I& get_lo_(const uint32_t& col) {
		return reinterpret_cast<I*>(m_data.data() + col)[0];
	}
	inline I& get_hi_(const uint32_t& col) {
		return reinterpret_cast<I*>(m_data.data() + col)[1];
	}
	inline const I& get_lo(const uint32_t& col) const {
		return reinterpret_cast<const I*>(m_data.data() + col)[0];
	}
	inline const I& get_hi(const uint32_t& col) const {
		return reinterpret_cast<const I*>(m_data.data() + col)[1];
	}

	inline uint32_t add_col(I lo, I hi, D def) {
		const auto c_ix = add_col_(lo, hi);
		get_lo_(c_ix) = lo;
		get_hi_(c_ix) = hi;
		for (I i = lo; i <= hi; i++)
			get_(c_ix, i - lo) = def;
		return c_ix;
	}
	Table_(I lo, I hi) {
		m_data.reserve(get_bytes(200));
		const auto c_ix = add_col_(lo, hi);
		get_lo_(c_ix) = lo, get_hi_(c_ix) = hi;
		for (I i = lo; i <= hi; i++)
			get_(c_ix, i - lo) = i - lo;
	}

	inline const D& get(const uint32_t& col, const I& row) const {
		return get_(col, row - get_lo(col));
	}

	inline const D& update(const I& i, const D& v, const uint32_t&, const I&) {
		assert(i >= get_lo(last_added) && i <= get_hi(last_added));
		D& res = get_(last_added, i - get_lo(last_added));
		res = std::min(res, v);
		return res;
	}

	inline void prune(const uint32_t& col_ix, const D& res, I& lo, I& hi) {
		while (lo <= hi && get(col_ix, lo) >= res) lo++;
		while (hi >= lo && get(col_ix, hi) >= res) { assert(hi > 0); hi--; }
	}
};
template<typename D, typename I, bool traceback>
struct Table;
template<typename D, typename I>
struct Table<D, I, false> : public Table_<D, I, false> {
	Table(I lo, I hi) : Table_<D, I, false>(lo, hi) {}
};
template<typename D, typename I>
struct Table<D, I, true> : public Table_<D, I, true> {
	using Base = Table_<D, I, true>;
	using Base::last_added;
	using Base::get_hi;
	using Base::get_lo;
	using Base::m_data;
	using Base::get_;
	using TB_type = std::pair<uint32_t, I>;
	Table(I lo, I hi) : Base(lo, hi) {}

	inline TB_type& get_pred(uint32_t col, I row) {
		const auto s = get_hi(col) + 1 - get_lo(col);
		return reinterpret_cast<TB_type*>(m_data.data() + col + pad(2 * sizeof(I)) + pad(s * sizeof(D)))[row - get_lo(col)];
	}
	inline const TB_type& get_pred(uint32_t col, I row) const {
		const auto s = get_hi(col) + 1 - get_lo(col);
		return reinterpret_cast<const TB_type*>(m_data.data() + col + pad(2 * sizeof(I)) + pad(s * sizeof(D)))[row - get_lo(col)];
	}

	inline const D& update(const I& i, const D& v, const uint32_t& pre_col, const I& pre_i) {
		assert(i >= get_lo(last_added) && i <= get_hi(last_added));
		D& res = get_(last_added, i - get_lo(last_added));
		if (v < res) {
			res = v;
			get_pred(last_added, i) = std::make_pair(pre_col, pre_i);
		}
		return res;
	}
};



// traceback = 0: no traceback, returns type D
// traceback = 1: traceback, returns a (reversed) operation string (cigar without RLE)
// traceback = 2: traceback, returns <error, distance from j to read_e in opt. alignment, rev. operation string>
template<uint8_t traceback = 0,
	typename D, // index type for error count
	typename PC, // iterator type for eds
	typename RP, // iterator type for read
	typename JI, // jump index
	typename R = std::conditional_t< // return type
		traceback == 0,
		D,
		std::conditional_t<traceback == 1, std::string, std::tuple<D, size_t, std::string>>>,
	typename Tab = Table<D, size_t, traceback != 0>>
R align(const JI& jump_index, PC j, const PC j_end, const RP& read_b, const RP& read_e, D max_err) {
	static_assert(std::is_unsigned_v<D>);
	static_assert(traceback <= 2);
	
	const size_t read_size = std::distance(read_b, read_e);
	assert(read_size > 0);
	if (read_size == 1) {
		if constexpr (traceback == 0)
			return 0;
		else if constexpr (traceback == 1)
			return "";
		else
			return R(0, 0, "");
	}

	assert(matches(*read_b, *j));

	size_t lo = 0, hi = std::min(static_cast<size_t>(read_size - 1u), (size_t)max_err); // inclusively
	Tab columns(lo, hi); // initial column
	
	D res = max_err + 1;
	size_t res_col;
	PC res_pos;
	size_t res_row; // TODO: if there is a solution, shouldn't this always be read_size - 1
	if (hi + 1 == read_size) { // initial col. already contains valid solution
		res = columns.get(0, hi);
		res_col = 0, res_row = hi, res_pos = j;
		hi--;
	}

	align<traceback, D, PC, RP, R, Tab, JI>(
		jump_index,
		std::next(j), j_end, read_b, max_err,
		read_size,
		columns,
		lo, hi,
		res, res_col, res_pos, res_row,
		0); // last_col (points to col. before cur.)
	
	if constexpr (traceback == 0)
		return res;
	else {
		assert(res <= max_err);

		// generate CIGAR string

		std::string cigar;
		for (D tmp_res = res; res_col > 0 || res_row > 0; ) {
			const auto[next_col, next_row] = columns.get_pred(res_col, res_row);
			const auto next_res = columns.get(next_col, next_row);
			if (next_col == res_col) {
				assert(next_row + 1 == res_row);
				assert(next_res + 1 == tmp_res);
				cigar += 'I';
			} else if (next_row == res_row) {
				if (next_res != tmp_res) // may be equal when merging alternatives
					cigar += 'D';
			} else {
				assert(next_row + 1 == res_row);
				assert(next_res <= tmp_res && next_res + 1 >= tmp_res);
				cigar += (next_res != tmp_res ? 'X' : '=');
			}
			res_col = next_col, res_row = next_row, tmp_res = next_res;
		}

		if constexpr (traceback == 1)
			return cigar;
		else
			return R(res, std::distance(j, res_pos), std::move(cigar));
	}
}

template<uint8_t traceback,
	typename D, // index type for error count
	typename PC, // iterator type for eds
	typename RP, // iterator type for read
	typename R, // return type
	typename Tab = Table<D, size_t, traceback != 0>,
	typename JI> // jump index
void align(const JI& jump_index,
		PC j, const PC j_end, const RP& read_b, D max_err,
		const size_t read_size,
		Tab& columns,
		size_t lo, size_t hi,
		D& res, size_t& res_col, PC& res_pos, size_t& res_row,
		size_t last_col) {

	constexpr size_t INVALID = std::numeric_limits<size_t>::max();

	size_t var_lo = 0, var_hi = 0;//=0 does nothing except it gets rid of a compiler warning
	size_t var_start = INVALID;
	std::vector<size_t> variants;

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

		const auto ed_c = *j;
		switch (ed_c) {
			case '#':
			{
				for (auto it : jump_index(j)) {
					align<traceback, D, PC, RP, R, Tab, JI>(
						jump_index,
						it, j_end, read_b, max_err,
						read_size,
						columns,
						lo, hi,
						res, res_col, res_pos, res_row,
						last_col);
					columns.prune(last_col, res, lo, hi);
					if (lo > hi) break;
				}
				return;
			}
			case '(':
			{
				assert(var_start == INVALID);
				assert(variants.empty());
				var_lo = read_size, var_hi = 0;
				var_start = last_col;
				break;
			}
			case ')':
			{
				if (var_start == INVALID) {
					// did not open variant
					//  -> just continue
					continue;
				}

				// terminate variant ending here
				if (lo <= hi) {
					variants.emplace_back(last_col);
					var_lo = std::min(var_lo, lo);
					var_hi = std::max(var_hi, hi);
				} else if (variants.empty()) {
					// no useful variant
					return;
				}

				// merge
				lo = var_lo, hi = var_hi;
				const auto c_ix = columns.add_col(lo, hi, max_err + 1);
				for (auto jj : variants) {
					for (size_t i = std::max(lo, columns.get_lo(jj));
						 i <= std::min(hi, columns.get_hi(jj));
						 i++) {
						columns.update(i, columns.get(jj, i), jj, i);
					}
				}
				columns.prune(c_ix, res, lo, hi);
				if (lo > hi) return;

				variants.clear();
				var_start = INVALID;

				last_col = c_ix;
				break;
			}
			case '|':
			{
				if (var_start == INVALID) {
					// never opened an alternative -> skip to end
					skip_to_var_end();
					if (j != j_end) {
						j++;
						assert (j == j_end || *j == ')');
					} else return;
					break;
				}
				if (lo <= hi) {
					variants.emplace_back(last_col);
					assert(lo >= columns.get_lo(last_col));
					assert(hi <= columns.get_hi(last_col));
					var_lo = std::min(var_lo, lo);
					var_hi = std::max(var_hi, hi);
				}

				// start new alt
				last_col = var_start;
				lo = columns.get_lo(last_col);
				hi = columns.get_hi(last_col);
				columns.prune(last_col, res, lo, hi);
				if (lo > hi) return; // no other alternative can beat optimum
				break;
			}
			default:
			{
				assert(lo <= hi);
				// ed_c is not a syntax symbol
				assert(ed_c == 'N' || ed_c=='A' || ed_c=='C' || ed_c=='G' || ed_c=='T');
				if (hi+1 < read_size) hi++;
				const auto c_ix = columns.add_col(lo, hi, max_err + 1);
				for (size_t i = lo; i <= hi; i++) {
					if (i > columns.get_lo(last_col))
						columns.update(i, columns.get(last_col, i-1) + !matches(read_b[i], ed_c), last_col, i-1);
					if (i <= columns.get_hi(last_col))
						columns.update(i, columns.get(last_col, i) + 1, last_col, i);
					if (i > lo)
						columns.update(i, columns.get(c_ix, i-1) + 1, c_ix, i-1);
				}
				
				if (lo < read_size && read_size <= hi+1) {
					assert(columns.get(c_ix, read_size - 1) <= res);
					res = columns.get(c_ix, read_size - 1);
					if constexpr (traceback != 0) {
						res_col = c_ix;
						res_row = read_size - 1;
						if constexpr (traceback == 2) res_pos = j;
					}
				}
				
				// prune
				columns.prune(c_ix, res, lo, hi);

				last_col = c_ix;

				if (lo > hi) {
					if (var_start == INVALID) {
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
	constexpr eds_reverse_iterator<It>  operator+(difference_type n) { return eds_reverse_iterator<It>(It::operator+(n)); }
	constexpr eds_reverse_iterator<It>& operator--() { It::operator--(); return *this; } // prefix
	constexpr eds_reverse_iterator<It>  operator--(int x) { return eds_reverse_iterator<It>(It::operator--(x)); } // postfix
	constexpr eds_reverse_iterator<It>& operator-=(difference_type n) { It::operator+=(n); return *this; }
	constexpr eds_reverse_iterator<It>  operator-(difference_type n) { return eds_reverse_iterator<It>(It::operator-(n)); }
};
template<typename It>
bool operator==(const eds_reverse_iterator<It>& lhs, const eds_reverse_iterator<It>& rhs) { return static_cast<const It&>(lhs) == static_cast<const It&>(rhs); }
template<typename It>
bool operator!=(const eds_reverse_iterator<It>& lhs, const eds_reverse_iterator<It>& rhs) { return static_cast<const It&>(lhs) != static_cast<const It&>(rhs); }

} // namespace edsm_align

// when traceback==true, it is assumed that the optimal alignment has exactly
// (at most?) max_err errors
template<bool traceback = false,
	typename D, // error type
	typename PC, // index for eds
	typename R = std::conditional_t<traceback, std::pair<std::string, size_t>, D>>
R align(std::string_view eds, const adjacency& adj, PC j, const Read& read, size_t i, D max_err) {
	using namespace edsm_align;
	if (read.size() == 0)
		return R{};
	
	const auto err_b = align<2*traceback, D>(
		[&](auto it) {
			assert(*it == '#');
			const auto i = std::distance(it, eds_reverse_iterator(eds.crend())) - 1;
			assert(eds[i] == '#');
			std::vector<eds_reverse_iterator<std::string_view::const_reverse_iterator>> res;
			for (auto p : adj(i, adjacency::BACKWARD)) {
				auto p_it = eds_reverse_iterator(eds.crbegin() + eds.size() - p);
				assert(*(p_it-1) == '#');
				res.emplace_back(p_it);
			}
			return res;
		},
		eds_reverse_iterator(eds.crbegin() + eds.size() - 1 - j),
		eds_reverse_iterator(eds.crend()),
		read.crbegin() + read.size() - 1u - i,
		read.crend(),
		max_err);

	if constexpr (!traceback)
		if (err_b > max_err)
			return err_b;

	const auto err_f = align<traceback, D>(
		[&](auto it) {
			assert(*it == '#');
			const auto i = std::distance(eds.cbegin(), it);
			assert(eds[i] == '#');
			std::vector<std::string_view::const_iterator> res;
			for (auto p : adj(i, adjacency::FORWARD)) {
				auto p_it = eds.cbegin() + p + 1;
				assert(*(p_it-1) == '#');
				res.emplace_back(p_it);
			}
			return res;
		},
		eds.cbegin() + j,
		eds.cend(),
		read.cbegin() + i,
		read.cend(),
		[&] { 
			if constexpr (traceback)
				return max_err - std::get<0>(err_b);
			else
				return max_err - err_b;
		}());
	if constexpr (!traceback)
		return err_b + err_f;
	else {
		auto cigar = RLE()
			.append(std::get<2>(err_b).cbegin(), std::get<2>(err_b).cend())
			.append('=')
			.append(err_f.crbegin(), err_f.crend())
			.finalize();
		assert(j >= std::get<1>(err_b));
		return std::make_pair(std::move(cigar), j - std::get<1>(err_b));
	}
}
