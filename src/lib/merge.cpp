/**
 * @author Thomas Buechler, Jannik Olbrich
 */
#pragma once

#include "arraywrapper.cpp"

/**
 * @author Jannik Olbrich
 *
 *@brief Given a vector A and a vector of boundaries in A, such that the subvectors of A (seperated by the boundaries) are already sorted.
 * Only neighbouring subvectors are merged. This procedure optimises the order
 * in which the subvectors are merged under this constraint in such a way that the number
 * of comparisons is minimised.
 *
 * O(m^3 + k) where m is the number of subvectors and k is the minimum number of
 * operations required for merging subvectors using the inplace_merge method.
 * Really only useful when m is small and the original subvectors vary
 * significantly in their length.
 */
template<typename Iterator>
void opt_merge_nb(Iterator A, std::vector<uint32_t>& a) {
	if (a.size() <= 2)
		return;

	const auto n = a.size() - 1;
	using MemT = std::pair<size_t, uint32_t>; // <num ops, total size>
	ArrayWrapper<MemT, 2> mem(n, n+1);
	ArrayWrapper<uint32_t, 2> split_mem(n, n+1);
	const auto comp = [](const MemT& a, const MemT& b) -> MemT {
		const auto& [o1, s1] = a;
		const auto& [o2, s2] = b;
		return std::make_pair(o1 + o2 + (s1 + s2 - 1), s1 + s2);
	};

	for (uint32_t i = 0; i < n; i++) {
		mem(i, i+1) = std::make_pair(0, a[i+1] - a[i]);
	}
	for (uint32_t len = 2; len <= n; len++) {
		for (uint32_t i = 0; i + len <= n; i++) {
			const uint32_t e = i + len;

			uint32_t opt_split = i + 1;
			auto opt_ops = comp(mem(i, opt_split), mem(opt_split, e));
			for (uint32_t split = i + 2; split < e; split++) {
				const auto ops = comp(mem(i, split), mem(split, e));
				if (ops < opt_ops)
					opt_ops = ops, opt_split = split;
			}
			mem(i, e) = opt_ops;
			split_mem(i, e) = opt_split;
		}
	}

	const std::function<void(uint32_t, uint32_t)> mmerge = [&](uint32_t i, uint32_t e) {
		if (i + 1 >= e) return;

		const auto split = split_mem(i, e);
		mmerge(i, split);
		mmerge(split, e);
		inplace_merge(A + a[i], A + a[split], A + a[e]);
	};

	mmerge(0, n);
}

/**
 * @author Jannik Olbrich
 * @brief Given a vector A and a vector of boundaries in A, such that the subvectors of A (seperated by the boundaries) are already sorted.
 * This procedure optimises the order * in which the subvectors are merged
 * in such a way that the number of comparisons is minimised.
 * NOTE: This procedure allocates much more memory than opt_merge_nb and may
 * hence be substantially slower. 
 */
template<typename Iterator, typename T = std::remove_cvref_t<decltype(*std::declval<Iterator>())>>
void opt_merge(Iterator A, std::vector<uint32_t>& a) {
	if (a.size() <= 2)
		return;

	const auto cmp = [](const auto& a, const auto& b) {
		return a.size() > b.size();
	};
	std::priority_queue<std::vector<T>, std::vector<std::vector<T>>, decltype(cmp)> q(cmp);
	for (size_t i = 0; i + 1 < a.size(); i++)
		q.emplace(std::vector(A + a[i], A + a[i+1]));
	while (q.size() > 1) {
		const auto s1 = std::move(q.top()); q.pop();
		const auto s2 = std::move(q.top()); q.pop();
		std::vector<T> res(s1.size() + s2.size());
		std::merge(s1.begin(), s1.end(),
				   s2.begin(), s2.end(),
				   res.begin());
		q.emplace(std::move(res));
	}
	const auto res = std::move(q.top());
	for (auto& x : res)
		*(A++) = std::move(x);
}

/**
 *@brief Given a vector A and a vector of boundaries in A, such that the subvectors of A (seperated by the boundaries) are already sorted.
 * This method sorts A by iteratively merging the subvectors
 */
template <class Iterator>
void multi_merge( Iterator A, std::vector<uint32_t> & a) {
	opt_merge_nb(A, a);
	/*
	while(a.size() > 2){
		//merge interval A[a_2i..a_2i+1] with A[a_2i+1..a_2i+2]
		for(unsigned int i = 0; (2*i+2) < a.size(); i++) {
			inplace_merge(A+a[2*i], A+a[2*i+1], A+a[2*i+2]);
		}
		//update boundaries
		for(unsigned int i = 1; (2*i) < a.size(); i++)		a[i] = a[2*i];
		if(a.size() % 2 == 0) a[a.size()/2] = a[a.size()-1];
		a.resize(a.size()/2 + 1);
	}
	*/
}


/**
 * SOURCE: https://cw.fel.cvut.cz/old/_media/courses/b4m35pag/lab6_slides_advanced_openmp.pdf
 * 
 */
template <typename T>
void mergeSortRecursive(std::vector<T>& v, size_t left, size_t right) {
	if (left < right) {
		if (right-left >= 100000 ) {
		unsigned long mid = (left+right)/2;
		#pragma omp taskgroup
		{
			#pragma omp task shared(v) untied if(right-left >= (1<<14))
			mergeSortRecursive(v, left, mid);
			#pragma omp task shared(v) untied if(right-left >= (1<<14))
			mergeSortRecursive(v, mid+1, right);
			#pragma omp taskyield
		}
		inplace_merge(v.begin()+left, v.begin()+mid+1, v.begin()+right+1);
		} else {
			sort(v.begin()+left, v.begin()+right+1);
		}
	}
}
template <typename T>
void mergeSortParallel(std::vector<T>& v) {
	#pragma omp parallel
	#pragma omp single
	mergeSortRecursive(v, 0, v.size()-1);
}
