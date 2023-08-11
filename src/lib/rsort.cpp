/**
 * @brief radix sort
 * @author Jannik Olbrich
 */
template<typename It, typename P>
std::vector<uint32_t> find_scores(It begin, It end, P key) {
	const auto n = std::distance(begin, end);
	std::vector<uint32_t> res(n);

	std::vector<uint32_t> st;
	auto r = res.begin();
	for (auto it = begin; it != end; it++, r++) {
		const uint32_t k = key(*it);

		const auto p = std::lower_bound(st.begin(), st.end(), k);

		*r = std::distance(st.begin(), p);

		if (p == st.end()) {
			st.emplace_back(k);
		} else if (*p > k) {
			*p = k;
		}
	}


	return res;
}

template<bool reverse = false, typename F, typename V, typename T = std::remove_cvref_t<decltype(std::declval<F>()(std::declval<V>()))>>
void rsort(V* const begin, V* const end, F f) {
	static_assert(std::is_unsigned_v<T>, "rsort: type must be unsigned");
	uint32_t count[256];
#define ibyte(x) ((f(x) >> shift) & 0xFF)

	const uint32_t n = std::distance(begin, end);
	const auto buffer = std::make_unique<V[]>(n);
	
	V* out = buffer.get();
	V* in = begin;
	for (uint32_t byte = 0; byte + (reverse ? 1u : 0u) < sizeof(T); byte++) {
		const uint32_t shift = byte * 8;
		memset(count, 0, sizeof(count));
		for (uint32_t i = 0; i < n; i++)
			count[ ibyte(in[i]) ]++;

		for (uint32_t i = 0, s = 0; i < 256; i++) {
			s += count[i];
			count[i] = s - count[i];
		}

		for (uint32_t i = 0; i < n; i++) {
			const auto p = count[ ibyte(in[i]) ]++;
			out[p] = std::move(in[i]);
		}

		std::swap(in, out);
	}

	if constexpr (reverse) {
		const uint32_t shift = (sizeof(T) - 1) * 8;
		memset(count, 0, sizeof(count));
		for (uint32_t i = 0; i < n; i++)
			count[ ibyte(in[i]) ]++;

		for (uint32_t i = 256, s = 0; i-- > 0; ) {
			s += count[i];
			count[i] = s - count[i];
		}
		
		for (uint32_t i = n; i-- > 0; ) {
			const auto p = count[ ibyte(in[i]) ]++;
			out[p] = std::move(in[i]);
		}

		std::swap(in, out);
	}

	if (in != begin)
		for (uint32_t i = 0; i < n; i++)
			begin[i] = std::move(in[i]);
#undef ibyte
}
