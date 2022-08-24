#pragma once

template<typename It>
struct it_range {
	using T = std::remove_cvref_t<decltype(*std::declval<It>())>;
	It b, e;
	it_range() = default;
	it_range(const It& b, const It& e) : b(b), e(e) {}
	const It& begin() const { return b; }
	const It& end() const { return e; }
	size_t size() const { return std::distance(begin(), end()); }
	T& operator[](size_t i) { return *(begin() + i); }
	const T& operator[](size_t i) const { return *(begin() + i); }
	std::vector<T> to_vector() const {
		return std::vector<T>(begin(), end());
	}
};

template<typename It>
std::ostream& operator<<(std::ostream& out, const it_range<It>& r) {
	out << "r<";
	size_t n = 0;
	for (const auto& x : r)
		out << (n++ > 0 ? ", " : "") << x;
	return out << ">";
}

template<typename A, typename B>
std::ostream& operator<<(std::ostream& out, const std::pair<A,B>& rhs) {
	return out << "<<" << rhs.first << ", " << rhs.second << ">>";
}

template<typename A, typename B>
std::ostream& operator<<(std::ostream& out, const std::vector<std::pair<A,B>>& rhs) {
	return out << it_range(rhs.begin(), rhs.end());
}

template<typename...T>
std::ostream& operator<<(std::ostream& out, const std::tuple<T...>& rhs) {
	std::apply([&out](T const&... t) {
		size_t n{0};
		((out << (n++ > 0 ? ", " : "<<") << t), ...);
		out << ">>";
	}, rhs);
	return out;
}
