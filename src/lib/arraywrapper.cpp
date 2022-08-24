#pragma once

// Just because
// src/lib/merge.cpp:39:57: sorry, unimplemented: capture of variably-modified type ‘uint32_t [n][(((long unsigned int) n) + 1)]’ {aka ‘unsigned int [n][(((long unsigned int)n) + 1)]’} that is not an N3639 array of runtime bound
//    39 |  const std::function<void(uint32_t, uint32_t)> mmerge = [&](uint32_t i, uint32_t e) {
template<typename T, size_t D>
struct ArrayWrapper {
public:
	static_assert(D > 0, "Cannot declare array with 0 dimensions");
	template<typename ... dims_t>
	ArrayWrapper(dims_t ... dims)
			: sizes{dims...}
			, data(compute_strides({dims...}))
			{
		static_assert(sizeof...(dims) == D, "Wrong number of dimensions in constructor");
	}
	template<size_t d>
	constexpr const size_t& size() {
		static_assert(d < D, "Invalid dimension index in size()");
		return std::get<d>(sizes);
	};
	template<typename ... idx>
	T& operator()(idx ... i) {
		static_assert(sizeof...(i) == D, "Wrong number of indices in element access");
		return get({i...});
	};
	template<typename ... idx>
	const T& operator()(idx ... i) const {
		static_assert(sizeof...(i) == D, "Wrong number of indices in element access");
		return get({i...});
	};
protected:
	std::array<const size_t, D> sizes;
	std::array<size_t, D-1> strides;
	std::vector<T> data;
	T& get(const std::array<size_t, D>& ix) const {
		const size_t i = std::inner_product(strides.begin(), strides.end(), ix.begin(), ix.back());
		if (i >= data.size()) {
			std::cerr << "i = " << i << " " << data.size() << std::endl;
			exit(1);
		}
		return const_cast<T&>( data[ i ] );
	}
	size_t compute_strides(const std::array<size_t, D>& dims) {
		if constexpr (D == 1) {
			return dims[0];
		}
		strides[D-2] = dims[D-1];
		for (ptrdiff_t i = static_cast<ptrdiff_t>(D) - 3; i >= 0; --i) {
			strides[i] = strides[i+1] * dims[i+1];
		}
		return dims[0] * strides[0];
	}
};
