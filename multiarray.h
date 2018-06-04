#ifndef __MULTIARRAY_H__
#define __MULTIARRAY_H__

#include <vector>

template <class T>
class MultiArray {
	public:
		MultiArray(
		const std::vector<T> &a,
		const std::size_t &n_rows,
		const std::size_t &n_cols): a(a), n_rows(n_rows), n_cols(n_cols) {}

		std::vector<T> a;
		std::size_t n_rows;
		std::size_t n_cols;

		inline const T operator() (const std::size_t &i, const std::size_t &j) {
			assert(i >= 0 && i < n_rows && j >= 0 && j < n_cols);
            return a[i * n_cols + j];
        };
	private:
};

#endif

