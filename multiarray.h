#ifndef __MULTIARRAY_H__
#define __MULTIARRAY_H__

#include <vector>

template <class T>
class MultiArray {
	public:
		MultiArray(
		std::vector<T> &a,
		const std::size_t &n_rows,
		const std::size_t &n_cols): a_(a), n_rows_(n_rows), n_cols_(n_cols) {}

		std::vector<T> a_;
		std::size_t n_rows_;
		std::size_t n_cols_;

		inline const T operator() (const std::size_t &i, const std::size_t &j) {
			assert(i >= 0 && i < n_rows_ && j >= 0 && j < n_cols_);
            return a_[i * n_cols_ + j];
        }

        inline void update (const std::size_t &i, const std::size_t &j, const T &t) {
            a_[i * n_cols_ + j] = t;
        }
};

#endif

