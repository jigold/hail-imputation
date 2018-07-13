#ifndef __MULTIARRAY_H__
#define __MULTIARRAY_H__

#include <vector>

template <class T>
class MultiArray {
	public:
		MultiArray(
		std::vector<T> &a_,
		std::size_t n_rows_,
		std::size_t n_cols_): a(a_), n_rows(n_rows_), n_cols(n_cols_) {}

		MultiArray(std::size_t n_rows_, std::size_t n_cols_): a(), n_rows(n_rows_), n_cols(n_cols_) {
			a.resize(n_rows * n_cols);
		}

		MultiArray() = default;

		std::vector<T> a;
		std::size_t n_rows;
		std::size_t n_cols;

		inline T& operator() (std::size_t i, std::size_t j) {
			assert(i >= 0 && i < n_rows && j >= 0 && j < n_cols);
            return a[i * n_cols + j];
        }
		
		inline T operator() (const std::size_t i, const std::size_t j) const {
			assert(i >= 0 && i < n_rows && j >= 0 && j < n_cols);
            return a[i * n_cols + j];
        }

        inline void resize (const std::size_t &i, const std::size_t &j) {
            a.resize(i * j);
        }

        bool same (const MultiArray<T> &other, double tolerance) const;

		std::string to_string() const;
};

#endif

