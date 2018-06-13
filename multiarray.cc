#include <vector>
#include <cmath>
#include "multiarray.h"

template <class T>
const bool MultiArray<T>::same (const MultiArray<T> &other, const double &tolerance) {
	assert(n_rows == other.n_rows && n_cols == other.n_cols);
	bool same = true;
	for (auto i = 0; i < n_rows; ++i){
		for (auto j = 0; j < n_cols; ++j) {
			double x = operator()(i, j);
			double y = other(i, j);
			bool tmp = std::abs(x - y) <= tolerance;
			if (!tmp) {
				printf("(%i, %i): (%f, %f)\n", i, j, x, y);
			}
			same &= tmp;
		}
    }
    return same;
}