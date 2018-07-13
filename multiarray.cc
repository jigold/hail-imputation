#include <vector>
#include <cmath>
#include <string>
#include "multiarray.h"

template <class T>
bool MultiArray<T>::same (const MultiArray<T> &other, double tolerance) const {
	assert(n_rows == other.n_rows && n_cols == other.n_cols);
	bool same = true;
	for (auto i = 0; i < n_rows; ++i) {
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

template <class T>
std::string MultiArray<T>::to_string() const {
	std::string str;
	for (auto i = 0; i < n_rows; ++i) {		
		for (auto j = 0; j < n_cols; ++j) {
			if (j == 0) {
				str += std::to_string(operator()(i, j));
			} else {
				str += " " + std::to_string(operator()(i, j));
			}
		}
		str += "\n";
	}

	return str;
}
