#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#include "io_plink.h"
#include "methods.h"
#include "multiarray.h"


TEST_CASE("call rate") {
	PLINKReader pr {"data/example1"};
	CHECK(call_rate(pr) == 1.0);
}

TEST_CASE("multiarray") {
	std::vector<int> a {0, 1, 2, 3, 4, 5, 6, 7};
	MultiArray<int> ma {a, 4, 2};
	int idx = 0;
	for (auto i = 0; i < 4; ++i) {
		for (auto j = 0; j < 2; ++j) {
			CHECK(ma(i, j) == idx);
			ma.update(i, j, idx - 1);
			CHECK(ma(i, j) == idx - 1);
			++idx;
		}
	}
}