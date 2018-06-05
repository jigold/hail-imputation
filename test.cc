#define DOCTEST_CONFIG_IMPLEMENT
#include <vector>
#include <functional>
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

TEST_CASE("forward_algorithm") {
	int n = 5;
	std::vector<int> obs {0, 2, 1, 2, 2, 1};
	std::function<double(int const&)> start_prob = [n](int i) { return 1.0 / n; };
	std::function<double(int const&, int const&)> trans_prob = [n](int i, int j) { return 1.0 / n; };
	std::function<double(int const&, int const&)> emission_prob = [](int i, int j) { return 1.0 / 3.0; };

	double likelihood = forward_algorithm(n, obs, start_prob, trans_prob, emission_prob);
	CHECK(doctest::Approx(likelihood) == 0.001372);

	likelihood = backward_algorithm(n, obs, start_prob, trans_prob, emission_prob);
	CHECK(doctest::Approx(likelihood) == 0.001372);
}
