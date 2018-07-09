#define DOCTEST_CONFIG_IMPLEMENT
#include <vector>
#include <functional>
#include <cmath>
#include "doctest.h"
#include "io_plink.h"
#include "methods.h"
#include "multiarray.h"
#include "multiarray.cc"
#include "variant.h"
#include "li_stephens.h"

TEST_CASE("variant") {
	Variant v1 {"1", 1, "C", "A", 0};
	Variant v2 {"1", 1, "C", "A", 0};
	Variant v3 {"2", 2, "T", "G", 0};
	Variant v4 {"1", 5, "C", "A", 0};
    Variant v5 {"1", 2, "T", "C", 0};
    Variant v6 {"1", 2, "T", "A", 0};


	CHECK(v1 == v2);
	CHECK(v1 != v3);
	CHECK(v1 < v3);
	CHECK(!(v1 < v2));

	std::vector<Variant> variants {v1, v3, v6, v4, v5};
    std::sort(variants.begin(), variants.end());
    CHECK(variants == std::vector<Variant> {v1, v6, v5, v4, v3});
}

TEST_CASE("plink loader") {
	PLINKReader pr {"data/example1"};

	std::vector<Variant> variants {
		Variant("1", 1, "C", "A", 0),
		Variant("1", 2, "C", "A", 0),
		Variant("1", 3, "C", "A", 0),
		Variant("1", 4, "C", "A", 0),
		Variant("1", 5, "C", "A", 0)
	};

	std::vector<std::string> samples {
		"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"
	};

	std::vector<int> data_expected {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 2, 0, 0, 0, 2, 2, 2, 0, 0,
		0, 0, 0, 2, 2, 2, 2, 0, 0, 2,
		2, 0, 2, 0, 0, 2, 0, 2, 0, 0,
		0, 0, 0, 2, 0, 0, 2, 0, 0, 2
	};

	std::vector<int> data;
	for (std::size_t i = 0; i < pr.n_variants; ++i) {
		for (std::size_t j = 0; j < pr.n_samples; ++j) {
			data.push_back(pr(i, j));
		}
	}

	CHECK(pr.n_variants == 5);
	CHECK(pr.n_samples == 10);
	CHECK(pr.variants == variants);
	CHECK(pr.samples == samples);
	CHECK(data == data_expected);
}


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
	CHECK(ma.same(ma, 1e-4));
}

TEST_CASE("forward_algorithm + backward_algorithm") {
	int n_states = 5;
	int t = 6;
	int n_possible_obs = 3;
	std::function<double(int const&)> start_prob = [n_states](int i) { return 1.0 / n_states; };
	std::function<double(int const&, int const&)> trans_prob = [n_states](int i, int j) { return 1.0 / n_states; };
	std::function<double(int const&, int const&)> emission_prob = [](int i, int t) { return 1.0 / 3.0; };

	MultiArray<double> probs(n_states, t);

	double likelihood = forward_algorithm(probs, start_prob, trans_prob, emission_prob);
	CHECK(doctest::Approx(likelihood) == 1.0 / pow(n_possible_obs, t));

	likelihood = backward_algorithm(probs, start_prob, trans_prob, emission_prob);
	CHECK(doctest::Approx(likelihood) == 1.0 / pow(n_possible_obs, t));
}

TEST_CASE("li_stephens") {
	PLINKReader reference {"data/example1"};
	PLINKReader sample {"data/example2"};
	LSModel ls {reference, sample};

	CHECK(ls.probs.n_rows == reference.n_samples);
    CHECK(ls.probs.n_cols == reference.n_variants);

	auto g = 0.01;
	auto theta = 0.2;
	auto c = 0.4;

	auto f_likelihood = ls.forward_algorithm(theta, c, g);
	std::vector<double> forward_probs_expected {
        0.099, 0.00152837, 0.000458004, 0.034414951, 0.000483335,
        0.099, 0.151308597, 0.001684301, 0.000357666, 0.000204497,
        0.099, 0.00152837, 0.000458004, 0.034414951, 0.000483335,
        0.099, 0.00152837, 0.045342383, 0.000715108, 0.020534974,
        0.099, 0.00152837, 0.045342383, 0.000715108, 0.000207424,
        0.099, 0.151308597, 0.166745765, 0.169198406, 0.001586849,
        0.099, 0.151308597, 0.166745765, 0.001709075, 0.021340627,
        0.099, 0.151308597, 0.001684301, 0.035408918, 0.000491473,
        0.099, 0.00152837, 0.000458004, 0.000347626, 0.000204415,
        0.099, 0.00152837, 0.045342383, 0.000715108, 0.020534974
	};

	CHECK(ls.probs.same(MultiArray<double> {forward_probs_expected, ls.probs.n_rows, ls.probs.n_cols}, 1e-4));
	CHECK(doctest::Approx(f_likelihood) == 0.066071902);

	auto b_likelihood = ls.backward_algorithm(theta, c, g);
	CHECK(doctest::Approx(b_likelihood) == 0.066071902);
}

// optimization -- transpose probs matrix
