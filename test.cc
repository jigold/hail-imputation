#define DOCTEST_CONFIG_IMPLEMENT
#include <vector>
#include <set>
#include <functional>
#include <cmath>
#include "doctest.h"
#include "io_plink.h"
#include "methods.h"
#include "multiarray.h"
#include "multiarray.cc"
#include "variant.h"
#include "li_stephens.h"
#include "site.h"

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
	CHECK(v4 > v1);

	std::vector<Variant> variants {v1, v3, v6, v4, v5};
    std::sort(variants.begin(), variants.end());
    CHECK(variants == std::vector<Variant> {v1, v6, v5, v4, v3});
}

TEST_CASE("site") {
	Site s1 {Variant {"1", 1, "C", "A", 0}, 2};
	Site s2 {Variant {"1", 1, "A", "T", 0}, 0};
	Site s3 {Variant {"1", 1, "C", "A", 0}, 1};

	std::vector<Site> sites {s1, s2, s3};
	std::sort(sites.begin(), sites.end());
	CHECK(sites == std::vector<Site> {s2, s3, s1});
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

TEST_CASE("zip sites") {
	PLINKReader pr1 {"data/example1"};
	PLINKReader pr3 {"data/example3"};

	auto zipped_result = zip_sites(pr1, pr3);
	auto zipped = zipped_result.zipped_sites;

	CHECK(zipped.size() == 7);
	CHECK(zipped_result.n_both == 2);
	CHECK(zipped_result.n_only1 == 3);
	CHECK(zipped_result.n_only2 == 2);
	CHECK(zipped_result.contig_boundaries == std::set<std::size_t> {0, 5});

	CHECK(*zipped[0].s1 == Site { Variant {"1", 1, "C", "A", 0}, 0 });
	CHECK(*zipped[0].s2 == Site { Variant {"1", 1, "C", "A", 0}, 0 });
	CHECK(*zipped[1].s1 == Site { Variant {"1", 2, "C", "A", 0}, 1 });
	CHECK(zipped[1].s2 == nullptr);
	CHECK(*zipped[2].s1 == Site { Variant {"1", 3, "C", "A", 0}, 2 });
	CHECK(zipped[2].s2 == nullptr);
	CHECK(*zipped[3].s1 == Site { Variant {"1", 4, "C", "A", 0}, 3 });
	CHECK(zipped[3].s2 == nullptr);
	CHECK(*zipped[4].s1 == Site { Variant {"1", 5, "C", "A", 0}, 4 });
	CHECK(*zipped[4].s2 == Site { Variant {"1", 5, "C", "A", 0}, 2 });
	CHECK(zipped[5].s1 == nullptr);
	CHECK(*zipped[5].s2 == Site { Variant {"5", 2, "C", "A", 0}, 1 });
	CHECK(zipped[6].s1 == nullptr);
    CHECK(*zipped[6].s2 == Site { Variant {"5", 4, "C", "A", 0}, 3 });
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
			ma(i, j) = idx - 1;
			CHECK(ma(i, j) == idx - 1);
			++idx;
		}
	}
	CHECK(ma.same(ma, 1e-4));
}

TEST_CASE("li_stephens") {
	PLINKReader reference {"data/example1"};
	PLINKReader sample {"data/example2"};
	LSModel ls {reference, sample};

	CHECK(ls.alpha.n_rows == 5);
    CHECK(ls.alpha.n_cols == 10);

	auto g = 0.01;
	auto theta = 0.2;
	std::vector<double> c(10, 0.4);

	auto f_likelihood = ls.forward_pass(theta, c, g);

	std::vector<double> forward_probs_expected {
	0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198,
	0.006113479, 0.605234387, 0.006113479, 0.006113479, 0.006113479, 0.605234387, 0.605234387, 0.605234387, 0.006113479, 0.006113479,
	0.003664031, 0.013474405, 0.003664031, 0.362739067, 0.362739067, 1.333966117, 1.333966117, 0.013474405, 0.003664031, 0.362739067,
	0.550639216, 0.005722653, 0.550639216, 0.011441728, 0.011441728, 2.707174489, 0.027345197, 0.566542685, 0.005562012, 0.011441728,
	0.015466719, 0.00654392, 0.015466719, 0.657119154, 0.006637567, 0.050779153, 0.682900059, 0.015727132, 0.006541289, 0.657119154
	};

	CHECK(ls.alpha.same(MultiArray<double> {forward_probs_expected, ls.alpha.n_rows, ls.alpha.n_cols}, 1e-4));
	CHECK(doctest::Approx(f_likelihood) == 2.11430086);

	auto b_likelihood = ls.backward_pass(theta, c, g);
	CHECK(doctest::Approx(b_likelihood) == 2.11430086);
}

TEST_CASE("li_stephens_uneven") {
	PLINKReader reference {"data/example1"};
	PLINKReader sample {"data/example3"};
	LSModel ls {reference, sample};

	CHECK(ls.alpha.n_rows == 2);
    CHECK(ls.alpha.n_cols == 10);

	auto g = 0.01;
	auto theta = 0.2;
	std::vector<double> c(10, 0.4);

	auto f_likelihood = ls.forward_pass(theta, c, g);

	std::vector<double> forward_probs_expected {
		0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198,
		0.010501972, 0.010501972, 0.010501972, 1.039695219, 0.010501972, 0.010501972, 1.039695219, 0.010501972, 0.010501972, 1.039695219
	};

	CHECK(ls.alpha.same(MultiArray<double> {forward_probs_expected, ls.alpha.n_rows, ls.alpha.n_cols}, 1e-4));
	CHECK(doctest::Approx(f_likelihood) == 3.19259946);

	auto b_likelihood = ls.backward_pass(theta, c, g);
	CHECK(doctest::Approx(b_likelihood) == 3.19259946);
}

TEST_CASE("li_stevens_multiple_contigs") {
	PLINKReader reference {"data/example3"};
	PLINKReader sample {"data/example4"};

	LSModel ls {reference, sample};

    CHECK(ls.alpha.n_rows == 2);
    CHECK(ls.alpha.n_cols == 1);

    auto g = 0.01;
    auto theta = 0.2;
    std::vector<double> c(1, 0.4);

    auto f_likelihood = ls.forward_pass(theta, c, g);

    std::vector<double> forward_probs_expected {
        1.98, 
		1.98
    };

	CHECK(ls.alpha.same(MultiArray<double> {forward_probs_expected, ls.alpha.n_rows, ls.alpha.n_cols}, 1e-4));
	CHECK(doctest::Approx(f_likelihood) == 1.98);

	auto b_likelihood = ls.backward_pass(theta, c, g);
	CHECK(doctest::Approx(b_likelihood) == 1.98);
}

TEST_CASE("li_stephens_identical_haplotype") {
	PLINKReader reference {"data/example5"};
	PLINKReader sample {"data/example5-chimera"};

    auto g = 0.01;
    auto theta = 0.2;
    std::vector<double> c(10, 0.4);

	LSModel ls {reference, sample, 0};
	auto f_likelihood = ls.forward_pass(theta, c, g);
	auto b_likelihood = ls.backward_pass(theta, c, g);
	CHECK(doctest::Approx(f_likelihood) == b_likelihood);
	ls.compute_gamma(f_likelihood);

	printf("%s\n", ls.gamma.to_string().c_str());
}

TEST_CASE("li_stephens_expectation_maximization") {
	PLINKReader reference {"data/example5"};
	PLINKReader sample {"data/example5-chimera"};

    auto g = 0.01;
    auto theta = 0.2;
    std::vector<double> c(10, 0.4);

	LSModel ls {reference, sample};
	EMResult emr = ls.em(theta, c, g, 0.001, 10);	
	for (auto i = 0; i < ls.n_states; ++i) {
		if (i == 0 || i == 1) {
			CHECK(emr.c[i] >= 0.4);
			CHECK(emr.c[i] <= 0.5);
		} else {
			CHECK(emr.c[i] < 0.05);
		}
	}
}