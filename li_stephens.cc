#include <math.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "io_plink.h"
#include "multiarray.h"
#include "li_stephens.h"

void
LSModel::normalize_row(MultiArray<double> ma, std::size_t row_idx) {
	double sum = 0.0;
	for (auto col_idx = 0; col_idx < n_states; ++col_idx) {
		sum += ma(row_idx, col_idx);
	}
	for (auto col_idx = 0; col_idx < n_states; ++col_idx) {
		ma(row_idx, col_idx) /= sum;
	}
}

double
LSModel::forward_pass(double theta, const std::vector<double> &c, double g) {
	assert(c.size() == n_states);

	std::size_t t = 0;
	std::size_t ct = 0;
	std::size_t last_ref_v_idx = 0;
	bool is_contig_boundary = true;
	auto zipped = zipped_result.zipped_sites.begin();
	auto cb = zipped_result.contig_boundaries;

	while (zipped != zipped_result.zipped_sites.end()) {
		if (zipped->s1 != nullptr && zipped->s2 != nullptr) {
			auto ref_idx = zipped->s1->idx();
			auto sample_idx = zipped->s2->idx();
			if (is_contig_boundary) {
			    for (auto i = 0; i < n_states; ++i) {
                    alpha(t, i) = emission_prob(i, ref_idx, sample_idx, g) / n_states;
                }
				is_contig_boundary = false;
			} else {
		        auto jump_prob = 0.0;
                auto dist = distance(ref_idx, last_ref_v_idx);

                for (auto i = 0; i < n_states; ++i) {
                    jump_prob += alpha(t - 1, i);
                }
                jump_prob *= (1.0 - exp(-1 * theta * dist));

                for (auto i = 0; i < n_states; ++i) {
                    alpha(t, i) = (alpha(t - 1, i) * (exp(-1 * theta * dist)) + jump_prob * c[i]) * emission_prob(i, ref_idx, sample_idx, g);
                }
			}
			++t;
			last_ref_v_idx = ref_idx;
		}
		++zipped;
		++ct;
		is_contig_boundary |= (cb.count(ct) != 0);
	}

	assert(t == n_obs);

    auto total = 0.0;
    for (auto i = 0; i < n_states; ++i) {
        total += alpha(n_obs - 1, i);
    }
    return total;
}

double
LSModel::backward_pass(double theta, const std::vector<double> &c, double g) {
	assert(c.size() == n_states);
	
	std::size_t t = n_obs - 1;
	std::size_t ct = zipped_result.zipped_sites.size() - 1;
	std::size_t last_ref_v_idx = 0;
	bool is_contig_boundary = true;

	auto zipped = zipped_result.zipped_sites.rbegin();
	auto cb = zipped_result.contig_boundaries;

	while (zipped != zipped_result.zipped_sites.rend()) {
		if (zipped->s1 != nullptr && zipped->s2 != nullptr) {
			auto ref_idx = zipped->s1->idx();
			auto sample_idx = zipped->s2->idx();

			if (is_contig_boundary) {
			    for (auto i = 0; i < n_states; ++i) {
					beta(t, i) = emission_prob(i, ref_idx, sample_idx, g);
                }
				is_contig_boundary = false;
			} else {
		        auto jump_prob = 0.0;
                auto dist = distance(last_ref_v_idx, ref_idx);

                for (auto i = 0; i < n_states; ++i) {
                    jump_prob += beta(t + 1, i);
                }
                jump_prob *= (1.0 - exp(-1 * theta * dist));

                for (auto i = 0; i < n_states; ++i) {
					beta(t, i) = (beta(t + 1, i) * (exp(-1 * theta * dist)) + jump_prob * c[i]) * emission_prob(i, ref_idx, sample_idx, g);
                }
			}
			--t;
			last_ref_v_idx = ref_idx;
		}
		++zipped;
		is_contig_boundary |= (cb.count(ct) != 0);
		--ct;
	}

	assert(t == -1);

	auto total = 0.0;
	for (auto i = 0; i < n_states; ++i) {
		total += beta(0, i) / n_states;
	}
	return total;
}

void
LSModel::compute_gamma(double p_obs) {
	for (auto t = 0; t < n_obs; ++t) {
		for (auto i = 0; i < n_states; ++i) {
			gamma(t, i) = alpha(t, i) * beta(t, i) / p_obs;
		}
	}
}

EMResult
LSModel::em(double theta, std::vector<double> c, double g, double tolerance, std::size_t max_iterations) {
	std::vector<double> c_new(n_states);
	std::size_t	iter_no = 0;
	while (iter_no < max_iterations) {
		forward_pass(theta, c, g);

		double row_sum = 0.0;
		for (auto t = 0; t < n_obs; ++t) {
			for (auto i = 0; i < n_states; ++i) {
				row_sum += alpha(t, i);
			}	
			for (auto i = 0; i < n_states; ++i) {
				c_new[i] += alpha(t, i) / row_sum;
			}		
		}

		double c_new_sum = 0.0;
		for (auto i = 0; i < n_states; ++i) {
			c_new_sum += c_new[i];
		}

		bool within_tolerance = true;
		printf("iteration #%zu\n", iter_no);
		for (auto i = 0; i < n_states; ++i) {
			c_new[i] /= c_new_sum;
			printf("idx %i: (%f, %f)\n", i, c[i], c_new[i]);
			if (std::abs(c_new[i] - c[i]) > tolerance) {
				within_tolerance = false;
			}
		}
		printf("\n");

		c = c_new;

		if (within_tolerance) {
			break;
		}

		++iter_no;
	}

	EMResult emr {theta, c, g, iter_no};

	return emr;
}