#include <math.h>
#include <iostream>
#include "io_plink.h"
#include "multiarray.h"
#include "li_stephens.h"

double
LSModel::forward_pass(const double &theta, const double &c, const double &g) {
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
                    alpha(i, t) = emission_prob(i, ref_idx, sample_idx, g) / n_states;
                }
				is_contig_boundary = false;
			} else {
		        auto jump_prob = 0.0;
                auto dist = distance(ref_idx, last_ref_v_idx);

                for (auto i = 0; i < n_states; ++i) {
                    jump_prob += alpha(i, t - 1);
                }
                jump_prob *= (1.0 - exp(-1 * theta * dist));

                for (auto i = 0; i < n_states; ++i) {
                    alpha(i, t) = (alpha(i, t - 1) * (exp(-1 * theta * dist)) + jump_prob * c) * emission_prob(i, ref_idx, sample_idx, g);
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
        total += alpha(i, n_obs - 1);
    }
    return total;
}

double
LSModel::backward_pass(const double &theta, const double &c, const double &g) {
	std::size_t t = n_obs - 1;
	std::size_t ct = zipped_result.zipped_sites.size() - 1;
	std::size_t last_ref_v_idx = 0;
	bool is_contig_boundary = true;

	auto zipped = zipped_result.zipped_sites.rbegin();
	auto cb = zipped_result.contig_boundaries;

	for(auto iter=cb.begin(); iter!=cb.end();++iter) {
		std::cout << *iter << std::endl;
	}

	while (zipped != zipped_result.zipped_sites.rend()) {
		if (zipped->s1 != nullptr && zipped->s2 != nullptr) {
			auto ref_idx = zipped->s1->idx();
			auto sample_idx = zipped->s2->idx();

			if (is_contig_boundary) {
			    for (auto i = 0; i < n_states; ++i) {
					beta(i, t) = emission_prob(i, ref_idx, sample_idx, g);
                }
				is_contig_boundary = false;
			} else {
		        auto jump_prob = 0.0;
                auto dist = distance(last_ref_v_idx, ref_idx);

                for (auto i = 0; i < n_states; ++i) {
                    jump_prob += beta(i, t + 1);
                }
                jump_prob *= (1.0 - exp(-1 * theta * dist));

                for (auto i = 0; i < n_states; ++i) {
					beta(i, t) = (beta(i, t + 1) * (exp(-1 * theta * dist)) + jump_prob * c) * emission_prob(i, ref_idx, sample_idx, g);
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
		total += beta(i, 0) / n_states;
	}
	return total;
}

void
LSModel::compute_gamma(const double &p_obs) {
	for (auto i = 0; i < n_states; ++i) {
		for (auto t = 0; t < n_obs; ++t) {
			gamma(i, t) = alpha(i, t) * beta(i, t) / p_obs;
		}
	}
}
