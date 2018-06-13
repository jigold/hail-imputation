#include <vector>
#include <functional>
#include <math.h>
#include "io_plink.h"
#include "multiarray.h"

const double
call_rate(PLINKReader &pr) {
	long total = 0;
	long n_missing = 0;
	for (auto i = 0; i < pr.n_variants; ++i) {
		for (auto j = 0; j < pr.n_samples; ++j) {
			auto gt = pr(i, j);
			if (gt == -1) {
				++n_missing;
			}
			++total;
		}
	}
	if (total != 0)
		return (double) (total - n_missing) / total;
	else {
        printf("divide by zero.\n");
        exit(1);
	}
}

const double
forward_algorithm(
	MultiArray<double> &probs,
	std::function<double(int)> start_prob,
	std::function<double(int, int)> transition_prob,
	std::function<double(int, int)> emission_prob) {

    auto n_states = probs.n_rows;
    auto n_obs = probs.n_cols;

	// initialize
    for (auto i = 0; i < n_states; ++i) {
        probs.update(i, 0, start_prob(i) * emission_prob(i, 0));
    }

    // recursion
    for (auto t = 1; t < n_obs; ++t) {
        for (auto j = 0; j < n_states; ++j) {
	        probs.update(j, t, 0);
            for (auto i = 0; i < n_states; ++i) {
                probs.update(j, t, probs(j, t) + probs(i, t - 1) * transition_prob(i, j) * emission_prob(j, t));
            }
        }
    }

    // termination
    double total = 0.0;
    for (auto i = 0; i < n_states; ++i) {
        total += probs(i, n_obs - 1);
    }

    return total;
}

const double
backward_algorithm(
	MultiArray<double> &probs,
	std::function<double(int)> start_prob,
	std::function<double(int, int)> transition_prob,
	std::function<double(int, int)> emission_prob) {

    auto n_states = probs.n_rows;
    auto n_obs = probs.n_cols;

	// initialize
	for (auto i = 0; i < n_states; ++i) {
	    probs.update(i, n_obs - 1, emission_prob(i, n_obs - 1));
	}

	// recursion
	for (long t = n_obs - 2; t >= 0; --t) {
	    for (auto i = 0; i < n_states; ++i) {
	        probs.update(i, t, 0);
	        for (auto j = 0; j < n_states; ++j) {
	            probs.update(i, t, probs(i, t) + probs(j, t + 1) * transition_prob(i, j) * emission_prob(i, t));
	        }
	    }
	}

	// termination
	double total = 0.0;
	for (auto i = 0; i < n_states; ++i) {
	    total += probs(i, 0) * start_prob(i);
	}

	return total;
}
