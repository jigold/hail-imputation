#include <vector>
#include <functional>
#include "io_plink.h"
#include "multiarray.h"

const double
call_rate(PLINKReader &pr) {
	long total = 0;
	long n_missing = 0;
	for (auto i = 0; i < pr.n_variants; ++i) {
		for (auto j = 0; j < pr.n_samples; ++j) {
			int gt = pr(i, j);
			if (gt == 1) {
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
forward_algorithm(const std::size_t &n_states,
			 const std::vector<int> &observations,
			 std::function<double(int)> start_prob,
			 std::function<double(int, int)> transition_prob,
			 std::function<double(int, int)> emission_prob) {
			    const std::size_t T = observations.size();
			    MultiArray<double> probs(n_states, T);

				// initialize
			    for (auto i = 0; i < n_states; ++i) {
			        probs.update(i, 0, start_prob(i) * emission_prob(i, observations[0]));
			    }

			    // recursion
			    for (auto t = 1; t < T; ++t) {
			        auto obs = observations[t];
			        for (auto j = 0; j < n_states; ++j) {
			            for (auto i = 0; i < n_states; ++i) {
			                probs.update(j, t, probs(j, t) + probs(i, t - 1) * transition_prob(i, j) * emission_prob(j, obs));
			            }
			        }
			    }

			    // termination
			    double total = 0.0;
			    for (auto i = 0; i < n_states; ++i) {
			        total += probs(i, T - 1);
			    }
			    return total;
			 }
