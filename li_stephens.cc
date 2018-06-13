#include <functional>
#include <math.h>
#include "io_plink.h"
#include "multiarray.h"
#include "li_stephens.h"

double
LSModel::forward_algorithm(double &theta, double &c, double &g) {
	// initialize
    for (auto i = 0; i < n_states; ++i) {
        probs.update(i, 0, (1.0 / n_states) * emission_prob(i, 0, g));
    }

    // recursion
    for (auto t = 1; t < n_obs; ++t) {
        auto jump_prob = 0.0;
        auto dist = distance(t, t - 1);

        for (auto i = 0; i < n_states; ++i) {
            jump_prob += probs(i, t - 1);
        }
        jump_prob = jump_prob * (1.0 - exp(-1 * theta * dist));

        for (auto i = 0; i < n_states; ++i) {
            probs.update(i, t, (probs(i, t - 1) * (exp(-1 * theta * dist)) + jump_prob * c) * emission_prob(i, t, g));
        }
    }

    // termination
    auto total = 0.0;
    for (auto i = 0; i < n_states; ++i) {
        total += probs(i, n_obs - 1);
    }
    assert(total >= 0.0 && total <= 1.0);
    return total;
}


double
LSModel::backward_algorithm(double &theta, double &c, double &g) {
	// initialize
    for (auto i = 0; i < n_states; ++i) {
        probs.update(i, n_obs - 1, emission_prob(i, n_obs - 1, g));
    }

    // recursion
    for (long t = n_obs - 2; t >= 0; --t) {
        auto jump_prob = 0.0;
        auto dist = distance(t, t + 1);

        for (auto i = 0; i < n_states; ++i) {
            jump_prob += probs(i, t + 1);
        }
        jump_prob = jump_prob * (1.0 - exp(-1 * theta * dist));

        for (auto i = 0; i < n_states; ++i) {
            probs.update(i, t, (probs(i, t + 1) * (exp(-1 * theta * dist)) + jump_prob * c) * emission_prob(i, t, g));
        }
    }

    // termination
    auto total = 0.0;
    for (auto i = 0; i < n_states; ++i) {
        total += probs(i, 0) / n_states;
    }
    assert(total >= 0.0 && total <= 1.0);
    return total;
}


// TODO
// 3. EM from reference to estimate parameters
// 4. Map from variant to idx
// 5. Handle missing variants