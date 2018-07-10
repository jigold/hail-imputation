#ifndef __LI_STEPHENS_H__
#define __LI_STEPHENS_H__

#include <functional>
#include "multiarray.h"
#include "site.h"


class LSModel {
	public:
		LSModel(PLINKReader &reference, PLINKReader &sample): LSModel(reference, sample, 0) {}

		LSModel(const PLINKReader &reference,
			const PLINKReader &sample, const std::size_t &s_idx): reference(reference), sample(sample),
			s_idx(s_idx) {
				zipped_result = zip_sites(reference, sample);
				n_states = reference.n_samples;
                n_obs = zipped_result.n_both;
                probs = MultiArray<double> {n_states, n_obs};
			}

		LSModel(const PLINKReader &reference, const PLINKReader &sample,
			const std::size_t &s_idx, MultiArray<double> &probs): reference(reference), sample(sample), s_idx(s_idx), probs(probs) {
			zipped_result = zip_sites(reference, sample);
            n_states = reference.n_samples;
            n_obs = zipped_result.n_both;
            assert(probs.n_rows == n_states && probs.n_cols == n_obs);
        }

		// FIXME: add destructor

		PLINKReader reference;
		PLINKReader sample;
		std::size_t n_states;
		std::size_t n_obs;
		std::size_t s_idx;
		MultiArray<double> probs;
		ZippedResult zipped_result;

		void set_sample_idx(std::size_t &i) { s_idx = i; }

		double forward_pass(double &theta, double &c, double &g);

//		double forward_algorithm(double &theta, double &c, double &g);
//		double backward_algorithm(double &theta, double &c, double &g);

	private:
		inline double emission_prob(const std::size_t &i, const std::size_t &reference_v_idx, const std::size_t &sample_v_idx, double const&) const;
		std::size_t distance(const std::size_t &i, const std::size_t &j) const;
};

inline
double
LSModel::emission_prob(const std::size_t &i, const std::size_t &reference_v_idx, const std::size_t &sample_v_idx, const double &g) const {
    auto gt_i = reference(reference_v_idx, i);
    auto gt_sample = sample(sample_v_idx, s_idx);

    if (gt_i == gt_sample) {
        return (1 - g);
    } else {
        return g;
    }

// FIXME: This isn't correct.
//    if (gt_i == -1 || gt_sample == -1) {
//        return 1.0;
//    } else if (gt_i == gt_sample) {
//        return 2 * (1 - g);
//    } else {
//        return 2 * g;
//    }
};

inline
std::size_t
LSModel::distance(const std::size_t &i, const std::size_t &j) const {
    return reference.distance(i, j);
};

#endif