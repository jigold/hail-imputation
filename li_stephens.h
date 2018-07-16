#ifndef __LI_STEPHENS_H__
#define __LI_STEPHENS_H__

#include <vector>
#include "multiarray.h"
#include "site.h"

struct EMResult {
	double theta;
	std::vector<double> c; // should this be a pointer?
	double g;
	std::size_t n_iterations;
};

class LSModel {
	public:
		LSModel(PLINKReader &reference, PLINKReader &sample): LSModel(reference, sample, 0) {}

		LSModel(const PLINKReader &reference,
			const PLINKReader &sample, const std::size_t &s_idx): reference(reference), sample(sample),
			s_idx(s_idx) {
				zipped_result = zip_sites(reference, sample);
				n_states = reference.n_samples;
                n_obs = zipped_result.n_both;
                alpha = MultiArray<double> {n_obs, n_states};
				beta = MultiArray<double> {n_obs, n_states};
				gamma = MultiArray<double> {n_obs, n_states};
			}

		PLINKReader reference;
		PLINKReader sample;
		std::size_t n_states;
		std::size_t n_obs;
		std::size_t s_idx;
		MultiArray<double> alpha;
		MultiArray<double> beta;
		MultiArray<double> gamma;
		ZippedResult zipped_result;

		void set_sample_idx(std::size_t &i) { s_idx = i; }

		double forward_pass(double theta, std::vector<double> c, double g);
		double backward_pass(double theta, std::vector<double> c, double g);
		void compute_gamma(double p_obs);
		void normalize_row(MultiArray<double> ma, std::size_t row_idx);
		EMResult em(double theta, std::vector<double> c, double g, double tolerance, std::size_t max_iterations);

	private:
		inline double emission_prob(std::size_t i, std::size_t reference_v_idx, std::size_t sample_v_idx, double g) const;
		std::size_t distance(std::size_t i, std::size_t j) const;
};

inline
double
LSModel::emission_prob(std::size_t i, std::size_t reference_v_idx, std::size_t sample_v_idx, double g) const {
    auto gt_i = reference(reference_v_idx, i);
    auto gt_sample = sample(sample_v_idx, s_idx);

   if (gt_i == -1 || gt_sample == -1) {
       return 1.0;
   } else if (gt_i == gt_sample) {
       return 2 * (1 - g);
   } else {
       return 2 * g;
   }
};

inline
std::size_t
LSModel::distance(std::size_t i, std::size_t j) const {
    return reference.distance(i, j);
};

#endif