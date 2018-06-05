#ifndef __METHODS_H__
#define __METHODS_H__

#include <functional>

const double call_rate(PLINKReader &pr);

const double
forward_algorithm(const std::size_t &n_states,
			const std::vector<int> &observations,
			std::function<double(int)> start_prob,
			std::function<double(int, int)> transition_prob,
			std::function<double(int, int)> emission_prob);

const double
backward_algorithm(const std::size_t &n_states,
			const std::vector<int> &observations,
			std::function<double(int)> start_prob,
			std::function<double(int, int)> transition_prob,
			std::function<double(int, int)> emission_prob);
#endif