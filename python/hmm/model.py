import random
from .multiarray import MultiArray2
from .utils import partition


class Model(object):
    def __init__(self, start_probs, trans_probs, emission_probs):
        assert(isinstance(trans_probs, MultiArray2))
        assert(isinstance(emission_probs, MultiArray2))
        assert(isinstance(start_probs, list))

        assert(trans_probs.n_rows == trans_probs.n_cols)
        assert(trans_probs.n_rows == emission_probs.n_rows)
        assert(trans_probs.n_rows == len(start_probs))

        self.start_probs = start_probs
        self.trans_probs = trans_probs
        self.emission_probs = emission_probs

        self.n_states = trans_probs.n_rows

    @classmethod
    def gen(cls, seed=None):
        if seed is not None:
            random.seed(seed)

        n_states = random.choice([1, 2, 4, 8])
        n_obs = random.choice([1, 2, 4, 8])
        t = random.choice([1, 2, 4, 5])

        start_probs = partition(1.0, n_states)

        trans_probs = []
        for i in range(n_states):
            trans_probs.extend(partition(1.0, n_states))
        trans_probs = MultiArray2(n_states, n_states, trans_probs)

        emission_probs = []
        for i in range(n_states):
            emission_probs.extend(partition(1.0, n_obs))
        emission_probs = MultiArray2(n_states, n_obs, emission_probs)

        observations = []
        for i in range(t):
            observations.append(random.choice(list(range(n_obs))))

        return Model(start_probs, trans_probs, emission_probs), observations

    def forward_algorithm(self, observations):
        assert(isinstance(observations, list))
        n_obs = len(observations)
        probs = MultiArray2(self.n_states, n_obs, [0] * self.n_states * n_obs)

        # initialize
        obs_0 = observations[0]
        for i in range(self.n_states):
            probs.update(i, 0, self.start_probs[i] * self.emission_probs[i, obs_0])

        # recursion
        for t in range(1, n_obs):
            obs = observations[t]
            for j in range(self.n_states):  # current state
                for i in range(self.n_states):  # prev state
                    probs.update(j, t, probs[j, t] + probs[i, t - 1] * self.trans_probs[i, j] * self.emission_probs[j, obs])

        return probs

    def forward_algorithm_likelihood(self, observations, probs=None):
        if probs is None:
            probs = self.forward_algorithm(observations)

        n_obs = len(observations)
        final_prob = 0
        for i in range(self.n_states):
            final_prob += probs[i, n_obs - 1]

        return final_prob

    def backward_algorithm(self, observations):
        assert(isinstance(observations, list))

        n_states = self.n_states
        n_obs = len(observations)

        probs = MultiArray2(n_states, n_obs, [0] * n_states * n_obs)

        # initialize
        obs_last = observations[n_obs - 1]
        for i in range(self.n_states):
            probs.update(i, n_obs - 1, self.emission_probs[i, obs_last])

        # recursion
        for t in reversed(range(n_obs - 1)):
            obs = observations[t]
            for i in range(n_states):  # current state
                for j in range(n_states):  # prev state
                    probs.update(i, t, probs[i, t] + probs[j, t + 1] * self.trans_probs[i, j] * self.emission_probs[i, obs])

        return probs

    def backward_algorithm_likelihood(self, observations, probs=None):
        if probs is None:
            probs = self.backward_algorithm(observations)

        final_prob = 0
        for i in range(self.n_states):
            final_prob += probs[i, 0] * self.start_probs[i]

        return final_prob

    def viterbi(self, observations):
        assert(isinstance(observations, list))

        n_states = self.n_states
        n_obs = len(observations)

        probs = MultiArray2(n_states, n_obs, [0] * n_states * n_obs)
        backpointers = MultiArray2(n_states, n_obs, [0] * n_states * n_obs)

        # initialize
        obs_0 = observations[0]
        for i in range(n_states):
            probs.update(i, 0, self.start_probs[i] * self.emission_probs[i, obs_0])

        # recursion
        for t in range(1, n_obs):
            obs = observations[t]
            for j in range(n_states):  # current state
                max_score = None
                bp_score = None
                bp = None
                for i in range(n_states):  # prev state
                    score = probs[i, t - 1] * self.trans_probs[i, j] * self.emission_probs[j, obs]
                    score2 = probs[i, t - 1] * self.trans_probs[i, j] * self.emission_probs[j, obs]
                    if max_score is None or score > max_score:
                        max_score = score
                    if (bp_score is None and bp is None) or score2 > bp_score:
                        bp = i
                        bp_score = score2
                probs.update(j, t, max_score)
                backpointers.update(j, t, bp)

        # termination
        max_score = None
        argmax = None
        for i in range(n_states):
            score = probs[i, n_obs - 1]
            if (max_score is None and argmax is None) or (score > max_score):
                max_score = score
                argmax = i

        backtrace = [0] * n_obs
        state = argmax
        backtrace[n_obs - 1] = state
        for t in reversed(range(1, n_obs)):
            state = backpointers[state, t]
            backtrace[t - 1] = state

        return backtrace, max_score

    def baum_welch(self, observations):
        forward = self.forward_algorithm(observations)
        backward = self.backward_algorithm(observations)
        prob_observations = self.forward_algorithm_likelihood(observations,
                                                              forward)

        estimated_transitions = self.trans_probs
        n_obs = len(observations)

        def prob_state(t, i, j):
            return (forward[i, t] * estimated_transitions[i, j] *
                    self.emission_probs[j, observations[t + 1]] *
                    backward[j, t + 1]) / prob_observations

        num = MultiArray2(self.n_states, self.n_states, [0] * self.n_states * self.n_states)
        denom = [0] * self.n_states

        for i in range(self.n_states):
            for j in range(self.n_states):
                for t in range(n_obs - 1):
                    num.update(i, j, num[i, j] + prob_state(t, i, j))
                    for k in range(self.n_states):
                        denom[i] += prob_state(t, i, k)






