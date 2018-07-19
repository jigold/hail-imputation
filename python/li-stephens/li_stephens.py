import numpy as np
import math

class LSModel(object):
    def __init__(self, reference, sample, theta, c, g, sample_idx=0):
        """
        Assumptions:
        1. reference and sample are same length (same number of variants)
        2. length of c is the same length as reference and sample
        3. distance between variants is 1
        """
        self._reference = reference # numpy matrix
        self._sample = sample # numpy array
        self._theta = theta
        self._c = c
        self._g = g
        self._sample_idx = sample_idx

        assert(self._c.size == self._reference.shape[1])
        assert(self._reference.shape[0] == self._sample.size)

        self._n_states = self._reference.shape[1]
        self._n_obs = self._reference.shape[0]

        self._alpha = np.zeros(self._reference.shape)
        self._beta = np.zeros(self._reference.shape)
        self._gamma = np.zeros(self._reference.shape)

    def emission(self, ref_v_idx, sample_v_idx, state):
        g_ref = self._reference[ref_v_idx, state]
        g_sample = self._sample[sample_v_idx, self._sample_idx]

        if g_ref == None or g_sample == None:
            return 1
        elif g_ref == g_sample:
            return 2 * (1 - self._g)
        else:
            return 2 * self._g

    def transition_matrix(self):
        transition_probs = np.zeros((self._n_states, self._n_states))
        for i in range(self._n_states):
            for j in range(self._n_states):
                if i == j:
                    transition_probs[i, j] = math.exp(-self._theta) + self._c[i] * (1 - math.exp(-self._theta))
                else:
                    transition_probs[i, j] = self._c[i] * (1 - math.exp(-self._theta))
        
        return transition_probs

    def forward_pass(self):
        jump_prob = 1 - math.exp(-1 * self._theta * 1) # hardcoded distance as 1 for balding nichols generated examples

        # initialize
        emission = np.array([self.emission(0, 0, i) for i in range(self._n_states)])
        self._alpha[0, :] = emission / self._n_states

        # recursion
        for t in range(1, self._n_obs):
            jump = jump_prob * self._alpha[t - 1, :].sum()
            emission = np.array([self.emission(t, t, i) for i in range(self._n_states)])
            self._alpha[t, :] = emission * (self._alpha[t - 1, :] * (1 - jump_prob) + jump * self._c)
        
        # terminate
        return self._alpha[self._n_obs - 1, :].sum()

    def backward_pass(self):
        jump_prob = 1 - math.exp(-1 * self._theta * 1) # hardcoded distance as 1 for balding nichols generated examples

        # initialize
        self._beta[self._n_obs - 1, :] = 1

        # recursion
        for t in reversed(range(self._n_obs - 1)):
            jump = jump_prob * self._beta[t + 1, :].sum()
            emission = np.array([self.emission(t + 1, t + 1, i) for i in range(self._n_states)])
            self._beta[t, :] = emission * (self._beta[t + 1, :] * (1 - jump_prob) + jump * self._c)

        # terminate
        emission = np.array([self.emission(0, 0, i) for i in range(self._n_states)])
        return (self._beta[0, :] * emission / self._n_states).sum()

    
class HMM(object):
    def __init__(self, observations, transition, emission, initial):
        self._observations = observations # numpy array
        self._transition = transition
        self._emission = emission
        self._initial = initial

        self._n_states = self._transition.shape[0]
        self._n_obs = self._observations.shape[0]

        assert(self._initial.size == self._n_states)

        shape = (self._n_obs, self._n_states)
        self._alpha = np.zeros(shape)
        self._beta = np.zeros(shape)
        self._gamma = np.zeros(shape)

    def normalize(self, a, row_idx):
        a[row_idx, :] = a[row_idx, :] / a[row_idx, :].sum()

    def forward_pass(self):
        # initialize
        e = np.array([self._emission[i, self._observations[0]] for i in range(self._n_states)])
        self._alpha[0, :] = self._initial * e
        self.normalize(self._alpha, 0)

        # recursion
        for t in range(1, self._n_obs):
            for j in range(self._n_states):
                total = 0.0
                for i in range(self._n_states):
                    total += self._alpha[t - 1, i] * self._transition[i, j]
                self._alpha[t, j] = total * self._emission[j, self._observations[t]]
            self.normalize(self._alpha, t)
        
        # terminate
        observation_likelihood = 0.0
        for i in range(self._n_states):
            observation_likelihood += self._alpha[self._n_obs - 1, i]

        print(self._alpha)
        print(observation_likelihood)
        return observation_likelihood
