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

    
