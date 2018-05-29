import unittest
import itertools
from .model import Model

count = 10


class TestHMM(unittest.TestCase):
    def test_forward_same_as_backward_and_exhaustive(self):
        for i in range(count):
            model, observations = Model.gen(seed=i)

            total_prob = 0
            for state_seq in itertools.product(range(model.n_states), repeat=len(observations)):
                prob = (model.start_probs[state_seq[0]] * model.emission_probs[state_seq[0], observations[0]])
                for t in range(1, len(observations)):
                    prob *= (model.trans_probs[state_seq[t - 1], state_seq[t]] * model.emission_probs[state_seq[t], observations[t]])
                total_prob += prob

            self.assertAlmostEqual(model.forward_algorithm_likelihood(observations), total_prob)
            self.assertAlmostEqual(model.backward_algorithm_likelihood(observations), total_prob)

    def test_viterbi(self):
        for i in range(count):
            model, observations = Model.gen(seed=i)

            # exhaustive search
            result = (None, float('-inf'))
            for state_seq in itertools.product(range(model.n_states), repeat=len(observations)):
                score = (model.start_probs[state_seq[0]] * model.emission_probs[state_seq[0], observations[0]])
                for t in range(1, len(observations)):
                    score *= (model.trans_probs[state_seq[t - 1], state_seq[t]] * model.emission_probs[state_seq[t], observations[t]])
                if score > result[1]:
                    result = (state_seq, score)

            v_backtrace, v_maxscore = model.viterbi(observations)

            self.assertEqual(list(result[0]), v_backtrace)
            self.assertAlmostEqual(result[1], v_maxscore)
