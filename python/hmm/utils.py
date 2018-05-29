import random


def partition(n, k):
    result = []
    cum_sum = 0
    for i in range(k - 1):
        p = random.uniform(0, n - cum_sum)
        result.append(p)
        cum_sum += p
    result.append(1 - cum_sum)
    random.shuffle(result)
    return result
