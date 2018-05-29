class MultiArray2(object):
    def __init__(self, n_rows, n_cols, a):
        assert (n_rows >= 0 and n_cols >= 0 and len(a) == n_rows * n_cols)
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.a = a

    def __getitem__(self, indices):
        if isinstance(indices, slice):
            raise TypeError("don't support slices")
        return self.apply(indices[0], indices[1])

    def update(self, i, j, t):
        assert (0 <= i < self.n_rows)
        assert (0 <= j < self.n_cols)
        self.a[i * self.n_cols + j] = t

    def apply(self, i, j):
        assert (0 <= i < self.n_rows)
        assert (0 <= j < self.n_cols)
        return self.a[i * self.n_cols + j]

    def array(self):
        return self.a
