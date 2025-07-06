import numpy as np

def reshape_cov(raw, n):
    arr = raw.flatten()
    if arr.size == n * n:
        return arr.reshape((n, n))
    if arr.size == n * n + 1:
        return arr[:-1].reshape((n, n))
    if arr.size == n * n + n:
        return arr[:-n].reshape((n, n))
    raise ValueError(f"Cannot reshape array of size {arr.size} into {n}×{n}")
