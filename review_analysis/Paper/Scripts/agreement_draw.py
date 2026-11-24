"""Reproduce the paper's uniform within-person sample with R's seed and indexing."""

import numpy as np


class RSample:
    """R's Mersenne Twister initialisation and integer rejection sampler."""

    def __init__(self, seed):
        word = seed
        for _ in range(50):
            word = (69069 * word + 1) & 0xFFFFFFFF
        words = []
        for _ in range(625):
            word = (69069 * word + 1) & 0xFFFFFFFF
            words.append(word)
        self.rng = np.random.RandomState()
        self.rng.set_state(('MT19937', np.array(words[1:], dtype=np.uint32), 624, 0, 0))

    def draw(self, n):
        bits = (int(n) - 1).bit_length()
        while True:
            word = 0
            for _ in range(0, bits + 1, 16):
                high = int(self.rng.randint(0, 2**32, dtype=np.uint32)) >> 16
                word = (word << 16) + high
            index = word & ((1 << bits) - 1)
            if index < n:
                return index


def select_one(frame):
    rng = RSample(52001)
    groups = frame.groupby('patient_id', sort=True).indices.values()
    positions = [rows[rng.draw(len(rows))] for rows in groups]
    selected = frame.iloc[positions].copy()
    if (
        len(selected) != selected.patient_id.nunique()
        or len(selected) != frame.patient_id.nunique()
    ):
        raise ValueError('Each eligible person must contribute exactly once.')
    return selected
