from gng import GrowingNeuralGas
from random import shuffle

class Network:
    def __init__(self, alphabet):
        self.alphabet = alphabet

    def evaluate(self, sequences, length, indices):
        score = 0
        for i in range(length):
            counts = dict(zip(self.alphabet, [0] * len(self.alphabet)))
            for seq,start in zip(sequences, indices):
                counts[seq[start+i]] += 1
            score += max(counts.values())
        return score

class GNGNetwork:
    def __init__(self, alphabet, size=20, verbose=False):
        self.gng = GrowingNeuralGas(size, feature_length=len(alphabet), verbose=verbose)
        self.alphabet = alphabet
        self.gng.lock(True)

    def train(self, columns, iterations, verbose=False):
        if verbose: print("Training GNG network (iterations: %d)..." % iterations)
        self.gng.lock(False)
        indices = list(range(len(columns)))
        for _ in range(iterations):
            shuffle(indices)
            if verbose: print("Iter %4d" % _)
            for i in indices:
                self.gng.feedforward(columns[i])
        self.gng.lock(True)

    def evaluate(self, sequences, length, indices):
        score = 0
        for i in range(length):
            counts = dict(zip(self.alphabet, [0] * len(self.alphabet)))
            for seq,start in zip(sequences, indices):
                counts[seq[start+i]] += 1.0
            score += self.evaluate_column([x / len(sequences) for x in counts.values()])
        return score

    def evaluate_column(self, column):
        output = self.gng.feedforward(column)
        active = self.gng.active_neurons
        return max([o for o,a in zip(output,active) if a])
