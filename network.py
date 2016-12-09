from gng import GrowingNeuralGas
from random import shuffle
from alphabet import handle_special_letters
import pickle

def gen_columns(alphabet, sequences, length, indices):
    for i in range(length):
        counts = dict(zip(alphabet, [0] * len(alphabet)))
        for seq,start in zip(sequences, indices):
            try:
                counts[seq[start+i]] += 1
            except KeyError:
                handle_special_letters(alphabet, counts, seq[start+i])
        yield tuple(counts[x] / len(sequences) for x in alphabet)

class Network:
    def __init__(self, alphabet):
        self.alphabet = alphabet

    def evaluate(self, sequences, length, indices):
        return sum(max(column)
            for column in gen_columns(self.alphabet, sequences, length, indices))

class GNGNetwork:
    def __init__(self, alphabet, size=20, verbose=False):
        self.gng = GrowingNeuralGas(size, feature_length=len(alphabet), verbose=verbose)
        self.alphabet = alphabet
        self.run_iterations = 0
        self.gng.lock(True)

    def save(self, filename):
        pickle.dump(self, open( filename, "wb" ) )

    @staticmethod
    def load(filename):
        return pickle.load( open( filename, "rb" ) )

    def train(self, columns, iterations, fraction=1.0, verbose=False):
        if verbose: print("Training GNG network (iterations: %d)..." % iterations)
        self.gng.lock(False)
        indices = list(range(len(columns)))
        for _ in range(iterations):
            self.run_iterations += 1
            shuffle(indices)
            if verbose: print("Iter %4d" % _)
            for i in indices[:int(fraction * len(indices))]:
                self.gng.feedforward(columns[i])
        self.gng.lock(True)

    def evaluate(self, sequences, length, indices):
        return sum(self.evaluate_column(column)
            for column in gen_columns(self.alphabet, sequences, length, indices))

    def evaluate_column(self, column, verbose=False):
        output = self.gng.feedforward(column)
        active = self.gng.active_neurons
        output = [o for o,a in zip(output,active) if a]
        result = max(output) ** 2
        if verbose:
            print(" ".join("%.4f" % x if x > 0.0 else "      " for x in column))
            print(" ".join("%.4f" % x for x in output))
            #print("Max output: %.4f" % result)
            print("")
        return result
