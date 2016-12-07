from gng import GrowingNeuralGas

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

    def train(self, columns, iterations):
        for _ in range(iterations):
            for column in columns:
                self.gng.feedforward(column)
        self.gng.lock()

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
        return max(output)
