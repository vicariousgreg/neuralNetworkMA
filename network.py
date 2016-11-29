from gng import GrowingNeuralGas

class Network:
    def evaluate(self, sequences, length, indices):
        score = 0
        for i in range(length):
            counts = dict((("A", 0), ("C", 0), ("T", 0), ("G", 0)))
            for seq,start in zip(sequences, indices):
                counts[seq[start+i]] += 1
            score += max(counts.values())
        return score

class GNGNetwork:
    def __init__(self, verbose=False):
        self.gng = GrowingNeuralGas(20, feature_length=4, verbose=verbose)

    def train(self, columns, iterations):
        for _ in range(iterations):
            for column in columns:
                self.gng.feedforward(column)
        self.gng.lock()

    def evaluate(self, sequences, length, indices):
        score = 0
        for i in range(length):
            counts = dict((("A", 0.0), ("C", 0.0), ("T", 0.0), ("G", 0.0)))
            for seq,start in zip(sequences, indices):
                counts[seq[start+i]] += 1.0
            score += self.evaluate_column([x / len(sequences) for x in counts.values()])
        return score

    def evaluate_column(self, column):
        output = self.gng.feedforward(column)
        return max(output)
