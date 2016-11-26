class Network:
    def evaluate(self, sequences, length, indices):
        score = 0
        for i in xrange(length):
            counts = dict((("A", 0), ("C", 0), ("T", 0), ("G", 0)))
            for seq,start in zip(sequences, indices):
                counts[seq[start+i]] += 1
            score += max(counts.values())
        return score
