import random

class Entry:
    def __init__(self, sequences, seq_ids, indices, length):
        if len(seq_ids) != len(indices):
            raise ValueError
        self.seq_ids = seq_ids
        self.indices = indices
        self.length = length

        self.columns = []
        for offset in range(length):
            counts = dict((("A", 0.0), ("C", 0.0), ("T", 0.0), ("G", 0.0)))
            for seq_id,start in zip(seq_ids, indices):
                counts[sequences[seq_id][start+offset]] += 1.0
            self.columns.append([x / len(sequences) for x in counts.values()])

class Dataset:
    def __init__(self):
        self.sequences = []
        self.entries = []

    def add_sequence(self, sequence):
        self.sequences.append(sequence)
        return len(self.sequences)-1

    def add_alignment(self, seq_ids, indices, length):
        # Check sequence IDs
        if any(seq_id >= len(self.sequences) for seq_id in seq_ids):
            raise ValueError
        # Check indices
        if any(i + length >= len(self.sequences[seq_id])
               for seq_id, i in zip(seq_ids, indices)):
            raise ValueError
        self.entries.append(Entry(self.sequences, seq_ids, indices, length))

    def get_columns(self):
        return [entry.columns[i]
                for entry in self.entries
                for i in range(len(entry.columns))]

def gen_sequences(num_sequences, length, seed_length=0, mut_rate=None):
    print("Generating %d sequences of length %d..." % (num_sequences, length))
    seqs = [[random.choice("ATCG") for _ in range(length)]
            for _ in range(num_sequences)]

    indices = []
    if seed_length > 0:
        seed = [random.choice("ATCG") for _ in range(seed_length)]
        print("Planting seed... (%s)" % "".join(seed))
        for seq in seqs:
            start = random.randint(0, length-seed_length-1)
            indices.append(start)
            for i in range(seed_length):
                if mut_rate and random.random() < mut_rate:
                    seq[start+i] = random.choice("ATCG")
                else:
                    seq[start+i] = seed[i]
    return ["".join(seq) for seq in seqs], indices

def generate_dataset():
    dataset = Dataset()

    num_sequences = 100
    sequence_length = 100
    sub_length = 10
    sequences,indices = gen_sequences(num_sequences, sequence_length, seed_length=sub_length, mut_rate=0.25)

    seq_ids = [dataset.add_sequence(sequence) for sequence in sequences]
    dataset.add_alignment(seq_ids, indices, sub_length)

    return dataset
