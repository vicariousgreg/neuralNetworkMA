import random
from alphabet import get_nucleotides, get_amino_acids
from readXml import FilterSequences
from os import listdir

class SequenceCluster:
    def __init__(self, alphabet, sequences, alignments):
        #print(sequences)
        #print(alignments)
        if any(len(indices) != len(sequences) for length,indices in alignments):
            raise ValueError("Invalid index list length!")

        alignments = [(l,i) for l,i in alignments if not any(x is None for x in i)]

        self.sequences = sequences
        self.alignments = alignments

        self.columns = []
        for length,indices in alignments:
            #print(length, indices)
            for offset in range(length):
                counts = dict(zip(alphabet, [0] * len(alphabet)))
                for sequence,start in zip(sequences, indices):
                    if start is None: continue
                    if start+offset >= len(sequence):
                        raise ValueError("Invalid index start(%d) offset(%d) into len(%d)!" % (start, offset, len(sequence)))
                    try:
                        counts[sequence[start+offset]] += 1.0
                    except KeyError:
                        raise ValueError("Invalid letter %s in sequence!" % sequence[start+offset])
                self.columns.append([x / len(sequences) for x in counts.values()])

class Dataset:
    def __init__(self, alphabet, directory):
        self.alphabet = alphabet
        self.sequence_clusters = []
        self.columns = []
        for f in listdir(directory):
            path = "%s/%s" % (directory, f)
            #print(path)
            filtered = FilterSequences(path)
            #print(filtered.AlignmentEntry)
            try:
                cluster = SequenceCluster(alphabet,
                    filtered.sequences, filtered.AlignmentEntry)
                self.sequence_clusters.append(cluster)
                self.columns += cluster.columns
            except ValueError as e:
                print(e)
                print("Invalid data in %s.  Skipping..." % path)

    def get_columns(self):
        return self.columns

def gen_sequences(alphabet, num_sequences, length, seed_length=0, mut_rate=None):
    print("Generating %d sequences of length %d..." % (num_sequences, length))
    seqs = [[random.choice(alphabet) for _ in range(length)]
            for _ in range(num_sequences)]

    indices = []
    if seed_length > 0:
        seed = [random.choice(alphabet) for _ in range(seed_length)]
        print("Planting seed... (%s)" % "".join(seed))
        for seq in seqs:
            start = random.randint(0, length-seed_length-1)
            indices.append(start)
            for i in range(seed_length):
                if mut_rate and random.random() < mut_rate:
                    seq[start+i] = random.choice(alphabet)
                else:
                    seq[start+i] = seed[i]
    return ["".join(seq) for seq in seqs], indices

def generate_dataset():
    dataset = Dataset(get_nucleotides())

    num_sequences = 100
    sequence_length = 100
    sub_length = 10
    sequences,indices = gen_sequences(get_nucleotides(), num_sequences, sequence_length, seed_length=sub_length, mut_rate=0.25)

    seq_ids = [dataset.add_sequence(sequence) for sequence in sequences]
    dataset.add_alignment(seq_ids, indices, sub_length)

    return dataset

def test():
    dataset = Dataset(get_amino_acids(), "data")
    print("Sequence_clusters: %d" % len(dataset.sequence_clusters))
    print("Columns: %d" % len(dataset.get_columns()))

    s = 0
    for cluster in dataset.sequence_clusters:
        for length,indices in cluster.alignments:
            s += length
    print(s)
