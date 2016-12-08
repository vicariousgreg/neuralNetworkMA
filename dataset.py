import random
from alphabet import get_nucleotides, get_amino_acids, handle_special_letters
from readXml import FilterSequences
from os import listdir

FACTOR = 1.0

def randomCategory(probs):
    r = random.random() # range: [0,1)
    total = 0           # range: [0,1]
    for i,prob in enumerate(probs):
        total += prob
        if total>r:
            return i
    raise Exception('distribution not normalized: {probs}'.format(probs=probs))

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
                        handle_special_letters(alphabet, counts, sequence[start+offset])
                        #raise ValueError("Invalid letter %s in sequence!" % sequence[start+offset])
                self.columns.append([FACTOR * x / len(sequences) for x in counts.values()])

    def get_unaligned_columns(self, alphabet):
        columns = []
        i = 0
        done = False
        while not done:
            counts = dict(zip(alphabet, [0] * len(alphabet)))
            for sequence in self.sequences:
                if i >= len(sequence):
                    done = True
                    break
                try:
                    counts[sequence[i]] += 1.0
                except KeyError:
                    handle_special_letters(alphabet, counts, sequence[i])
            columns.append([FACTOR * x / len(self.sequences) for x in counts.values()])
            i += 1
        return columns

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
            if len(filtered.sequences) <= 20: continue

            try:
                cluster = SequenceCluster(alphabet,
                    filtered.sequences, filtered.AlignmentEntry)
                if len(cluster.columns) > 0:
                    self.sequence_clusters.append(cluster)
                    self.columns += cluster.columns
            except ValueError as e:
                print(e)
                print("Invalid data in %s.  Skipping..." % path)
        self.calc_statistics()
        self.generate_random_columns(len(self.columns))

    def get_columns(self):
        return self.columns

    def get_random_columns(self):
        return self.random_columns

    def generate_random_columns(self, count, num_seq=None):
        if num_seq is None: num_seq = int(self.average_seq_per_alignment)

        self.random_columns = []
        for _ in range(count):
            counts = dict(zip(self.alphabet, [0] * len(self.alphabet)))
            for i in range(num_seq):
                counts[self.alphabet[
                    randomCategory(self.overall_letter_distribution)]] += 1.0
            self.random_columns.append([FACTOR * x / num_seq for x in counts.values()])

    def calc_statistics(self):
        # Alphabet distribution overall
        counts = dict(zip(self.alphabet, [0] * len(self.alphabet)))
        for cluster in self.sequence_clusters:
            for sequence in cluster.sequences:
                for i in range(len(sequence)):
                    try:
                        counts[sequence[i]] += 1.0
                    except KeyError:
                        handle_special_letters(self.alphabet, counts, sequence[i])
        total = sum(counts.values())
        self.overall_letter_distribution = \
            [(counts[letter] / total) for letter in self.alphabet]

        # Alphabet distribution in alignments
        counts = dict(zip(self.alphabet, [0] * len(self.alphabet)))
        for cluster in self.sequence_clusters:
            for length,indices in cluster.alignments:
                for sequence,start in zip(cluster.sequences, indices):
                    if start is not None:
                        for offset in range(length):
                            try:
                                counts[sequence[start+offset]] += 1.0
                            except KeyError:
                                handle_special_letters(self.alphabet, counts, sequence[start+offset])
        total = sum(counts.values())
        self.alignment_letter_distribution = \
            [(counts[letter] / total) for letter in self.alphabet]

        # Number of alignments
        self.num_alignments = \
            sum(len(cluster.alignments)
            for cluster in self.sequence_clusters)

        # Number of sequences and clusters
        self.num_sequences = sum(len(cluster.sequences) for cluster in self.sequence_clusters)

        # Average number of sequences per alignment
        self.average_seq_per_alignment = \
            (sum(len([i for i in indices if i is not None])
            for cluster in self.sequence_clusters
            for length,indices in cluster.alignments) / self.num_alignments)

        # Average column
        self.average_column = [0.0] * len(self.alphabet)
        for column in self.columns:
            for i in range(len(column)):
                self.average_column[i] += column[i]
        for i in range(len(self.average_column)):
            self.average_column[i] /= len(self.columns)

    def print_statistics(self):
        # Alphabet distribution overall
        print("Overall letter distribution:")
        print(" ".join("%.4f" % x for x in self.overall_letter_distribution))

        # Alphabet distribution in alignments
        print("Alignment letter distribution:")
        print(" ".join("%.4f" % x for x in self.alignment_letter_distribution))

        # Number of alignments
        print("Total alignment count: %d" % self.num_alignments)

        # Number of columns
        print("Total column count: %d" % len(self.columns))

        # Number of sequences and clusters
        print("Total sequences: %d" % self.num_sequences)
        print("Total sequence clusters: %d" % len(self.sequence_clusters))

        # Average number of sequences per alignment
        print("Average sequences per alignment: %f" % self.average_seq_per_alignment)

        # Average Column
        print("Average column:")
        print(" ".join("%.4f" % x for x in self.average_column))

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

def generate_random_columns(alphabet, num_columns):
    columns = []
    for _ in range(num_columns):
        column = [random.random() for x in range(len(alphabet))]
        total = sum(x for x in column)
        columns.append([FACTOR * x / total for x in column])
    return columns

def test():
    dataset = Dataset(get_amino_acids(), "data")
    print("Sequence_clusters: %d" % len(dataset.sequence_clusters))
    print("Columns: %d" % len(dataset.get_columns()))

    s = 0
    for cluster in dataset.sequence_clusters:
        for length,indices in cluster.alignments:
            s += length
    print(s)
