import random
from galg import GeneticAlgorithm, Population, Selection, SelMethod
from network import Network, GNGNetwork
from dataset import Dataset, generate_dataset
from alphabet import get_nucleotides, get_amino_acids

def gen_sequences(alphabet, num_sequences, length, seed_length=0):
    print("Generating %d sequences of length %d..." % (num_sequences, length))
    seqs = [[random.choice(alphabet) for _ in range(length)]
            for _ in range(num_sequences)]

    if seed_length > 0:
        seed = [random.choice(alphabet) for _ in range(seed_length)]
        print("Planting seed... (%s)" % "".join(seed))
        for seq in seqs:
            start = random.randint(0, length-seed_length-1)
            for i in range(seed_length):
                seq[start+i] = seed[i]
    return ["".join(seq) for seq in seqs]

def print_sequences(sequences):
    print("Sequences:")
    for seq in sequences: print("".join(seq))
    print("")

def print_subsequences(sequences, length, indices):
    print("Subsequences:")
    for seq,i in zip(sequences, indices):
        print(seq[i:i+length])
    print("")

def ga_main():
    # Sequences
    num_sequences = 100
    sequence_length = 100
    sub_length = 10
    sequences = gen_sequences(get_nucleotides(), num_sequences, sequence_length, seed_length=sub_length)
    print_sequences(sequences)

    # Population
    pop_size = 500
    num_random = 50
    pop = Population(Network(get_nucleotides()), sequences, sub_length, pop_size, num_random)

    # Genetic Algorithm
    #selection = Selection(SelMethod.truncation, 10)
    # RouletteWheel Selection
    #selection = Selection(SelMethod.roulettewheel, 0)
    # Stochastic Universal Selection
    #selection = Selection(SelMethod.stochastic, 0)
    # Ranking selection (quick version)
    #selection = Selection(SelMethod.ranking_quick, 0)
    # Ranking selection (slow version)
    #selection = Selection(SelMethod.ranking_slow, 0)
    # Tournament Selection
    selection = Selection(SelMethod.tournament, 10)
    mutation_rate = 0.05
    galg = GeneticAlgorithm(selection, mutation_rate)

    iterations = 100
    best = galg.run(pop, num_iterations=iterations, verbose=True)

    print_subsequences(sequences, sub_length, best)

def gng_main():
    dataset = Dataset(get_amino_acids(), "data")
    columns = dataset.get_columns()[:1000]

    network = GNGNetwork(get_amino_acids(), size=20, verbose=False)
    print("Columns: %d" % len(columns))
    total_score = 0
    for column in columns:
        #print(" ".join("%.4f" % x if x > 0.0 else "      " for x in column))
        #print(column, network.evaluate_column(column))
        total_score += network.evaluate_column(column)
    print("Pretraining total score: %f" % (total_score / len(columns)))
    network.train(columns, 10)

    print("GNG Nodes:")
    for loc in network.gng.locations: print(loc)

    print("Columns:")
    total_score = 0
    for column in columns:
        #print(column, network.evaluate_column(column))
        total_score += network.evaluate_column(column)
    print("Posttraining total score: %f" % (total_score / len(columns)))

if __name__ == "__main__":
    #ga_main()
    gng_main()
