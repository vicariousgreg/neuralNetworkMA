import random
from galg import GeneticAlgorithm, Population, Selection, SelMethod
from network import Network, GNGNetwork
from dataset import Dataset
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

def evaluate(network, columns, verbose=False):
    total_score = 0
    count = 0
    for column in columns:
        count += 1
        total_score += network.evaluate_column(column, verbose=verbose)
    return total_score / count

def gng_main():
    # Parameters
    max_columns = 50000
    network_size = 200
    iterations = 100
    training_fraction = 0.05

    # Generate dataset and print statistics
    training_dataset, test_dataset = \
        Dataset.create_split_dataset(get_amino_acids(), "data")
    print("Training set:")
    training_dataset.print_statistics()
    print("")
    print("Test set:")
    test_dataset.print_statistics()
    print("")

    # Pull columns
    training_columns = training_dataset.get_columns()[:max_columns]
    training_random_columns = training_dataset.get_random_columns()[:max_columns]
    test_columns = test_dataset.get_columns()[:max_columns]
    test_random_columns = test_dataset.get_random_columns()[:max_columns]

    # Create network
    network = GNGNetwork(get_amino_acids(), size=network_size, verbose=False)

    def print_scores():
        print("Scores: Training[ %.6f / %.6f ]  Test[ %.6f / %.6f ]" % \
            (evaluate(network, training_columns),
             evaluate(network, training_random_columns),
             evaluate(network, test_columns),
             evaluate(network, test_random_columns)))

    # Print pre-scores
    print("Columns: Training[ %6d / %6d ]  Test[ %6d / %6d]" % \
        (len(training_columns), len(training_random_columns),
         len(test_columns),      len(test_random_columns)))
    print_scores()

    # Perform training
    for _ in range(iterations):
        network.train(training_columns, 1, fraction=training_fraction, verbose=False)
        print_scores()

    # Print resulting network and post-scores
    network.gng.print_nodes()
    print_scores()

if __name__ == "__main__":
    #ga_main()
    gng_main()
