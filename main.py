import random
from galg import GeneticAlgorithm, Population, Selection, SelMethod
from network import Network, GNGNetwork
from dataset import Dataset, SequenceCluster
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

def gen_cluster(num_sequences, sequence_length, sub_length):
    sequences = gen_sequences(get_amino_acids(), num_sequences, sequence_length, seed_length=sub_length)
    print_sequences(sequences)
    return SequenceCluster(get_amino_acids(), sequences, [])

def print_sequences(sequences):
    print("Sequences:")
    for seq in sequences: print("".join(seq))
    print("")

def print_subsequences(sequences, length, indices):
    print("\nSubsequences:")
    for seq,i in zip(sequences, indices):
        print(seq[i:i+length])
    print("")

def ga_main(sequence_cluster, length, network):
    # Population
    pop_size = 500
    num_random = 50
    pop = Population(network, sequence_cluster.sequences, length, pop_size, num_random)

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

    print_subsequences(sequence_cluster.sequences, length, best)

    print("Resulting indices:")
    print(" ".join("%3d" % i for i in best))

    def calc_diff(i_a, i_b):
        return sum(((x-y) ** 2) for x,y in zip(i_a, i_b)) ** 0.5

    print("\nActual alignment indices:")
    comparisons = \
        [(length, calc_diff(indices, best),
         " ".join("%3d" % i for i in indices))
         for length,indices in sequence_cluster.alignments]

    for length,diff,indices in sorted(comparisons, key=lambda x:x[1]):
        print("len: %3d diff: %11.6f\n  %s" % \
            (length, diff, indices))

def evaluate(network, columns, verbose=False):
    total_score = 0
    count = 0
    for column in columns:
        count += 1
        total_score += network.evaluate_column(column, verbose=verbose)
    return total_score / count

def gng_main():
    # Parameters
    max_columns = 1000
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

    def print_scores(iteration):
        print("Scores (%4d): Training[ %.6f / %.6f ]  Test[ %.6f / %.6f ]" % \
            (iteration,
             evaluate(network, training_columns),
             evaluate(network, training_random_columns),
             evaluate(network, test_columns),
             evaluate(network, test_random_columns)))

    # Print pre-scores
    print("Columns: Training[ %6d / %6d ]  Test[ %6d / %6d]" % \
        (len(training_columns), len(training_random_columns),
         len(test_columns),      len(test_random_columns)))
    print_scores(0)

    # Perform training
    for _ in range(iterations):
        network.train(training_columns, 1, fraction=training_fraction, verbose=False)
        print_scores(_)

    # Print resulting network and post-scores
    network.gng.print_nodes()
    print_scores(_+1)

    return training_dataset, test_dataset, network

if __name__ == "__main__":
    train, test, network = gng_main()

    # Create sequence cluster
    sub_length=10
    #cluster = gen_cluster(num_sequences=100, sequence_length=100, sub_length=sub_length)
    cluster = test.sequence_clusters[0]

    #ga_main(cluster, sub_length, Network(get_amino_acids()))
    ga_main(cluster, sub_length, network)
