import random
from sys import argv
from time import time
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

def ga_main(sequence_cluster, search_length, network):
    # Population
    pop_size = 500
    num_random = 50
    pop = Population(network, sequence_cluster.sequences, search_length, pop_size, num_random)

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

    iterations = 25
    best = galg.run(pop, num_iterations=iterations, verbose=True)

    print_subsequences(sequence_cluster.sequences, search_length, best)

    print("Resulting indices (score: %6.4f):" % \
        network.evaluate(sequence_cluster.sequences, search_length, best))
    print(" ".join("%3d" % i for i in best))

    def calc_diff(i_a, i_b):
        return sum(((x-y) ** 2) for x,y in zip(i_a, i_b)) ** 0.5

    print("\nActual alignment indices:")
    comparisons = \
        [(length,
         calc_diff(indices, best),
         network.evaluate(
            sequence_cluster.sequences, length, indices),
         " ".join("%3d" % i for i in indices))
         for length,indices in sequence_cluster.alignments]

    for length,diff,score,indices in sorted(comparisons, key=lambda x:x[1]):
        print("len: %3d score: %11.6f: diff: %11.6f\n  %s" % \
            (length, score, diff, indices))

def evaluate(network, columns, verbose=False):
    total_score = 0
    count = 0
    for column in columns:
        count += 1
        total_score += network.evaluate_column(column, verbose=verbose)
    return total_score / count

def gng_main(training_dataset, test_dataset, network_size, mean, time_limit):
    # Parameters
    max_columns = 50000 # Use all
    training_fraction = 0.05

    # Print statistics
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
    network = GNGNetwork(get_amino_acids(), size=network_size, mean=mean, verbose=False)

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

    # Perform training
    start = time()
    iterations = 0
    while time() - start < time_limit:
        iterations += 1
        network.train(training_columns, 1, fraction=training_fraction, verbose=False)
        #print_scores(iterations)

    # Print resulting network and post-scores
    network.gng.print_nodes()
    #print_scores(iterations+1)

    return network

def build_dataset(fraction=0.7):
    return Dataset.create_split_dataset(
        get_amino_acids(), "data", training_fraction=fraction)

def evaluate_network(training_dataset, test_dataset, network):
    # Pull columns
    training_columns = training_dataset.get_columns()
    training_random_columns = training_dataset.get_random_columns()
    test_columns = test_dataset.get_columns()
    test_random_columns = test_dataset.get_random_columns()

    tr_score, tr_rand_score, te_score, te_rand_score = \
        (evaluate(network, training_columns),
         evaluate(network, training_random_columns),
         evaluate(network, test_columns),
         evaluate(network, test_random_columns))

    print("Scores: Training[ %.6f / %.6f (ratio %10.2f)]  Test[ %.6f / %.6f (ratio %10.2f)]" % \
        (tr_score, tr_rand_score, tr_score / tr_rand_score,
         te_score, te_rand_score, te_score / te_rand_score))

if __name__ == "__main__":
    try:
        train = Dataset.load("dataset/train.dataset")
        test = Dataset.load("dataset/test.dataset")
        print("Loaded datasets from ./datasets/")
    except FileNotFoundError:
        print("Creating training and test sets ...")
        train,test = build_dataset()
        print(" ... saving to ./datasets/")
        train.save("dataset/train.dataset")
        test.save("dataset/test.dataset")

    mean = float(argv[1])
    for size in [20, 100, 200, 1000]:
        path = "networks/size%d-mean%0.3f.network" % (size, mean)
        try:
            network = GNGNetwork.load(path, mean)
            print("Loaded network (%d nodes, mean %f) from %s" % (size, mean, path))
        except FileNotFoundError:
            # Train for 3 hours
            minutes = 60 * 3
            print("Training network (%d nodes, mean %f) for %d minutes..." % (size, mean, minutes))
            network = gng_main(train,test,size, mean, 60 * minutes)
            print("... saving to %s" % path)
            network.save(path)
        print("Evaluating network...")
        print("Training iterations: %d" % network.run_iterations)
        evaluate_network(train, test, network)
        print("")

    '''
    network = gng_main(train,test,200,30)

    # Create sequence cluster
    sub_length=10
    #cluster = gen_cluster(num_sequences=100, sequence_length=100, sub_length=sub_length)
    cluster = test.sequence_clusters[0]

    #ga_main(cluster, sub_length, Network(get_amino_acids()))
    ga_main(cluster, sub_length, network)
    '''
