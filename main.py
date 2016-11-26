import random
from galg import GeneticAlgorithm, Population
from network import Network

def gen_sequences(num_sequences, length, seed_length=0):
    print("Generating %d sequences of length %d..." % (num_sequences, length))
    seqs = [[random.choice("ATCG") for _ in xrange(length)]
            for _ in xrange(num_sequences)]

    if seed_length > 0:
        seed = [random.choice("ATCG") for _ in xrange(seed_length)]
        print("Planting seed... (%s)" % "".join(seed))
        for seq in seqs:
            start = random.randint(0, length-seed_length)
            for i in xrange(seed_length):
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

def main():
    # Sequences
    num_sequences = 100
    sequence_length = 100
    sub_length = 10
    sequences = gen_sequences(num_sequences, sequence_length, seed_length=sub_length)
    print_sequences(sequences)

    # Population
    pop_size = 500
    num_random = 50
    pop = Population(Network(), sequences, sub_length, pop_size, num_random)

    # Genetic Algorithm
    tournament_size = 10
    mutation_rate = 0.05
    galg = GeneticAlgorithm(tournament_size, mutation_rate)

    iterations = 100
    best = galg.run(pop, num_iterations=iterations, verbose=True)

    print_subsequences(sequences, sub_length, best)

if __name__ == "__main__":
    main()
