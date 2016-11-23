import random
import copy
import time

class Individual:
    def __init__(self, chromosome):
        self.chromosome = chromosome
        self.fitness = None

    def copy(self):
        return copy.deepcopy(self)

    @staticmethod
    def random(max_gene_values):
        return Individual(
            [int(random.randint(0, x))
             for x in max_gene_values])

    def mutate(self, rate, max_gene_values):
        for i in xrange(len(self.chromosome)):
            if random.random() < rate:
                # Reset fitness and mutate
                self.fitness = None
                self.chromosome[i] = int(random.randint(0, max_gene_values[i]))

    def crossover(self, other):
        return Individual(
            [random.choice(genes) for genes in
                zip(self.chromosome, other.chromosome)])

class Population:
    def __init__(self, network, sequences, chromosome_length, size, num_random=0):
        # Check number of new random individuals per generation
        if num_random < 0 or num_random > size:
            raise ValueError("Number of new random individuals " + \
                "must be between zero and population size!")

        # Determine maximum indices for each sequence
        # If any sequence is shorter than the chromosome_length, raise an error
        self.max_indices = [len(seq) - chromosome_length for seq in sequences]
        if any(index < 0 for index in self.max_indices):
            raise ValueError("Not all sequences are long enough for alignment!")

        # Set attributes
        self.sequences = sequences
        self.chromosome_length = chromosome_length
        self.network = network
        self.num_random = num_random

        # Generate individuals
        self.individuals = \
            [Individual.random(self.max_indices)
            for _ in xrange(size)]

        # Calculate initial fitnesses
        self.calc_fitness()

    def calc_fitness(self):
        # Use network to evaluate fitness
        for individual in self.individuals:
            # Only calculate fitness if it is None to avoid recalculation
            if individual.fitness is None:
                individual.fitness = self.network.evaluate(
                    self.sequences, self.chromosome_length, individual.chromosome)

    def mutate(self, rate):
        for individual in self.individuals:
            individual.mutate(rate, self.max_indices)

    def tournament(self, size):
        # Return the most fit from a random sample
        return max(
            (self.individuals[
                random.randint(0, len(self.individuals)-1)]
                for _ in xrange(size)),
            key = lambda indiv: indiv.fitness)

    def generation(self, tournament_size, elitism=True):
        # Add new random individuals.
        new_individuals = \
            [Individual.random(self.max_indices)
             for _ in xrange(self.num_random)]

        # Keep best if elitism is true.
        if elitism: new_individuals.append(self.get_best())

        # Perform tournaments to fill population
        while len(new_individuals) < len(self.individuals):
            # Run tournaments to determine parents
            a = self.tournament(tournament_size)
            b = self.tournament(tournament_size)

            # Ensure unique parents
            while a == b: b = self.tournament(tournament_size)

            # Crossover
            new_individuals.append(a.crossover(b))

        # Set the new population
        self.individuals = new_individuals

    def get_best(self):
        return max(self.individuals, key=lambda indiv:indiv.fitness)

class GeneticAlgorithm:
    def __init__(self, tournament_size=5, mutation_rate=0.005):
        self.tournament_size = tournament_size
        self.mutation_rate = mutation_rate

    def run(self, population, num_iterations=None, fitness_threshold=None,
            time_limit=None, verbose=False):
        iterations = 0
        best_of_best = population.get_best().copy()

        ##################################
        ### Set up stopping conditions ###
        ##################################
        # Iteration cap checker.
        if num_iterations: iterations_done = lambda: iterations >= num_iterations
        else: iterations_done = lambda: False

        # Fitness threshold checker.
        if fitness_threshold:
            fitness_done = lambda: best_of_best.fitness >= fitness_threshold
        else:
            fitness_done = lambda: False

        # Time checker.
        if time_limit:
            time_limit += time.time()
            time_done = lambda: time.time() > time_limit
        else:
            time_done = lambda: False

        # OR done checker.
        def is_done(): return iterations_done() or fitness_done() or time_done()

        ######################
        ### Run iterations ###
        ######################
        while not is_done():
            iterations += 1
            if verbose:
                print("Iteration %d: (best fitness %f)" \
                    % (iterations, best_of_best.fitness))

            # New generation
            population.generation(self.tournament_size, elitism=True)

            # Mutate
            population.mutate(self.mutation_rate)

            # Recalculate fitness
            population.calc_fitness()

            # Find new best
            best = population.get_best()
            if best.fitness > best_of_best.fitness:
                best_of_best = best.copy()

        return best_of_best.chromosome



class Network:
    def evaluate(self, sequences, length, indices):
        score = 0
        for i in xrange(length):
            counts = dict((("A", 0), ("C", 0), ("T", 0), ("G", 0)))
            for seq,start in zip(sequences, indices):
                counts[seq[start+i]] += 1
            score += max(counts.values())
        return score



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
