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
