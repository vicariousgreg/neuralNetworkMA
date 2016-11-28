import random
import copy
import time
import math
from enum import  Enum


class SelMethod(Enum):
    truncation = 1
    roulettewheel = 2
    stochastic = 3
    ranking_quick = 4
    ranking_slow = 5
    tournament = 6


class Selection:
    def __init__(self, method, size=0):
        self.method = method
        self.cnt = 0
        self.size = size
        self.threshold = 0.0
        self.rank_val = (0, 0)

    def reset(self, population):
        self.cnt = 0
        if self.method is SelMethod.stochastic:
            self.threshold = random.random()\
                              * (1 / len(population.individuals))
        if population is not None and\
            self.method is SelMethod.ranking_quick or\
            self.method is SelMethod.ranking_slow:
            self.prep_ranking(population)

    def prep_ranking(self, population):
        i = 0
        pop_size = len(population.individuals)
        for individual in population.individuals:
            if individual.fitness is not None:
                break
            else:
                i += 1
        median = int((i + pop_size - 1) / 2)
        pressure = population.individuals[pop_size - 1].fitness\
                    / population.individuals[median].fitness
        if pressure <= 0:
            raise ValueError("pressure value is below or equal to zero")
        a = (2 * pop_size - pressure * (pop_size + 1))\
            / (pop_size * (pop_size - 1))
        b = (2 * (pressure - 1)) / (pop_size * (pop_size - 1))
        self.rank_val = (a, b)
        if self.method is SelMethod.ranking_slow:
            population.calc_proportion(self.rank_val)


class Individual:
    def __init__(self, chromosome):
        self.chromosome = chromosome
        self.fitness = None
        self.prop = 0.0

    def copy(self):
        return copy.deepcopy(self)

    @staticmethod
    def random(max_gene_values):
        return Individual(
            [int(random.randint(0, x))
             for x in max_gene_values])

    def mutate(self, rate, max_gene_values):
        for i in range(len(self.chromosome)):
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
            for _ in range(size)]

        # Calculate initial fitnesses
        self.calc_fitness()

        # Calculate initial proportion
        self.calc_proportion()

    def calc_fitness(self):
        # Use network to evaluate fitness
        for individual in self.individuals:
            # Only calculate fitness if it is None to avoid recalculation
            if individual.fitness is None:
                individual.fitness = self.network.evaluate(
                    self.sequences, self.chromosome_length, individual.chromosome)
       
    def calc_proportion(self, rank_val=None):
        # Calculate proportion of an individual
        if rank_val is None:
            sum_fitness = 0
            for individual in self.individuals:
                if individual.fitness is not None:
                    sum_fitness += individual.fitness
            if sum_fitness <= 0:
                raise ValueError("The sum of fitnesses is zero")
            for individual in self.individuals:
                individual.prop = individual.fitness / sum_fitness
        else:
            a, b = rank_val
            i = 1
            for individual in self.individuals:
                individual.fitness = a + b * i
                i += 1
        
    def sort(self, reverse=True):
        # sort the individuals 
        self.individuals = sorted(self.individuals, key=lambda indiv: indiv.fitness, reverse=reverse)

    def mutate(self, rate):
        for individual in self.individuals:
            individual.mutate(rate, self.max_indices)

    def truncation(self, select):
        # Use truncation method for selection
        # "size" = # of individuals in the ordered generation
        cnt = (select.cnt % len(self.individuals)) % select.size
        return self.individuals[cnt - 1]

    def proportional(self, threshold):
        # Use proportional selection
        sum_prop = 0.0
        for individual in self.individuals:
            sum_prop += individual.prop
            if sum_prop >= threshold:
                return individual
        raise ValueError("Proportional selection should not reach here")
        return None

    def roulettewheel(self, select=None):
        # Use RouletteWheel selection
        return self.proportional(random.random())

    def stochastic(self, select=None):
        # Use Stochastic Universal selection
        if select.cnt > 1:
            select.threshold += (1 / len(self.individuals))
        if select.threshold > 1.0:
            select.threshold -= 1.0
        return self.proportional(select.threshold)

    def ranking_quick(self, select):
        # Use (Linear) Ranking selection
        # Note: this version might fail if the proportion of the median is
        # the same as that of the fittest one
        a, b = select.rank_val
        if b == 0 :
            raise ValueError("b value becomes zero during ranking select")
        r = random.random()
        in_root = (2 * a + b) ** 2 + 8 * b * r
        k = (-1 * (2 * a + b) + math.sqrt(in_root)) / (2 * b)
        return self.individuals[math.ceil(k-1)]

    def ranking_slow(self, select=None):
        # Use Linear Ranking selection
        return self.proportional(random.random())

    def tournament(self, select):
        # Return the most fit from a random sample
        return max(
            (self.individuals[
                random.randint(0, len(self.individuals)-1)]
                for _ in range(select.size)),
            key = lambda indiv: indiv.fitness)

    def pre_selection(self, select):
        if select.method is SelMethod.truncation:
            self.sort()
        elif select.method is SelMethod.ranking_quick or\
              select.method is SelMethod.ranking_slow:
            self.sort(reverse=False)

    def selection(self, select):
        funcs = {1 : self.truncation,
                 2 : self.roulettewheel,
                 3 : self.stochastic,
                 4 : self.ranking_quick,
                 5 : self.ranking_slow,
                 6 : self.tournament,
                }
        select.cnt += 1
        return funcs[select.method.value](select)

    def generation(self, select, elitism=True):
        # Add new random individuals.
        new_individuals = \
            [Individual.random(self.max_indices)
             for _ in range(self.num_random)]

        # Keep best if elitism is true.
        if elitism: new_individuals.append(self.get_best())

        # Initialize selection
        self.pre_selection(select)
        select.reset(self)
        while len(new_individuals) < len(self.individuals):
            # Run selection to determine parents
            a = self.selection(select)
            b = self.selection(select)

            # Ensure unique parents
            while a == b: b = self.selection(select)

            # Crossover
            new_individuals.append(a.crossover(b))

        # Set the new population
        self.individuals = new_individuals

    def get_best(self):
        return max(self.individuals, key=lambda indiv:indiv.fitness)

class GeneticAlgorithm:
    def __init__(self, selection, mutation_rate=0.005):
        self.selection = selection
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
            population.generation(self.selection, elitism=True)

            # Mutate
            population.mutate(self.mutation_rate)

            # Recalculate fitness and proportion
            population.calc_fitness()

            # Recalculate proportion
            population.calc_proportion()

            # Find new best
            best = population.get_best()
            if best.fitness > best_of_best.fitness:
                best_of_best = best.copy()

        return best_of_best.chromosome
  
