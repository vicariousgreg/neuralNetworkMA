import random
import copy
from numpy import array, zeros, ones, matrix, dot, multiply, add, subtract, copyto
from numpy.random import rand
from math import exp

def distance(a,b):
    return sum((a[i] - b[i])**2 for i in range(len(a)))

def rbf(distance, mean):
    return exp(-distance / (2 * (mean**2)))

def move(location, destination, factor):
    return [location[i] + ((destination[i]-location[i])*factor)
                for i in range(len(location))]

class GrowingNeuralGas:
    def __init__(self, size, seed_size=2, es=0.05, en=0.0005, beta=0.9999,
                       max_age=25, error_threshold=25.0, feature_length=13, verbose=False):
        self.size = size
        self.output = zeros(self.size)
        self.error = [0.0 for _ in range(self.size)]
        self.verbose = verbose
        self.locked = False

        ##################
        # GNG attributes #
        ##################
        # Neuron motion factors
        self.es = es
        self.en = en
        # Max edge age
        self.max_age = max_age
        # Error coefficient
        self.beta = beta
        # Error threshold for adding new nodes
        self.error_threshold = error_threshold

        self.means = zeros(self.size)
        self.distances = zeros(self.size)
        self.active_neurons = [False] * self.size
        self.num_active_neurons = 0

        # Neuron locations
        self.locations = zeros((self.size, feature_length))

        # Edge age matrix
        # Be sure to update symmetrically
        # -1 means no edge, otherwise value is age
        self.edges = zeros((self.size, self.size))
        self.edges.fill(-1)

        # Create seed neurons.
        for i in range(seed_size):
            self.add_neuron()
            # Randomize locations
            for j in range(feature_length):
                self.locations[i][j] = random.random()

        # Add edges
        for i in range(seed_size):
            for j in range(seed_size):
                if i != j:
                    self.set_edge_age(i, j, 0)

        # Calculate means
        for i in range(seed_size):
            self.recalculate_mean(i)

    def lock(self, value=True):
        self.locked = value

    def insertion_criteria(self):
        return self.num_active_neurons < self.size and \
            any(error > self.error_threshold for error in self.error)

    def recalculate_mean(self, index):
        # Set RBF mean of moved neurons to avg dist to neighbors
        # Calculate average distance to neighbors
        total_distance = 0.0
        count = 0
        for n in range(self.size):
            if self.edges[index][n] >= 0:
                count += 1
                total_distance += distance(self.locations[index], self.locations[n]) ** 0.5
        self.means[index] = total_distance / count

    def remove_neuron(self, index):
        self.active_neurons[index] = False
        self.num_active_neurons -= 1
        if self.verbose: print("removed neuron (size %d)" % self.num_active_neurons)

    def add_neuron(self):
        n = self.active_neurons.index(False)
        self.active_neurons[n] = True
        self.num_active_neurons += 1

        if self.verbose: print("added neuron (size %d)" % self.num_active_neurons)
        if self.verbose: print(max(self.error))
        return n

    def set_edge_age(self, a, b, age):
        if self.verbose:
            if age == -1: print("Removing edge %d %d" % (a,b))
            elif self.edges[a][b] == -1: print("Adding edge %d %d" % (a,b))
        self.edges[a][b] = age
        self.edges[b][a] = age

    def feedforward(self, values):
        # Reset output.
        self.output.fill(0.0)

        # Calculate distances and activate neurons.
        for i in range(self.size):
            if not self.active_neurons[i]: continue
            d = distance(self.locations[i], values)
            self.distances[i] = d
            self.output[i] = rbf(d, self.means[i])

        # If unlocked, perform additional computation for learning.
        if not self.locked: 
            # Get the closest neurons that are active
            closest = list(self.distances.argsort())
            while not self.active_neurons[closest[0]]: del closest[0]
            while not self.active_neurons[closest[1]]: del closest[1]

            # Add error to closest neuron
            self.error[closest[0]] += self.distances[closest[0]]

            # Adjust weights
            self.adjust_weights(values, closest[0], closest[1])

        return self.output

    def adjust_weights(self, last_input, closest_neuron, second_closest_neuron):
        # Move closest neuron closer to input vector by es
        self.locations[closest_neuron] = \
            move(self.locations[closest_neuron], last_input, self.es)
        # Move s's neighbors closer by en
        changed = [closest_neuron]
        for i in range(self.size):
            if self.edges[closest_neuron][i] >= 0:
                changed.append(i)
                self.locations[i] = \
                    move(self.locations[i], last_input, self.en)

        # Set RBF mean of moved neurons to avg dist to neighbors
        for i in changed:
            self.recalculate_mean(i)

        # Increment s's edge ages
        # Expire edges if they are too old
        for i in range(self.size):
            age = self.edges[closest_neuron][i]
            if age >= 0:
                new_age = age+1
                if new_age > self.max_age:
                    self.set_edge_age(i, closest_neuron, -1)
                else:
                    self.set_edge_age(i, closest_neuron, new_age)

        # Create edge b/w closest and second closest (set age to 0)
        self.set_edge_age(closest_neuron, second_closest_neuron, 0)

        # Remove neurons that don't have edges
        for i in range(self.size):
            if self.active_neurons[i]:
                if all(self.edges[i][n] < 0 for n in range(self.size)):
                    self.remove_neuron(i)

        # Insert new neuron
        if self.insertion_criteria():
            # Find neuron with largest error
            max_error = sorted(self.error)[-1]
            u = self.error.index(max_error)

            worst_neighbor = 0
            worst_error = 0
            for i in range(self.size):
                if self.edges[u][i] >= 0:
                    err = self.error[i]
                    if err > worst_error:
                        worst_neighbor = i
                        worst_error = err

            # Insert new node
            n = self.add_neuron()

            # Set location to midpoint
            self.locations[n] = move(self.locations[u], self.locations[worst_neighbor], 0.5)

            # Set up edges
            self.set_edge_age(u, worst_neighbor, -1)
            self.set_edge_age(u, n, 0)
            self.set_edge_age(n, worst_neighbor, 0)

            # Recalculate means
            self.recalculate_mean(u)
            self.recalculate_mean(n)
            self.recalculate_mean(worst_neighbor)

            # Recalculate errors
            self.error[u] /= 2
            self.error[worst_neighbor] /= 2
            self.error[n] = (self.error[u] + self.error[worst_neighbor]) / 2

        # Decrement all errors
        for i in range(self.size):
            self.error[i] *= self.beta
