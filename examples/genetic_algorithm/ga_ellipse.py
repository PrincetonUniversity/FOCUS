import sys
sys.path.append('../../python')
from mpi4py import MPI
from famuspy import FAMUS
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import time
import numpy as np

# MPI_INIT
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# run FOCUS
if rank==0:
    print("##### Begin FOCUS run with {:d} CPUs. #####".format(size))
    master = True
else:
    master = False
test = FAMUS(comm=comm, extension='ga_ellipse', verbose=True)
test.initialize()
self = test

niter=1000; npop=100; CXPB=0.5; MUTPB=0.2;cross_rate=0.5; mutate_rate=0.2
import random
from deap import base
from deap import creator
from deap import tools
from deap import algorithms
# initialization
# set size, master uses full-size
master = False
if self.globals.myid == self.globals.master:
    master = True
IND_SIZE = self.globals.ncoils_total
# create types
creator.create("FitnessMin", base.Fitness, weights=(-1.0, ))
creator.create("Individual", list, fitness=creator.FitnessMin)

lower = -3
upper = 3
toolbox = base.Toolbox()
toolbox.register("attr_int", random.randint, lower, upper)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_int, n=IND_SIZE)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
# evaluation function
def evaluate(individual):
    self.comm.Barrier()
    dof = np.ravel(self.comm.allgather(individual))
    return self.func(dof), np.count_nonzero(dof)/len(dof)
def evaluate2(individual, quit=False):
    if master:
        dof = np.array(individual)
    else:
        dof = np.empty(len(individual), dtype=int)
    while True:
        self.comm.Barrier()
        self.comm.Bcast(dof, root=self.globals.master)
        quit = self.comm.bcast(quit, root=self.globals.master)
        fitness = self.func(dof)
        if master:
            return fitness,
        if quit:
            return fitness,

toolbox.register("mate", tools.cxUniformPartialyMatched, indpb=cross_rate)
toolbox.register("mutate", tools.mutUniformInt, low=lower, up=upper, indpb=mutate_rate)
toolbox.register("select", tools.selTournament, tournsize=10)
toolbox.register("evaluate", evaluate2)
# Statics
stats = tools.Statistics(key=lambda ind: ind.fitness.values)
stats.register("avg", np.mean)
stats.register("std", np.std)
stats.register("min", np.min)
stats.register("max", np.max)       
# initialization, mutation, crossover
pop = toolbox.population(n=npop)
evaluate2(pop[0])
if master:
    pop, logbook = algorithms.eaSimple(pop, toolbox, cxpb=CXPB, mutpb=MUTPB, ngen=niter, stats=stats, verbose=master)
    fits = [abs(np.sum(np.array(ind.fitness.values)*np.array(ind.fitness.weights))) for ind in pop]
    best = np.array(pop[np.argmin(fits)])
    print(evaluate2(best, quit=True))
else:
    best = np.empty(len(pop[0]), dtype=int)
self.comm.Bcast(best, root=self.globals.master)
self.unpacking(best)
self.save()
