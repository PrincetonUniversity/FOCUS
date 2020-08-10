import sys
import os.path
from mpi4py import MPI
import famus
import numpy as np
import time

class FAMUS(object):
    def __init__(self, comm=MPI.COMM_WORLD, extension=None, verbose=True):
        """Initialize FAMUS class
        Args:
            comm (MPI communicator): the comunicator assigned to famus. Default: MPI.COMM_WORLD
            extension (str): extention for FOCUS runs, identical command line argument. Default: None
            verbose (bool): if False, suppress the screen output and redirect to a file "tmp.focus_output". Default: True
        Returns:
            None
        """
        # Initialize MPI related settings
        self.comm = comm
        fcomm = self.comm.py2f()
        famus.globals.mpi_comm_famus = fcomm
        famus.globals.myid = comm.Get_rank()
        famus.globals.ncpu = comm.Get_size()
        # Read input namelist
        if extension is not None:
            inputFile = extension+'.input'
            assert os.path.exists(inputFile), "File not existed, please check again!"
            famus.globals.ext = extension
            famus.globals.inputfile = inputFile
            famus.read_namelist(inputFile)
        else:
            # assign a temp name and will not read the input file
            famus.globals.ext = 'focus_tmp'
        # check verbose status
        self.verbose = verbose
        self.dumb = False
        self.time = time.time()
        self.initialized = False
        if not self.verbose:
            famus.mute(1)
            self.dumb = True
        self.globals = famus.globals
        self.packdof = famus.packdof
        self.unpacking = famus.unpacking
        self.allocdata = famus.allocdata
        return

    def initialize(self, **kwargs):
        """prepare runs;

        """
        famus.check_input()
        famus.fousurf()
        famus.rdcoils()
        famus.packdof(famus.globals.xdof)
        self.allocdata(0)
        self.initialized = True
        return

    def run(self, verbose=False, **kwargs):
        """FAMUS run
        Args:
            verbose (bool): if False, suppress the screen output and 
                            redirect to a file "tmp.focus_output". Default: True
        Returns:
            None        
        """
        # re-check verbose status
        self.verbose = verbose
        if self.verbose:
            # cancel redirect
            if self.dumb:
                famus.mute(0) 
        else:
            # redirect
            if not self.dumb:
                famus.mute(1)
                self.dumb = True 
        # standard run
        if not self.initialized:
            self.initialize()
        famus.solvers()
        return

    def save(self, **kwargs):
        famus.saving()
        return

    def func(self, x):
        """FOCUS basic optimization functions
        """
        n = len(x)
        return famus.myvalue(x, n)
    
    def grad(self, x):
        """Gradient
        """
        n = len(x)
        g = np.zeros(n)
        famus.mygrad(g, x, n)
        return g
    
    def callback(self, x):
        famus.unpacking(x)
        famus.costfun(0)
        famus.output(time.time() - self.time)
        return famus.globals.chi

    def ga_opt(self, niter=100, npop=100, CXPB=0.5, MUTPB=0.2, cross_rate=0.5, mutate_rate=0.2, **kwargs):
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
        IND_SIZE = self.globals.ncoils
        # create types
        creator.create("FitnessMin", base.Fitness, weights=(-1E3, -1.0))
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
            while not quit:
                self.comm.Barrier()
                fitness = self.func(dof)
                if master:
                    return fitness
        toolbox.register("mate", tools.cxUniformPartialyMatched, indpb=cross_rate)
        toolbox.register("mutate", tools.mutUniformInt, low=lower, up=upper, indpb=mutate_rate)
        toolbox.register("select", tools.selTournament, tournsize=3)
        toolbox.register("evaluate", evaluate)
        # Statics
        stats = tools.Statistics(key=lambda ind: ind.fitness.values)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)       
        # initialization, mutation, crossover
        pop = toolbox.population(n=npop)
        pop, logbook = algorithms.eaSimple(pop, toolbox, cxpb=CXPB, mutpb=MUTPB, ngen=niter, 
                                   stats=stats, verbose=master)
        # for i in range(niter):
        #     # Select the next generation individuals
        #     offspring = toolbox.select(pop, len(pop))
        #     # Clone the selected individuals
        #     offspring = list(map(toolbox.clone, offspring))

        #     # Apply crossover on the offspring
        #     for child1, child2 in zip(offspring[::2], offspring[1::2]):
        #         if random.random() < CXPB:
        #             toolbox.mate(child1, child2)
        #             del child1.fitness.values
        #             del child2.fitness.values

        #     # Apply mutation on the offspring
        #     for mutant in offspring:
        #         if random.random() < MUTPB:
        #             toolbox.mutate(mutant)
        #             del mutant.fitness.values

        #     # Evaluate the individuals with an invalid fitness
        #     invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        #     fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        #     for ind, fit in zip(invalid_ind, fitnesses):
        #         ind.fitness.values = fit

        #     # The population is entirely replaced by the offspring
        #     pop[:] = offspring
        #     # print status
        #     fits = [ind.fitness.values[0] for ind in pop]
        #     if master:
        #         print(i, np.min(fits), np.max(fits), np.mean(fits))
        fits = [abs(np.sum(np.array(ind.fitness.values)*np.array(ind.fitness.weights))) for ind in pop]
        best = np.argmin(fits)
        dof = np.ravel(self.comm.allgather(best))
        self.unpacking(dof)
        self.save()
        return 

if __name__== "__main__":
    ext = sys.argv[1]
    if '.input' in ext:
        ind = ext.index
        ext = ext[:ind]
    print('Begin to run FAMUS from python with input file at ', ext+'.input')
    comm = MPI.COMM_WORLD
    test = FAMUS(comm=comm, extension=ext, verbose=True)
    #print('after __init__', test.globals.inputfile)
    test.initialize()
    test.ga_opt()
    famus.diagnos()
    test.save()