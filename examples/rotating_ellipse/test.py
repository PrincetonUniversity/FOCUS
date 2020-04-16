import sys
sys.path.append('../../python')
from mpi4py import MPI
from focuspy import FOCUSpy
import focus
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
test = FOCUSpy(comm=comm, extension='ellipse', verbose=True)

maxiter = 100
focus.globals.cg_maxiter = maxiter
# standard execution
#test.run(verbose=True)
#sys.exit()

# customize optimizers
test.prepare() 
x0 = focus.globals.xdof
x_copy = np.copy(x0)
# if master:
#     print(np.min(x0), np.mean(x0), np.max(x0))

collections = ('newton-cg', 'cg', 'bfgs', 'l-bfgs-b', 'tnc', 'SLSQP', 'nelder-mead', 'powell') #'dogleg', 'trust-ncg', 'trust-constr', 'trust-krylov', 'trust-exact',
#collections = ('cg', 'bfgs', 'l-bfgs-b')
if master:
    f, ax = plt.subplots(2)
for method in collections:
    if master:
        print('----------------------{:}-----------------------'.format(method))
    focus.globals.iout = 0
    x0 = np.copy(x_copy)
    test.time = time.time()
    test.callback(x0)
    res = minimize(test.func, x0, method=method, jac=test.grad, callback=test.callback, options={'maxiter':maxiter})
    test.callback(res.x)
    if master:
        ax[0].semilogy(focus.globals.evolution[:,1], label=method)
        ax[1].semilogy(focus.globals.evolution[:,0], focus.globals.evolution[:,1], label=method)

if master:
    ax[0].set_ylabel('target function')
    ax[0].set_xlabel('iteration')
    ax[1].set_ylabel('target function')
    ax[1].set_xlabel('time [s]')
    ax[0].legend(loc=1)
    ax[1].legend(loc=1)
    #plt.show()
    plt.savefig('scipy_optimizer_comp.png')
sys.exit()

