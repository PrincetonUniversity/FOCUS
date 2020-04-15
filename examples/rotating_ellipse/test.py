import sys
sys.path.append('/Users/czhu/Documents/Code/FOCUS/python')
from mpi4py import MPI
from focuspy import FOCUSpy

# MPI_INIT
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# run FOCUS
if rank==0:
    print("##### Begin FOCUS run with {:d} CPUs. #####".format(size))
test = FOCUSpy(comm=comm, extension='ellipse', verbose=True)
test.run(verbose=True)
sys.exit()

