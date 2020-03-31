from mpi4py import MPI
from focuspy import focus

# initial settings
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

color = rank % 2
if color==0:
    key = +rank
else :
    key = -rank

sub_comm = comm.Split(color, rank)
sub_rank = sub_comm.Get_rank()
sub_size = sub_comm.Get_size()

for i in range(size):
    comm.Barrier()
    if rank == i:
        print('Global: rank={:d} / {:d};  Local: rank={:d} / {:d}. color={:d}, key={:d}'.format(
            rank, size, sub_rank, sub_size, color, key))
'''
try:
    fcomm = comm.py2f()
    focus(fcomm)
except:
    print('Got errors!')
    pass
'''
path = '../examples/rotating_ellipse/ellipse'
if color == 0:
    scomm = sub_comm.py2f()
    focus(path, scomm)

exit
