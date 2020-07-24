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
        fcomm = comm.py2f()
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
        return

    def initialize(self, **kwargs):
        """prepare runs;

        """
        famus.check_input()
        famus.fousurf()
        famus.rdcoils()
        famus.packdof(famus.globals.xdof)
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
    test.run(verbose=True)
    famus.diagnos()
    test.save()