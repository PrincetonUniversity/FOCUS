import sys
import os.path
from mpi4py import MPI
import focus
import numpy as np
import time

class FOCUSpy(object):
    def __init__(self, comm=MPI.COMM_WORLD, extension=None, verbose=True):
        """Initialize FOCUSpy class
        Args:
            comm (MPI communicator): the comunicator assigned to FOCUS. Default: MPI.COMM_WORLD
            extension (str): extention for FOCUS runs, identical command line argument. Default: None
            verbose (bool): if False, suppress the screen output and redirect to a file "tmp.focus_output". Default: True
        Returns:
            None
        """
        # Initialize MPI related settings
        fcomm = comm.py2f()
        focus.globals.mpi_comm_focus = fcomm
        focus.globals.myid = comm.Get_rank()
        focus.globals.ncpu = comm.Get_size()
        # Read input namelist
        if extension is not None:
            inputFile = extension+'.input'
            assert os.path.exists(inputFile), "File not existed, please check again!"
            focus.globals.ext = extension
            focus.globals.inputfile = inputFile
            focus.read_namelist(inputFile)
        else:
            # assign a temp name and will not read the input file
            focus.globals.ext = 'focus_tmp'
        # check verbose status
        self.verbose = verbose
        self.dumb = False
        self.time = time.time()
        self.globals = focus.globals
        if not self.verbose:
            focus.mute(1)
            self.dumb = True
        return

    def prepare(self, **kwargs):
        """prepare runs;

        """
        focus.check_input()
        focus.surface()
        focus.rdcoils()
        focus.packdof(focus.globals.xdof)
        focus.allocdata(abs(focus.globals.case_optimize))
        focus.diagnos()
        if focus.globals.isnormweight:
            focus.normweight()
        return

    def run(self, verbose=False, **kwargs):
        """standard run
        Args:
            verbose (bool): if False, suppress the screen output and redirect to a file "tmp.focus_output". Default: True
        Returns:
            None        
        """
        # re-check verbose status
        self.verbose = verbose
        if self.verbose:
            # cancel redirect
            if self.dumb:
                focus.mute(0) 
        else:
            # redirect
            if not self.dumb:
                focus.mute(1)
                self.dumb = True 
        # standard run
        focus.focus()
        return

    def save(self, **kwargs):
        famus.saving()
        return

    def func(self, x):
        """FOCUS basic optimization functions
        """
        n = len(x)
        return focus.myvalue(x, n)
    
    def grad(self, x):
        """Gradient
        """
        n = len(x)
        g = np.zeros(n)
        focus.mygrad(g, x, n)
        return g
    
    def callback(self, x):
        focus.unpacking(x)
        focus.costfun(0)
        focus.output(time.time() - self.time)
        return focus.globals.chi

if __name__== "__main__":
    ext = sys.argv[1]
    if '.input' in ext:
        ind = ext.index
        ext = ext[:ind]
    print('Begin to run FAMUS from python with input file at ', ext+'.input')
    comm = MPI.COMM_WORLD
    test = FOCUSpy(comm=comm, extension=ext, verbose=True)
    test.run(verbose=True)
