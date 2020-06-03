#!/usr/bin/env python

from scipy.optimize import Bounds, minimize
import numpy
import ruamel_yaml as yaml

import jax.numpy as np
from jax.config import config
config.update("jax_enable_x64", True)
from jax import grad, jit
from jax import jacfwd, jacrev
from jax.ops import index, index_add, index_update

import sys
import os
from CoilPy import FourSurf, Coil

methods = ['nelder-mead', 'powell', 'cg', 'bfgs', 'newton-cg',
           'l-bfgs-b', 'tnc', 'cobyla', 'slsqp', 'trust-constr',
           'dogleg', 'trust-ncg', 'trust-exact', 'trust-krylov']

# read input parameters
def read_parameters():    
    # get commondline args
    name = sys.argv[-1]
    print("read input parameters from {:}".format(name))
    assert os.path.exists(name), "File not exist!"
    # read yaml into dicts
    with open(name, 'r') as stream:
        params = yaml.safe_load(stream)
    # check and set default values
    params.setdefault('plasma_file', None)
    params.setdefault('coil_file', None)
    params.setdefault('npol', 64)
    params.setdefault('ntor', 64)
    params.setdefault('optimizer', 'l-bfgs-b')
    params.setdefault('random_perturbation', False)
    params.setdefault('output_key', 'coil_trans')
    params.setdefault('final_coils', 'coils.'+params['output_key'])
    params.setdefault('plot_coils', False)
    params.setdefault('maxiter', 50)
    params.setdefault('coil', [{'x0':[0,0,0,0,0,0],
                           'vary':[True, True, True, True, True, True],
                           'lbound':[None, None, None, -3.141592653, -3.141592653, -3.141592653],
                           'ubound':[None, None, None,  3.141592653,  3.141592653,  3.141592653]}])
    return params

# read and parse plasma boundary
def read_plasma_boundary(filename, npol=64, ntor=64, theta0=0, theta1=2*np.pi, zeta0=0, zeta1=2*np.pi):
    plasma = FourSurf.read_focus_input(filename)
    _theta = np.linspace(theta0, theta1, npol, endpoint=False)
    _zeta = np.linspace(zeta0, zeta1, ntor, endpoint=False)
    _tv, _zv = np.meshgrid(_theta, _zeta, indexing='ij')
    plasma_data = {}
    plasma_data['npol'] = npol
    plasma_data['ntor'] = ntor
    plasma_data['x'], plasma_data['y'], plasma_data['z'], plasma_data['n'] = plasma.xyz(_tv, _zv, normal=True)
    plasma_data['nn'] = np.linalg.norm(plasma_data['n'], axis=1)
    plasma_data['n'] = plasma_data['n'] / plasma_data['nn'][:, np.newaxis]
    # read plasma Bn
    with open(filename, 'r') as f:
        line = f.readline() #skip one line
        line = f.readline()
        num = int(line.split()[0]) #harmonics number
        nfp = int(line.split()[1]) #number of field periodicity
        nbn = int(line.split()[2]) #number of Bn harmonics
        if nbn > 0:
            line = f.readline() #skip one line
            line = f.readline() #skip one line
            for i in range(num):
                line = f.readline()
            line = f.readline() #skip one line
            line = f.readline() #skip one line
            xm = []; xn = []; bnc = []; bns = []
            for i in range(nbn):
                line = f.readline()
                line_list = line.split()
                n = int(line_list[0])
                m = int(line_list[1])
                xm.append(m)
                xn.append(n)
                bnc.append(float(line_list[2]))
                bns.append(float(line_list[3]))
            _mtnz = numpy.matmul( numpy.reshape(xm, (-1,1)), numpy.reshape(_tv, (1,-1)) ) \
                  - numpy.matmul( numpy.reshape(xn, (-1,1)), numpy.reshape(_zv, (1,-1)) ) 
            _cos = numpy.cos(_mtnz)
            _sin = numpy.sin(_mtnz)
            plasma_data['plas_bn'] = numpy.ravel(numpy.matmul( numpy.reshape(bnc, (1,-1)), _cos ) \
                                            + numpy.matmul( numpy.reshape(bns, (1,-1)), _sin))
        else:
            plasma_data['plas_bn'] = numpy.zeros_like(plasma_data['nn'])
    return plasma_data

# read coil data
class RigidCoil(Coil):
    @classmethod
    def read_makegrid(cls, filename):
        return super().read_makegrid(filename)

    # initialize DOF settings from input params
    def initialize_dof(self, params):
        lncoils = len(params['coil'])
        vary = []
        dof = []
        x0 = []
        lb = []
        ub = []
        for i in range(self.num):
            if i < lncoils:
                j = i
            else:
                j = 0
            vary.append(params['coil'][j]['vary'])
            x0.append(params['coil'][j]['x0'])
            ub.append(params['coil'][j]['ubound'])
            lb.append(params['coil'][j]['lbound'])
        self.vary = numpy.ravel(numpy.array(vary, dtype=bool))
        self.x0 = numpy.ravel(numpy.array(x0, dtype=float))
        self.lb = numpy.ravel(lb)
        self.ub = numpy.ravel(ub)
        self.dof = numpy.copy(self.x0)
        self.dim = numpy.count_nonzero(self.vary)
        print("Number of free parameters: {:d}".format(self.dim))
        # in case lb and ub are provided in full list
        if len(self.lb) != self.dim:
            print("Lower bounds are clipped according to 'vary'.")
            self.lb = self.lb[self.vary]
        if len(self.ub) != self.dim:
            print("Upper bounds are clipped according to 'vary'.")
            self.ub = self.ub[self.vary]
        return

    # pack dof into 1D array
    def pack_dof(self):
        return self.dof[self.vary]

    # unpack dof
    def unpack(self, x):
        self.dof = numpy.copy(self.x0)
        self.dof[self.vary] = x
        return self.dof

    # pack xyz & currents into a dict
    def pack_coils(self):
        coilxyz = []
        currents = []
        for icoil in list(self):
            coilxyz.append(numpy.array([icoil.x, icoil.y, icoil.z]).T)
            currents.append(icoil.I)
        coil = {}
        coil['xyz'] = coilxyz
        coil['currents'] = currents
        return coil

    # overide coil data
    def update(self, x):
        ldof = numpy.reshape(self.unpack(x), (-1, 6))
        for inum, icoil in enumerate(list(self)):
            sx, sy, sz, alpha, beta, gamma = ldof[inum, :]
            xyz = np.asarray([icoil.x, icoil.y, icoil.z]).T
            new_xyz = rigid_trans(xyz, sx, sy, sz, alpha, beta, gamma)
            icoil.x = new_xyz[:, 0]
            icoil.y = new_xyz[:, 1]
            icoil.z = new_xyz[:, 2]            
        return 

# Obtain rotation matrix
def rotation_matrix(alpha=0.0, beta=0.0, gamma=0.0):
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    cb = np.cos(beta)
    sb = np.sin(beta)
    cc = np.cos(gamma)
    sc = np.sin(gamma)
    return np.array([[ca*cb, ca*sb*sc-sa*cc, ca*sb*cc+sa*sc],
                     [sa*cb, sa*sb*sc+ca*cc, sa*sb*cc-ca*sc],
                     [-sb  , cb*sc         , cb*cc         ]])

# Obtain new curvs after rigid transformation
def rigid_trans(xyz, sx=.0, sy=.0, sz=.0, alpha=0.0, beta=0.0, gamma=0.0):
    assert (xyz.shape)[1] == 3, "The dimension of xyz should be Nx3." 
    rot = rotation_matrix(alpha, beta, gamma)
    oxyz = np.mean(xyz, axis=0)
    new_xyz = np.matmul(xyz-oxyz, rot.T) + oxyz
    return  new_xyz + np.array([sx, sy, sz])

# Biot-Savart Law
@jit
def bfield(xyz, I, pos, **kwargs):
    """Calculate B field at an arbitrary point

    Arguments:
        pos {array-like} -- Cartesian coordinates for the evaluation point

    Returns:
        ndarray -- the calculated magnetic field vector
    """        
    u0_d_4pi = 1.0E-7
    assert (pos.shape)[1]  == 3
    assert (xyz.shape)[1] == 3, "The dimension of xyz should be Nx3." 
    Rvec = pos[:,np.newaxis,:] - xyz[np.newaxis,:,:]
    assert (Rvec.shape)[-1] == 3
    RR = np.linalg.norm(Rvec, axis=2)
    Riv = Rvec[:, :-1, :]
    Rfv = Rvec[:, 1:, :]
    Ri = RR[:,:-1]
    Rf = RR[:,1:]
    B = np.sum(np.cross(Riv, Rfv)*((Ri+Rf)/((Ri*Rf)*(Ri*Rf+np.sum(Riv*Rfv, axis=2))))[:, :, np.newaxis], axis=1)\
        *u0_d_4pi*I
    return B

# objective function
# both rotate and shift
@jit
def chi2b(dof, coil, plasma):
    dof = dof.reshape((-1, 6))
    assert (dof.shape)[0] == len(coil['xyz'])
    # calculate B
    B = np.zeros_like(plasma['n'])
    pos = np.asarray([plasma['x'], plasma['y'], plasma['z']]).T
    # unpack coils
    coilxyz = coil['xyz']
    currents = coil['currents']
    for inum in range(len(coilxyz)):
        sx, sy, sz, alpha, beta, gamma = dof[inum, :]
        nxyz = rigid_trans(coilxyz[inum], sx, sy, sz, alpha, beta, gamma)   
        B = index_add(B, index[:, 0:3], bfield(nxyz, currents[inum], pos))
    Bn = np.sum(B*plasma['n'], axis=1) - plasma['plas_bn']
    Bn_integral = np.sum(Bn*Bn*plasma['nn'])*2*np.pi**2/(plasma['npol']*plasma['ntor'])
    return Bn_integral
    #return np.max(np.abs(Bn))

@jit
def grad_chi2b(dof, coil, plasma):
    return grad(chi2b, argnums=[0])(dof, coil, plasma)

# wrapper for functions
def func(x, *args):
    dof = np.array(coils.unpack(x))
    return float(chi2b(dof, *args))

def jac(x, *args):
    dof = np.array(coils.unpack(x))
    return numpy.ravel(grad_chi2b(dof, *args))[coils.vary]

def output(x):
    global it, gcoil, gplasma
    it += 1
    # plot_coils(dof, it/params['maxiter'])
    print('  iter: {:4d}, function value: {:12.5E}'.format(it, func(x, gcoil, gplasma)))
    return

def getBn(dof, coil, plasma, filename='Bn.png'):
    dof = dof.reshape((-1, 6))
    # calculate B
    B = np.zeros_like(plasma['n'])
    pos = np.asarray([plasma['x'], plasma['y'], plasma['z']]).T
    # unpack coils
    coilxyz = coil['xyz']
    currents = coil['currents']
    for inum in range(len(coilxyz)):
        sx, sy, sz, alpha, beta, gamma = dof[inum, :]
        nxyz = rigid_trans(coilxyz[inum], sx, sy, sz, alpha, beta, gamma)   
        B = index_add(B, index[:, 0:3], bfield(nxyz, currents[inum], pos))
    Bn = np.sum(B*plasma['n'], axis=1) - plasma['plas_bn']
    return Bn

# compare two matries
def compare_Bn(init, final, filename='Bn_comparison.png', **kwargs):
    # used on cluster with X display
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(figsize=[14, 5], ncols=3, sharex=True, sharey=True)
    axes = numpy.atleast_1d(axes)
    kwargs['aspect'] = kwargs.get('aspect', 'auto')
    kwargs['cmap'] = kwargs.get('cmap','RdBu_r') # colormap
    kwargs['origin'] = kwargs.get('origin', 'lower') 
    kwargs['extent'] = kwargs.get('extent', [0, 2*numpy.pi, 0, 2*numpy.pi]) # extent
    vmin = numpy.min([init, final])
    vmax = numpy.max([init, final])
    title = ['initial Bn', 'optimized Bn', '|init| - |final|']
    data = [init, final, numpy.abs(init)-numpy.abs(final)]
    for i in range(3):
        ax = axes[i]
        pcm = ax.imshow(data[i], **kwargs)
        ax.set_ylabel(r'$\theta$', fontsize=14)
        ax.set_xlabel(r'$\phi$', fontsize=14)
        ax.set_title(title[i], fontsize=14)
        fig.colorbar(pcm, ax=ax) 
    plt.savefig(filename)
    return

# Main program
if __name__ == '__main__':
    params = read_parameters()
    gplasma = read_plasma_boundary(params['plasma_file'], npol=params['npol'], ntor=params['ntor'])
    coils = RigidCoil.read_makegrid(params['coil_file'])
    gcoil = coils.pack_coils()
    coils.initialize_dof(params)

    # optimize
    mydof = coils.pack_dof()
    if params['random_perturbation']:
        print('Initial random perturbation applied.')
        mydof += numpy.random.rand(len(mydof))/5
    print('Initial chi2b = {:12.5E} T^2m^2, |gradient| = {:12.5E}.'.format(func(mydof, gcoil, gplasma), np.linalg.norm(jac(mydof, gcoil, gplasma))))
    init_Bn = getBn(coils.unpack(mydof), gcoil, gplasma)
    print('Initial residual Bn: max={:12.5E}, L2-norm={:12.5E}'.format(numpy.max(numpy.abs(init_Bn)), numpy.linalg.norm(init_Bn)))
    bounds = Bounds(coils.lb, coils.ub)
    it = 0
    if params['optimizer'] in methods:
        res = minimize(func, mydof, args=(gcoil, gplasma), bounds=bounds, method=params['optimizer'], jac=jac, 
            callback=output, options={'maxiter':params['maxiter'], 'gtol': 1E-8})
        print('Optimization status: \n', res)
        xf = res.x
    else:
        xf = mydof
    final_Bn = getBn(coils.unpack(xf), gcoil, gplasma)
    print('Final residual Bn: max={:12.5E}, L2-norm={:12.5E}'.format(numpy.max(numpy.abs(final_Bn)), numpy.linalg.norm(final_Bn)))

    # update and plot coils
    if params['plot_coils']:
        from mayavi import mlab
        mlab.options.offscreen = True
        fig = mlab.figure(size=[1600, 1200], bgcolor=(1,1,1), fgcolor=(0,0,0))
        shape2d = (params['npol'], params['ntor'])
        coils.plot(engine='mayavi', color=(0.5, 0.5, 0.5))
        coils.update(xf)
        coils.plot(engine='mayavi', color=(0, 0, 1))
        mlab.mesh(numpy.reshape(gplasma['x'], shape2d), numpy.reshape(gplasma['y'], shape2d), numpy.reshape(gplasma['z'], shape2d),
                scalars=numpy.reshape(final_Bn, shape2d), colormap='coolwarm')
        mlab.colorbar()
        mlab.savefig(params['output_key']+'_coils.png')
    else:
        coils.update(xf)
    coils.save_makegrid(params['final_coils'])

    # plot Bn comparison
    dim = (params['npol'], params['ntor'])
    compare_Bn(init_Bn.reshape(dim), final_Bn.reshape(dim), params['output_key']+'_Bn.png')