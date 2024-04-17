import os
from ase.optimize import BFGS
from ase.atoms import Atoms as ase_Atoms
from xtb.ase.calculator import XTB
from ase.calculators.emt import EMT
from ase.io import read as ase_read
DEFAULT_MAX_STEPS = 100_000_000
class ASEOptimizer(object):
    def __init__(self, atoms=None, filename=None, engine='xtb', engine_params={'method': 'GFN2-xTB'}, optimizer='bfgs'):
        assert atoms is not None or filename is not None, "Either atoms or filename must be provided."

        self.filename = filename
        if filename is not None:
            self.fn = os.path.basename(filename).split('.')[0]
        else:
            self.fn = 'tmp'
        self.engine = engine.lower()
        self.optimizer = optimizer.lower()
        if atoms is not None and type(atoms) is ase_Atoms:
            if self.engine == 'xtb':
                atoms.calc = XTB(method=engine_params['method'])
            else:
                atoms.calc = EMT()
            self.atoms = atoms
        else:
            atoms = ase_read(filename)
            if self.engine == 'xtb':
                atoms.calc = XTB(method=engine_params['method'])
            else:
                atoms.calc = EMT()
            self.atoms = atoms
    def set_optimizer(self, optimizer_params={'restart':None,'logfile':'-','trajectory':None,'append_trajectory':False,
                                            'maxstep':None,'master':None,'alpha':None}):
        if self.optimizer == 'bfgs':
            if optimizer_params['trajectory'] is None:
                optimizer_params['trajectory'] = f'{self.fn}.traj'
            opt = BFGS(self.atoms, logfile=optimizer_params['logfile'], trajectory=optimizer_params['trajectory'],
                    maxstep=optimizer_params['maxstep'], master=optimizer_params['master'], alpha=optimizer_params['alpha'])
            return opt
        
    def run(self, fmax=0.05, steps=DEFAULT_MAX_STEPS):
        try:
            opt = self.set_optimizer()
        except:
            raise ValueError("Optimizer not set.")
        opt.run(fmax=fmax,steps=steps)
    
    def save_xyz(self, filename=None):
        if filename is None:
            filename = f'{self.fn}.xyz'
        '''
        N
        information
        Atom  X   Y   Z   Fx  Fy  Fz  charge
        ...
        '''
        self.atoms.write(filename)