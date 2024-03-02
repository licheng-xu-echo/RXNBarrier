import os
from rdkit import Chem
from rdkit.Chem import AllChem
from subprocess import run,PIPE
from .utils import gen_xtb_sh

HA2KCAL = 627.503
CWD = os.getcwd()

def mol2xyz(mol,xyz_file):
    '''
    Save the RDKit molecule object Mol as a .xyz molecular file.

    Parameters
    ----------
    mol : TYPE
        RDKit Mol object.
    xyz_file : string
        .xyz file path.

    Returns
    -------
    None.

    '''
    symbols = [at.GetSymbol() for at in mol.GetAtoms()]
    positions = mol.GetConformer().GetPositions()
    xyz_inform = [f'{len(symbols)}','From RDKit']
    for sym,xyz in zip(symbols,positions):
        xyz_inform.append(f'{sym:5s} {xyz[0]:15f} {xyz[1]:15f} {xyz[2]:15f}')
    with open(xyz_file,'w') as fw:
        fw.writelines('\n'.join(xyz_inform))

def mol2gjf(mol,gjf_file,solvent='',chrg=0,mult=1,g16_param={'method':'b3lyp','basis':'def2svp',
                                    'nproc':'8','mem':'8GB'},engine='g16'):
    '''
    Save the RDKit molecule object Mol as a .gjf molecular file.

    Parameters
    ----------
    mol : TYPE
        RDKit Mol object.
    gjf_file : string
        .gjf file path.
    g16_param: Dict
        G16 parameters
    Returns
    -------
    None.

    '''
    symbols = [at.GetSymbol() for at in mol.GetAtoms()]
    positions = mol.GetConformer().GetPositions()
    g16_kwd = f'#p opt freq {g16_param["method"]}/{g16_param["basis"]}'
    if solvent != '':
        g16_kwd += f' scrf=(smd,solvent={solvent})'
    if engine == 'xtb':
        g16_kwd += " external='./xtb.sh'"
    gjf_inform = [f'%nproc={g16_param["nproc"]}',f'%mem={g16_param["mem"]}',g16_kwd,'','Generate by RDKit','',f'{chrg} {mult}']
    for sym,xyz in zip(symbols,positions):
        gjf_inform.append(f'{sym:5s} {xyz[0]:15f} {xyz[1]:15f} {xyz[2]:15f}')
    gjf_inform.append('')
    gjf_inform.append('')
    with open(gjf_file,'w') as fw:
        fw.writelines('\n'.join(gjf_inform))

class Mole():
    def __init__(self,inchi=None,smiles=None,working_dir='.',fn='temp',chrg=0,mult=1,optimizer='G16',engine='G16',solvent='',
                 xtb_param={'gfn':2,'thread':8},g16_param={'method':'b3lyp','basis':'def2svp','nproc':8,'mem':'8GB'}):
        '''
        Molecule object, used for generating and optimizing three-dimensional
        molecular structures as well as for energy calculations.

        Parameters
        ----------
        inchi : string
            InChI string of the molecule. InChI and SMILES cannot both be None.
        smiles : string
            SMILES string of the molecule. InChI and SMILES cannot both be None.
        working_dir : string
            The working directory for quantum chemical calculations. The default 
            is the current directory.
        fn : string
            The filename for the three-dimensional molecular structure generated 
            during the process. The default is 'temp'.
        chrg : int
            The charge of molecule. The default value is 0.
        mult : int
            The spin multiplicity of molecule. The default value is 1.
        optimizer : string
            The optimizer for molecular structure, with options including G16 and 
            geomeTRIC. The default is G16.
        engine : string
            The quantum chemistry computation engine, with options for the
            calculation engine being G16 and xTB. The default is G16.
        solvent : string
            The solvent setting.
        xtb_param : Dict
            Parameters for xTB methods
        g16_param : Dict
            Parameters for G16 methods
        Returns
        -------
        Molecule object.
        '''
        self.inchi = inchi
        self.smiles = smiles
        self.working_dir = working_dir
        self.fn = f'{fn}.xyz' if optimizer.lower() == 'geometric' else f'{fn}.gjf'
        self.chrg = chrg
        self.mult = mult
        self.optimizer = optimizer.lower()
        self.engine = engine.lower()
        self.solvent = solvent
        self.xtb_param = xtb_param
        self.g16_param = g16_param
        os.chdir(self.working_dir)
    
    def gen_3d_geom(self,ff_opt=True):
        '''
        Utilize RDKit to parse InChI or SMILES containing molecular 2D structural information, 
        thereby constructing a 2D molecular graph. Subsequently, invoke ETKDG to convert the 
        molecular graph into a 3D structure, and perform a simple optimization using the UFF force field.

        Parameters
        ----------
        ff_opt : bool
            Whether to optimize the molecule using a molecular force field.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        if self.inchi is not None:
            mol = Chem.MolFromInchi(self.inchi)
            print(f'[INFO] {self.inchi}')
        elif self.smiles is not None:
            mol = Chem.MolFromSmiles(self.smiles)
            print(f'[INFO] SMILES: {self.smiles}')
        else:
            raise ValueError('No input structure')
        mol = AllChem.AddHs(mol)
        flag = AllChem.EmbedMolecule(mol)
        if ff_opt:
            AllChem.UFFOptimizeMolecule(mol)
        if flag == -1:
            raise ValueError('3D geometry generation failed')
        else:
            if self.optimizer.lower() == 'geometric':
                mol2xyz(mol,self.fn)
            elif self.optimizer.lower() == 'g16':
                mol2gjf(mol,gjf_file=self.fn,solvent=self.solvent,chrg=self.chrg,mult=self.mult,g16_param=self.g16_param,engine=self.engine)
            self.rdkit_mol = mol

    def opt_freq(self):
        '''
        Invoke quantum chemical computation software to further optimize the molecular structure, 
        and conduct frequency analysis upon completion of the optimization. This allows for the 
        calculation of a series of thermodynamic parameters, including Gibbs free energy.

        Returns
        -------
        None.

        '''
        if self.engine.lower() == 'xtb' and self.optimizer == 'geometric':
            cmd = "geometric-optimize --engine ase --ase-class=xtb.ase.calculator.XTB "+\
                "--ase-kwargs='{\"method\":\"GFN%d-xTB\", \"charge\":%d, \"mult\":%d, \"solvent\":\"%s\"}' --hessian first+last "%(self.xtb_param['gfn'],self.chrg,self.mult,self.solvent)+self.fn
        elif self.engine.lower() == 'g16' and self.optimizer == 'g16':
            cmd = f"g16 {self.fn} {self.fn[:-4]}.log"
        elif self.engine.lower() == 'xtb' and self.optimizer == 'g16':
            gen_xtb_sh(xtb_sh=f'./xtb.sh',thread_num=self.xtb_param['thread'],gfn=self.xtb_param['gfn'])
            cmd = f"g16 {self.fn} {self.fn[:-4]}.log"
        result = run(cmd,stdout=PIPE,stderr=PIPE,universal_newlines=True,
                     cwd=None,shell=True,executable='/bin/bash',check=False)

    def check_log(self):
        
        '''
        Check whether the quantum calculation process has concluded normally. 

        Returns
        -------
        None.

        '''
        
        with open(f'{self.fn[:-4]}.log','r') as fr:
            lines = fr.readlines()
        if self.optimizer.lower() == 'geometric':
            for line in reversed(lines):
                if 'Converged! =D' in line:
                    return True
            return False
        elif self.optimizer.lower() == 'g16':
            if ' Normal termination of Gaussian 16' in lines[-1]:
                return True
            else:
                return False

    def read_free_energy(self):
        '''
        Read Gibbs free energy from log file

        Returns
        -------
        float
            Gibbs free energy in kcal/mol.

        '''
        with open(f'{self.fn[:-4]}.log','r') as fr:
            lines = fr.readlines()
        if self.optimizer.lower() == 'geometric':
            for line in reversed(lines):
                if 'Gibbs Free Energy' in line:
                    gibbs = eval(line.strip().split()[-5])
                    break
        elif self.optimizer.lower() == 'g16':
            for line in reversed(lines):
                if 'Sum of electronic and thermal Free Energies=' in line:
                    gibbs = eval(line.strip().split()[-1])
                    break
        return gibbs * HA2KCAL # convert Hartree to kcal/mol
    
    def run(self):
        '''
        The core function of the end-to-end computational workflow that transforms molecular strings 
        into Gibbs free energy. The workflow comprises:

        1. Invoking the gen_3d_geom method to parse InChI or SMILES to generate a 2D molecular graph 
        and obtain the 3D structure of the molecule.

        2. Calling the opt_freq method to optimize the 3D molecular structure obtained in the first 
        step, and then performing frequency analysis on the optimized structure to acquire the 
        thermodynamic information of the molecule.

        3. Utilizing the read_free_energy method to read the Gibbs free energy information 
        from the log files produced in the second step.

        Returns
        -------
        float
            Gibbs free energy.

        '''
        print('[INFO] Generate 3D geometry...')
        self.gen_3d_geom()
        print('[INFO] Optimize 3D geometry and frequency calculation...')
        self.opt_freq()
        calc_done = self.check_log()
        if calc_done:
            self.gibbs = self.read_free_energy()
            print(f'[INFO] Gibbs free energy: {self.gibbs:.2f} kcal/mol')
        else:
            self.gibbs = None
            print(f'[ERROR] Calculation is failed')
        os.chdir(CWD)
        print('[INFO] Done\n')
        return self.gibbs