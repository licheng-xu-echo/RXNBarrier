import os
from rdkit import Chem
from rdkit.Chem import AllChem
from subprocess import run,PIPE
def mol2xyz(mol,xyz_file):
    '''
    将RDKit的分子对象Mol保存成.xyz分子文件

    Parameters
    ----------
    mol : TYPE
        RDKit的分子对象.
    xyz_file : string
        .xyz文件路径.

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

def get_opted_geom(xyz_trj):
    '''
    获取geomeTRIC优化轨迹文件的最后一帧，也即优化好的结构，保存在后缀是_opted.xyz的文件里

    Parameters
    ----------
    xyz_trj : string
        geomeTRIC优化轨迹文件.

    Returns
    -------
    string
        优化结构保存的文件名.

    '''
    with open(xyz_trj) as fr:
        lines = fr.readlines()
    natom = int(lines[0].strip())
    opted_geom_lst = lines[-(natom+2):]
    with open(xyz_trj[:-4]+'_opted.xyz','w') as fw:
        fw.writelines(''.join(opted_geom_lst))
    return xyz_trj[:-4]+'_opted.xyz'
HA2KJ = 2625.5002
CWD = os.getcwd()
class Mole():
    def __init__(self,inchi=None,smiles=None,working_dir='.',fn='temp.xyz',
                 chrg=0,mult=1,optimizer='G16',engine='G16',method={'method':'b3lyp',
                                                                    'basis':'defsvp',
                                                                    'solvent':''}):
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
            during the process. The default is 'temp.xyz'.
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
        method : dict
            Parameters of calculation methods, including choices for computational 
            accuracy and solvent model settings."
        Returns
        -------
        Molecule object.
        '''
        self.inchi = inchi
        self.smiles = smiles
        self.working_dir = working_dir
        self.fn = fn
        self.chrg = chrg
        self.mult = mult
        self.optimizer = optimizer
        self.engine = engine
        self.method = method
        os.chdir(self.working_dir)
    def gen_3d_geom(self):
        '''
        使用RDKit解析包含分子2D结构信息的InChI或者SMILES，从而构建2D分子图，
        然后调用ETKDG将分子图转化为3D结构，并用UFF力场进行简单优化

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
        AllChem.UFFOptimizeMolecule(mol)
        if flag == -1:
            raise ValueError('3D geometry generation failed')
        else:
            mol2xyz(mol,self.fn)
            self.rdkit_mol = mol

    def optimize(self):
        '''
        调用geomeTRIC优化器与xTB方法对分子结构进行优化

        Returns
        -------
        None.

        '''
        cmd = "geometric-optimize --engine ase --ase-class=xtb.ase.calculator.XTB "+\
            "--ase-kwargs='{\"method\":\"GFN%d-xTB\", \"charge\":%d, \"mult\":%d, \"solvent\":\"%s\"}' "%(self.method,self.chrg,self.mult,self.solvent)+\
            self.fn
        result = run(cmd,stdout=PIPE,stderr=PIPE,universal_newlines=True,
                     cwd=None,shell=True,executable='/bin/bash',check=False)
        opted_fn = get_opted_geom(self.fn[:-4]+'_optim.xyz')
        self.opted_fn = opted_fn
    def opt_freq(self):
        '''
        调用geomeTRIC优化器与xTB方法对分子结构进行优化，同时在优化结束后进行频率分析，
        生成包括Gibbs自由能在内的一系列热力学参数

        Returns
        -------
        None.

        '''
        cmd = "geometric-optimize --engine ase --ase-class=xtb.ase.calculator.XTB "+\
            "--ase-kwargs='{\"method\":\"GFN%d-xTB\", \"charge\":%d, \"mult\":%d, \"solvent\":\"%s\"}' --hessian first+last "%(self.method,self.chrg,self.mult,self.solvent)+\
            self.fn
        result = run(cmd,stdout=PIPE,stderr=PIPE,universal_newlines=True,
                     cwd=None,shell=True,executable='/bin/bash',check=False)
        opted_fn = get_opted_geom(self.fn[:-4]+'_optim.xyz')
        self.opted_fn = opted_fn

    def write_xcontrol(self):
        '''
        生成控制xTB计算的中间文件xcontrol

        Returns
        -------
        None.

        '''
        inform = ['$thermo',f'    temp={self.temperature}']
        with open('xcontrol','w') as fw:
            fw.writelines('\n'.join(inform))
    def read_free_energy_xtb(self,hess_out):
        '''
        读取xTB的hessian计算输出文件中的Gibbs自由能信息

        Parameters
        ----------
        hess_out : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        with open(hess_out,'r') as fr:
            lines = fr.readlines()
        for line in reversed(lines):
            if 'TOTAL FREE ENERGY' in line:
                gibbs = eval(line.strip().split()[-3])
                break
        return gibbs * HA2KJ  ## Hartree to kJ/mol
    
    def calc_free_energy_xtb(self):
        '''
        调用xTB进行hessian计算，从而获得热力学信息，并读取其中的Gibbs自由能

        Returns
        -------
        gibbs : TYPE
            DESCRIPTION.

        '''
        self.write_xcontrol()
        cmd = f"xtb {self.opted_fn} --gfn {self.method} --hess --gbsa {self.solvent} --input xcontrol > {self.opted_fn[:-4]}_hess.out"
        result = run(cmd,stdout=PIPE,stderr=PIPE,universal_newlines=True,
                     cwd=None,shell=True,executable='/bin/bash',check=False)
        gibbs = self.read_free_energy(f'{self.opted_fn[:-4]}_hess.out')
        return gibbs
    def remove_temp_files(self):
        '''
        清除计算过程中的临时文件

        Returns
        -------
        None.

        '''
        temp_files = ['g98.out','hessian','vibspectrum','wbo','xcontrol','xtbrestart','xtbtopo.mol','charges']
        for file in temp_files:
            if os.path.exists(file):
                os.remove(file)
    def read_free_energy(self):
        '''
        从geomeTRIC的日志文件中读取Gibbs自由能

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        with open(f'{self.fn[:-4]}.log','r') as fr:
            lines = fr.readlines()
        for line in reversed(lines):
            if 'Gibbs Free Energy' in line:
                gibbs = eval(line.strip().split()[-5])
                break
        return gibbs * HA2KJ # convert Hartree to kJ/mol
    def run(self):
        '''
        执行分子字符串到Gibbs自由能的端到端计算工作流的核心函数。工作流包括：
        1. 调用gen_3d_geom方法解析InChI或SMILES生成2D分子图，并得到分子的3D结构
        2. 调用opt_fre方法，对第一步中获得的分子3D结构进行优化，并对优化得到的结构
           进行频率分析，获得分子的热力学信息
        3. 调用read_free_energy方法从第二步的日志文件中读取Gibbs自由能信息

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        # 获得初始化3D几何结构
        print('[INFO] Generate 3D geometry...')
        self.gen_3d_geom()
        # 优化分子结构并进行频率计算
        print('[INFO] Optimize 3D geometry and frequency calculation...')
        self.opt_freq()
        # 读取分子吉布斯自由能
        self.gibbs = self.read_free_energy()
        print(f'[INFO] Gibbs free energy: {self.gibbs:.2f} kJ/mol')
        os.chdir(CWD)
        print('[INFO] Done\n')
        return self.gibbs