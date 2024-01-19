import os,argparse,time
import pandas as pd
from mol import Mole
from rdkit import Chem
from rdkit.Chem import AllChem
from subprocess import run,PIPE
HA2KJ = 2625.5002   # 能量单位从Hartree到kJ/mol的转化系数
CWD = os.getcwd()

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


def main():
    '''
    主函数，执行反应Gibbs自由能计算，包含：
    1. 解析命令行输入的参数；
    2. 计算反应物的Gibbs自由能之和；
    3. 计算产物的Gibbs自由能之和；
    4. 基于反应物与产物的Gibbs自由能计算反应的Gibbs自由能

    Returns
    -------
    None.

    '''
    # 获取命令行中输入的参数
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--working_dir',type=str,default='./working_dir/',help='working directory')
    parser.add_argument('--rct_inf_lst',type=str,default="[('Cc(c(F)cc1)c2c1cc([C@@H](NC(OC(C)(C)C)=O)C)c(Cl)n2',1),('CN(C(CN1CCNCC1)=O)C',1)]",help='reactant information')
    parser.add_argument('--pdt_inf_lst',type=str,default="[('Cc(c(F)cc1)c2c1cc([C@@H](NC(OC(C)(C)C)=O)C)c(N3CCN(CC3)CC(N(C)C)=O)n2',1),('[H]Cl',1)]",help='product information')
    parser.add_argument('--string_type',type=str,default='smiles',help='string type')
    parser.add_argument('--temperature',type=float,default=298.15,help='temperature')
    parser.add_argument('--solvent',type=str,default='water',help='solvent')

    args = parser.parse_args()
    working_dir = args.working_dir
    rct_inf_lst = eval(args.rct_inf_lst)
    pdt_inf_lst = eval(args.pdt_inf_lst)
    string_type = args.string_type.lower()
    temperature = args.temperature
    solvent = args.solvent
    assert string_type in ['smiles','inchi'], "only support SMILES and InChI string type"
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    rct_gibbs_sum = 0
    pdt_gibbs_sum = 0
    rct_gibbs_lst = []
    pdt_gibbs_lst = []
    # 计算反应物的Gibbs自由能
    for i,rct_inf in enumerate(rct_inf_lst):
        if string_type == 'smiles':
            mol = Mole(smiles=rct_inf[0],fn=f'rct_{i}.xyz',working_dir=working_dir,solvent=solvent,temperature=temperature)
        else:
            mol = Mole(inchi=rct_inf[0],fn=f'rct_{i}.xyz',working_dir=working_dir,solvent=solvent,temperature=temperature)
        rct_gibbs_sum += mol.run() * rct_inf[1]
        rct_gibbs_lst.append(mol.gibbs)
    # 计算产物的Gibbs自由能
    for i,pdt_inf in enumerate(pdt_inf_lst):
        if string_type == 'smiles':
            mol = Mole(smiles=pdt_inf[0],fn=f'pdt_{i}.xyz',working_dir=working_dir,solvent=solvent,temperature=temperature)
        else:
            mol = Mole(inchi=pdt_inf[0],fn=f'pdt_{i}.xyz',working_dir=working_dir,solvent=solvent,temperature=temperature)
        pdt_gibbs_sum += mol.run() * pdt_inf[1]
        pdt_gibbs_lst.append(mol.gibbs)
    # 计算反应的Gibbs自由能
    delta_gibbs = pdt_gibbs_sum - rct_gibbs_sum
    print(f'[INFO] Delta gibbs = {delta_gibbs:.2f} kJ/mol')
    max_len = max([len(rct_inf_lst),len(pdt_inf_lst)])
    gibbs_inf_dict = {'Reactant':[item[0] for item in rct_inf_lst]+['']*(max_len-len(rct_inf_lst)),
                      'Reactant Gibbs Free Energy (kJ/mol)': rct_gibbs_lst+['']*(max_len-len(rct_inf_lst)),
                      'Product':[item[0] for item in pdt_inf_lst]+['']*(max_len-len(pdt_inf_lst)),
                      'Product Gibbs Free Energy (kJ/mol)': pdt_gibbs_lst+['']*(max_len-len(pdt_inf_lst)),
                      'Delta Gibbs Free Energy (kJ/mol)':[delta_gibbs]+['']*(max_len-1)}
    gibbs_inf_df = pd.DataFrame.from_dict(gibbs_inf_dict)
    # 保存结果信息
    gibbs_inf_df.to_csv(f'{working_dir}/gibbs_inf.csv',index=False)
    end_time = time.time()
    print(f'[INFO] Results save into gibbs_inf.csv')
    print(f'[INFO] Time duration: {end_time-start_time:.2f} s')
if __name__ == '__main__':
    main()