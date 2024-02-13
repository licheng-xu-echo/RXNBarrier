import os,argparse,time
import pandas as pd
from mol import Mole

def main():
    '''
    Main function for executing the calculation of the reaction Gibbs free energy, which includes:

    1. Parsing parameters input from the command line;

    2. Calculating the sum of Gibbs free energies of the reactants;

    3. Calculating the sum of Gibbs free energies of the products;

    4. Determining the reaction Gibbs free energy based on the Gibbs free energies of both reactants and products.

    Returns
    -------
    None.

    '''
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--working_dir',type=str,default='./working_dir/',help='working directory')
    parser.add_argument('--rct_inf_lst',type=str,default="[('Cc(c(F)cc1)c2c1cc([C@@H](NC(OC(C)(C)C)=O)C)c(Cl)n2',1),('CN(C(CN1CCNCC1)=O)C',1)]",help='reactant information')
    parser.add_argument('--pdt_inf_lst',type=str,default="[('Cc(c(F)cc1)c2c1cc([C@@H](NC(OC(C)(C)C)=O)C)c(N3CCN(CC3)CC(N(C)C)=O)n2',1),('[H]Cl',1)]",help='product information')
    parser.add_argument('--string_type',type=str,default='smiles',help='string type')
    parser.add_argument('--solvent',type=str,default='',help='solvent')
    parser.add_argument('--engine',type=str,default='g16',help='quantum chemical calculation engine')
    parser.add_argument('--optimizer',type=str,default='g16',help='optimizer for geometry optimization')

    args = parser.parse_args()
    working_dir = args.working_dir
    rct_inf_lst = eval(args.rct_inf_lst)
    pdt_inf_lst = eval(args.pdt_inf_lst)
    string_type = args.string_type.lower()
    engine = args.engine.lower()
    optimizer = args.optimizer.lower()

    solvent = args.solvent
    assert string_type in ['smiles','inchi'], "only support SMILES and InChI string type"
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    rct_gibbs_sum = 0
    pdt_gibbs_sum = 0
    rct_gibbs_lst = []
    pdt_gibbs_lst = []
    for i,rct_inf in enumerate(rct_inf_lst):
        if string_type == 'smiles':
            mol = Mole(smiles=rct_inf[0],fn=f'rct_{i}.xyz',working_dir=working_dir,solvent=solvent,engine=engine,optimizer=optimizer)
        else:
            mol = Mole(inchi=rct_inf[0],fn=f'rct_{i}.xyz',working_dir=working_dir,solvent=solvent,engine=engine,optimizer=optimizer)
        rct_gibbs_sum += mol.run() * rct_inf[1]
        rct_gibbs_lst.append(mol.gibbs)
    for i,pdt_inf in enumerate(pdt_inf_lst):
        if string_type == 'smiles':
            mol = Mole(smiles=pdt_inf[0],fn=f'pdt_{i}.xyz',working_dir=working_dir,solvent=solvent,engine=engine,optimizer=optimizer)
        else:
            mol = Mole(inchi=pdt_inf[0],fn=f'pdt_{i}.xyz',working_dir=working_dir,solvent=solvent,engine=engine,optimizer=optimizer)
        pdt_gibbs_sum += mol.run() * pdt_inf[1]
        pdt_gibbs_lst.append(mol.gibbs)
    delta_gibbs = pdt_gibbs_sum - rct_gibbs_sum
    print(f'[INFO] Delta gibbs = {delta_gibbs:.2f} kcal/mol')
    max_len = max([len(rct_inf_lst),len(pdt_inf_lst)])
    gibbs_inf_dict = {'Reactant':[item[0] for item in rct_inf_lst]+['']*(max_len-len(rct_inf_lst)),
                      'Reactant Gibbs Free Energy (kcal/mol)': rct_gibbs_lst+['']*(max_len-len(rct_inf_lst)),
                      'Product':[item[0] for item in pdt_inf_lst]+['']*(max_len-len(pdt_inf_lst)),
                      'Product Gibbs Free Energy (kcal/mol)': pdt_gibbs_lst+['']*(max_len-len(pdt_inf_lst)),
                      'Delta Gibbs Free Energy (kcal/mol)':[delta_gibbs]+['']*(max_len-1)}
    gibbs_inf_df = pd.DataFrame.from_dict(gibbs_inf_dict)
    gibbs_inf_df.to_csv(f'{working_dir}/gibbs_inf.csv',index=False)
    end_time = time.time()
    print(f'[INFO] Results save into gibbs_inf.csv')
    print(f'[INFO] Time duration: {end_time-start_time:.2f} s')

if __name__ == '__main__':
    main()