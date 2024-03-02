import argparse,time,os
from rxnb.mol import Mole


def main():
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--working_dir',type=str,default='./working_dir/',help='working directory')
    parser.add_argument('--string',type=str,default='c1ccccc1',help='molecule string')
    parser.add_argument('--string_type',type=str,default='smiles',help='string type')
    parser.add_argument('--solvent',type=str,default='',help='solvent')
    parser.add_argument('--engine',type=str,default='g16',help='quantum chemical calculation engine')
    parser.add_argument('--optimizer',type=str,default='g16',help='optimizer for geometry optimization')
    parser.add_argument('--threads_num',type=int,default=8,help='threads number')
    parser.add_argument('--memory',type=str,default='8GB',help='memory')   

    args = parser.parse_args()
    working_dir = args.working_dir
    string = args.string
    string_type = args.string_type.lower()
    engine = args.engine.lower()
    optimizer = args.optimizer.lower()
    threads_num = args.threads_num
    mem = args.memory
    solvent = args.solvent
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    xtb_param = {'gfn':2,'thread':threads_num}
    g16_param = {'method':'b3lyp','basis':'def2svp','nproc':threads_num,'mem':mem}
    if string_type == 'smiles':
        mol = Mole(smiles=string,working_dir=working_dir,optimizer=optimizer,engine=engine,solvent=solvent,xtb_param=xtb_param,g16_param=g16_param)
    elif string_type == 'inchi':
        mol = Mole(inchi=string,working_dir=working_dir,optimizer=optimizer,engine=engine,solvent=solvent,xtb_param=xtb_param,g16_param=g16_param)
    gibbs = mol.run()
    end_time = time.time()
    print(f'[INFO] Gibbs = {gibbs:.2f} kcal/mol, time duration: {end_time-start_time:.2f} s')
if __name__ == '__main__':
    main()