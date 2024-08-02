import os

def write_ts_gjf(symbols,positions,file_path,fixed_inf_lst=[[11,18],[18,19]],
                 params={'method':'b3lyp','basis':'def2svp','nproc':'8','mem':'8GB'},
                 chrg=0,mult=1):
    part1_inf = [f'%nproc={params["nproc"]}',
                 f'%mem={params["mem"]}',
                 f'%chk={os.path.basename(file_path)[:-4]}.chk',
                 f"#p opt=modredundant freq external='./xtb.sh'",
                 '',
                 'Generate by RDKit',
                 '',
                 f'{chrg} {mult}']
    
    geom_inf = []
    for sym,xyz in zip(symbols,positions):
        geom_inf.append(f'{sym:5s} {xyz[0]:15f} {xyz[1]:15f} {xyz[2]:15f}')
    
    fiexed_inf = []

    for fixed_inf in fixed_inf_lst:
        fiexed_inf.append(f'B {fixed_inf[0]} {fixed_inf[1]} F')

    
    link2_inf = ['--link1--',
                 f'%oldchk={os.path.basename(file_path)[:-4]}.chk',
                 f'%mem={params["mem"]}',
                 f'%nproc={params["nproc"]}',
                 f"#p opt=(ts,readfc,noeigen,nofreeze,notrust,maxcycle=250) freq geom=allcheck guess=tcheck external='./xtb.sh'",
                 '',
                 '']
    
    tot_inf = part1_inf + geom_inf + [''] + fiexed_inf + [''] + link2_inf + [''] + ['']
    
    with open(file_path,'w') as fw:
        fw.writelines('\n'.join(tot_inf))