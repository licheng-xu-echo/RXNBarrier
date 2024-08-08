import numpy as np
from .read import read_xyz
from subprocess import run,PIPE
from .write import write_mod_gjf,write_gjf
from molop import AutoParser,molopconfig
molopconfig.quiet()

def perform_xtbts_calc(symbols,init_coords,charge,multiplicity,fixed_inf_lst,
                    mod_opt_gjf_file="1-modopt.gjf",mod_opt_log_file="1-modopt.log",
                    mod_freq_gjf_file="2-modfreq.gjf",mod_freq_log_file="2-modfreq.log",
                    ts_opt_gjf_file="3-tsopt.gjf",ts_opt_log_file="3-tsopt.log",
                    ts_freq_gjf_file="4-tsfreq.gjf",ts_freq_log_file="4-tsfreq.log"):
    # Step 1
    print("[INFO] (1/4) Start performing constrained optimization...")
    write_mod_gjf(symbols,init_coords,mod_opt_gjf_file,
                fixed_inf_lst=fixed_inf_lst,
                    kwd_line="#p opt=(modredundant,nomicro) external='./xtb.sh'",
                    chrg=charge,mult=multiplicity)           ## 生成固定优化的计算文件

    mod_opt_cmd = f'g16 {mod_opt_gjf_file} {mod_opt_log_file}'
    run(mod_opt_cmd,stdout=PIPE,stderr=PIPE,universal_newlines=True,cwd=None,shell=True,executable='/bin/bash',check=False) ## 执行固定优化
    print("[INFO] (1/4) Constrained optimization finished")
    
    mod_opted_mol = AutoParser(mod_opt_log_file)[0][-1]
    if not mod_opted_mol.is_optimized:
        print("[ERROR] Constrained optimization failed")
        return
    
    # Step 2
    print("[INFO] (2/4) Start performing frequency calculation...")
    mod_opted_mol.to_XYZ_file('./modopted.xyz')
    symbols,mod_opted_coords = read_xyz('./modopted.xyz')
    write_gjf(symbols,mod_opted_coords,mod_freq_gjf_file,
                kwd_line="#p freq external='./xtb.sh'",
                chrg=charge,mult=multiplicity)   ## 生成频率的计算文件
    mod_freq_cmd = f'g16 {mod_freq_gjf_file} {mod_freq_log_file}'
    run(mod_freq_cmd,stdout=PIPE,stderr=PIPE,universal_newlines=True,cwd=None,shell=True,executable='/bin/bash',check=False)  ## 执行频率计算
    print("[INFO] (2/4) Frequency calculation finished")

    mod_freq_mol = AutoParser(mod_freq_log_file)[0][-1]
    if mod_freq_mol.is_error:
        print("[ERROR] Frequency calculation failed")
        return
    mod_freqs = np.array(mod_freq_mol.vibrations.frequencies)
    mod_first_freq = mod_freqs[0]
    mod_img_freq_num = np.sum(mod_freqs < 0)
    print(f"[INFO] The first frequency is {mod_first_freq} cm-1, the number of imaginary frequencies is {mod_img_freq_num}")

    # Step 3
    print("[INFO] (3/4) Start performing TS optimization...")
    write_gjf(symbols,mod_opted_coords,ts_opt_gjf_file,
                kwd_line="#p opt(ts,calcfc,noeigen,nomicro) external='./xtb.sh'",
                chrg=charge,mult=multiplicity)   ## 生成TS优化的计算文件
    ts_opt_cmd = f'g16 {ts_opt_gjf_file} {ts_opt_log_file}'
    run(ts_opt_cmd,stdout=PIPE,stderr=PIPE,universal_newlines=True,cwd=None,shell=True,executable='/bin/bash',check=False)  ## 执行TS优化计算
    print("[INFO] (3/4) TS optimization finished")
    ts_opted_mol = AutoParser(ts_opt_log_file)[0][-1]
    if not ts_opted_mol.is_optimized:
        print("[ERROR] TS optimization failed")
        return
    # Step 4
    print("[INFO] (4/4) Start performing frequency calculation...")
    ts_opted_mol.to_XYZ_file('./tsopted.xyz')
    symbols,ts_opted_coords = read_xyz('./tsopted.xyz')
    write_gjf(symbols,ts_opted_coords,ts_freq_gjf_file,
                kwd_line="#p freq external='./xtb.sh'",
                chrg=charge,mult=multiplicity)   ## 生成频率的计算文件
    ts_freq_cmd = f'g16 {ts_freq_gjf_file} {ts_freq_log_file}'
    run(ts_freq_cmd,stdout=PIPE,stderr=PIPE,universal_newlines=True,cwd=None,shell=True,executable='/bin/bash',check=False)  ## 执行频率计算
    print("[INFO] (4/4) Frequency calculation finished")

    ts_freq_mol = AutoParser(ts_freq_log_file)[0][-1]
    if ts_freq_mol.is_error:
        print("[ERROR] Frequency calculation failed")
        return
    
    ts_freqs = np.array(ts_freq_mol.vibrations.frequencies)
    ts_first_freq = ts_freqs[0]
    ts_img_freq_num = np.sum(ts_freqs < 0)
    print(f"[INFO] The first frequency is {ts_first_freq} cm-1, the number of imaginary frequencies is {ts_img_freq_num}")
