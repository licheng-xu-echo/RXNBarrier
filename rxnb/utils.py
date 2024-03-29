import os

def gen_xtb_sh(xtb_sh='./xtb.sh',thread_num=28,gfn=0):
    inf = ['#!/bin/bash',
           '',
           '#This script was written by Dr. Tian Lu at Beijing Kein Research Center for Natural Sciences (www.keinsci.com)',
           '#Contact: sobereva@sina.com',
           'read atoms derivs charge spin < $2',
           '',
           '#Create temporary .xyz file',
           '#the element index should be replaced with element name, and the coordinate should be convert to Angstrom',
           'echo "Generating mol.tmp"',
           'cat >> mol.tmp <<EOF',
           '$atoms',
           '',
           '$(sed -n 2,$(($atoms+1))p < $2 | cut -c 1-72)',
           'EOF',
           '',
           'echo "Generating mol.xyz via genxyz"',
           'genxyz',
           'rm -f mol.tmp',
           '',
           'rm -f charges energy xtbrestart gradient hessian xtbout',
           'rm -f hessian xtb_normalmodes g98_canmode.out g98.out wbo xtbhess.coord',
           '',
           f'export OMP_NUM_THREADS={thread_num}',
           f'export MKL_NUM_THREADS={thread_num}',
           'uhf=`echo "$spin-1" | bc` #nalpha-nbeta',
           'if [ $derivs == "2" ] ; then',
           '    echo "Running: xtb mol.xyz --chrg $charge --uhf $uhf --gfn%d --hess --grad > xtbout"'%gfn,
           '    xtb mol.xyz --chrg $charge --uhf $uhf --gfn%d --hess --grad > xtbout'%gfn,
           'elif [ $derivs == "1" ] ; then',
           '    echo "Running: xtb mol.xyz --chrg $charge --uhf $uhf --gfn%d --grad > xtbout"'%gfn,
           '    xtb mol.xyz --chrg $charge --uhf $uhf --gfn%d --grad > xtbout'%gfn,
           'fi',
           'echo "xtb running finished!"',
           '',
           'echo "Extracting data from xtb outputs via extderi"',
           'extderi $3 $atoms $derivs',
           '',
           'rm -f charges energy xtbrestart gradient hessian xtbout mol.xyz tmpxx vibspectrum',
           'rm -f hessian xtb_normalmodes g98_canmode.out g98.out wbo xtbhess.coord .tmpxtbmodef',
           'rm -f .engrad *.engrad xtbtopo.mol xtbhess.xyz'
           ]
    with open(xtb_sh,'w') as fw:
        fw.writelines('\n'.join(inf))
    os.system(f'chmod +x {xtb_sh}')