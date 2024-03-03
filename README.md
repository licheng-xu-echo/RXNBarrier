# RXNBarrier

## Introduction

This repository hosts a suite of scripts designed for comprehensive end-to-end reaction barrier calculations.

The code in this repository can currently be used to complete the end-to-end calculation of reaction thermodynamic free energy barrier. In the near feature, we will also add the ability to perform calculations of reaction kinetic free energy barrier.

The quantum calculation engine within can be selected to use either G16[<sup>1</sup>](#refer-anchor-1) or xTB[<sup>2</sup>](#refer-anchor-2), and the structure optimizer can be chosen from either G16 or geomeTRIC[<sup>3</sup>](#refer-anchor-3). Thus, users can set the calculation engine and optimizer based on their own conditions and requirements.

| Engine | Optimizer | Features                                 |
| ------ | --------- | ---------------------------------------- |
| xTB    | geomeTRIC | free                                     |
| xTB    | G16       | fast[<sup>4</sup>](#refer-anchor-4) |
| G16    | G16       | accurate                                 |

## Dependencies

The following software packages are required to run the code in this repository:

- Python 3.8.18
- pandas 2.0.3
- RDKit 2022.9.5
- xtb 6.5.0
- xtb-python 22.1
- geometric 1.0.1
- ase 3.22.1
- Gaussian 16 (optional)

You can run commands below to install these packages
```
pip install rdkit-pypi ase pandas geometric
conda install xtb-python -c conda-forge
```

## Installation

Execute the command within the project directory ('RXNBarrier' folder), to install `rxnb` as a package:

```
pip install .
```

## Usage

### Thermodynamic free energy barrier calculations
For reaction thermodynamic free energy barrier calculations, users can run the following command Eg. 1 to Eg. 4 in `scripts` folder:

#### Example 1

Use **xTB** as the calculation engine and **G16** as the structure optimizer, input string is **SMILES**:

```
python thermo_dG_calc.py --working_dir ./eg_1 --string_type smiles --rct_inf_lst "[('Cc(c(F)cc1)c2c1cc([C@@H](NC(OC(C)(C)C)=O)C)c(Cl)n2',1),('CN(C(CN1CCNCC1)=O)C',1)]" --pdt_inf_lst "[('Cc(c(F)cc1)c2c1cc([C@@H](NC(OC(C)(C)C)=O)C)c(N3CCN(CC3)CC(N(C)C)=O)n2',1),('[H]Cl',1)]" --engine xtb --optimizer g16 --solvent water --threads_num 8 --memory 8GB
```

#### Example 2

Use **G16** as the calculation engine and **G16** as the structure optimizer, input string is **SMILES**:

```
python thermo_dG_calc.py --working_dir ./eg_2 --string_type smiles --rct_inf_lst "[('C',1),('ClCl',1)]" --pdt_inf_lst "[('CCl',1),('Cl',1)]" --engine g16 --optimizer g16 --solvent water --threads_num 8 --memory 8GB
```

#### Example 3

Use **xTB** as the calculation engine and **geomeTRIC** as the structure optimizer, input string is **InChI**:

```
python thermo_dG_calc.py --working_dir ./eg_3 --string_type inchi --rct_inf_lst "[('InChI=1S/C14H12O2/c15-13-8-12(9-14(16)10-13)7-6-11-4-2-1-3-5-11/h1-10,15-16H',1),('InChI=1S/O2/c1-2',1/2)]" --pdt_inf_lst "[('InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H',1)]" --engine xtb --optimizer geometric --threads_num 8
```

#### Example 4

To calculate the Gibbs free energy of a single molecule, execute the following command:

```
python G_calc.py --working_dir ./eg_4 --string c1ccccc1 --string_type smiles --engine xtb --optimizer g16 --solvent water --threads_num 8 --memory 8GB
```

### G16 log file parser

#### Example 5

```python
from rxnb.analyse import G16Log_Parser
log_file = './example/test-ts.log'
g16parser = G16Log_Parser(log_file,types=['sp','freq','opt'])
# Check if the stationary point is valid
print(g16parser.is_valid_statpt(key_bond=[1,23],tsk_type='ts',th=1.0))
# Read single point energy
print(g16parser.read_energy())
# Read the optimized geometry
atom_types, coords = g16parser.read_opted_geom()
# Read the thermal correction information
TCG, TCH = g16parser.read_therm_inf()
```

<div id="refer-anchor-1"></div>

- [1] https://gaussian.com/ (accessed Feb 13, 2024)

<div id="refer-anchor-2"></div>

- [2] C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme WIREs Comput. Mol. Sci., 2020, 11, e01493. DOI: 10.1002/wcms.1493

<div id="refer-anchor-3"></div>

- [3] Wang, L.-P.; Song, C.C. "Geometry optimization made simple with translation and rotation coordinates", J. Chem, Phys. 2016, 144, 214108. http://dx.doi.org/10.1063/1.4952956

<div id="refer-anchor-4"></div>

- [4] Lu T., gau_xtb: A Gaussian interface for xtb code, http://sobereva.com/soft/gau_xtb (accessed Feb 13, 2024)
