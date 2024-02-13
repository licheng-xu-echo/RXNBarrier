# RXNBarrier

## Introduction 
This repository hosts a suite of scripts designed for comprehensive end-to-end reaction barrier calculations.

The code in this repository can currently be used to complete the end-to-end calculation of reaction thermodynamic free energy barrier. In the near feature, we will also add the ability to perform calculations of reaction kinetic free energy barrier.

The quantum calculation engine within can be selected to use either G16 or xTB, and the structure optimizer can be chosen from either G16 or geomeTRIC. Thus, users can set the calculation engine and optimizer based on their own conditions and requirements.

| engine | optimizer | requirement |
|---------|---------|---------|
| xTB | geomeTRIC | free and fast |
| xTB | G16 | fast and accurate |
| G16 | G16 | high accurate |

## Dependencies
- Python 3.8.18
- pandas 2.0.3
- RDKit 2022.9.5
- xtb 6.5.0
- xtb-python 22.1
- geometric 1.0.1
- ase 3.22.1
- Gaussian 16 (optional)