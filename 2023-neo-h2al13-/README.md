# Nuclear-Electronic Orbital Quantum Dynamics of Plasmon-Driven H<sub>2</sub> Photodissociation

The folders here contain all necessary files to generate figures in the following paper:

- Li, T. E., Hammes-Schiffer, S. Nuclear-Electronic Orbital Quantum Dynamics of Plasmon-Driven H<sub>2</sub> Photodissociation. [J. Am. Chem. Soc. **2023**](https://doi.org/10.1021/jacs.3c04927)

## File structure

  1. **plotting/** There are several python scripts (ending with **.py**) in this folder. Enter this folder, try running each of them (by **python xxx.py**) to obtain all the figures in the manuscript and SI.

  2. **run_pbe_6-31gd/** Folder for simulations with PBE/6-31G(d) for both RT-NEO-TDDFT and conventional RT-Ehrenfest calculations.
     - **neo_density_dynamics/** Cube files for plottings in Fig. 3 and 4 of the manuscript

## Caution

The simulations are performed in a developer version of Q-Chem. The current release version of Q-Chem does not support BO-RT-NEO simulations. 
