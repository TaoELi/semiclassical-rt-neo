# Born-Oppenheimer Real-Time Nuclear-Electronic Orbital (BO-RT-NEO) Dynamics

The folders here contain all necessary files to generate figures in the following paper:

- Li, T. E., Hammes-Schiffer, S. Electronic Born-Oppenheimer Approximation in Nuclear-Electronic Orbital Dynamics. [arXiv:2301.04076. **2023**](https://arxiv.org/abs/2301.04076).

## File structure

  - **plotting/** There are several python scripts (ending with **.py**) in this folder. Enter this folder, try running each of them (by *python xxx.py*) to obtain all the figures in the manuscript and SI.

  - **benchmark_RTNEO/**: Input files and output data for the HCN molecule outside the cavity.

  - **benchmark_incav/**: Input files and output data for the HCN molecule inside the cavity.

  - **proton_transfer/**: Input files and output data for the MDA molecule proton transfer.

## Caution

The simulations are performed in a developer version of Q-Chem. The current release version of Q-Chem does not support BO-RT-NEO simulations. 
