# Born-Oppenheimer Real-Time Nuclear-Electronic Orbital (BO-RT-NEO) Dynamics

The folders here contain all necessary files to generate figures in the following paper:

- Li, T. E., Hammes-Schiffer, S. Electronic Born-Oppenheimer Approximation in Nuclear-Electronic Orbital Dynamics. [J. Chem. Phys. **2023**, 158, 114118](https://doi.org/10.1063/5.0142007).

## File structure

  1. **plotting/** There are several python scripts (ending with **.py**) in this folder. Enter this folder, try running each of them (by **python xxx.py**) to obtain all the figures in the manuscript and SI.

  2. **benchmark_RTNEO/**: Input files and output data for the HCN molecule outside the cavity.

     * lr_tddft/ : Linear-response NEO-TDDFT calculation of HCN in the frequency domain.

     * rt_delta_proton/ : Real-time NEO-TDDFT calculation of HCN.
        * dt_xx/ : simulation when time step is xx.
            * HCN.in : Q-Chem input file.
            * mu_n.txt : protonic dipole moment as a function of time; each line contains [time step, time, mux, muy, muz] in atomic units.

     * rtBO_delta_proton/ : Born-Oppenheimer real-time NEO-TDDFT calculation of HCN; the file structure is similar as the above rt_delta_proton/.

  3. **benchmark_incav/**: Input files and output data for the HCN molecule inside the cavity. The file structure is the same as the above **benchmark_RTNEO/**.

  4. **proton_transfer/**: Input files and output data for the MDA molecule proton transfer. 

     * neo_boeh_3cen/ : BO-NEO-Ehrenfest simulation using three proton basis centers with dt = 4.2 au.
        * config.xyz : standard xyz molecular geometry as a function of time
        * config.txt : molecular geometry as a function of time; each line contains the atomic coordinates of the molecule in the order of [x1, y1, z1, x2, y2, ...]
        * e_tot.txt : total energy of the molecular system as a function of time; each line contains [time step, energy]
        * den_p_xxx_0.cube : standard cube file for the protonic density at time step xxx (only for neo_boeh_3cen/).

     * neo_boeh_3cen_smalldt/ : BO-NEO-Ehrenfest simulation using three proton basis centers with dt = 0.4 au.

     * neo_eh_3cen/ : NEO-Ehrenfest simulation using three proton basis centers with dt = 0.4 au.

     * eh_cl/ : Conventional Ehrenfest simulation with classical proton with dt = 0.4 au.

     * boeh_cl/ : Conventional BOMD simulation with classical proton with dt = 4.2 au.
    

## Caution

The simulations are performed in a developer version of Q-Chem. The current release version of Q-Chem does not support BO-RT-NEO simulations. 
