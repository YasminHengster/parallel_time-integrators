RIDC_simple.cpp
===========

To solve a simple ODE, use the RIDC_simple.cpp script.

In the script the following has to be defined:
RHS of the ODE func
Number of processes M
End Time T 
Initial condition y0
Number of timesteps N

To use the script for a different number of processes/ corrections, the corresponding matrix S and the vector beta has to be defined. 
For M=2-6 these are already written down and need to be commented/ uncommented (ll 35-49)

To compile the code on cirrus, the compile_now executable will load all the necessary modules (gcc and mpt) 

With the run_job.slurm file the program can be executed where the number of processes has to be defined.

The current version is for M=4 processes.



RIDC_particle.cpp
===========

To solve the particle interactions ODE, use the RIDC_particle.cpp script.

In the script the following has to be defined:

fx_posi, fv_posi, fx_nega, fv_nega: RHS of the ODE function for ions and electrons
N_posi: Number of ions
N_nega: Number of electrons
q_posi: Charge of ions
q_nega: Charge of electrons
m_posi: Mass of ions
m_nega: Mass of electrons
d: Distance between particles
M: Number of processes 
T: End Time 
y0: Initial condition 
N: Number of timesteps

Ohter settings are similar to RIDC_simple.cpp. Prediction results at time T will be saved in a txt file.



