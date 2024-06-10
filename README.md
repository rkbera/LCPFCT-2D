# LCPFCT-2D
This is a 2D multi-fluid plasma code, based on the LCPFCT scheme
- It solves a coupled set of relativistic fluid-Maxwell equations - suitable for laser/beam-plasma interaction system
- A Test Particle Simulation (TPS) module is also added to acquire some knowledge of the kinetic properties of particles
- The code is written in Fortran 77. It is a serial code. We are working on parallelization

# DEPENDENCY 
-  "gfortran" is only needed to run this code. 


# EXECUTION
- save the code "LCPFCT-2D" and go to inside the directory "$cd LCPFCT-2D"
- run "$ make all", it will create the executable "lcpfct2D"
  
- To clean the executable and object files, run "make clean" on the same directory
  
# RUN SIMULATION
- create another directory where you would like to run the simulation, e.g. create "run_simulation" directory
- save "input_file" in the directory "run_simulation" directory. A test "input_file" for the beam-plasma system is given inside the "LCPFCT-2D" code
- run "$ /path/lcpfct2D  input_file"
- The data will be saved in two folders "PLASMA_DATA" and "TEST_PARTICLE_DATA"
- PLASMA_DATA saves the space-time data foe plasma and beam quantities
- TEST_PARTICLE_DATA saves the space-time data for test particles

# FOR Plotting data
- Python-based plotting files are also given inside the LCPFCT-2D/plot_routines/. 


  
