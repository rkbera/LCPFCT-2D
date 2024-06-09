# LCPFCT_2D
This is a 2D multi-fluid code based on the LCPFCT scheme
- It solves a coupled set of relativistic fluid-Maxwell equations - suitable for laser/beam-plasma interaction system
- A Test Particle Simulation (TPS) module is also added to acquire some knowledge of the kinetic properties of particles
- The code is written in Fortran 77. It is a serial code. We are working on parallelization

# DEPENDENCY 
-  "gfortran" is only needed to run this code. 


# EXECUTION
- save the code- "LCPFCT_2D" and go to inside the directory "$cd LCPFCT_2D"
- run "$ make all", it will create the executable "lcpfct2D"
- To clean the executable, run "make clean" on the same directory
  
# RUN SIMULATION
- create another directory "run_simulation"
- save "input_file" in the directory "run_simulation" directory. A test "input_file" for beam-plasma system is given inside the "LCPFCT_2D" code
- then run "$ /path-to-the-executable/lcpfct2D input_file"
- It will create data at different time steps. The data will be saved in two folders "PLASMA_DATA" and "TEST_PARTICLE_DATA"


  
