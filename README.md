# object-trajectory-control
MATLAB Simulation for controlling magnetic object in human tissue by external magnetic field
#### trajectory_sim.m is the simulation code
#### traj_S & traj_circle is the trajectory generation code(for circle and "S" shape)
#### traj_carrot is the path following control code (close loop)
#### cyl_Mag is the permanent magnet model 
in cyl_Mag (input: object position,orien of magnet,Diameter of external magnet, iteration(i),D: distance from magnet to the centoid of space(0,0))
           (output: magnetic field and gradient (axial&radial direction (z,p)), magnetic force of object in x,y direction)
