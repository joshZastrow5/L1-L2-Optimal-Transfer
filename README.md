# L1-L2 Optimal Transfer Project

This work studies the process in determining the minimum-fuel optimal transfer trajectory
between two halo orbits under the Circular Restricted Three-Body Problem with a low-thrust
capable satellite. First the shooting method is explored with the co-state vector initialized from
previous works, and then the steepest descent philosophy is used to determine trajectories
without knowledge of co-states that converge to solutions that inform the shooting method for proper convergence. A continuation method is
considered to transition minimum-energy solutions to minimum-fuel solutions and vice-versa.

## Script Summaries

 - L1_to_L2_transfer_full.m <br>
    > This script presents the process from transfering from an L1 Halo orbit to an L2 Halo orbit. It begins with an initial guess for the total thrust and thrust direction trajectory, simple and constant values. A very bad guess! We first use the method of Steepest Descent which involves guess the entire time history of the control input. The script displays the iterations using Steepest Descent, and the final result is also not very good! But that's okay, because then it uses that result as the initial guess for the Shooting Method, which converges very well and quite quickly. Then the script converts the minimum-energy to a minimum fuel solution using the continuation method outlined in the final writeup. 

- L2_to_L1_transfer_full.m <br>
    > This script presents the process from transfer from an L2 Halo orbit to and L1 Halo orbit. It performs the same processes as L1_to_L2_transfer_full.m, simply switching the starting and final orbits. 

- L2_to_L1_minFuel_transferTests.m
    > This script was a confirmation of the functions for propagating the cislunar continous thrust simulation. A lot of this work, including the orbit and transfer parameters, is developed from references [1] and [4]. As always, the first step in doing anything is to redo what others have done. Both of those references presented their final optimal co-states, and this script verifies that my code can use those to obtain the correct trajectories. 

- 
