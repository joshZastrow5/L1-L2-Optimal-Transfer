# L1-L2 Optimal Transfer Project

This work studies the process in determining the minimum-fuel optimal transfer trajectory
between two halo orbits under the Circular Restricted Three-Body Problem with a low-thrust
capable satellite. First the shooting method is explored with the co-state vector initialized from
previous works, and then the steepest descent philosophy is used to determine trajectories
without knowledge of co-states that converge to solutions that inform the shooting method for proper convergence. A continuation method is
considered to transition minimum-energy solutions to minimum-fuel solutions and vice-versa.

## Script Summaries

All scripts can be run independently of each other, just hit run in MATLAB. <br> 
I reccomend reading the [final write-up](Report/Optimal_Control_Final_Report.pdf) for the project to get aquainted with what these scripts do.

### Main Scripts
 - L1_to_L2_transfer_full.m <br>
    > This script presents the process from transfering from an L1 Halo orbit to an L2 Halo orbit. It begins with an initial guess for the total thrust and thrust direction trajectory, simple and constant values. A very bad guess! We first use the method of Steepest Descent which involves guess the entire time history of the control input. The script displays the iterations using Steepest Descent, and the final result is also not very good! But that's okay, because then it uses that result as the initial guess for the Shooting Method, which converges very well and quite quickly. Then the script converts the minimum-energy to a minimum fuel solution using the continuation method outlined in the final writeup. 

- L2_to_L1_transfer_full.m <br>
    > This script presents the process from transfer from an L2 Halo orbit to and L1 Halo orbit. It performs the same processes as L1_to_L2_transfer_full.m, simply switching the starting and final orbits. 

### Testing and Development Scripts
- L2_to_L1_minFuel_transferTests.m
    > This script was a confirmation of the functions for propagating the cislunar continous thrust simulation. A lot of this work, including the orbit and transfer parameters, is developed from references [1] and [4] in the final writeup (see Report\\). As always, the first step in doing anything is to redo what others have done. Both of those references presented their final optimal co-states, and this script verifies that my code can use those to obtain the correct trajectories. 

- minFuel_convertTo_minEnergy.m
    > One is often concerned with finding the minimum fuel trajectory. This is done usually by finding the minimum energy case and using a continuation method to convert to the minimum fuel. However, in my case this script tests the continuation method by going from minimum fuel (known solutions given in references [1] and [4]) to the minimum energy solutions. The script outputs the images showcasing the trajectories associated with the continuation. 

- steepest_descent.m
    > This script was used for testing the Steepest Descent method of generating the transfer trajectory. It allows for "convergence" given an initially terrible guess, something the Shooting Method has issues with. However, the "convergence" is not great. Which is why the results from this method are used to initialize the

- STM_troubleshooting
    > Oh like you get everything right on the first try? The state transition matrix for this system is additionally complicated due to the switching function that determines if the thrust is on or off. When the switching curve is at zero it creates a discontinuity that needs to be corrected for (all of this is outlined in the final writeup, you should check it out if you haven't). Anyway, the first many attempts at the Shooting Method did not work, but I used this to figure out where and what the issue was (it was an indexing issue of all thing, all that math and I messed up the order in which matrices were multiplied). 

- xf_continuation
    > This was the beginning of something I thought was neat but ultimately didn't lead anywhere without having to reformulate the problem. I was curious if a transfer from apolune to apolune was the best choice, or if hopping onto the target orbit somewhere else would be more beneficial. So I used the continuation method to move the final state further away from the original final position. It works for small deviations but I found that the final time being held constant bottlnecks anything too far away. I would need to work on the free-time problem (which would would be a good idea anyway). 
    
