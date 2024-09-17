# Goal and principle

We developed a MATLAB-based program called TrackParticle, designed to quantify vesicle trajectories along dendrites for independent analysis and to derive global behavioral patterns. 
The program processes data from a track datasheet by first measuring the "dendrite backbone," a curvilinear axis extending from the cell body to the dendritic tip. 
Each vesicle track is then projected onto this backbone to assess track properties. 

Subsequently, tracks are segmented into subcategories based on their movement: immobile, forward, and backward. Then, each subcategory can be quantify in terms of mean length, speed and duration.
The program facilitates visualization of individual dendrites and provides statistical analyses for entire datasets. 
Additionally, it supports the export of curvilinear speed data for forward and backward movements. 
The program is accessible at: https://github.com/SebastienSchaub/TrackParticleDendrite.

# How it works 
To run, the program requires :

- a datasheet(s) with sequentially all vesicle tracks obtained from other tool. The program will concatenate all datasheet in a folder as a single condition
- pixel size and timelapse unit,
- the dendrite direction to properly assign forward and backward movement,
 
 ## References

- TrackParticleDendrite has been published in ***to be announced***

- contact : **sebastien.schaub@imev-mer.fr**
