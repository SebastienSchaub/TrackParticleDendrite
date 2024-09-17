# Goal and principle

We developped in Matlab a Trackparticle, a program quantifying vesicle tracks along a dendrite to analyze them independently and extract a global behavior. 
From a track datasheet, the program will measure the "dendrite backbone", a curvilinear axis oriented from the cell body to the tip.
Then each track is projected in the backbone to measure the track properties.
Then each track is segmented in subpart based on their movement : immobile, forward and backward. 
The program provides visualisation of both the dendrite and the statistics for a dendrite set. 
Currently export of the the curvilinear speed in forward/backward.
The program is available here : https://github.com/SebastienSchaub/TrackParticleDendrite

# How it works 
To run, the program requires :

- a datasheet with sequentially all vesicle tracks obtained from other tool,
- pixel size and timelapse unit,
- the dendrite direction to properly assign forward and backward movement,
 
 ## References

- TrackParticleDendrite has been published in ***to be announced***

- contact : **sebastien.schaub@imev-mer.fr**
