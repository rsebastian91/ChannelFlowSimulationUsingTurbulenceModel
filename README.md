# ChannelFlowSimulationUsingTurbulenceModel
 
A simple FORTRAN code to perform Reynolds Averaged Navier Stokes simulation of a channel flow using a mixing-length turbulence model. 

Flow conditions in the simulation match the experiments of EL Telbany and A. J. Reynolds [1], and Anne Gilliot-Ottavy [2]. 

The simulation is performed with the finite volume method with non-uniform distribution of the mesh points. An algebraic Eddy-viscosity model using mixing length concept is used.

The simulation uses an iterative method for converging the solution. The results obtained from the simulation is then cross checked with the experimental results obtained for verifying the model.



![image](https://user-images.githubusercontent.com/84919918/161615086-e8824529-22fc-4ffb-a44d-2b6d44dab80e.png)

Fig. 1: Schematic of the flow.


![image](https://user-images.githubusercontent.com/84919918/161615385-e2ad3440-1933-4dec-b5de-9fb42fb33080.png)

Fig. 2: Pseudo-code of the code.


**References**

[1] EL Telbany and A. J. Reynolds, “Velocity distributions in plane turbulent channel flows”

[2] Anne Gilliot-Ottavy, “Characterization of turbulent flow using hot wire anemometer measurement on Poiseuille and Coutte-Poiseuille for the validation of turbulence models”
