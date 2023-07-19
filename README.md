# Simulation manager
This Python module manages the simulations produced by the code presented [here](https://github.com/LoreCip/LQG-WENO-Reconstruction).

The minimal structure of the simulation folder must be:

```
📦id#_m#_dx#_xMax#_tf#_r0#_a0#
 ┣ 📂outputs
 ┃ ┗ 📜output.h5
 ┗ 📜ParameterFile.par
 ```

where in place of # the value of the respective parameter is expected. The single simulation is loaded via a simple command:
```python
import simRead

simulation = simRead.Sim('path/to/id#_m#_dx#_xMax#_tf#_r0#_a0#')
```
The simulation object is a class object containing the relevant information.

It is also possible to collect many simulations is a single folder like 
```
📦simulations
 ┣📦id#_m#_dx#_xMax#_tf#_r0#_a0#
 ┃ ┣ 📂outputs
 ┃ ┃ ┗ 📜output.h5
 ┃ ┗ 📜ParameterFile.par
 ┗📦id##_m##_dx##_xMax##_tf##_r0##_a0##
   ┣ 📂outputs
   ┃ ┗ 📜output.h5
   ┗ 📜ParameterFile.par
```
and load the entire set using the command
```python
import simRead

simulation = simRead.Sims('path/to/simulations')
```

## Dependencies list

This package has been developed employing the following packages:

- numpy==1.25.1
- scipy==1.11.1
- matplotlib==3.7.2
- h5py==3.9.0
- p_tqdm==1.4.0

