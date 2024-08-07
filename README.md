# Long-term G-functions for bidirectional aquifer thermal energy storage system operation
[Trained KNN Models](https://drive.google.com/drive/folders/16ggPohiPyBKsgfQ8kYWGHRXsoD9BObes?usp=sharing) |
[COMSOL Benchmark Data](https://drive.google.com/drive/folders/1PxxW9ovcw0Zc_Yc7MPGFkF-T5g5LbLkR?usp=sharing)
## Requirements
- Anaconda (with Python installed)
```
pip install -r requirements.txt
```
- Matlab

Curve Fitting Toolbox

Symbolic Math Toolbox

## COMSOL Model
The benchmark COMSOL model can be checked from [the previous work](https://github.com/Kecheng-Chen/1d-hlm).

## Getting Started
* Download KNN models and place them in ``model/``.
* Download benchmark data and place them in ``validation_data/``.
* type ``matlab`` in Anaconda console.
* Open ``main.m``, change case study id, and run.
* Results will be saved in ``results/``.

## To Do
* Build the software as a fully Python package.

## Acknowledgments
* This repository and code have been developed and maintained by [Kecheng Chen](https://geomechanics.berkeley.edu/people/kecheng-chen/).

## References
<a id="1">[1]</a> 
Chen, K., Sun, X., Soga, K., Nico, P. S., and Dobson, P. F. (2024). 
Machine-learning-assisted long-term G functions for bidirectional aquifer thermal energy storage system operation. 
Energy, 301, 131638.
