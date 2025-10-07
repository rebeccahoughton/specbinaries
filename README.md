## Estimating the orbital parameters of spectroscopic binaries
===========================================================

### Overview

This repository contains code to simulate and estimate selection effects in samples of spectroscopic binary (SB) stars. The approach is to (a) generate a large synthetic population of binary systems with astrophysical priors, (b) produce "fake observations" of these binaries by calculating the radial velocity they would have at a random position in their orbit, and then (c) compare simulated detections to the actual observed sample to quantify selection biases. Through this method, we can constrain the true binary fraction of stars that match the properies of our observed sample. 

For example: "We observe 100 stars, and conclude that 20 of them are in binary systems. This gives us an observed binary fraction of 20%. However, when we run a binary population synthesis model, we find that out of a fake test sample of 1000 binaries, we would only be able to actually identify 650 of them as binary systems (i.e. for the other 350, the radial velocity shift between epochs is too small for us to see with our observing setup). Therefore, our true binary fraction is 20% * (1000/650) = ~31%."

Population synthesis also allows us to test out different
* Semi-major axis distributions,
* Eccentricity distributions, 
* Mass-ratio distributions between the primary and secondary stars. 


### Directory setup
--------------

All of this software is written in python. The directory is structured as follows:

#### Data
* TableRV_nonSB_SB1.csv - the time between epochs, shift in radial velocities, and uncertainties.
* TableSB_MassAge.csv - The classification, masses, and ages of the objects in the sample. 

#### Code
* rv_binary_pop.py - python file containing the binary population synthesis model. 
* makebinariesrv.py - functions used in rv_binary_pop.py (required for use, don't change)
* binary_tree.py - create a 2D binary tree that can be used to statistically compare two 2D distributions (**Note: not needed for this work!**)
* btree_plot.py - plotting file for the binary tree (**Note: not needed for this work!**)


### How to use
--------------

The main file for population synthesis is rv_binary_pop.py. You will also need the data files saved in the same directory, as well as makebinariesrv.py. The best way to use these files is therefore to either git clone the repository into your own workspace, or **download this reposity as a zip folder** and extract the contents. You can then edit the rv_binary_pop.py file locally and run the code using 
```
python rv_binary_pop.py
```
from the command line on your machine. 

### Dependencies
----------------

To run this code, you will need python installed. You will also need the modules \texttt{numpy, matplotlib, pandas, astropy} and \texttt{scipy}. If you do not have any of these installed, you can use

```
pip install astropy
pip install numpy
```
etc. 




