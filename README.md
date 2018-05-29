# TLEED

Tensor LEED (Low Energy Electron Diffraction) Optimization Project
 
 This project uses a Tensor LEED physics simulation code as a black-box objective function for testing derivative-free 
 optimization codes. The simulation codes consist of a Tensor LEED (TLEED) program which should be considered a high-fidelity simulation of an
 electron diffraction experiment that computes theoretical IV curves that are in turn compared to experimental data. The second
 simulation code is called KLEED and is a kinematic approximation to TLEED.  This means that KLEED is can be considered a low-fidelity
 physcics based approximation to TLEED. We use it as a surrogate for TLEED in the optimization codes.
 
 Both codes take as input (x,y,z) coordinates of 14 atoms representing a material that is being used to represent the surface that we 
 wish to consider.  Each of the atoms also has a "chemical identity", which is represented by an integer. 
 
 If we optimize with respect to all the variables, this represents a so-called Mixed Variable Problem (MVP).


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

There are two sets of black box simulators (KLEED and TLEED) that need to be built.  The code
for everything these are in the src/ directory. The convention is to call the executables from the main project directory.

Each of the executables assumes that the required input data sets are in a directory called tleed_data.  In addition you need to have
work directories for each of KLEED/TLEED called respectively kwork000/twork000.  Again, these should already be there but if not you
will have to create them and populate them with the required input data files.


### Prerequisites

In order to run this you will need the following:
 ```
 Fortran compiler (e.g. gfortran, ifort, etc.)
 Matlab
 DFO optimization code (e.g. nomadm, bfo, etc.)
```

### Installing

There is a makefile included that should work with a Linux-type environment.
I have tested this out on Mac OS X 10.13 (High Sierra) but it should also work in earlier versions.
I have also tested this out on Linux (ubuntu), but not lately.

To build the executables go into the src directory and type
```
make 
```
This should create 4 executables: kleedfcn.exe, tleedfcn.exe, runkleed, runtleed.
cd back into the main directory and create symbolic links to the 4 executables if they do not already exist.
```
cd ..
ln -is src/tleedfcn.exe .
ln -is src/kleedfcn.exe .
ln -is src/runtleed.exe .
ln -is src/runkleed.exe .
```
## Running the tests

The two files - runtleed.exe/runkleed.exe - are standalone KLEED/TLEED executables. They can be used to test the black box functions.
From the project directory type 
```
./runkleed.exe
```
or
```
./runtleed.exe
```
If all goes well you should get some output to the screen for several test input values.

## Running the optimization code

There are two black-box executables that are called from Matlab: tleedfcn.exe/kleedfcn.exe. In the project directory there should
also be corresponding mfiles to call these black-box functions - tleedfcn.m/kleedfcn.m and tleedfcn2.m/kleedfcn2.m.
The only difference between these two sets is that the first set opimize with respect to the continuous variables only, while the
second set of functions optimizes with respect to both the continuous and categorical variables at once.

In addition, there are 2 matlab scripts to run a series of optimization runs over 10 trial initial points.
Beware that these can take a long time to run (especially if you're using TLEED.

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Juan C. Meza** - *Initial work* - [TLEED](https://github.com/jcmeza/TLEED)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Mark Abramson
* Aran Garcia-Leuke
* Zhengji Zhao
