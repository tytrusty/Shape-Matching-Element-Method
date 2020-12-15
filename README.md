# Shape Matching Virtual Element Method

For continuum mechanics simulation of boundary only representations!

vem_sim_2d is a minimal example that *should* only depend on gptoolbox... So see if that runs!

### External Dependencies  ###
1. [GPToolbox](https://github.com/alecjacobson/gptoolbox) 
2. OpenMP (OPTIONAL)

### Included Submodules (Installed Automatically) ###
1. Libigl https://github.com/libigl/libigl
2. Eigen >= 3.2 (uses the libigl Eigen install)

> **To get started:** Clone this repository and all its [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) dependencies using:
> 
>     git clone --recursive https://github.com/tytrusty/virtual-element-method.git

### MATLAB Setup
To let MATLAB recognize the scripts, we need to add this folder and all the project's subfolders to the MATLAB path. Assuming your installation directory `/usr/local/virtual-element-method/`, then you could issue the following command in the MATLAB command prompt:

    addpath(genpath('/usr/local/virtual-element-method/'))
    savepath

### C++ MEX Compilation
**TODO** -- until I add instruction for this you cannot run vem_nurbs.m
Costly functions are implemented in C++ and executed from Matlab from the C++ MEX API. Compiling these is **required** for the NURBs simulation example. The following steps may be used to compile the C++ files:

    cd ${SOURCE_DIRECTORY}/matlab
    mkdir build
    cd build
    cmake ..
    make all
    
