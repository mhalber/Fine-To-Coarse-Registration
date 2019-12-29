This repository is no longer maintained

# Scan Registration

This is the code repository for the code for scan registration project, 
specifically the ["Fine-to-Coarse Global Registration of RGB-D Scans"][1]. This
repository contains four main parts:
- fetregister - main application implementing algorithm described in ["Fine-to-Coarse Global Registration of RGB-D Scans"][1]
- fetbenchmark - simple application for benchmarking results of registration algorithms (for data visit the [project website][2]). 
- gaps - general computer graphics library by Thomas Funkhouser
- basics - visualization and utility helper library
- sun3dsfm - modified version of [code from Xiao et al.][3]. For our purpose we only perform pairwise alignment.

## Building
Provided that you have dependencies installed (see below), you should be able to simply type:
~~~~
build.sh
~~~~
This should build all libraries and programs and put them into *bin/* folder. Moreover each of the **fetregister**, **fetbenchmark** and **gaps** folders contains a makefile, so they can be build separately by running: 
~~~~
make -C <folder_name>
~~~~
You might need to edit those, if your dependencies are not installed in /usr/local/lib. Also note that by default the code builds in window-less mode. To enable visualization and window support in *fetregister* run:
~~~~
make -C fetregister USE_WINDOW=1
~~~~

This code has been tested on both MacOs (Sierra) and Linux (Springdatle Linux 7). We have not tested the code on Windows, however we expect compiling this code using Cygwin should not be hard.

### Dependencies
- [glfw3](http://www.glfw.org/docs/latest/) (windowing)
- [glew](http://glew.sourceforge.net) (opengl extension loading)

## Running registration
For description on how to run registration algorithm, please go to [fetregister](fetregister/Readme.md).

## Running benchmark
For description on how to run benchmark, please go to [fetbenchmark](fetregister/Readme.md).

## Citing
If you use this code or benchmark, please cite our work:
~~~~
@article{HalberF2C17,
  title     = {Fine-To-Coarse Global Registration of RGB-D Scans},
  author    = {Maciej Halber and Thomas Funkhouser},
  booktitle = {CVPR}, 
  year      = {2017} 
}
~~~~


[1]:  http://scanregistration.cs.princeton.edu/assets/HalberF2C_CVPR2017.pdf "Fine-to-Coarse Global Registration of RGB-D Scans"
[2]:  http://scanregistration.cs.princeton.edu "Scan Registation Project websites"
[3]:  https://github.com/PrincetonVision/SUN3Dsfm "Sun3dsfm"
