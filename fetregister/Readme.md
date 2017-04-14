# fetregister

Fetregister is main application in this code-base, implementing the algorithm described in ["Fine-to-Coarse Global Registration of RGB-D Scans"][1].

## Compiling
It should be possible to compile this code by simply running make in current folder. This requires gaps to be compiled, and dependencies ([glfw3](http://www.glfw.org/docs/latest/) and [glew](http://glew.sourceforge.net)) to be installed. Note that by default, the code builds in a window-less mode. If you would like to use visualization features, build this program using following command: 
~~~~ 
make USE_WINDOW=1
~~~~

## Running 
To run the code you need 2 files for any RGB-D sequence: the configuration(.conf) file and feature(.fcb) file. Their formats are described at [our project website](http://scanregistration.cs.princeton.edu). In short .conf file will store additional information about sequence like number of images and poses of each individual frame, while feature file stores pre-processed information from RGB-D frames. 

### Extended SUN3D Dataset
You simply need to download the scenes from [the project website](http://scanregistration.cs.princeton.edu). Extract the contents and navigate to one of the sequences folders. Once there, simply run:
~~~~ 
fetregister fet/5_0.05.fcb <configuration_name.conf> -v
~~~~
To give a concrete example, if you navigate to harvard_c11/hv_c11_2 you should run:
~~~~ 
fetregister fet/5_0.05.fcb hv_c11_2_initial.conf -v
~~~~

### Your own sequence 
To run fine-to-coarse registration on your own sequence you must first output your images as a sequence of color and depth images, as well as have a depth camera intrinsics. Then you need to generate a configuration file, describing this sequence - minimally you
will need a file that looks like this:
~~~~ 
dataset <your_dataset>
n_images <number_of_images>
image_directory <color_dir_name>
depth_directory <depth_dir_name>
intrinsics <depth_camera_intrinsics_filename>
matches matches/pairwise_matches.txt

scan <depth_name1> <color_name1> 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1
scan <depth_name2> <color_name2> 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1
....
~~~~

For a more detailed description of this format, please refer to [our project website](http://scanregistration.cs.princeton.edu)

Now you will need to run [conf2fet](../gaps/apps/conf2fet/Readme.md) to generate the feature file, as well as our modified version of [sun3dsfm](../sun3dsfm/Readme.md) pairwise_matches.txt storing the initial transformations. Follow the links to see details for both.

Once both processes are complete you should be able to run fetregister in sequence's root directory as:
~~~~ 
fetregister <your_feature_file.fcb> <your_configuration_name.conf> -v
~~~~


[1]:  http://scanregistration.cs.princeton.edu/assets/HalberF2C_CVPR2017.pdf "Fine-to-Coarse Global Registration of RGB-D Scans"