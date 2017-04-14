# conf2fet

conf2fet utility processes provided configuration file and produces feature(.fcb) files. 

## Compiling
This utility should be compiled while compiling gaps. To recompile, you should be able to simply run **make** command in this folder.

## Running 
Once configuration file is generated, you can generate fet file running following command:
~~~~ 
conf2fet <input_configuration_filename> <output_feature_filename>
~~~~
While conf2file has multiple parameters, two most important ones are -load_every_kth_image and -min_feature_spacing. As described in [our paper][1], for our experiments we used n=5 and m=0.05 for these repesctive arguments. Thus the command in the whole is:
~~~~ 
conf2fet <input_configuration_filename> <output_feature_filename> -load_every_kth_image n -min_feature_spacing m -v
~~~~
Here we also added -v for more verbose printing.

[1]:  https://arxiv.org/pdf/1607.08539v3.pdf "Fine-to-Coarse Global Registration of RGB-D Scans"