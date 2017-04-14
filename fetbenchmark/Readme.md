# fetregister

Fetregister is main application in this code-base, implementing the algorithm described in ["Fine-to-Coarse Global Registration of RGB-D Scans"][1].

## Compiling
It should be possible to compile this code by simply running make in current folder.

## Running 
To run the code you need the configuration(.conf) file which stores your computed poses. The format of configuration file is described at [our project website](http://scanregistration.cs.princeton.edu).
Simple statistic( RMSE, Mean and Standard Deviation) can be obtained by simply running 
~~~~ 
fetbenchmark <configuration_name.conf>
~~~~
To get more detailed, per correspondence errors, run:
~~~~ 
fetregister <configuration_name.conf> -v
~~~~

[1]:  http://scanregistration.cs.princeton.edu/assets/HalberF2C_CVPR2017.pdf "Fine-to-Coarse Global Registration of RGB-D Scans"