# HFrame

///////////////////////////////////////////////////////////////////////////////
//
//
//
//                 HFrame:A Tool to Detect Stave Core Thermal Flaws
//                 Author: Jie Yu, Editor: Will Heidorn
//
//
//
///////////////////////////////////////////////////////////////////////////////

This code reads output from the seqToRoot program, which can be found here
https://jie-yu.web.cern.ch/jie-yu/work/2016_StaveQC/seqToRoot_package.tgz 

The seqToRoot tool converts *.seq files to *.root files. The HFrame Tool 
assumes that you have taken 200 frames of data where the stave core has fluid 
running through it at a temperature of either +50 or -40 C.

To install and test the program...

1. From HFrame directory in command line

$ make

This will create the library that will be used. Next you need to link that
library. Type

$ export DYLD_LIBRARY_PATH=$PWD/bin:$DYLD_LIBRARY_PATH

If this does not work instead replace $PWD with you exact directory address.
There are three example programs that have been written. Each of them
analyzes the data given in the folder roo/LongStaveNo5_L_80lens_m55/ 
If you wish to change the analyzed data create a new folder under roo and also
update the input.txt file so that they are the same. There are 3 example
macros. Each can be run by going to the directory and running the run.sh file.

$ cd examples/A1_rooToFrames
$ ./run.sh

Each of the three macros creates new files and requires the previous macro
to have been run.

A1_rooToFrames takes 200 frames and averages them. It then cuts out the 
unnecessary background data.

B1_frameToPipe fits a double peak function perpendicular to the pipes and
finds the maximum/minimum at each point along the stave. These values are then
used to create an average temperature along the stave cooling pipes.

C1_pipeToPeaks uses the TSpectrum tool to find peaks along the stave. This
often finds extra stuff and is not very reliable, but it may be updated in the
future.
 

