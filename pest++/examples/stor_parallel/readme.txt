BACKGROUND
-----------
This directory contains files which run a simple storage problem 
using PEST++ in parallel mode using the GENIE parallel run manager.
This is the same example problem whose serial solution was published
in Appendix 6 of "PEST++, a Parameter ESTimation code optimized for large 
environmental models. U.S. Geological Survey Techniques and Methods,
book 7, section C5". 


HOW TO RUN THIS EXAMPLE
-----------------------
To run the this example:
  1) make sure that PEST++, gman and gslave are installed and specified 
     in your system PATH environment variable
  2) from this directory run the windows batch file: example.bat
  3) this will start the GENIE gman run manager, three GENIE gslaves and
     PEST++.  Once launched, this example will pop up a number of
     command windows.  PEST++, gman and each of the gslaves will start in 
     its own window and the gslaves will open additional windows to perform 
     the model runs.  Once the run is complete, all of the windows except 
     the main PEST++ window should close automatically.


HOW THIS EXAMPLE IS STRUCTURED
------------------------------
The PEST++ input files and the files necessary to run the simple storage 
model are contained in the directory ./template.  When run, the script
run_example.bat performs the following tasks:
  1) the template directory and all of its content are copied to the
     directories master, slave1, slave2 and slave3.
  2) gman is started in this directory (ie the root directory of the
     example problem)
  3) the gslaves are started in the slave directories created in step 1
  4) PEST++ is started in the master directory created in step 1
     
     
To tell PEST++ to use the GENIE parallel run manager, the following 
line has been added to the storage5.pst PEST++ control file:
    ++gman_socket(127.0.0.1:24772)
    
This line instructs PEST++ to run in parallel mode and tells it to look for
the GENIE gman run manager at IP address 127.0.0.1 and port 24772.
The address 127.0.0.1 refers to the local computer and is often referred to
as the loop back address.  Using this IP address allows this example to
be run on any windows computer without modification.

The run_example.bat file can be easily changed by modifying the follow
parameters at the top of the script:
   set nslaves=3            -> sets the number of slaves
   set gman_ip=127.0.0.1    -> sets the IP address for gman
   set gman_port=24772      -> sets the port for gman
   set pestpp_file=storage5 -> sets the name of the PEST++ control file
   
   
   
