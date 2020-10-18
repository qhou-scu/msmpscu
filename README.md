# msmpscu
 A GPU-based package to conduct multiple-scal simulations the evolution of material structure.
 
  Platform requirements:  
 
    operation platform: Linux workstation with at least one GPU card in Nvidia Kx or later serials.  
    compiler:  PGI Fortran community Edition version.
    Cuda: 8.0 or higher.
 
 Installation requirements:  
 This package dose not provide a single excutable, but consists of many libraries used in MD, or MC applications. The libraries and the applications are built seperately. 
 
 To correctlly build the msmpscu libraries, the envirenment varibale should be set in file .bashrc. Run command "resetenv" as follows:

    . .resetenv
    
 to setup the enviremental variables, or manually change the .bashrc file as follows.
 
 By default, the installed cuda is assumed to be 8.0 with the caplbility 3.0. 
 If the other version, for example, cuda10.1 is installed, add the variable to .bashrc file

    export CUDAV=10.1
    export CUDAC=7.0

If the path of msmpscu sources is SOMEWHERE in youe home path, add the following statement in your .bashrc file:

    MSMPSCUSOR=$HOME/SOMEWHERE/msmpscu; export MSMPSCUSOR
    PATH=$MSMPSCUSOR:$PATH; export PATH
 
It is recommanded to output the intermediate files (.o, .mod) and libraries (lib*.a) to a location beyond the localtion of the mdpscu sources.
For example, if the your workspace is SOMEWHERE in your home path,  add the following statement in your .bashrc file:

    WORKSPACE=$HOME/SOMEWHERE; export WORKSPACE
    PATH=$WORKSPACE/applications:$PATH; export PATH

Then, the intermediate files and libraries will be output to 	$HOME/SOMEWHERE/LIB, and the excutables to be build will be output to $HOME/SOMEWHERE/applications.
Run the enirenment variables set, run command "bidlib" to build the libraries (see bidlib.readme). Run command "gapp" to build an specific application(see gapp.readme).  

  
  
