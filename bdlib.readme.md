# bdlib
bdlib is to build the static libraries for applications based on this packaged. 
The LIB and MSMLIB directory contain the sources for some mathematical and data operations either on CPU and GPU.  
The MDLIB directory contains sources for MD simulations. The MCLIB directories contains sources for MC simulations. Both MDLIB and MCLIB are dependent of LIB and MSMLIB.
MDLIB and MCLIB are independent, and thus can be build seperately.

    Run "bidlib MDLIB", to build the static libraries for MD applications
    
    Run "bidlib MCLIB", to build the static libraies for MC applications
    
    Run "bidlib all",   to build all the static libraies. 
   


