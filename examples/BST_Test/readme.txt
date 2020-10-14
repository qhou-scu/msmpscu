This folder conatins a number of input files for testing the
application: Cluster_Growth

to run the example:
On linux:  
   run Commandline: Cluster_Growth_GPU.exe   "SampleSetup.dat" 0 1
   
In VisualStudio(windows): 
   open the project Cluster_Growth
   open the property page of the poject
   set  "Working Directory" to be the present directory
   set  "Application Arguments" to be:   
        "SetupSampe_He_Cluster"      0 1 ,   for running He cluster growth, 
        "SetupSampe_H_Cluster"       0 1 ,   for running H  cluster growth, 
        "SetupSampe_HHe_Cluster"     0 1 ,   for running H-He  cluster growth, 
        "SetupSampe_Vac_Cluster.dat" 0 1 ,   for running vacancy cluster growth, 
   choice the menu: "Debug"-"Start Without Debugging"
      
         

 