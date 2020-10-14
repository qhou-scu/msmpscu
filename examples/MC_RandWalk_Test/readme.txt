This folder conatins a number of input files for testing MC_RandWalker programs

1. Walk-1D-MSD-Test1.exe
   this program is used to calculate MSD and spread function of random-walk with the 
   wating time distribution (WTD) is exponetial or pow-law
   
   to run Walk-1D-MSD-Test1.exe example:
   cat Walk-1D-MSD-Test1.stdin  | Walk-1D-MSD-Test1.exe     
   cat Walk-1D-MSD-Test1.stdin1 | Walk-1D-MSD-Test1.exe     
      
         

2. Walk-1D.exe
   this program is an advance of Walk-1D-MSD-Test1,  used to calculate the random walk process, 
   with or without continueous inserting of particles. The particles are assumed without ineracting 
   between them.

   Running  Walk-1D.exe could be commandline, for example:
     Walk-1D.exe -CF Walk_Surf1D.setup1 -WF Walker_Exp.setup1 
   or 
     cat Walk-1D.in |  Walk-1D.exe
   where, Walk_Surf1D.setup1 is a control file, 
          Walker_Exp.setup1  is a desciption file of the walker. 
     