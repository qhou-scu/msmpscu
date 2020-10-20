gapp is to build excutables for specific applications of of this packaged.
Check the following points before running this command.
1. bdlib have been run for building the basic libraiies needed for MD or MC applications.   
2. For the convenience of access the application sources, make linkes to the directories where the main programs of applications in


    ln -s $MSMPSCUSOR/MDLIB/sor/Applications -T $WORKSPACE/mdapps
    
    ln -s $MSMPSCUSOR/MDLIB/sor/Analysis     -T $WORKSPACE/mdanaly
    
    ln -s $MSMPSCUSOR/MCLIB/sor/Applications -T $WORKSPACE/mcapps. 
   

for the variables $MSMPSCUSOR and $WORKSPACE, see README.md 


3. enter to the directories and run gapp. For example:


    cd $WORKSPACE/mdapp
    gapp SUBSTRATE_GPU
    
to generate the excutable SUBSTRATE.exe, in the dierectory  $WORKSPACE/applications.
    



