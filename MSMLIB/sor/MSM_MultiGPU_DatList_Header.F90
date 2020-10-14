  !---------------------------------------------------------
     
          type::LIISTNAME
                character(len=32)      ::Tag = "" 
                DATATYPE,       pointer::TheData
                type(LIISTNAME),pointer::Next =>null() 
          end type LIISTNAME

          private:: Add_0,     &
                    Add_1 
          public::  AddDataToList
          interface AddDataToList
                    module procedure Add_0
                    module procedure Add_1    
          end interface AddDataToList
          
          private:: Get_0,     &
                    Get_1      
          public::  GetDataFromList
          interface GetDataFromList
                    module procedure Get_0
                    module procedure Get_1    
          end interface GetDataFromList 

          private:: Numberof_0
          public::  NumberofDataFromList
          interface NumberofDataFromList
                    module procedure Numberof_0
          end interface NumberofDataFromList 

          private:: Clear_0
          public::  ClearDataFromList
          interface ClearDataFromList
                    module procedure Clear_0
          end interface ClearDataFromList 

