program miniExpress
use Express_Module
use Express_LOGICAL_Module
use Express_Ext_Module
use INTER_IO_module
use STRING_OP_module
implicit none

integer*4, dimension(4)::Switch 
integer*4::hPipeI=0, hPipeO=0, ERR, I
character*80::para(4)
character*32,dimension(3)::V
real*8,dimension(3,1)::VV
logical, dimension(1)::log

!**** to initialize the IO files and handles
	Switch(1) = ISTR4('iF')
	Switch(2) = ISTR4('oF')
	Switch(3) = ISTR4('iP')
	Switch(4) = ISTR4('oP')
    para(1)='express.tab'
	para(2)='express.dat'
	para(3)=''
	para(4)=''

	if( GET_COMD_LINE(4,switch,para) .gt. 0) then
	   stop 'Syntax error in command line'
      end if

	if(len_trim(para(3)) .gt. 0) then
	   read(para(3),*)hPipeI
      end if
	if(len_trim(para(4)) .gt. 0) then
	   read(para(4),*)hPipeO
      end if

    call Get_SubProc_form_File(Para(1)) 
    call Begin_Calculation()
    call Putout_results_to_File(Para(2)) 

    V=(/'A','B','C'/)
	VV(1:3,1)=(/10.0,20.0,30.0/)

    ERR = GetLogExpress('A+B<B+C',log,V,VV)
    ERR = GetLogExpress('A*B>B*C',log,V,VV)
    ERR = GetLogExpress('B/A>C/B',log,V,VV)

  
stop 'miniExpress'
end program miniExpress
