 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to simulate cluster growth.
 !                  The growth of a cluster is simulated by place a new atom in the neigbor region of
 !                  the cluster. The Wigner-Seitz cells occupied by the cluster are firstly calcualted,
 !                  Then, the new atoms will be added to the WS sites that are a few WS site away from the
 !                  cluster. The module RefVoronoiVacancy_GPU to be used to identify the occupied WS cells.
 !
 !                  DEPENDENCE____________________________________________________________________________
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  by HOU Qing, 2018-01
 !

 !****************************************************************************************
  module Cascade_GPU
    use MD_TypeDef_Projectile

    contains
  end module Cascade_GPU

 !****************************************************************************************
 !****************************************************************************************
  Program Cascade_GPU_Main
  use MD_SimBoxArray_AppShell_16_GPU
  use Cascade_GPU
  implicit none
  integer::numprocs=1, procid=0

       call APPSHELL_AddRecord( PRERECORD =Initialize_Projectile)
       call APPSHELL_Main(numprocs, procid, INICONFIGPROC=Deposit_Projectile)

       stop
  end  program Cascade_GPU_Main
