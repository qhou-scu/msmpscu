  module EMBEDMENT_COMMON_2010
  !***  DESCRIPTION:
  !     This program is used to calculate the integeraction radius for small He clusters
  !     embedded in materials.
  !
  !
  !    Adopted by HOU Qing, Dec, 2010
  !

  !*** the application specific module ********************************
  use EMBEDMENT_TypeDef_CtrlParam_2010
  implicit none

      type(EMBEDCtrlParam), private::m_EMBEDCtrlParam
      character(len=11),parameter, private::mp_FTAGI="&AUXF_EMBED"

  contains

  !****************************************************************************************
  subroutine Initialize12_EMBEDMENT(SimBox,CtrlParam, dummy)
  !***  DESCRIPTION:  to load controal parameters and initialize the parameters
  !                   when the APPShell module MD_EM_TB_ForceTable_Shell_12_GPU to be used
  !     INPUT:        Simbox, the description of the simulation box
  !                   CtrlParam, the generic control parameters
   use MD_CONSTANTS
   use MD_SimboxArray
   use MD_TYPEDEF_SimMDCtrl
   implicit none
   !--- dummy variables and subroutines
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam
       optional       ::dummy
       external         dummy
       !--- local
       call Initialize14_EMBEDMENT(SimBox,CtrlParam)

    return
  end subroutine Initialize12_EMBEDMENT
  !****************************************************************************************

  subroutine Initialize14_EMBEDMENT(SimBox,CtrlParam)
  !***  DESCRIPTION:  to load controal parameters and initialize the parameters
  !                   when the APPShell module MD_SimBoxArray_ToolShell_14_GPU to be used
  !     INPUT:        Simbox, the description of the simulation box
  !                   CtrlParam, the generic control parameters
   use MD_CONSTANTS
   use MD_SimboxArray
   use MD_TYPEDEF_SimMDCtrl
   implicit none
   !--- dummy variables and subroutines
       type(SimMDBox)  ::SimBox
       type(SimMDCtrl) ::CtrlParam
       !--- local
       integer::I


      !*** to load control parameters specific for embeddment calculations
           do I=1, size(CtrlParam%f_tag)
               if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
                  exit
               end if
           enddo
           if(I.gt. size(CtrlParam%f_tag)) goto 100
           if(len_trim(CtrlParam%f_others(I)).le.0) goto 100

           call Initialize_EMBEDCtrlParam(CtrlParam%f_others(I),m_EMBEDCtrlParam)
           call Print_Parameter_EMBEDCtrlParam(6, m_EMBEDCtrlParam)


           return
    100    write(*,*)"MDPSCU error: control file for embedment of atoms  is not given"
           write(*,*)"              check the setup file for keyword ", mp_FTAGI
           stop

    return
  end subroutine Initialize14_EMBEDMENT
  !****************************************************************************************

  !****************************************************************************************
  subroutine IniConfig_EMBEDMENT(SimBox, CtrlParam, RESTART)
  !***  PORPOSE: to initialize the configure of the box by emdedding an atom
  !              in a box
  !     INPUT:  SimBox0, the original substarte
  !             SimBox,  the box array to be created
  !             CtrlParam, the control parameter
  !             RESTART, indictor to indicating if restarting a new session
   use MD_CONSTANTS
   use MD_Globle_Variables
   use MD_SimboxArray
   use MD_TYPEDEF_SimMDCtrl
   implicit none
   !--- dummy variables and subroutines
       type(SimMDBox), dimension(:)::SimBox
       type(SimMDCtrl)             ::CtrlParam
       integer::RESTART

  !--- Local variables
  integer::I,K,K1
  real(KINDDF)::RR


           RR = SimBox(1)%RR
           if(RESTART.eq. 0) then
           !--- now we embeded an atom to the substrate
                do K=1,m_EMBEDCtrlParam%nType
                   do K1=1, m_EMBEDCtrlParam%PNPRT(K)
                      select case(m_EMBEDCtrlParam%EmStyle(K))
                             case(EM_STYLE_RANDOM)
                                 call AddOneAtom_SimBoxArray(SimBox, m_EMBEDCtrlParam%aType(K),&
                                                                     DISMI0=RR*m_EMBEDCtrlParam%EMMinSep(K))
                             case(EM_STYLE_CLUSTER)
                                 call AppendOneAtom_SimBoxArray(SimBox, m_EMBEDCtrlParam%aType(K),       &
                                                                        CTYPE0 = m_EMBEDCtrlParam%CTYPE(K),     &
                                                                        DISMI0=RR*m_EMBEDCtrlParam%EMMinSep(K), &
                                                                        INR0=RR*m_EMBEDCtrlParam%EMINR(K),      &
                                                                        OUTR0=RR*m_EMBEDCtrlParam%EMOUTR(K))
                              case(EM_STYLE_REPLACE)
                                   call ReplaceOneAtom_SimBoxArray(SimBox, m_EMBEDCtrlParam%aType(K), &
                                                                                    OLDTY=m_EMBEDCtrlParam%cType(K) )

                              case(EM_STYLE_VACANCY)
                                   call DeleteOneAtom_SimBoxArray(SimBox, m_EMBEDCtrlParam%aType(K) )
                        end select
                     end do
                 end do
           else
                 do I=1, RESTART
                 do K=1,m_EMBEDCtrlParam%nType
                    do K1=1, m_EMBEDCtrlParam%PNPRT(K)
                        select case(m_EMBEDCtrlParam%EmStyle(K))
                               case(EM_STYLE_RANDOM,EM_STYLE_CLUSTER)
                                   call AddOneAtom_SimBoxArray(SimBox, m_EMBEDCtrlParam%aType(K))
                        end select
                    end do
                 end do
                 end do
           end if
       return
  end subroutine IniConfig_EMBEDMENT
  !****************************************************************************************


  end module EMBEDMENT_COMMON_2010
