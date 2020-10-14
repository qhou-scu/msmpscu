module CALLBACK_MOD

character(len=256), private::cBuff, fname, curPath
integer, parameter, private::slen=32
character(len=slen), private::str

contains

!/****************************************************************************
!
!    FUNCTION: MainWndProc(HWND, unsigned, WORD, LONG)
!
!    PURPOSE:  Processes messages
!
!    MESSAGES:
!
!        WM_COMMAND    - application menu (About dialog box)
!        WM_CREATE     - create window and objects
!        WM_PAINT      - update window, draw objects
!        WM_DESTROY    - destroy window
!
!    COMMENTS:
!
!        Handles to the objects you will use are obtained when the WM_CREATE
!        message is received, and deleted when the WM_DESTROY message is
!        received.  The actual drawing is done whenever a WM_PAINT message is
!        received.
!
!****************************************************************************/

integer*4 function MainWndProc(hWnd, message, wParam, lParam)
!MS$ ATTRIBUTES STDCALL, ALIAS : '_MainWndProc@16' :: MainWndProc
use WndFig

implicit none
integer*4::hWnd, message, wParam, lParam, hThread, ThreadID
type (T_CLIENTCREATESTRUCT)::clientcreate
type(T_SECURITY_ATTRIBUTES)::lpThreadAttributes 
integer*4::hActiveChild, hDlg
integer*4::TOOL1H =40
type(T_POINT)::POINT


select case (message) 

    case (WM_CREATE)

		call InitSciFig(hWnd)
	      !logistat = GetClientRect(hWnd, RC)
          !hWndJobList = CreateWindowEx(WS_EX_TOPMOST, "LISTBOX"C, ""C,               &
          !        IOR(WS_VISIBLE , IOR(WS_THICKFRAME, IOR(WS_CHILD,  IOR(WS_CLIPSIBLINGS,   &
	      !        IOR(LBS_STANDARD,LBS_OWNERDRAWFIXED))))),                                 &
          !        0, 32, 150, RC%BOTTOM-RC%TOP,                                              &
          !        hWnd, NULL, hInst, NULL)

		 STAT = SetWindowLong(hWnd, 0, NULL)
         clientcreate%hWindowMenu  = hMenuMain
         clientcreate%idFirstChild = 1
         hWndMainClient = CreateWindow("MDICLIENT"C, ""C,         &
                            IOR(WS_CHILD , IOR(WS_CLIPCHILDREN,   &
                                WS_VISIBLE)), 0, 0, 0, 0,         &
                        hWnd, NULL, hInst, LOC(clientcreate))

         !hWndTool1   = CreateDialogParam(hInst, loc("DLG_TOOL1"c), hWnd,  &
		 !                          loc(proDlgTool1), 0)

         !hWndJoblist = CreateDialogParam(hInst, loc("DLG_JOBLIST"c), hWnd,  &
		 !                          loc(proDlgJoblist), 0)


        MainWndProc = 0

          
    case (WM_COMMAND)
        select case (INT4(LoWord(wParam)))

            case (ID_MENU_NewWindow)
				  hWnd    = InitWndFig("")
				  MainWndProc = 0
                  return
				   
            case (ID_MENU_EXIT)
                STAT = SendMessage(hWnd,WM_CLOSE, wParam, lParam)
                 MainWndProc = 0
                 return

           case (ID_MENU_Tile)
                LOGISTAT = SendMessage(hWndMainClient, WM_MDITILE, 0, 0)
                MainWndProc = 0
                return

            case (ID_MENU_Cascade)
                LOGISTAT = SendMessage(hWndMainClient,WM_MDICASCADE, 0, 0)
                MainWndProc = 0

            case (ID_MENU_ArrangeIcon)
                LOGISTAT= SendMessage(hWndMainClient, WM_MDIICONARRANGE,&
				                       0, 0)
                 MainWndProc = 0
            case (ID_MENU_ABOUT)
                  !stat = Dlg_About()
                  MainWndProc = 0


!*******************

            case (ID_MENU_NewFig, ID_MENU_OPENFig, ID_MENU_SELECTFIG, &
			      ID_MENU_AddCurve, ID_MENU_SelCurve, ID_MENU_PRINT)
                hActiveChild = SendMessage(hWndMainClient,            &
                                            WM_MDIGETACTIVE, 0, 0)
                if(hActiveChild .NE. 0) then
                   logistat = SendMessage(hActiveChild, WM_COMMAND,   &
                                    wParam, lParam)
                end if
                MainWndProc = 0

            case DEFAULT
                MainWndProc = DefFrameProc(hwnd,  hWndMainClient,  &
                                message, wParam, lParam)
		end select


!******
	  case (WM_DESTROY) 
	    call End_SciFig()
        call PostQuitMessage(0)

      case DEFAULT
		MainWndProc = DefFrameProc(hwnd, hWndMainClient, message, wParam, lParam)

        return
end select
MainWndProc = 0
return
end function MainWndProc


end module CALLBACK_MOD
