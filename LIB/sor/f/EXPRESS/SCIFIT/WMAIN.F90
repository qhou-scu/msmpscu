 
!/******************************************************************************\
!*       This is a part of the Microsoft Source Code Samples. 
!*       Copyright (C) 1993 Microsoft Corporation.
!*       All rights reserved. 
!*       This source code is only intended as a supplement to 
!*       Microsoft Development Tools and/or WinHelp documentation.
!*       See these sources for detailed information regarding the 
!*       Microsoft samples programs.
!\******************************************************************************/

!/****************************************************************************
!
!    PROGRAM: Output.f90
!
!    PURPOSE: Output template for Windows applications
!
!    FUNCTIONS:
!
!        WinMain() - calls initialization function, processes message loop
!        InitApplication() - initializes window data and registers window
!        InitInstance() - saves instance handle and creates main window
!        MainWndProc() - processes messages
!        About() - processes messages for "About" dialog box
!
!****************************************************************************/

interface
 integer(4) function WinMain (hInstance, hPrevInstance, lpszCmdLine, nCmdShow)
 !MS$ATTRIBUTES STDCALL, ALIAS : '_WinMain@16' :: WinMain
 integer(4) hInstance
 integer(4) hPrevInstance
 integer(4) lpszCmdLine
 integer(4) nCmdShow
 end function WinMain
end interface
end

!/****************************************************************************
!
!    FUNCTION: WinMain(HANDLE, HANDLE, LPSTR, int)
!
!    PURPOSE: calls initialization function, processes message loop
!
!****************************************************************************/

integer(4) function WinMain( hInstance, hPrevInstance, lpCmdLine,  nCmdShow )
!MS$ ATTRIBUTES STDCALL, ALIAS : '_WinMain@16' :: WinMain
use msfwina

implicit none

interface
 integer*4 function InitApplication (hInstance)
 !MS$ ATTRIBUTES STDCALL, ALIAS : '_InitApplication@4' :: InitApplication
 integer             hInstance
 end  function InitApplication
end interface
interface 
 integer*4 function InitInstance (hWnd,hInstance, nCmdShow)
 !MS$ ATTRIBUTES STDCALL, ALIAS : '_InitInstance@12' :: InitInstance
 integer*4   hInstance, nCmdShow
 integer*4   hWnd
 end function InitInstance
end interface

integer        hWnd,hInstance, hPrevInstance, lpCmdLine, nCmdShow
type (T_MSG)     mesg
logical(4)::bret
integer(4)::ret

lpcmdline = lpcmdline
   
if (hPrevInstance == 0) then
    if (InitApplication(hInstance) == 0) then
        WinMain = 0
        return
    end if
end if

if (InitInstance(hWnd,hInstance, nCmdShow) == 0) then
    WinMain = 0
    return
end if
do while (GetMessage(mesg, NULL, 0, 0))
    bret = TranslateMessage(mesg)
    ret = DispatchMessage(mesg)
end do
WinMain = mesg%wParam
return
end 


!/****************************************************************************
!
!    FUNCTION: InitApplication(HANDLE)
!
!    PURPOSE: Initializes window data and registers window class
!
!****************************************************************************/
integer*4 function InitApplication (hInstance)
!MS$ ATTRIBUTES STDCALL, ALIAS : '_InitApplication@4' :: InitApplication
use msfwina
use WndFig, only:FigWndProc
use CALLBACK_MOD, only:MainWndProc


implicit none
integer             hInstance
type (T_WNDCLASS)     wc

wc%style = IOR(CS_HREDRAW,IOR(CS_VREDRAW,CS_OWNDC))
wc%lpfnWndProc      = LOC(MainWndProc)
wc%cbClsExtra       = 0
wc%cbWndExtra       = 4
wc%hInstance = hInstance
!wc%hIcon = LoadIcon(hInstance, LOC("PANDAICO1"C))
wc%hIcon            = LoadIcon(NULL, IDI_APPLICATION)
!wc%hCursor = LoadCursor(NULL, IDC_ARROW)
wc%hbrBackground = GetStockObject(WHITE_BRUSH) 
wc%hbrBackground = (COLOR_APPWORKSPACE) 
wc%lpszMenuName =  LOC("MENUMAIN"C)
wc%lpszClassName = LOC("SciFit"C)

if( RegisterClass(wc)==0) then
    InitApplication = 0
	return
endif

wc%lpfnWndProc      = LOC(FigWndProc)
wc%hIcon            = LoadIcon(NULL, IDI_APPLICATION)
wc%hbrBackground = GetStockObject(WHITE_BRUSH) 
wc%lpszMenuName     = LOC("MENUMAIN"C)
wc%lpszClassName    = LOC("FigWnd"C)

if (RegisterClass(wc) == 0) then
    InitApplication = 0
    return
end if


return
end                                                                                                                                                                       


!/****************************************************************************
!                                      
!    FUNCTION:  InitInstance(HANDLE, int)
!
!    PURPOSE:  Saves instance handle and creates main window
!
!****************************************************************************/

integer*4 function InitInstance (hWnd,hInstance, nCmdShow)
!MS$ ATTRIBUTES STDCALL, ALIAS : '_InitInstance@12' :: InitInstance
use wghDef

implicit none

integer*4   hInstance, nCmdShow
integer*4   hWnd
logical(4)  bret
integer*4:: ret

hInst = hInstance
!STYLE = IOR(INT(WS_OVERLAPPEDWINDOW),INT(WS_CLIPCHILDREN)) 
!STYLE = IOR(STYLE,INT((WS_VSCROLL))) 
!STYLE = IOR(STYLE,INT((WS_VSCROLL))) 

hMenuMain    = LoadMenu(hInst, LOC("MENUMAIN"C))
hMenuWindow  = GetSubMenu(hMenuMain, 2)


!hWnd  = CreateWindow(                                             &
!            "PANDA"C,                                             &
!            "Sample1 Application"C,                               &
!            STYLE,                                                &
!            CW_USEDEFAULT,                                        &
!            CW_USEDEFAULT,                                        &
!            CW_USEDEFAULT,                                        &
!            CW_USEDEFAULT,                                        &
!            NULL,                                                 &
!            NULL,                                                 &
!            hInstance,                                            &
!            NULL                                                  &
!       )
	      
hwnd = CreateWindowEx(0, "SciFit"C,   "SciFit"C,                    &
      IOR(WS_OVERLAPPED, IOR(WS_CAPTION, IOR(WS_BORDER,           &
      IOR(WS_THICKFRAME, IOR(WS_MAXIMIZEBOX, IOR(WS_MINIMIZEBOX,  &
      IOR(WS_CLIPCHILDREN, IOR(WS_VISIBLE, WS_SYSMENU)))))))),    &
      CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, &
      NULL, hMenuMain, hInstance, NULL)
		                                                         
hWndMain = hWnd
if (hWnd == 0) then
    InitInstance = 0
    return
end if
! Just to verify if SetWindowLong is working
! First call sets the window value to 4
! The return value after the second SetWindowLong call should be 4
ret = SetWindowLong(hWndMain, GWL_USERDATA, 4) 
ret = SetWindowLong(hWndMain, GWL_USERDATA, 0)
ret = SetFocus(hWndMain)    !/* set initial focus */

bret = ShowWindow(hWnd, SW_MAXIMIZE)
bret = UpdateWindow(hWnd)


InitInstance = 1
return
end 



