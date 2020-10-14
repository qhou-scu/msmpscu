module WndFig
use wghDef
use FIG_MODULES


!//The definition for a display handle
integer*4::hScrDC, hPrnDC
integer*4::hPen                !/* "---" pen handle   */
integer*4::hBrush              !/* brush handle       */

!//Number of Child windows
integer*4::nmFigWnd	 = 0

contains  

integer*4 function InitWndFig(Title) 
implicit none
integer*4::hWnd

character*(*)::Title
character*40::Til
character*32::Caption


	   nmFigWnd = nmFigWnd + 1
	   write(caption, 1) nmFigWnd
1      format('Untitled', I3)
       if(len_trim(Title) .gt. 1) Caption = Title 
	   Til = Caption//" "C
       mdicreate%szClass = LOC("FigWnd"C)
       mdicreate%szTitle = LOC(Til)
	   mdicreate%x       = CW_USEDEFAULT
	   mdicreate%y       = CW_USEDEFAULT
       mdicreate%cx      = CW_USEDEFAULT
       mdicreate%cy      = CW_USEDEFAULT
	   mdicreate%hOwner  = hInst
	   mdicreate%style   = IOR(WS_OVERLAPPED, IOR(WS_CAPTION,     &
	                       IOR(WS_BORDER,                         &
      IOR(WS_THICKFRAME, IOR(WS_MAXIMIZEBOX, IOR(WS_MINIMIZEBOX,  &
      IOR(WS_CLIPCHILDREN, IOR(WS_VISIBLE,                        &
	  IOR(WS_VSCROLL,IOR(WS_HSCROLL,WS_SYSMENU))))))))))
!      mdicreate%lParam  = hInfo

!      /*Create Child Window*/
       hWnd = SendMessage(hWndMainClient, WM_MDICREATE,        &
		                            INT4(0), LOC(mdicreate))
       if (hWnd == NULL) then
        LOGISTAT = MessageBox(hWndMain,"Failed in Creating Child Window"C,& 
                            "Error"C, MB_OK)
		return
       end if

  InitWndFig = hWnd
return
end function InitWndFig


integer*4 function FigWndProc(hWnd, message, wParam, lParam)
!MS$ ATTRIBUTES STDCALL, ALIAS : '_FigWndProc@16'::FigWndProc
use FIG_MODULES

implicit none
integer*4::hWnd, message, wParam, lParam, hDC
type (T_PAINTSTRUCT)::ps


!Testing data
real*4::xdat(100), ydat(100)
integer*4::I,N,p, numc
data p/0/, I/0/, numc/0/
save p, I
!

  select case (message) 
  	case (WM_MDIACTIVATE)

		FigWndProc = 0

    case (WM_SETFOCUS)
            logistat = SendMessage(GetParent(hWnd), WM_MDISETMENU,     &
                       hMenuMain, hMenuWindow)
            logistat = DrawMenuBar(Getparent(hWnd) )
		FigWndProc = 0
        return
    case (WM_PAINT)
		hDC = BeginPaint (hWnd, ps)
		call Show_SciFIG(hDC,hWnd)
		logistat = validateRect(hWnd, NULL_RECT)
		logistat = EndPaint(hWnd,  ps)
	 	FigWndProc = 0

    case (WM_COMMAND)
	     n=	INT4(LoWord(wParam))
	     select case(INT4(LoWord(wParam)))
            case (ID_MENU_NewFig)
			     I = I+1
				 if(I .ge. 5) then
		            call New_Fig(hWnd=hWnd, xtype = IOR(axS_Log, IOR(axS_Vis, axS_labed)), &	 
				    ytype = IOR(axS_Log, IOR(axS_Vis, axS_labed)))
				 else 
		            call New_Fig(hWnd=hWnd,  &
					       xtil = "Incident energy", &
                           ytil = "Reflection coefficients.....1234")
				 endif
                 logistat = InvalidateRect(hWnd, NULL_RECT, .FALSE.)
            case (ID_MENU_AddCurve)
			     numc = numc+1
				 if(numc .gt. 16) then
				    numc = numc+1
				    numc = numc+1
				    numc = numc+1
				    numc = numc+1
				    numc = numc+1
				    numc = numc-1
                 end if
			     P=p+1
			     if(I .le. 2) then
		    	    do n=1, 10
				      xdat(n) = P*3.14156*real(n-1)/9.0
				      ydat(n) = sin(xdat(n))*cos(xdat(n))
					enddo
   				call New_CURVE(DATX=xdat,DATY=ydat,DATN=10, &
				     Symbol=iChar(char(ichar('A')+P-1))	 )
                else
		    	   do n=1, 10
				    xdat(n) = .1*(real(P)*real(n)/9.0)
				    ydat(n) = xdat(n)**p
				  enddo
   			      call New_CURVE(DATX=xdat,DATY=ydat,DATN=10)
			    endif 
            case (ID_MENU_PRINT)
                 call Print_SciFig(hWnd)
         end select
	   
  end select

   FigWndProc = DefMDIChildProc(hwnd, message, wParam, lParam)

return
end  function FigWndProc


end module WndFig

