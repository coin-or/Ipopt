C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: subclass.f 531 2004-03-11 01:31:07Z andreasw $
!**********************************************************************
!
! Windows tricks with the windows that could not implemented by QWIN 
!          functions
!
!======================================================================
!
!     Module Subclass   
!

      module subclass

      use dflib
      use dfwin
      implicit none

      private
      public :: SubClassInit
       

      integer(4) FrameProc, ChildProc

      contains

      subroutine SubclassInit()

      integer i, j
      logical l

      FrameProc = SetWindowLong( GetHWndQQ(QWIN$FRAMEWINDOW), 
     1                     GWL_WNDPROC, loc(SubclassFrame) )
      j = GetWindowLong( GetHWndQQ(QWIN$FRAMEWINDOW), GWL_STYLE )
      j = ior( iand( j, not(WS_THICKFRAME) ), WS_BORDER )
      j = iand( j, not(WS_MAXIMIZEBOX) )
      i = SetWindowLong( GetHWndQQ(QWIN$FRAMEWINDOW), GWL_STYLE, j )    

      l = SetWindowText( GetHwndQQ(QWIN$FRAMEWINDOW), 
     1                "Visual Dynamic Optimization Outputs "C )

      ChildProc = SetWindowLong( GetHWndQQ(0), 
     1                     GWL_WNDPROC, loc(SubclassChild) )
      !j = GetWindowLong( GetHWndQQ(0), GWL_STYLE )
      !j = ior( iand( j, not(WS_CAPTION.or.WS_THICKFRAME.or.WS_SYSMENU) ), WS_BORDER )
      !i = SetWindowLong( GetHWndQQ(0), GWL_STYLE, j )      
      l = SetWindowText( GetHwndQQ(0), 
     1                "Message from IPOPT "C )

      l = MoveWindow( GetHWndQQ(0), -1, 600, 1100, 120, .TRUE. )

      !ChildProc = SetWindowLong( GetHWndQQ(CurveUnit), &
      !                    GWL_WNDPROC, loc(SubclassChild) )
      !j = GetWindowLong( GetHWndQQ(CurveUnit), GWL_STYLE )
      !j = ior( iand( j, not(WS_CAPTION.or.WS_THICKFRAME.or.WS_SYSMENU) ), WS_BORDER )
      !i = SetWindowLong( GetHWndQQ(CurveUnit), GWL_STYLE, j )    

      !l = MoveWindow( GetHWndQQ(CurveUnit),100,100,800,400, .TRUE. )

      end subroutine

      subroutine Subclass1_CW()
!**************************************************************
! to delete the scroll bars and to disable resize of that Child
! window
!          YDL 
!          Feb. 21, 2001
!**************************************************************
!
      end subroutine

      integer(4) function SubclassFrame( Hwnd, Msg, wParam, lParam )
      !DEC$ attributes stdcall :: SubclassFrame
      ! use pokerregistry
      automatic
      integer(4) Hwnd, Msg, wParam, lParam

      if( Msg == WM_CLOSE ) then
!     call SaveRegistry()     
      end if
      SubclassFrame = CallWindowProc( FrameProc, 
     1                              Hwnd, Msg, wParam, lParam )
      end function

      integer(4) function SubclassChild( Hwnd, Msg, wParam, lParam )
      !use layout
      !DEC$ attributes stdcall :: SubclassChild
      automatic
      integer(4) Hwnd, Msg, wParam, lParam

      SubclassChild = CallWindowProc( ChildProc, Hwnd, 
     1                                    Msg, wParam, lParam )
      end function

      end module
