! Copyright (C) 2009, Jonathan Hogg <jdh41.at.cantab.net>
! All Rights Reserved.
! This code is published under the Common Public License.
!
! $Id$
!
! Authors: Jonathan Hogg                           2009-07-29 

module ma77_iface_data
   use hsl_ma77_double
   implicit none
   ! NOTE: Global data. This is bad as it is not thread safe. We are limited in
   ! passing pointers around between C and Fortran - we cannot use bind(C) as
   ! these types contain allocatable components.
   type(ma77_keep), save :: keep
   type(ma77_info), save :: info77

contains 

subroutine set_control(control, icntl, rcntl)
   integer, dimension(7) :: icntl
   ! 1   print_level
   ! 2   bits
   ! 3   buffer_lpage
   ! 4   buffer_npage
   ! 5   file_size       ! should be long
   ! 6   maxstore        ! should be long
   ! 7   nemin
   double precision, dimension(3) :: rcntl
   ! 1   small
   ! 2   static
   ! 3   u
   type(ma77_control), intent(inout) :: control

   control%print_level = icntl(1)

   control%bits = icntl(2)
   control%buffer_lpage(:) = icntl(3)
   control%buffer_npage(:) = icntl(4)
   control%file_size = icntl(5)
   control%maxstore = icntl(6)

   control%nemin = icntl(7)

   control%small = rcntl(1)
   control%static = rcntl(2)
   control%u = rcntl(3)

end subroutine set_control
end module ma77_iface_data

!subroutine ma77_iface_set_control(print_level, bits, buffer_lpage, &
!      buffer_npage, file_size, maxstore, nemin, small, static, u)
!   use ma77_iface_data
!   implicit none
!
!   integer, intent(in) :: print_level
!   integer, intent(in) :: bits
!   integer, intent(in) :: buffer_lpage
!   integer, intent(in) :: buffer_npage
!   integer, intent(in) :: file_size       ! should be long
!   integer, intent(in) :: maxstore        ! should be long
!   integer, intent(in) :: nemin
!   double precision, intent(in) :: small
!   double precision, intent(in) :: static
!   double precision, intent(in) :: u
!
!   control%print_level = print_level
!
!   control%bits = bits
!   control%buffer_lpage(:) = buffer_lpage
!   control%buffer_npage(:) = buffer_npage
!   control%file_size = file_size
!   control%maxstore = maxstore
!
!   control%nemin = nemin
!
!   control%small = small
!   control%static = static
!   control%u = u
!
!   print *, "setup ", control%small, control%static, control%u
!end subroutine ma77_iface_set_control

subroutine ma77_iface_initStructure(ndim,nnz,ptr,row,order,icntl,rcntl,info)
   use ma77_iface_data
   implicit none

   integer, intent(in) :: ndim
   integer, intent(in) :: nnz
   integer, intent(in) :: ptr(ndim+1)
   integer, intent(in) :: row(nnz)
   integer, intent(inout) :: order(ndim)
   integer, dimension(7), intent(in) :: icntl
   double precision, dimension(3), intent(in) :: rcntl
   integer, intent(out) :: info

   !type(ma77_info) :: info77
   type(ma77_control) :: control
   character(len=100) :: filename(4) = &
      (/ "ma77_int  ", "ma77_real ", "ma77_work ", "ma77_delay" /)
   integer :: i, st

   ! Initialise control
   call set_control(control, icntl, rcntl)

   ! Initialise MA77
   call ma77_open(ndim, filename, keep, control, info77)

   ! Copy data into files
   do i = 1, ndim
      call ma77_input_vars(i, ptr(i+1)-ptr(i), row(ptr(i):ptr(i+1)-1), &
         keep, control, info77)
      if(info77%flag.lt.0) then
         info = info77%flag
         return
      endif
   end do

   ! Call analyse
   call ma77_analyse(order,keep,control,info77)

   info = info77%flag
end subroutine ma77_iface_initStructure

subroutine ma77_iface_factor(ndim, ptr, val, icntl, rcntl, numneg, info)
   use ma77_iface_data
   implicit none
   
   integer, intent(in) :: ndim
   integer, intent(in) :: ptr(ndim+1)
   double precision, intent(in) :: val(*)
   integer, dimension(7), intent(in) :: icntl
   double precision, dimension(3), intent(in) :: rcntl
   integer, intent(out) :: numneg
   integer, intent(out) :: info

   !type(ma77_info) :: info77
   type(ma77_control) :: control
   integer :: i

   ! Initialise control
   call set_control(control, icntl, rcntl)

   ! Initialise return values
   numneg = -1

   ! Copy values into files
   do i = 1, ndim
      call ma77_input_reals(i, ptr(i+1)-ptr(i), val(ptr(i):ptr(i+1)-1), &
         keep, control, info77)
      if(info77%flag.lt.0) then
         print *, "wtf!", info77%flag
         info = info77%flag
         return
      endif
   end do

   ! Call factorise
   call ma77_factor(.false., keep, control, info77)
   info = info77%flag
   numneg = info77%num_neg
end subroutine ma77_iface_factor

subroutine ma77_iface_solve(nrhs, lx, rhs, icntl, rcntl, info)
   use ma77_iface_data
   implicit none
   
   integer, intent(in) :: nrhs
   integer, intent(in) :: lx
   double precision, intent(inout) :: rhs(lx,*)
   integer, dimension(7), intent(in) :: icntl
   double precision, dimension(3), intent(in) :: rcntl
   integer, intent(out) :: info

   !type(ma77_info) :: info77
   type(ma77_control) :: control

   ! Initialise control
   call set_control(control, icntl, rcntl)
   
   call ma77_solve(nrhs, lx, rhs, keep, control, info77)
   info = info77%flag
end subroutine ma77_iface_solve

subroutine ma77_iface_finalise(icntl,rcntl)
   use ma77_iface_data
   implicit none

   integer, dimension(7), intent(in) :: icntl
   double precision, dimension(3), intent(in) :: rcntl

   !type(ma77_info) :: info77
   type(ma77_control) :: control

   ! Initialise control
   call set_control(control, icntl, rcntl)

   call ma77_finalise(keep, control, info77)
end subroutine ma77_iface_finalise
