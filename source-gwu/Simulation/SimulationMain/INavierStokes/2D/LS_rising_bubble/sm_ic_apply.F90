!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Kinematics/sm_pk_apply_entireBody
!!
!!
!! NAME
!!
!! 
!!
!!
!! SYNOPSIS
!!    apply the prescribed kinematics of RestSurface 1 to the entire body
!!  
!!
!!
!! DESCRIPTION
!!  
!!
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine sm_ic_apply(ibd,time)

  use SolidMechanics_data, only :  sm_BodyInfo, sm_structure
  use Driver_interface, only: Driver_abortFlash
  implicit none

  ! IO
  integer, intent(in) :: ibd
  real, intent(in)    :: time

  ! Internal variables
  type(sm_structure),  pointer :: body

  integer :: A, i, idx
  real, allocatable, dimension(:,:) :: v, vd, vdd
  character(len=6) :: libd

  ! Get the body
  body => sm_BodyInfo(ibd)

  allocate( v(body%max_dofs_per_node,1))
  allocate( vd(body%max_dofs_per_node,1))
  allocate( vdd(body%max_dofs_per_node,1))

  write(libd,'(I6.6)') ibd
  open(123, file='sm_body_ic.'//libd)
    read(123,*) 
    read(123,*) 
    read(123,*) v
    read(123,*) 
    read(123,*) vd
    read(123,*) 
    read(123,*) vdd
  close(123)

  ! Apply kinematics to dofs
  do A = 1,1
     do i = 1,body%max_dofs_per_node
        idx           = body%ID(i,A)
        body%qn(idx)  = v(i,A)
        body%qdn(idx) = vd(i,A)
        body%qddn(idx)= vdd(i,A)
     end do
  end do

  ! deallocate
  deallocate(v, vd, vdd)

end subroutine sm_ic_apply
