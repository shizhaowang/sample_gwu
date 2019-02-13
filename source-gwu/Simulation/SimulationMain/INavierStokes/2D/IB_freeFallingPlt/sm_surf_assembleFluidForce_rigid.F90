!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/sm_surf_assembleFluidForce_rigid
!!
!! NAME
!! 
!!
!! SYNOPSIS
!!
!!  
!! DESCRIPTION 
!! 
!!
!! ARGUMENTS 
!!
!!***

!#include "Flash.h"
!#include "SolidMechanics.h"
!#include "constants.h"
!#include "Flash_mpi.h"
!#define WRITEFORCE 1


subroutine sm_surf_assembleFluidForce_rigid(ibd)

  use SolidMechanics_data, only: sm_meshMe,sm_NumBodies,sm_bodyInfo, sm_structure
  use Driver_interface, only: Driver_abortFlash
  use sm_surf_interface, only : sm_surf_assembleFluidForce_toPoint
  use Driver_data, only : dr_simtime,dr_dt

  implicit none
#include "SolidMechanics.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: ibd
  
  ! Internal variables
  type(sm_structure),  pointer :: body
  real, dimension(NDIM) :: point, force_pres, force_visc, moment_pres, moment_visc
  real, allocatable, dimension(:) :: qn_master,qdn_master,qddn_master
  real :: xm(MDIM)

  integer :: maxdofs
  integer :: i, idx, i_dim, imaster

  character(len=6) :: str_ibd

  logical, save :: first_call = .true.

  real, allocatable, dimension(:) :: sendbuf
  integer :: ierr

  ! Get the body
  body => sm_BodyInfo(ibd)

  ! Extract info from Master:
  ! Reference x,y,z:
  imaster = body % borigin_node
  xm = 0.
  xm(1) = body % x(imaster)
  xm(2) = body % y(imaster)
#if NDIM == MDIM
  xm(3) = body % z(imaster)
#endif          

  force_pres  = 0.
  force_visc  = 0.
  moment_pres = 0.
  moment_visc = 0.

  if (first_call) then
    first_call = .false.
    ! Either Hs_pres, Hs_visc are zero or have been read from checkpoint file. 
    ! Defer.

  else

    ! Hs_pres and Hs_visc to zero:
    body%Hs_pres = 0.
    body%Hs_visc = 0.

    ! Actual state in qn, qdn, qddn for master Node:
    maxdofs   = body%max_dofs_per_node 
    allocate(qn_master(maxdofs),qdn_master(maxdofs),qddn_master(maxdofs))
    qn_master(1:maxdofs)  = body%  qn(body%ID(1:maxdofs,imaster)) 
    qdn_master(1:maxdofs) = body% qdn(body%ID(1:maxdofs,imaster))
    qddn_master(1:maxdofs)= body%qddn(body%ID(1:maxdofs,imaster))
 
    ! Body frame origin actual location:
    point(1:NDIM) = xm(1:NDIM) + qn_master(1:NDIM)

    ! Now Call the integration Routine:
    call sm_surf_assembleFluidForce_toPoint(ibd, point, force_pres, force_visc, moment_pres, moment_visc)


    ! debug , Shizhao
    !force_pres = 0.0
    !force_pres(1) = -0.5*0.1934
    !force_pres(2) = 4.0969*0.1934*0.75
    !force_visc = 0.0
    !moment_pres = 0.0162*0.5
    !moment_visc = 0.0 

    ! Add to whatever is in Hs:
    ! Forces:
    i_dim = 0
    do i = body%ix,body%ex
     i_dim = i_dim+1
     idx = body%ID(i,imaster)
     body%Hs_pres(idx) = body%Hs_pres(idx) + force_pres(i_dim) 
     body%Hs_visc(idx) = body%Hs_visc(idx) + force_visc(i_dim)
    end do

    ! Moments: 
    ! Remember, in 2D moment(CONSTANT_ONE) is the moment for theta variable 
    ! (rotation around axis normal to x-y plane): 
    i_dim = 0
    do i = body%iw,body%ew
     i_dim = i_dim+1
     idx = body%ID(i,imaster)
     body%Hs_pres(idx) = body%Hs_pres(idx) + moment_pres(i_dim) 
     body%Hs_visc(idx) = body%Hs_visc(idx) + moment_visc(i_dim)
    end do

    ! Here if this body is part of a metabody we will perform an all reduce sum for all 
    ! Pressure and viscous forces. 
    ! All processors managing bodies part of a metabody will run for their particular ibd
    ! through the routines sm_IntegX_advance, sm_assemble_ExtForce, sm_assembleFluidForces_rigid
    ! in this case. At this point they should have their contributions to body%Hs_pres, body%Hs_visc
    ! computed. An all reduce on the group is performed to gather the total corresponding
    ! fluid forces for the set:
    if (body%Metabody .gt. CONSTANT_ZERO) then
      allocate(sendbuf(Body%ndofs))

      ! Pressure forces:
      sendbuf(:) = body%Hs_pres(:)
      call MPI_allReduce(sendbuf, body%Hs_pres, Body%ndofs, FLASH_REAL, FLASH_SUM, body%mbcomm, ierr)

      ! Viscous forces:
      sendbuf(:) = body%Hs_visc(:)
      call MPI_allReduce(sendbuf, body%Hs_visc, Body%ndofs, FLASH_REAL, FLASH_SUM, body%mbcomm, ierr)

      deallocate(sendbuf)
    endif

  endif

  !if (sm_meshMe .eq. MASTER_PE) then
  !write(*,*) ' '
  !write(*,*) 'ibd=',ibd,'imet=',body%Metabody,body%Hs_pres(body%ID(body%ix:body%ex,imaster))
  !write(*,*) 'ibd=',ibd,'imet=',body%Metabody,body%Hs_visc(body%ID(body%ix:body%ex,imaster))
  !write(*,*) 'ibd=',ibd,'imet=',body%Metabody,body%Hs_pres(body%ID(body%iw:body%ew,imaster))
  !write(*,*) 'ibd=',ibd,'imet=',body%Metabody,body%Hs_visc(body%ID(body%iw:body%ew,imaster))
  !endif
  
  ! Add to whatever is in Hs:
  ! Forces:
  do i = body%ix,body%ex
     idx = body%ID(i,imaster)
     if (idx .le. Body%neq) then
        body%Hs(idx)  = body%Hs(idx) + body%Hs_pres(idx) + body%Hs_visc(idx)
     endif
  end do

  ! Moments: 
  do i = body%iw,body%ew
     idx = body%ID(i,imaster)
     if (idx .le. Body%neq) then
        body%Hs(idx)  = body%Hs(idx) + body%Hs_pres(idx) + body%Hs_visc(idx)
     endif
  end do

#ifdef DEBUG_SOLID
  write(*,*) 'Body, Proc =',ibd,sm_meshMe
  write(*,*) 'FPres=',force_pres
  write(*,*) 'FVisc=',force_visc
  write(*,*) 'Momt =',moment_pres + moment_visc
#endif

  return

end subroutine sm_surf_assembleFluidForce_rigid
