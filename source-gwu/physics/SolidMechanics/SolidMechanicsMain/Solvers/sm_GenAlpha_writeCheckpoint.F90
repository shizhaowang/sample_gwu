!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/sm_GenAlpha_writeCheckpoint.F90
!!
!! NAME
!!  stub
!!
!!
!! SYNOPSIS
!!
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!!
!!
!!***
#include "SolidMechanics.h"

subroutine sm_GenAlpha_writeCheckpoint(ibd, sm_checkpt_num)
  use Driver_interface,    only: Driver_abortFlash
  implicit none

  ! IO Varialbes
  integer, intent(in) :: ibd
  integer, intent(in) :: sm_checkpt_num

  call Driver_abortFlash('Stub')
  
  return

end subroutine sm_GenAlpha_writeCheckpoint
