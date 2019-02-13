! SolidMechanicsMain/Solvers/PredCorr
! Modified based on sm_PredCorr_advance.F90
! Shizhao Wang
! Aug 11, 2015

#include "SolidMechanics.h"

subroutine sm_PredCorr_advance_beam(ibd,restart)
  
  use SolidMechanics_data, only: sm_BodyInfo, sm_meshMe
  use sm_PredCorr_data, only: sm_PredCorr_info
  use sm_integinterface, only: sm_pc_compute_qddn, sm_predcorrInteg
  use sm_pk_interface, only: sm_pk_apply, sm_pk_updatekinematics_rigid, &
                             sm_pk_angvelconstraint_rigid
  use sm_assemble_interface, only: sm_assemble_mass, sm_assemble_IntForce, sm_assemble_ExtForce
  use Driver_interface, only: Driver_getSimTime
  use sm_beamData, only : sm_beamDh

  implicit none
#include "sm_integrator.h"
  ! IO Variables
  integer, intent(in) :: ibd
  logical, intent(in) :: restart

  ! Internal Variables
  real :: time
  integer, save :: first_call=1

  if(sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then

  ! Get the current DT
  call Driver_getDT(sm_PredCorr_info(ibd)%dt)

  call Driver_getSimTime(time)

  if(first_call == 1) then
    first_call = 0
  else
    call sm_beamLoad(ibd)
  endif
!  call sm_pc_compute_qddn_beam(ibd) 
!  call sm_predcorrInteg_beam(ibd,sm_PredCorr_info(ibd)%dt/sm_beamDh,time)

   call sm_EulerBeam(ibd, sm_PredCorr_info(ibd)%dt)

  endif

  return
end subroutine sm_PredCorr_advance_beam

