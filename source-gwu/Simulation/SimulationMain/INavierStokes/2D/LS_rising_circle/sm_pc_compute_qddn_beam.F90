! Modified based on sm_pc_compute_qddn_beam
! Shizhao Wang
! Aug 11, 2015

#include "SolidMechanics.h"

subroutine sm_pc_compute_qddn_beam(ibd)

  use sm_beamData

  implicit none
  integer :: ibd, i
  
  sm_beamGrad = 0.0
  ! Gradient
  do i = 3, sm_beamNp+2
    sm_beamGrad(i) = (sm_beam_qn(i-2)-4.0*sm_beam_qn(i-1) &
                   +6.0*sm_beam_qn(i)-4.0*sm_beam_qn(i+1) &
                   +sm_beam_qn(i+2))/(sm_beamDh**2)
  enddo

  sm_beam_qddn = (-sm_beamGrad*sm_beamModulus*sm_beamInertia + sm_beamPres*sm_beamDh*sm_beamDh)/sm_beamRhou/sm_beamArea ! to reduce the error. The time step will be dt/sm_beamDh

  return
end subroutine sm_pc_compute_qddn_beam

