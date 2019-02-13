! Data for Eulerian beam
! Shizhao Wang
! Aug 11, 2015

!
      module sm_beamData

       implicit none
  
#include "constants.h"

        integer, save :: sm_beamNp
        real, save :: sm_beamDt, sm_beamDh
        real, save :: sm_beamModulus, sm_beamInertia, sm_beamRhou, sm_beamArea
        real, save, allocatable, dimension(:) :: sm_beam_qn,sm_beam_qdn, &
            sm_beam_qddn, sm_beam_qi, sm_beam_qdi, sm_beam_qddi
        real, save, allocatable, dimension(:,:) :: sm_beam_qms, sm_beam_qdms, sm_beam_qddms
        real, save, allocatable, dimension(:) :: sm_beamPos, sm_beamShear, sm_beamPres
        real, save, allocatable, dimension(:) :: sm_beamGrad
        real, save, allocatable, dimension(:,:) :: sm_beamMat

      end module sm_beamData
