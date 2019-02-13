! Setup pressure load on beam
! Shizhao 
! Aug 04, 2015

#include "Flash.h"
#include "SolidMechanics.h"
#include "constants.h"

subroutine sm_beamLoad(ibd)

  use SolidMechanics_Data, only : sm_meshMe,sm_meshComm,sm_NumBodies,sm_BodyInfo
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_sbData, only : gr_sbBodyInfo
  use sm_beamData, only: sm_beamNp, sm_beamPres

  implicit none

  integer :: i,ibd
  integer :: ip1, ip2
  real :: p1, p2, p3, p4

!  allocate( sm_beamPres(sm_beamNP+4) ) 
!  sm_beamPres = 0.0d0
  
  !! Loop over bodies and write:
!  do ibd = 1,sm_NumBodies
     if (sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster)then
       do i = 4, sm_beamNP+1
         ip1 = i - 3
         p1 = gr_sbBodyInfo(ibd)%particles(PRES_PART_PROP,ip1)
         p2 = gr_sbBodyInfo(ibd)%particles(PRES_PART_PROP,ip1+1)
         ip2 = sm_beamNp*2 + 3 - i 
         p3 = gr_sbBodyInfo(ibd)%particles(PRES_PART_PROP,ip2)
         p4 = gr_sbBodyInfo(ibd)%particles(PRES_PART_PROP,ip2-1)
         sm_beamPres(i) = -0.5d0*(p1+p2) + 0.5*(p3+p4)
       enddo
       ip1 = sm_beamNP+2-3 
       p2 = gr_sbBodyInfo(ibd)%particles(PRES_PART_PROP,ip1)
       ip2 = sm_beamNP+1 
       p4 = gr_sbBodyInfo(ibd)%particles(PRES_PART_PROP,ip2)
       sm_beamPres(sm_beamNP+2) = -p2 + p4
     endif
!  enddo
  

  ! Testing
  sm_beamPres = 0.0d0
  write(*,*) 'Warning: beam load is forced to be 0.0'
  return

end subroutine sm_beamLoad
