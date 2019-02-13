!     
! File:   sm_pk_init.F90
! Author: tim
!
! 

#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_pk_init()
     
  use sm_pk_interface, only: sm_pk_ReadHDF5
  use sm_pk_data, only: sm_pk_timedelay
  use SolidMechanics_data, only : sm_MeshMe
  use SolidMechanics_data, only : sm_numBodies
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  
  implicit none
  integer :: ibd

  if(sm_meshMe .eq. MASTER_PE) write(*,'(A)') 'Prescribed Kinematics Init.'

  ! read in the file containing the parmeters for the kinematics functions
  call sm_pk_ReadHDF5('kinematics.input.h5')
  
  ! set the offset time per the flash.par option
  call RuntimeParameters_get("sm_pk_timedelay",sm_pk_timedelay)

  ! do other things as needed.
  do ibd=1,sm_NumBodies
    call sm_pk_readPara(ibd)
  enddo

  return

end subroutine sm_pk_init
