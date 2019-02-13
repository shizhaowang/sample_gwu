
#include "Flash.h"
#include "constants.h"

subroutine sm_pk_readPara(ibd)

  use SolidMechanics_data, only :  sm_BodyInfo, sm_structure
  use Driver_interface, only: Driver_abortFlash
  implicit none

  ! IO
  integer, intent(in) :: ibd

  ! Internal variables
  type(sm_structure),  pointer :: body

  real, allocatable, dimension(:,:) :: params
  integer, allocatable, dimension(:) :: nparams, rtype
  character(len=6) :: libd
  integer :: ircoord

  ! Get the body
  body => sm_BodyInfo(ibd)

  allocate ( params(Body%maxrestparams,Body%nrcoords) )
  allocate ( nparams(Body%nrcoords) )
  allocate ( rtype(Body%nrcoords) )

  write(libd,'(I6.6)') ibd
  open(123, file='sm_body_pk.'//libd)
    read(123,*) 
    do ircoord = 1, Body%nrcoords
      read(123,*) 
      read(123,*) 
      read(123,*) rtype(ircoord)
      read(123,*) 
      read(123,*) nparams(ircoord)
      read(123,*) 
      read(123,*) params(:,ircoord)
    enddo
  close(123)

  do ircoord = 1,Body%nrcoords
    write(*,*) 'parameters before:',body%restraints(ircoord)%param(1:Body%maxrestparams), ircoord
    body%restraints(ircoord)%restype = rtype(ircoord)
    body%restraints(ircoord)%nparam = nparams(ircoord)
    body%restraints(ircoord)%param(1:Body%maxrestparams)=params(1:Body%maxrestparams,ircoord)
    write(*,*) 'parameters after:',body%restraints(ircoord)%param(1:Body%maxrestparams), ircoord
  enddo

  ! deallocate
  deallocate(params)
  deallocate(rtype, nparams)

end subroutine sm_pk_readPara
