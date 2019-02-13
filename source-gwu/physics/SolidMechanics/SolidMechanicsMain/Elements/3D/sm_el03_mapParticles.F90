!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Elements/sm_el02_mapParticles
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

#include "Flash.h"
#include "SolidMechanics.h"
#include "constants.h"

subroutine sm_el03_mapParticles(body, e, ptelem,  &
                                xpos,ypos,zpos,   &
                                xvel,yvel,zvel,   &
                                xacc,yacc,zacc,   &
                                xnrm,ynrm,znrm,   &
                                areai, loc_num )
  use SolidMechanics_data, only: sm_structure
  use Driver_interface, only: Driver_abortFlash
  use sm_misc_interface,only : sm_crossprod
  implicit none
  
  ! IO Variables
  type(sm_structure)   :: body     ! entire body structure
  integer, intent(in)  :: e        ! element number
  integer, intent(in)  :: ptelem
  real, dimension(ptelem) :: xpos, ypos, zpos, xvel, yvel, zvel, xacc, yacc, zacc, xnrm, ynrm, znrm, areai
  integer, dimension(ptelem) :: loc_num

  ! Internal variables
  integer, parameter :: nen_e = FOUR_NODES
  integer :: a, i, ii, j, idx1, idx2, nEta, tmp
  real, dimension(nen_e,NDIM) :: u, ud, udd
  real, dimension(NDIM) :: norm,v12,v13, tmpv
  real :: del,delp,area,area_sub
  real :: ei,ej,N1,N2,N3
  real :: phi
  !
  ! Get absolute position of each node in the element
  !
  do a = 1,nen_e
     idx1 = body%ws_IEN(a,e)
     u(a,1) = body%x(idx1)
     u(a,2) = body%y(idx1)
     u(a,3) = body%z(idx1)
!     write(4100,*) u(a,:) 
     if (Body%BodyType.eq.BODYTYPE_RBC) u(a,1:NDIM)=0.;

     do i = 1,NDIM
        idx2     = body%ID(i,idx1)
        u(a,i)   = u(a,i) + body%qn(idx2)
        ud(a,i)  = body%qdn(idx2)
        udd(a,i) = body%qddn(idx2)
     end do
!     write(4200,*) u(a,:) 
  enddo

  ! nXi = nEta
  nEta  = body%ws_nXi(e)
  
  !Get the normals
  v12(1)=u(2,IAXIS)-u(1,IAXIS);
  v12(2)=u(2,JAXIS)-u(1,JAXIS);
  v12(3)=u(2,KAXIS)-u(1,KAXIS);
  v13(1)=u(3,IAXIS)-u(1,IAXIS);
  v13(2)=u(3,JAXIS)-u(1,JAXIS);
  v13(3)=u(3,KAXIS)-u(1,KAXIS);

  norm=sm_crossprod(v12,v13);
  
  area=0.5*sqrt(sum(norm*norm));
  norm=norm/(2.*area);

  v12(1)=u(3,IAXIS)-u(1,IAXIS);
  v12(2)=u(3,JAXIS)-u(1,JAXIS);
  v12(3)=u(3,KAXIS)-u(1,KAXIS);
  v13(1)=u(4,IAXIS)-u(1,IAXIS);
  v13(2)=u(4,JAXIS)-u(1,JAXIS);
  v13(3)=u(4,KAXIS)-u(1,KAXIS);

  tmpv=sm_crossprod(v12,v13);
  
  area= area + 0.5*sqrt(sum(tmpv*tmpv));

  area_sub = area/real(ptelem);
  del=1./real(nEta);

  ii = 1
  ! x,y,z positions of internal particles:
  xpos(ii) = 0.25d0*(u(1,IAXIS) + u(2,IAXIS) + u(3,IAXIS) + u(4,IAXIS));
  ypos(ii) = 0.25d0*(u(1,JAXIS) + u(2,JAXIS) + u(3,JAXIS) + u(4,JAXIS));
  zpos(ii) = 0.25d0*(u(1,KAXIS) + u(2,KAXIS) + u(3,KAXIS) + u(4,KAXIS));

  ! Vel
  xvel(ii) = 0.25d0*(ud(1,IAXIS) + ud(2,IAXIS) + ud(3,IAXIS) + ud(4,IAXIS));
  yvel(ii) = 0.25d0*(ud(1,JAXIS) + ud(2,JAXIS) + ud(3,JAXIS) + ud(4,JAXIS));
  zvel(ii) = 0.25d0*(ud(1,KAXIS) + ud(2,KAXIS) + ud(3,KAXIS) + ud(4,KAXIS));

  ! Acc
  xacc(ii) = 0.25d0*(udd(1,IAXIS) + udd(2,IAXIS) + udd(3,IAXIS) + udd(4,IAXIS));
  yacc(ii) = 0.25d0*(udd(1,JAXIS) + udd(2,JAXIS) + udd(3,JAXIS) + udd(4,JAXIS));
  zacc(ii) = 0.25d0*(udd(1,KAXIS) + udd(2,KAXIS) + udd(3,KAXIS) + udd(4,KAXIS));

  xnrm(ii) = norm(IAXIS) ;
  ynrm(ii) = norm(JAXIS) ;
  znrm(ii) = norm(KAXIS) ; 
  loc_num(ii)=ii ;
  areai(ii) = area_sub;


  return

end subroutine sm_el03_mapParticles
