! Modified based on  SUBROUTINE sm_predcorrInteg:
! Shizhao Wang
! Aug 11, 2015

!======================================================================

subroutine sm_predcorrInteg_beam(ibd,dt,t)

  use sm_PredCorr_data, only: sm_PredCorr_type, sm_PredCorr_info

  use SolidMechanics_data, only : sm_structure,sm_BodyInfo 

  use sm_integinterface, only : sm_pcem, sm_pcmem, sm_pc_AB2, sm_pc_AM2,   &
                                sm_pc_AB3, sm_pc_AM3, sm_pc_AB4, sm_pc_AM4

  use Driver_interface, ONLY : Driver_abortFlash

  use sm_beamData, only : sm_beam_qn, sm_beam_qdn, sm_beam_qddn, sm_beam_qi, sm_beam_qdi, sm_beam_qddi, &
                          sm_beam_qms, sm_beam_qdms, sm_beam_qddms, sm_beamNp, sm_beamPos

  implicit none
#include "constants.h"
#include "sm_integrator.h"

  integer, intent(in) :: ibd
  real, intent(in) :: dt, t

  real, parameter :: CONST_ONE = 1.

  real :: e1
  integer :: neq
  integer :: i

  type(sm_structure), pointer :: BodyInfo
  type(sm_PredCorr_type), pointer :: integ


  ! Assign Pointers
  BodyInfo => sm_BodyInfo(ibd)
  integ => sm_PredCorr_info(ibd)
  
  ! number of unrestrained dofs
  neq = sm_beamNp+4

  if (neq .eq. 0) return ! Only unrestrained dofs are time-advanced here. 

  select case(integ%pciter) ! EULER = 1, AB2 = 2, AB3 = 3, HMG = 4

  case(SM_PCEULER)

     if (integ%pcflag .eq. SM_PCPREDICTOR) then
        ! Y,YDOT(-1)                         Y, YDOT (0)
        sm_beam_qms(:,-1) = sm_beam_qn(:)
        sm_beam_qdms(:,-1) = sm_beam_qdn(:)
        sm_beam_qddms(:,-1) = sm_beam_qddn(:)

        ! shift the dt's
        integ%vardt(-1) = integ%vardt(0)
        integ%vardt( 0) = dt


! ...   Predict the Solution ( Euler Method )
!                        Y0            F0                  pY1
        ! Positions q:
        call sm_pcem(neq,sm_beam_qms(1:neq,-1),sm_beam_qdms(1:neq,-1),dt,sm_beam_qn(1:neq))
        !call sm_pcem(neq,sm_beam_qms(1:neq,-1),sm_beam_qdms(1:neq,-1),dt,sm_beam_qn(1:neq))
        ! Velocities qd:
        call sm_pcem(neq,sm_beam_qdms(1:neq,-1),sm_beam_qddms(1:neq,-1),dt,sm_beam_qdn(1:neq))        

        integ%pcflag     = SM_PCCORRECTOR
        integ%pcconvflag = SM_PCNOTCONVERGED
        integ%pcerr      = CONST_ONE !!! We are not ok with an error of 1. !!!
        
     else

! ...   Correct the Solution ( Modified Euler Method )
!                         Y0               F0              kF1
        ! Positions q:
        call sm_pcmem (neq,sm_beam_qms(1:neq,-1),sm_beam_qdms(1:neq,-1),sm_beam_qdn(1:neq), &
                     dt, sm_beam_qi(1:neq))
!                                   k+1Y1       
        ! Velocities qd:
        call sm_pcmem (neq,sm_beam_qdms(1:neq,-1),sm_beam_qddms(1:neq,-1),sm_beam_qddn(1:neq), &
                     dt, sm_beam_qdi(1:neq))
!                                   k+1Y1                  
!       Check for convergence
!
!       e1 = L1 norm of k+1Y1-kY1 
        e1 = max(maxval(abs(sm_beam_qi(1:neq)-sm_beam_qn(1:neq))), &
                 maxval(abs(sm_beam_qdi(1:neq)-sm_beam_qdn(1:neq)))) 
        

!       kY1            k+1Y1  
        sm_beam_qn(1:neq) = sm_beam_qi(1:neq)
        sm_beam_qdn(1:neq) = sm_beam_qdi(1:neq)

!       Load convergence error for the body:        
        integ%pcerr = e1

        if (e1 .le. integ%pcepsilon) then
           integ%pcflag      = SM_PCPREDICTOR
           integ%pcconvflag  = SM_PCCONVERGED
           if (integ%pcmethod .gt. SM_PCEULER) integ%pciter = SM_PCAB2
        endif

     endif
        
  case(SM_PCAB2)

     if (integ%pcflag .eq. SM_PCPREDICTOR) then

        ! Y,YDOT(-2)                         Y, YDOT (-1)
        sm_beam_qms(1:neq,-2) = sm_beam_qms(1:neq,-1) ! n-2 Previous steps pos, vel, accel
        sm_beam_qdms(1:neq,-2) = sm_beam_qdms(1:neq,-1)
        sm_beam_qddms(1:neq,-2) = sm_beam_qddms(1:neq,-1)

        ! Y,YDOT(-1)                         Y, YDOT (0)
        sm_beam_qms(1:neq,-1) = sm_beam_qn(1:neq) ! Previous steps pos, vel, accel
        sm_beam_qdms(1:neq,-1) = sm_beam_qdn(1:neq)
        sm_beam_qddms(1:neq,-1) = sm_beam_qddn(1:neq)

        ! shift the dt's
        integ%vardt(-2) = integ%vardt(-1)
        integ%vardt(-1) = integ%vardt( 0)
        integ%vardt( 0) = dt

        ! ...   Predict the Solution ( AB2 )
        ! Positions q:
        call sm_pc_AB2(neq,sm_beam_qms(1:neq,-1),sm_beam_qdms(1:neq,-1),   &
                       sm_beam_qdms(1:neq,-2),integ%vardt(-1:0),1,sm_beam_qn(1:neq))
        ! Velocities qd:
        call sm_pc_AB2(neq,sm_beam_qdms(1:neq,-1),sm_beam_qddms(1:neq,-1), &
                       sm_beam_qddms(1:neq,-2),integ%vardt(-1:0),1,sm_beam_qdn(1:neq))        

        integ%pcflag     = SM_PCCORRECTOR
        integ%pcconvflag = SM_PCNOTCONVERGED
        integ%pcerr      = CONST_ONE !!! We are not ok with an error of 1. !!!
        
     else

        ! ...   Correct the Solution ( Adams-Moulton 2 )
        ! Positions q:
        call sm_pc_AM2(neq,sm_beam_qms(1:neq,-1),sm_beam_qdms(1:neq,-1), &
                           sm_beam_qdms(1:neq,-2),sm_beam_qdn(1:neq),    &
                           integ%vardt(-1:0),1, sm_beam_qi(1:neq))
      
        ! Velocities qd:
        call sm_pc_AM2 (neq, sm_beam_qdms(1:neq,-1),  sm_beam_qddms(1:neq,-1), &
                             sm_beam_qddms(1:neq,-2), sm_beam_qddn(1:neq), &
                             integ%vardt(-1:0),1, sm_beam_qdi(1:neq))

        ! Check for convergence
        ! e1 = L1 norm of k+1Y1-kY1 
        e1 = max(maxval(abs(sm_beam_qi(1:neq)-sm_beam_qn(1:neq))), &
                 maxval(abs(sm_beam_qdi(1:neq)-sm_beam_qdn(1:neq)))) 
        

        ! kY1            k+1Y1  
        sm_beam_qn(1:neq) = sm_beam_qi(1:neq)
        sm_beam_qdn(1:neq) = sm_beam_qdi(1:neq)

        ! Load convergence error for the body:        
        integ%pcerr = e1

        if (e1 .le. integ%pcepsilon) then
           integ%pcflag      = SM_PCPREDICTOR
           integ%pcconvflag  = SM_PCCONVERGED
           if (integ%pcmethod .gt. SM_PCAB2) integ%pciter = SM_PCAB3
        endif

     endif

  case(SM_PCAB3)

     if (integ%pcflag .eq. SM_PCPREDICTOR) then

        ! Y,YDOT(-3)                         Y, YDOT (-2)
        sm_beam_qms(1:neq,-3) = sm_beam_qms(1:neq,-2) ! n-3 Previous steps pos, vel, accel
        sm_beam_qdms(1:neq,-3) = sm_beam_qdms(1:neq,-2)
        sm_beam_qddms(1:neq,-3) = sm_beam_qddms(1:neq,-2)

        ! Y,YDOT(-2)                         Y, YDOT (-1)
        sm_beam_qms(1:neq,-2) = sm_beam_qms(1:neq,-1) ! n-2 Previous steps pos, vel, accel
        sm_beam_qdms(1:neq,-2) = sm_beam_qdms(1:neq,-1)
        sm_beam_qddms(1:neq,-2) = sm_beam_qddms(1:neq,-1)

        ! Y,YDOT(-1)                         Y, YDOT (0)
        sm_beam_qms(1:neq,-1) = sm_beam_qn(1:neq) ! Previous steps pos, vel, accel
        sm_beam_qdms(1:neq,-1) = sm_beam_qdn(1:neq)
        sm_beam_qddms(1:neq,-1) = sm_beam_qddn(1:neq)

        ! shift the dt's
        integ%vardt(-3) = integ%vardt(-2)
        integ%vardt(-2) = integ%vardt(-1)
        integ%vardt(-1) = integ%vardt( 0)
        integ%vardt( 0) = dt

        ! ...   Predict the Solution ( AB3 )
        ! Positions q:
        call sm_pc_AB3(neq,sm_beam_qms(1:neq,-1),sm_beam_qdms(1:neq,-1),   &
                           sm_beam_qdms(1:neq,-2), sm_beam_qdms(1:neq,-3), &
                           integ%vardt(-2:0),2,sm_beam_qn(1:neq))
        ! Velocities qd:
        call sm_pc_AB3(neq, sm_beam_qdms(1:neq,-1),sm_beam_qddms(1:neq,-1),  &
                            sm_beam_qddms(1:neq,-2),sm_beam_qddms(1:neq,-3), &
                            integ%vardt(-2:0),2,sm_beam_qdn(1:neq))        

        integ%pcflag     = SM_PCCORRECTOR
        integ%pcconvflag = SM_PCNOTCONVERGED
        integ%pcerr      = CONST_ONE !!! We are not ok with an error of 1. !!!
        
     else

        ! ...   Correct the Solution ( Adams-Moulton 3 )
        ! Positions q:
        call sm_pc_AM3(neq, sm_beam_qms(1:neq,-1),sm_beam_qdms(1:neq,-1), &
                            sm_beam_qdms(1:neq,-2), sm_beam_qdms(1:neq,-3), &
                            sm_beam_qdn(1:neq),    &
                            integ%vardt(-2:0),2, sm_beam_qi(1:neq))
      
        ! Velocities qd:
        call sm_pc_AM3(neq, sm_beam_qdms(1:neq,-1),  sm_beam_qddms(1:neq,-1), &
                            sm_beam_qddms(1:neq,-2), sm_beam_qddms(1:neq,-3), &
                            sm_beam_qddn(1:neq), &
                            integ%vardt(-2:0),2, sm_beam_qdi(1:neq))

        ! Check for convergence
        ! e1 = L1 norm of k+1Y1-kY1 
        e1 = max(maxval(abs(sm_beam_qi(1:neq)-sm_beam_qn(1:neq))), &
                 maxval(abs(sm_beam_qdi(1:neq)-sm_beam_qdn(1:neq)))) 
        
        ! kY1            k+1Y1  
        sm_beam_qn(1:neq) = sm_beam_qi(1:neq)
        sm_beam_qdn(1:neq) = sm_beam_qdi(1:neq)

        ! Load convergence error for the body:        
        integ%pcerr = e1

        if (e1 .le. integ%pcepsilon) then
           integ%pcflag      = SM_PCPREDICTOR
           integ%pcconvflag  = SM_PCCONVERGED
           if (integ%pcmethod .gt. SM_PCAB3) integ%pciter = SM_PCAB4
        endif

     endif

  case(SM_PCAB4)        

     if (integ%pcflag .eq. SM_PCPREDICTOR) then

        ! Y,YDOT(-4)                         Y, YDOT (-3)
        sm_beam_qms(1:neq,-4) = sm_beam_qms(1:neq,-3) ! n-4 Previous steps pos, vel, accel
        sm_beam_qdms(1:neq,-4) = sm_beam_qdms(1:neq,-3)
        sm_beam_qddms(1:neq,-4) = sm_beam_qddms(1:neq,-3)

        ! Y,YDOT(-3)                         Y, YDOT (-2)
        sm_beam_qms(1:neq,-3) = sm_beam_qms(1:neq,-2) ! n-3 Previous steps pos, vel, accel
        sm_beam_qdms(1:neq,-3) = sm_beam_qdms(1:neq,-2)
        sm_beam_qddms(1:neq,-3) = sm_beam_qddms(1:neq,-2)

        ! Y,YDOT(-2)                         Y, YDOT (-1)
        sm_beam_qms(1:neq,-2) = sm_beam_qms(1:neq,-1) ! n-2 Previous steps pos, vel, accel
        sm_beam_qdms(1:neq,-2) = sm_beam_qdms(1:neq,-1)
        sm_beam_qddms(1:neq,-2) = sm_beam_qddms(1:neq,-1)

        ! Y,YDOT(-1)                         Y, YDOT (0)
        sm_beam_qms(1:neq,-1) = sm_beam_qn(1:neq) ! Previous steps pos, vel, accel
        sm_beam_qdms(1:neq,-1) = sm_beam_qdn(1:neq)
        sm_beam_qddms(1:neq,-1) = sm_beam_qddn(1:neq)

        integ%vardt(-4) = integ%vardt(-3)
        integ%vardt(-3) = integ%vardt(-2)
        integ%vardt(-2) = integ%vardt(-1)
        integ%vardt(-1) = integ%vardt( 0)
        integ%vardt( 0) = dt

        ! ...   Predict the Solution ( AB4 )
        ! Positions q:
        call sm_pc_AB4(neq, sm_beam_qms(1:neq,-1),sm_beam_qdms(1:neq,-1),   &
                            sm_beam_qdms(1:neq,-2), sm_beam_qdms(1:neq,-3), &
                            sm_beam_qdms(1:neq,-4), &
                            integ%vardt(-3:0),3,sm_beam_qn(1:neq))
        ! Velocities qd:
        call sm_pc_AB4(neq, sm_beam_qdms(1:neq,-1),sm_beam_qddms(1:neq,-1),  &
                            sm_beam_qddms(1:neq,-2),sm_beam_qddms(1:neq,-3), &
                            sm_beam_qddms(1:neq,-4), &
                            integ%vardt(-3:0),3,sm_beam_qdn(1:neq))  

        integ%pcflag     = SM_PCCORRECTOR
        integ%pcconvflag = SM_PCNOTCONVERGED
        integ%pcerr      = CONST_ONE !!! We are not ok with an error of 1. !!!
        
     else

        ! ...   Correct the Solution ( Adams-Moulton 4 )
        ! Positions q:
        call sm_pc_AM4(neq, sm_beam_qms(1:neq,-1),sm_beam_qdms(1:neq,-1), &
                            sm_beam_qdms(1:neq,-2), sm_beam_qdms(1:neq,-3), &
                            sm_beam_qdms(1:neq,-4), &
                            sm_beam_qdn(1:neq),     &
                            integ%vardt(-3:0),3, sm_beam_qi(1:neq))
      
        ! Velocities qd:
        call sm_pc_AM4(neq, sm_beam_qdms(1:neq,-1),  sm_beam_qddms(1:neq,-1), &
                            sm_beam_qddms(1:neq,-2), sm_beam_qddms(1:neq,-3), &
                            sm_beam_qddms(1:neq,-4), &
                            sm_beam_qddn(1:neq), &
                            integ%vardt(-3:0),3, sm_beam_qdi(1:neq))
        
        ! Check for convergence
        ! e1 = L1 norm of k+1Y1-kY1 
        e1 = max(maxval(abs(sm_beam_qi(1:neq)-sm_beam_qn(1:neq))), &
                 maxval(abs(sm_beam_qdi(1:neq)-sm_beam_qdn(1:neq)))) 
        
        ! kY1            k+1Y1  
        sm_beam_qn(1:neq) = sm_beam_qi(1:neq)
        sm_beam_qdn(1:neq) = sm_beam_qdi(1:neq)

        ! Load convergence error for the body:        
        integ%pcerr = e1

        !if (e1 .le. integ%pcepsilon) then
        !   integ%pcflag      = SM_PCPREDICTOR
        !   integ%pcconvflag  = SM_PCCONVERGED
           !!if (integ%pcmethod .gt. SM_PCAB4) integ%pciter = SM_PCAB5
        !endif

        ! Shizhao
        write(*,*) 'Warning: forced to exit FSI, error:', e1
           integ%pcflag      = SM_PCPREDICTOR
           integ%pcconvflag  = SM_PCCONVERGED

     endif
        

  case default
     
     call Driver_abortFlash("Error in sm_predcorrInteg: sm_pciter not equal 1-4.")
        
  end select

  !do i = 1, sm_beamNp+4
  !  sm_beam_qn(i) = 0.1*sin(2.0*PI*sm_beamPos(i)/4.0)*cos(2.0*PI*t*10.0)
  !enddo

  !integ%pcconvflag  = SM_PCCONVERGED

  call sm_beamSetBC()
  
  do i = 1, sm_beamNp 
    BodyInfo%yB(i) = sm_beam_qn(i+2) + 0.005 
    BodyInfo%yB(2*sm_beamNp+1-i) = sm_beam_qn(i+2)-0.005 
  enddo

  nullify(BodyInfo)
  nullify(integ)

  return

end subroutine sm_predcorrInteg_beam
