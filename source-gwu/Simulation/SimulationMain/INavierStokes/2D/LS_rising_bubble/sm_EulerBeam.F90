! Solve the linear Euler Beam equation
! Shizhao Wang
! Jul 22, 2015

! =====================================================      
      subroutine sm_EulerBeam(ibd, dt) 

        use SolidMechanics_data, only : sm_structure,sm_BodyInfo
        
        use sm_beamData
        implicit none
        integer :: ibd
        real :: dt
        integer :: i, neq
        integer, allocatable, dimension(:) :: indx
        real, allocatable, dimension(:) :: rhs
        real, allocatable, dimension(:,:) :: mat
        real :: d
        type(sm_structure), pointer :: BodyInfo

        sm_beam_qms(:,-4) = sm_beam_qms(:,-3)
        sm_beam_qms(:,-3) = sm_beam_qms(:,-2)
        sm_beam_qms(:,-2) = sm_beam_qms(:,-1)
        sm_beam_qms(:,-1) = sm_beam_qn(:)

        sm_beamDt = dt

        neq = sm_beamNp
        allocate( indx(sm_beamNp), rhs(sm_beamNp) )
        allocate( mat(neq, neq) )

        ! Matrix
        sm_beamMat = 0.0
        do i = 3, sm_beamNp-2
          sm_beamMat(i-2,i) =  1.0d0
          sm_beamMat(i-1,i) = -4.0d0
          sm_beamMat(i,  i) =  6.0d0
          sm_beamMat(i+1,i) = -4.0d0
          sm_beamMat(i+2,i) =  1.0d0
        enddo
        
        ! Modified matrix at boundry
        !call sm_beamSetMatBC(2,2)
        call sm_beamSetMatBC(1,3)
        

        ! Coefficient
        sm_beamMat = sm_beamMat*sm_beamModulus*sm_beamInertia &
                     /sm_beamRhou/sm_beamArea* & 
                     (sm_beamDt/sm_beamDh)*(sm_beamDt/sm_beamDh) &
                     /sm_beamDh/sm_beamDh

        ! Dynamic term
        do i = 1, sm_beamNp
          sm_beamMat(i,i) = sm_beamMat(i,i) + 1.0d0
        enddo

!        call sm_outputMat(neq, sm_beamMat, 'mat')

        ! Set the load
        do i = 1, sm_beamNp
          rhs(i) = sm_beamPres(i)/sm_beamRhou/sm_beamArea*sm_beamDt*sm_beamDt &
                 + 2.0d0*sm_beam_qms(i,-1) - sm_beam_qms(i,-2)
!          write(1000,*) i, rhs(i)  
        enddo

        mat = transpose(sm_beamMat)

!        call sm_outputMat(neq, mat, 'mat-tran')

        call sm_ludcmp(mat,neq,neq,indx,d)

!        call sm_outputMat(neq, mat, 'mat-lu')

        call sm_lubksb(mat,neq,neq,indx,rhs)

        do i = 1, sm_beamNp
          sm_beam_qn(i) = rhs(i)
!         write(2000,*) i, rhs(i)  
        enddo
        
        sm_beam_qdn(:) = (sm_beam_qn(:) - sm_beam_qms(:,-1))/sm_beamDt

        BodyInfo => sm_BodyInfo(ibd)
        do i = 1, sm_beamNp
          BodyInfo%yB(i) = sm_beam_qn(i) + 0.005
          BodyInfo%yB(2*sm_beamNp+1-i) = sm_beam_qn(i)-0.005
        enddo
                      
        deallocate( indx, rhs, mat )
        nullify(BodyInfo)
 
        return
      endsubroutine sm_EulerBeam  

! =====================================================      
! =====================================================      
