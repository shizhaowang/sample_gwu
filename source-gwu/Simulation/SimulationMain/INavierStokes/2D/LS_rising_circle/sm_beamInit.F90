! Initial beam parameters
! Shizhao Wang
! Aug 11, 2015

      subroutine sm_beamInit(fn)
        use sm_beamData
        implicit none
        character(len=*) :: fn
        integer :: i
        real :: pi, beta, a, x

        open(123,file=fn)
          read(123,*)
          read(123,*)
          read(123,*) sm_beamNp
          read(123,*)
          read(123,*) sm_beamDh
          read(123,*)
          read(123,*) sm_beamModulus
          read(123,*)
          read(123,*) sm_beamInertia
          read(123,*)
          read(123,*) sm_beamArea
          read(123,*)
          read(123,*) sm_beamRhou
        close(123)
        write(*,*) 'No. of points on beam', sm_beamNp
        write(*,*) 'length of beam element', sm_beamDh
        write(*,*) 'Young modulus', sm_beamModulus
        write(*,*) 'ineria of cross section', sm_beamInertia
        write(*,*) 'area of cross section', sm_beamArea
        write(*,*) 'densith of beam', sm_beamRhou
        write(*,*) ''

        allocate(sm_beam_qn(sm_beamNp),sm_beam_qdn(sm_beamNp),sm_beam_qddn(sm_beamNp)) 
        allocate(sm_beam_qi(sm_beamNp),sm_beam_qdi(sm_beamNp),sm_beam_qddi(sm_beamNp)) 
        allocate(sm_beam_qms(sm_beamNp,-4:-1),sm_beam_qdms(sm_beamNp,-4:-1),sm_beam_qddms(sm_beamNp,-4:-1)) 
        allocate(sm_beamPos(sm_beamNp))
        allocate(sm_beamShear(sm_beamNp))
        allocate(sm_beamPres(sm_beamNp))
        allocate(sm_beamGrad(sm_beamNp))
        allocate(sm_beamMat(sm_beamNp,sm_beamNp))

        sm_beam_qn = 0.0d0; sm_beam_qdn = 0.0d0; sm_beam_qddn = 0.0d0
        sm_beam_qi = 0.0d0; sm_beam_qdi = 0.0d0; sm_beam_qddi = 0.0d0
        sm_beam_qms = 0.0d0; sm_beam_qdms = 0.0d0; sm_beam_qddms = 0.0d0
        sm_beamShear = 0.0d0
        sm_beamPres = 0.0d0
        sm_beamGrad = 0.0d0
        sm_beamMat = 0.0d0

        do i = 1, sm_beamNp
          sm_beamPos(i) = (i-1)*sm_beamDh
        enddo

        ! Set the initial displacement
        pi = acos(-1.0d0)
        beta = 1.875d0
        !beta = 4.694d0
        !beta = 7.855d0
        !beta = 10.9955d0
!        a = 0.01/(cosh(beta) - cos(beta) + & 
!           ((cos(beta)+cosh(beta))*(sin(beta)-sinh(beta)))/(sin(beta)+sinh(beta)))

!        do i = 1, sm_beamNp
!          x = sm_beamPos(i)
!          sm_beam_qn(i) = a*(cosh(beta*x) - cos(beta*x) + &
!            ((cos(beta)+cosh(beta))*(sin(beta*x)-sinh(beta*x)))/(sin(beta)+sinh(beta)))
!!          sm_beam_qn(i) = 0.1d0*sin(2.0d0*pi*x/4.0d0)
!          enddo

        sm_beam_qms(:,-1) = sm_beam_qn(:)
        sm_beam_qms(:,-2) = sm_beam_qn(:)
        sm_beam_qms(:,-3) = sm_beam_qn(:)
        sm_beam_qms(:,-4) = sm_beam_qn(:)

        return

      endSubroutine sm_beamInit     
