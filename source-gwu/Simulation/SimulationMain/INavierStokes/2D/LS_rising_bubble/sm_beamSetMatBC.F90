! =====================================================      
      subroutine sm_beamSetMatBC(left,right) 
        
        use sm_beamData
        implicit none
        integer :: left, right
        integer :: i

        ! left,right == 1, fixed
        ! left,right == 2, hinge
        ! left,right == 3, free
        
        if(left == 1) then
        ! Boundary conditions
        ! Fixed at left, displacement
        sm_beamMat(1,1) = 1.0d0
        sm_beamPres(1)  = 0.0d0
        !sm_beamMat(1,1) =  6.0d0
        !sm_beamMat(2,1) = -8.0d0
        !sm_beamMat(3,1) =  2.0d0

        sm_beamMat(1,2) = -4.0d0
        sm_beamMat(2,2) =  7.0d0
        sm_beamMat(3,2) = -4.0d0
        sm_beamMat(4,2) =  1.0d0
        
        elseif(left ==2) then
        ! hinge at left
        sm_beamMat(1,1) =  6.0d0
        sm_beamMat(2,1) =  0.0d0
        sm_beamMat(3,1) =  0.0d0

        sm_beamMat(1,2) = -4.0d0
        sm_beamMat(2,2) =  5.0d0
        sm_beamMat(3,2) = -4.0d0
        sm_beamMat(4,2) =  1.0d0
        
        elseif(left == 3) then
          write(*,*) 'The free end at left is not supported'
          stop
        endif
       
        ! right end
        if (right == 1) then
          write(*,*) 'Fixed at right end is not supported'
          stop
        elseif( right == 2) then
          i = sm_beamNp      
          sm_beamMat(i,  i) =  6.0d0

          i = sm_beamNp - 1      
          sm_beamMat(i-2,i) =  1.0d0
          sm_beamMat(i-1,i) = -4.0d0
          sm_beamMat(i,  i) =  5.0d0
          sm_beamMat(i+1,i) = -4.0d0

        elseif( right == 3) then
          i = sm_beamNp      
          sm_beamMat(i-2,i) =  2.0d0
          sm_beamMat(i-1,i) = -4.0d0
          sm_beamMat(i,  i) =  2.0d0

          i = sm_beamNp - 1      
          sm_beamMat(i-2,i) =  1.0d0
          sm_beamMat(i-1,i) = -4.0d0
          sm_beamMat(i,  i) =  5.0d0
          sm_beamMat(i+1,i) = -2.0d0

        endif

        return
      endsubroutine sm_beamSetMatBC  

! =====================================================      
! =====================================================      
