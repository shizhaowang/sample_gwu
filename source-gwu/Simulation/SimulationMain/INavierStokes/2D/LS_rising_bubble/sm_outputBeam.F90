! output Euler beam information
       subroutine sm_outputBeam(fn,id)
         use sm_beamData
         implicit none
         character(len=*) :: fn
         integer :: id
         character(len=8) :: lid
         integer :: i
 
         write(lid,'(I8.8)') id
 
         open(123,file=fn//lid//'.dat')
           do i = 1, sm_beamNp
             !write(123,'(12f20.12)') sm_beamPos(i), sm_beamDisp(i), sm_beamLoad(i), sm_beamGrad(i) 
             write(123,'(12f20.12)') sm_beamPos(i), sm_beam_qn(i), sm_beamGrad(i), sm_beamPres(i)
           enddo
         close(123)
 
         return
       endSubroutine sm_outputBeam

