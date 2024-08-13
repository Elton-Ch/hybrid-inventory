!234567
      subroutine count_err

      use global_data

      implicit none

      ! Variables
      integer   :: i, count
      real*8    :: sum1, sum2, L1_err1, L1_err2, L2_err1, L2_err2
      real*8    :: temp1, temp2, temp_obs, temp_state

      sum1 = 0.0
      sum2 = 0.0
      do i=1,ndat
        sum1 = sum1 + abs( dev_dat(i,1) ) ! average error
        sum2 = sum2 + abs( dev_dat(i,2) )
      enddo
      L1_err1 = sum1/float(ndat)
      L1_err2 = sum2/float(ndat)
      write(*,*) 'avrage error: ', L1_err1, L1_err2

      sum1 = 0.0
      sum2 = 0.0
      do i=1,ndat
        sum1 = sum1 + dev_dat(i,1)**2 ! RMS error
        sum2 = sum2 + dev_dat(i,2)**2
      enddo
      L2_err1 = sqrt( sum1/float(ndat) )
      L2_err2 = sqrt( sum2/float(ndat) )
      write(*,*) 'RMS error: ', L2_err1, L2_err2

      sum1 = 0.0
      sum2 = 0.0
      do i=1,ndat
        sum1 = sum1 + dev_dat(i,1) ! sistematic error
        sum2 = sum2 + dev_dat(i,2)
      enddo
      write(*,*) 'systematic error: ', sum1/float(ndat), sum2/float(ndat)

!======= chi_square==============
      temp1 = 0.0
      do i=1,ndat
        chi_sq(i,1) = (dev_dat(i,1)**2)/cdat(i,i)
        temp1 = temp1 + chi_sq(i,1)
      enddo
 
      temp2 = 0.0
      do i=1,ndat
        chi_sq(i,2) = (dev_dat(i,2)**2)/cdat(i,i)
        temp2 = temp2 + chi_sq(i,2)
      enddo

      temp_obs = 0.0
      do i=1,ndat
        temp_obs = temp_obs + chi_sq(i,2)
      enddo

      temp_state = 0.0
      do i=1,nsrc
        temp_state = temp_state + chi_sq_st(i)
      enddo
!================================
      write(*,*) '-------------------------------'
      write(*,*) (temp_obs + temp_state)/float(ndat), temp1/float(ndat), temp2/float(ndat)
      write(*,*) '(temp_obs + temp_state)/float(ndat)', temp_obs,temp_state,float(ndat)    ! for check
      write(*,*) '-------------------------------'

      open(1, file=trim(dir_kalman)//'obs_errors/error_'//year_of_obs//mon_of_obs//'.dat')
      write(1,*) 'average error   ', L1_err1, L1_err2
      write(1,*) 'RMS error       ', L2_err1, L2_err2
      write(1,*) 'systematic error', sum1/float(ndat), sum2/float(ndat)
      write(1,*) (temp_obs + temp_state)/float(ndat), temp1/float(ndat), temp2/float(ndat)
      write(1,*) '(temp_obs + temp_state)/float(ndat)', temp_obs,temp_state,float(ndat)    ! for check
      close(1)

 end subroutine count_err
