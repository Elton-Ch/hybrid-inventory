!===================================================================
      program kalmfilt

      use global_data

      implicit none

! Variables

     !integer count_num, i, j    ! count_num is now defined in global data
      integer                              :: i,j
      real*8, dimension(:,:),allocatable   :: sim_obs
      real*8                               :: temp_obs, temp_state, temp1, temp2
      character*6                          :: buff

! Body of kalmfilt

      call getarg(1,buff)
      read(buff,'(i4.4)') year
      call getarg(2,buff)
      read(buff,'(i2.2)') mon
      call getarg(3,buff)
      read(buff,'(i2.2)') count_num
      call getarg(4,buff)
      read(buff,'(i2.2)') lag_window
      call getarg(5,buff)
      read(buff,'(i4.4)') st_year
      call getarg(6,buff)
      read(buff,'(i2.2)') st_mon

      call getarg(7,buff)
      read(buff,'(i2.2)') n_select

      call getarg(8,dir_obs)
      call getarg(9,list_obs)
      call getarg(10,dir_sf_init)
      call getarg(11,file_config)
      call getarg(12,dir_kalman)


      
! definition of the state vector length
      !if(count_num.lt.lag_window)then
      !  monperinv = count_num
      !else
      !  monperinv = lag_window
      !endif
      monperinv = lag_window
 
      allocate( ndat_mon(lag_window))
      ndat_mon = 0
      call read_ndat           ! added to use all obs in window  mi20190801
      ndat = sum(ndat_mon)  ! added to use all obs in window  mi20190801
      write(*,*) 'kalman', year, mon, monperinv

      call read_config

      nsrc = reg_num*monperinv
      nsrc = nsrc + number_of_offset ! offset

      write(mon_of_obs,'(i2.2)')  mon
      write(year_of_obs,'(i4.4)') year

     !write(*,*) nsrc, ndat
      write(*,*) nsrc, ndat ,(ndat_mon(i),i=1,monperinv)  ! mi20190801
      
      call build_global_data

!============================ cycle ============================================
      allocate( sim_obs(ndat,2), chi_sq(ndat,2) )
      sim_obs = 0.0
      chi_sq = 0.0

      sim_obs(:,1) = matmul(pdm,src)
      dev_dat(:,1) = dat(:) - sim_obs(:,1)
      do i=1,ndat
        chi_sq(i,1) = (dev_dat(i,1)**2)/cdat(i,i)
      enddo

  !
  ! if(nsrc.ge.ndat) then
  !   write(*,*) 'calculation using kalman_1: nsrc >= ndat'
  !   call kalman_1 
  ! else
      write(*,*) 'calculation using kalman_2: nsrc < ndat'
      call kalman_2
  ! endif

      sim_obs(:,2) = matmul(pdm,src)
      dev_dat(:,2) = dat - sim_obs(:,2)
      do i=1,ndat
        chi_sq(i,2) = (dev_dat(i,2)**2)/cdat(i,i)
      enddo 
 
      call output_data
      call count_err
!============================ cycle ============================================

      call build_result_fld
	  
      deallocate (sim_obs, chi_sq, ndat_mon)

      end program kalmfilt
