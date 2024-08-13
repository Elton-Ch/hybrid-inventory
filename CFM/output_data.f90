      subroutine output_data

      use global_data

      implicit none

     ! Variables
      integer                       :: i, k_old, k_new
      real*8                        :: work_old, work_new
                    
      integer                                :: ic, j,m,irec
      real                                   :: pri,post
      integer                                :: n_site
      integer                                :: iflg
      character*2                            :: temp_mon
      character*4                            :: temp_year
      integer                                :: iyear,imonth,iday,ihour,imin 
      real                                   :: temp_lon, temp_lat,temp_alt,temp_val,rsd
      character*3                            :: temp_ch3      ! for site with 13 characters
      character*10                           :: temp_date10   ! for date YYYYMMDDGG

      real*4,dimension( : ), allocatable     :: temp_pdm, temp_reg, temp_prsb, temp_bkrn
      real                                   :: bg,pri_local,post_local
      real                                   :: pri_work,post_work

     ! Body of output_data

     k_old = 0
     do i=1,ndat
       if(abs( dev_dat(i,1) ).gt.(2.*sqrt(cdat(i,i))) ) k_old = k_old + 1
     enddo
     work_old = maxval( abs( dev_dat(:,1) ) )

     k_new = 0
     do i=1,ndat
       if(abs( dev_dat(i,2) ).gt.(2.*sqrt(cdat(i,i))) ) k_new = k_new + 1
     enddo
     work_new = maxval( abs( dev_dat(:,2) ) )

     write(*,*) k_old, work_old
     write(*,*) k_new, work_new
     open(1, file=trim(dir_kalman)//'obs_profit/obs_'//year_of_obs//mon_of_obs//'.txt')
     write(1,*) ndat
     do i=1,ndat
!       write(1,'(i4,2x,i2,2x,i2,2x,i2,2x,i2,2x,4f15.4,2x,a13,2x,i3,3x,5f15.4)') obs_year(i), obs_month(i),obs_day(i),obs_hour(i),obs_min(i),&
!                                                         obs_lat(i), obs_lon(i), obs_alt(i), obs_val(i),st_name(i),obs_site(i),&
!                                                         dev_dat(i,1), dev_dat(i,2), sqrt(cdat(i,i)),chi_sq(i,1), chi_sq(i,2)
        write(1,'(a3,2x,a10,2x,6f15.10)') st_name(i),obs_date(i), obs_val(i),dev_dat(i,1), dev_dat(i,2), sqrt(cdat(i,i)),chi_sq(i,1), chi_sq(i,2)
 
     enddo
     close(1)

!======================================================
     
      do m=monperinv,1,-1 ! read obs in the previous and current months (old-->new)

        ic = mon - (m-1)
        if(ic.le.0) then
          write(temp_mon,'(i2.2)') ic + 12
          write(temp_year,'(i4.4)') year - 1
        else
          write(temp_mon,'(i2.2)') ic
          write(temp_year,'(i4.4)') year
        endif
  
        allocate( temp_pdm(nsrc), temp_reg(reg_num*monperinv), temp_prsb(prsb_fld), temp_bkrn(number_of_offset+bckgrnd_fld) )
  
        ! observations
        open(3, file=trim(dir_kalman)//'obs_data/obs_'//temp_year//temp_mon//'_tmp.txt')

        ! response function
        open(5, file=trim(dir_kalman)//'response_data/reg_resp/reg_'//temp_year//temp_mon//'.txt')
        read(5,*)

        ! pre-subtracted concentrations 
        open(6, file=trim(dir_kalman)//'response_data/prsb_resp/prsb_'//temp_year//temp_mon//'.txt')
        read(6,*)

        !background concentrations
        open(7, file=trim(dir_kalman)//'response_data/bkrn_resp/bkrn_'//temp_year//temp_mon//'.txt') 
        read(7,*)

        open(13, file=trim(dir_kalman)//'obs_data/obs_'//temp_year//temp_mon//'.txt')

        irec = 1
        do i=1,ndat_mon(m)
          read(3,*) temp_ch3,n_site,temp_date10,temp_val,iflg,rsd,bg,pri_local,pri,pri_work

        
      !=========build full pdm
          temp_pdm = 0.0
          read(5,*) temp_ch3,temp_date10,temp_reg
          read(6,*) temp_ch3,temp_date10,temp_prsb
          read(7,*) temp_ch3,temp_date10,temp_bkrn

          temp_pdm(1:reg_num*(monperinv-m+1)) = temp_reg(reg_num*(m-1)+1:reg_num*monperinv)
          do j=1,prsb_fld
            temp_bkrn(1) = temp_bkrn(1) + temp_prsb(j)
          enddo
    
!=======================

          post_local = dot_product(temp_pdm,src)
          post       = temp_bkrn(1) + post_local
          post_work  = temp_val - post
          write(13,330) temp_ch3,n_site,temp_date10,temp_val,iflg,rsd,bg,pri_local,pri,pri_work,post_local,post,post_work

          irec = irec + 1
        enddo
        close(5)
        close(6)
        close(7)
        close(3)
        close(13)  
        deallocate(temp_pdm, temp_reg, temp_prsb, temp_bkrn)
      enddo

 330  format(a3,2x,i3,2x,a10,2x,f15.10, 2x,i2,2x,f15.10,7f15.10)

     end subroutine output_data
