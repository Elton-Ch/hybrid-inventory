      subroutine build_global_data

      use global_data

      implicit none

     ! Variables
      integer                                :: i, ic, j, temp_ndat, jc
      integer                                :: irec, count, len, temp_int
      real                                   :: temp_real

      real*8                                 :: work, temp_val, temp_cval, temp_lon, temp_lat,temp_alt
      real*4,dimension( : ), allocatable     :: temp_pdm, temp_reg, temp_prsb, temp_bkrn

      character*13                           :: temp_ch3      ! for site with 3 characters mi20191212
      character*10                           :: temp_date10   ! for date YYYYMMDDHH mi20191212
      character*16                           :: temp_ch16     ! for time (year,month,day,hour,minumte)   
      character*2                            :: temp_mon
      character*4                            :: temp_year

      character*100                          :: filename

      integer                                :: nf,nf_site     ! to exclude difficult sites mi20140117
      character*13,dimension(:), allocatable :: site_ex        ! to read site_id with 13 characters mi20191002
      integer                                :: nf_rsd         ! total number of observation sites
      real,dimension(:), allocatable         :: rsd            ! to read site_id with 13 characters mi20191002

      integer                                :: iyear,imonth,iday,ihour,imin          
      integer,dimension(:), allocatable      :: tm_ex1,tm_ex2  ! to exclude difficult sites by time mi20140124

      integer                                :: m               ! mi20190801     
      integer                                :: n_site

      integer                                :: idummy,inv
      integer                                :: iflg
      real                                   :: local

      real*8,dimension(:,:), allocatable       :: csrc_temp

      allocate(site_ex(ndat))

    
      ! obs uncertainties & selection  ---------------  
      open(unit=1,file=trim(list_obs)//year_of_obs//mon_of_obs//'.txt')
      write(*,*) trim(list_obs)//year_of_obs//mon_of_obs//'.txt'
      read(1,*)
      read(1,*) nf_rsd   ! total number of sites
      write(*,*) 'nf_rsd = ', nf_rsd
      allocate(rsd(nf_rsd))
      nf_site = 0   ! number of sites to be excluded
      do nf=1,nf_rsd
        read(1,*) idummy,temp_ch3,inv,rsd(nf)
        write(*,'(i3,2x,a3,2x,i1,2x,f10.3)') idummy,temp_ch3,inv,rsd(nf)
        if(inv.ne.1) then
          nf_site = nf_site+1
          site_ex(nf_site) = temp_ch3
        endif
      enddo
      close(1)
      

      ! build state vector for prior (current month)
      allocate( src(nsrc), csrc(nsrc,nsrc), csrc_temp(nsrc,nsrc))
      src = 0.0; csrc = 0.0; csrc_temp=0.9
      open(1, file=trim(dir_sf_init)//'src_'//year_of_obs//mon_of_obs//'.txt')
      read(1,*)
      
      do j=1,reg_num ! add new part of state vector
        read(1,*) temp_val, temp_cval
        src( (monperinv-1)*reg_num + j) = temp_val
        csrc((monperinv-1)*reg_num + j,(monperinv-1)*reg_num + j) = (1.0*temp_cval)**2   ! prior uncertainty
      enddo
      close(1)

      ! build remaining part of state vector (previous months)
      
      if(monperinv.ge.2) then

        do i=monperinv,2,-1 
          ic = mon - (i-1)
          if(ic.le.0) then
            write(temp_mon,'(i2.2)') ic + 12
            write(temp_year,'(i4.4)') year - 1
          else
            write(temp_mon,'(i2.2)') ic
            write(temp_year,'(i4.4)') year
          endif
        
          open(1, file=trim(dir_sf_init)//'src_'//temp_year//temp_mon//'.dat') 
          read(1,*)
          do j=1,reg_num ! add new part of state vector
            read(1,*) temp_val, temp_cval
            src( (monperinv-i)*reg_num + j) = temp_val
            csrc((monperinv-i)*reg_num + j,(monperinv-i)*reg_num + j) = (1.0*temp_cval)**2   !  prior uncertainty
          enddo
          close(1)
          
        enddo
      endif


 !====== define real ndat ========== (remove the obs data if model-obs mismatch is too large)

      count = 0
      do m=monperinv,1,-1 ! read obs in current and previous months  mi20190731

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
        open(3, file=trim(dir_obs)//'obs_'//temp_year//temp_mon//'.txt')
        read(3,*)

        ! response function
        open(5, file=trim(dir_kalman)//'response_data/reg_resp/reg_'//temp_year//temp_mon//'.txt')
        read(5,*)

        ! pre-subtracted concentrations 
        open(6, file=trim(dir_kalman)//'response_data/prsb_resp/prsb_'//temp_year//temp_mon//'.txt')
        read(6,*)

        !background concentrations
        open(7, file=trim(dir_kalman)//'response_data/bkrn_resp/bkrn_'//temp_year//temp_mon//'.txt')
        read(7,*)

        irec = 1
        do i=1,ndat_mon(m)
          !write(*,*) trim(dir_obs)//'points_'//temp_year//temp_mon//'.txt'
          !read(3,*) iyear,imonth,iday,ihour,imin,temp_lat,temp_lon,temp_alt,temp_val,temp_ch13,n_site  ! to read site_id with 9 characters mi20160528
          read(3,*) n_site,temp_ch3,temp_date10,temp_val

!=========build full pdm
          temp_pdm  = 0.0
          read(5,*) temp_ch3,temp_date10,temp_reg
          read(6,*) temp_ch3,temp_date10,temp_prsb
          read(7,*) temp_ch3,temp_date10,temp_bkrn

          temp_pdm(1:reg_num*(monperinv-m+1)) = temp_reg(reg_num*(m-1)+1:reg_num*monperinv) 

          ! check matrix
          !do j=1,reg_num*(monperinv)
          !  write(*,*) j,temp_pdm(j),temp_reg(j)
          !enddo
          !write(*,*) i,temp_year,temp_mon,temp_bkrn(1), prsb_fld,temp_prsb(1),temp_prsb(2)
                 
          do j=1,prsb_fld
            temp_bkrn(1) = temp_bkrn(1) + temp_prsb(j)
          enddo
          
!=======================
          work = temp_val - temp_bkrn(1) - dot_product(temp_pdm,src)  ! obs-model mismatch
          !write(*,*) iyear,imonth,iday,temp_ch13,n_site,temp_val,temp_bkrn(1), dot_product(temp_pdm,src),work
          
          if(nf_site.gt.0) then
            do nf=1,nf_site                          
              if(temp_ch3.eq.site_ex(nf)) then    
                !write(*,'(a20,a13)') 'site to be excluded ',temp_ch13
                goto 888 ! sites to be excluded mi20140117
              endif
            enddo
          endif
      
       !
          if(n_select.gt.0) then
            if(abs(work).lt.rsd(n_site)*n_select) then    ! rejection criterion : unc*3 mi20131021
              count = count + 1
            endif
          else
            count = count + 1     ! no selection by model-data mismatch
          endif

 888      continue
          irec = irec + 1
          !pause
        enddo
        close(5)
        close(6)
        close(7)
        close(3)
        deallocate(temp_pdm, temp_reg, temp_prsb, temp_bkrn)
        
      enddo
    
    !----------------
      temp_ndat = ndat
      ndat      = count ! real number of observations after filtering

      write(*,*) nsrc, ndat,temp_ndat
      !pause

!======================================================
      !allocate( pdm(ndat,nsrc), dat(ndat), obs_val(ndat), cdat(ndat,ndat), dev_dat(ndat,2), obs_lon(ndat), obs_lat(ndat), obs_alt(ndat), st_name(ndat), obs_site(ndat))
      !allocate( obs_year(ndat), obs_month(ndat),obs_day(ndat),obs_hour(ndat),obs_min(ndat))

      allocate( pdm(ndat,nsrc), dat(ndat), obs_val(ndat), cdat(ndat,ndat), dev_dat(ndat,2), st_name(ndat))
      allocate( obs_date(ndat))
    
      pdm       = 0.0
      dat       = 0.0
      cdat      = 0.0
      obs_val   = 0.0
      dev_dat   = 0.0
 
      count = 1
      iflg  = 0     ! 1 data to be used, 0 : not to be used
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
        open(3, file=trim(dir_obs)//'obs_'//temp_year//temp_mon//'.txt')
        read(3,*)

        ! response function
        open(5, file=trim(dir_kalman)//'response_data/reg_resp/reg_'//temp_year//temp_mon//'.txt') 
        read(5,*)

        ! pre-subtracted concentrations 
        open(6, file=trim(dir_kalman)//'response_data/prsb_resp/prsb_'//temp_year//temp_mon//'.txt') 
        read(6,*)

        !background concentrations
        open(7, file=trim(dir_kalman)//'response_data/bkrn_resp/bkrn_'//temp_year//temp_mon//'.txt')
        read(7,*)

        open(13, file=trim(dir_kalman)//'obs_data/obs_'//temp_year//temp_mon//'_tmp.txt') 

        irec = 1
        do i=1,ndat_mon(m)
          read(3,*)  n_site,temp_ch3,temp_date10,temp_val
          write(*,*) n_site,temp_ch3,temp_date10,temp_val
        
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
          local = dot_product(temp_pdm,src)
          work  = temp_val - (temp_bkrn(1) + local)


          write(*,330) temp_ch3,n_site,temp_date10,temp_val,iflg,rsd(n_site),temp_bkrn(1),local,temp_bkrn(1)+local,work 
	    

          if(nf_site.gt.0) then
            do nf=1,nf_site                          
              if(temp_ch3.eq.site_ex(nf)) then 
                  iflg=-1
                  write(13,330) temp_ch3,temp_date10,temp_val,iflg,rsd(n_site),temp_bkrn(1),local,temp_bkrn(1)+local,work
                  goto 889 ! sites to be excluded mi20140117
              endif
            enddo
          endif

          if(n_select.eq.0) then
              dat(count)        = temp_val - temp_bkrn(1)
              obs_val(count)    = temp_val
              cdat(count,count) = (1.0*rsd(n_site))**2
              st_name(count)    = temp_ch3
              obs_date(count)   = temp_date10
              pdm(count,:)      = temp_pdm(:)

              iflg              = 1
              write(13,330) temp_ch3,n_site,temp_date10,temp_val,iflg,rsd(n_site),temp_bkrn(1),local,temp_bkrn(1)+local,work

              count             = count + 1
              goto 889 
          endif

          if((n_select.gt.0).and.(abs(work).lt.rsd(n_site)*n_select) ) then   ! rejection criterion 
              !write(*,*) count,work,temp_val,temp_bkrn(1)
              dat(count)        = temp_val - temp_bkrn(1)
              obs_val(count)    = temp_val
              cdat(count,count) = (1.0*rsd(n_site))**2
              st_name(count)    = temp_ch3
              obs_date(count)   = temp_date10
              pdm(count,:)      = temp_pdm(:)

              iflg              = 1
              write(13,330) temp_ch3,n_site,temp_date10,temp_val,iflg,rsd(n_site),temp_bkrn(1),local,temp_bkrn(1)+local,work

              count             = count + 1
              goto 889
          else
            write(13,330) temp_ch3,n_site,temp_date10,temp_val,iflg,rsd(n_site),temp_bkrn(1),local,temp_bkrn(1)+local,work
          endif

 889      continue
          irec = irec + 1
          iflg = 0

        enddo
        close(5)
        close(6)
        close(7)
        close(3)
        close(13)  
        deallocate(temp_pdm, temp_reg, temp_prsb, temp_bkrn)
      enddo

      write(*,*) nsrc, ndat,count-1
 330  format(a3,2x,i3,2x,a10,2x,f15.10, 2x,i2,5(2x,f15.8))

      end subroutine build_global_data
