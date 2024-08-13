      subroutine read_ndat
 
      use global_data

      implicit none
      integer :: m,ic
      integer :: temp_mon,temp_year
      integer  :: iyear,imon,idat
 
      open(unit=1,file=trim(dir_obs)//'ndat_month.txt',status='old')
      write(*,*) trim(dir_obs)//'ndat_month.txt'
  
      write(*,*) year,mon
      do while (.not.eof(1))
        read(1,*) iyear,imon,idat
        !write(*,*) iyear,imon,idat
        if((iyear.eq.year).and.(imon.eq.mon)) then
          ndat_mon(1) = idat 
          write(*,*) 1,iyear,imon,idat
          goto 10
        endif
      enddo
 10   continue
  
      if(monperinv-1.le.0) goto 30

      do m=1,monperinv-1 ! read obs in current and previous months  mi20190731

        ic = mon - m
        if(ic.le.0) then
          temp_mon  = ic + 12
          temp_year = year - 1
        else
          temp_mon  = ic
          temp_year = year
        endif

        rewind(1)
        do while (.not.eof(1))
          read(1,*) iyear,imon,idat
          !write(*,*) iyear,imon,idat
          if((iyear.eq.temp_year).and.(imon.eq.temp_mon)) then 
            ndat_mon(1+m) = idat 
            write(*,*) 1+m,iyear,imon,idat
            goto 20
          endif
        enddo
 20   enddo
 30   close(1)
      end
