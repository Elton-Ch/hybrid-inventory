      subroutine build_result_fld

      use global_data

      implicit none

      ! Variables
      integer         :: i, ic, j
      integer         :: irec
      character*2     :: temp_mon
      character*4     :: temp_year

      ! Body of build_fld_from_EOF
      do i=monperinv,1,-1
        ic = mon - (i-1)
        if(ic.le.0)then
          write(temp_mon,'(i2.2)') ic + 12
          write(temp_year,'(i4.4)') year - 1
        else
          write(temp_mon,'(i2.2)') ic
          write(temp_year,'(i4.4)') year
        endif
        open(1, file=trim(dir_kalman)//'state_vector/src_'//temp_year//temp_mon//'.dat') !", status='replace')" was deleted. by TSMI 2013feb13
        write(1,*) reg_num
        do j=1,reg_num
          write(1,*) src( (monperinv-i)*reg_num + j  ), sqrt(csrc( (monperinv-i)*reg_num + j, (monperinv-i)*reg_num + j ))
        enddo
        close(1)
      enddo

      write(temp_mon,'(i2.2)') mon
      write(temp_year,'(i4.4)') year

      open(1, file=trim(dir_kalman)//'B_matrices/csrc_'//temp_year//temp_mon//'.bin', form='unformatted', access='direct', convert='big_endian', recl=2*monperinv*reg_num)
      irec = 1
      do j=1,monperinv*reg_num
        write(1,rec=irec) (csrc(j,i),i=1,monperinv*reg_num)
        irec = irec + 1
      enddo
      close(1)

      end subroutine build_result_fld
