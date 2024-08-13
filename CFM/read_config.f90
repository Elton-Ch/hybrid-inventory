       subroutine read_config

       use global_data

       implicit none

       integer i
       character*2  temp_ch2
       character*13 temp_ch13
       character*14 temp_ch14
       character*17 temp_ch17
       character*18 temp_ch18
       character*21 temp_ch21

       open(1, file=trim(file_config),status='old')

       read(1,*) temp_ch2
       do while(temp_ch2.ne.'#0')
         read(1,*) temp_ch2
       enddo

!--------------------------
      read(1,*)
      read(1,*) temp_ch18, reg_num
      read(1,*)
      read(1,*) temp_ch13, land_reg
      read(1,*)
      read(1,*) temp_ch14, ocn_reg
      read(1,*)
      read(1,*) temp_ch21, prsb_fld
      read(1,*)
      read(1,*) temp_ch18, bckgrnd_fld
      read(1,*)
      read(1,*) temp_ch17, number_of_offset	  
!--------------------------

      close(1)

      end subroutine read_config
