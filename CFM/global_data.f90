 module global_data

     real*8,       dimension( : , : ), allocatable:: csrc             ! [ nsrc , nsrc ]
     real*8,       dimension( : , : ), allocatable:: cdat             ! [ ndat , ndat ]
     real*8,       dimension( : , : ), allocatable:: pdm              ! [ ndat , nsrc ]
     real*8,       dimension(   :   ), allocatable:: src              ! [    nsrc     ]
     real*8,       dimension(   :   ), allocatable:: dat              ! [    ndat     ]
     real*8,       dimension(   :   ), allocatable:: obs_val          ! [    ndat     ]
     real*8,       dimension( : , : ), allocatable:: dev_dat          ! [ ndat ,   2  ]
     real*8,       dimension( : , : ), allocatable:: chi_sq           ! [ ndat ,   2  ]
     real*8,       dimension(   :   ), allocatable:: chi_sq_st        ! [ ndat ,   2  ]
     real*8,       dimension(   :   ), allocatable:: obs_lon          ! [    ndat     ]
     real*8,       dimension(   :   ), allocatable:: obs_lat          ! [    ndat     ]
     real*8,       dimension(   :   ), allocatable:: obs_alt          ! [    ndat     ]
     integer,      dimension(   :   ), allocatable:: obs_year         ! [    ndat     ]
     integer,      dimension(   :   ), allocatable:: obs_month        ! [    ndat     ]
     integer,      dimension(   :   ), allocatable:: obs_day          ! [    ndat     ]
     integer,      dimension(   :   ), allocatable:: obs_hour         ! [    ndat     ]
     integer,      dimension(   :   ), allocatable:: obs_min          ! [    ndat     ]
     character*3,  dimension(   :   ), allocatable:: st_name          ! [    ndat     ]   !site name
     integer,      dimension(   :   ), allocatable:: obs_site         ! [    ndat     ]   !site number
     character*10, dimension(   :   ), allocatable:: obs_date         ! [    ndat     ]   !site number

     integer            :: nsrc, ndat
     integer            :: monperinv, mon, year, lag_window
     integer            :: st_year, st_mon
     integer            :: reg_num, prsb_fld, number_of_offset, bckgrnd_fld
     integer            :: land_reg, ocn_reg
     integer, parameter :: lon = 360, lat = 180
     character*2        :: mon_of_obs
     character*4        :: year_of_obs
     character*200      :: dir_obs
     character*200      :: list_obs
     character*200      :: dir_kalman
     character*200      :: dir_sf_init
     character*200      :: file_config

     real*8             :: use_prsb
     integer            :: n_select

     integer            :: count_num                 ! mi201907
     integer,dimension(:),allocatable :: ndat_mon    ! [log_window]   mi20190801

 end module global_data
