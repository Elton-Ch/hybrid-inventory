 subroutine kalman_1 ! if nsrc > ndat

 use global_data

 implicit none

! Variables

 integer i, k(ndat)

 real*8, dimension( : , : ), allocatable:: KH, icsrc, HBHO, KALM, BH
 real*8, dimension(   :   ), allocatable:: isrc

!    KH( : , : )                         ! [ nsrc , nsrc ]
! icsrc( : , : )                         ! [ nsrc , nsrc ]
!  HBHO( : , : )                         ! [ ndat , ndat ]
!  KALM( : , : )                         ! [ nsrc , ndat ]
!    BH( : , : )                         ! [ nsrc , ndat ]
!  isrc(   :   )                         ! [    nsrc     ]

! Body of kalman_1

!   K = B*H'*(H*B*H' + O)^-1
!   x = x0 + K*(y - H*x0)
!   B = (I - K*H)*B
!----------------------------------------------------                !
 allocate( HBHO(ndat,ndat) , BH(nsrc,ndat) )   !  

   BH = matmul( csrc , transpose(pdm) )        !   calculate (H*B*H' + O)^-1
   HBHO = matmul( pdm , BH ) + cdat            !
!----------------------------------------------------
 call mtxchi( HBHO , ndat )                    !
!----------------------------------------------------                 !
 allocate( KALM(nsrc,ndat) )                   !
                                               !   calculate K = B*H'*(H*B*H' + O)^-1
   KALM = matmul( BH , HBHO )                  !
 deallocate( HBHO , BH )                       !
!----------------------------------------------------
   write(*,*) 'calculate new src'              !
 allocate( isrc(nsrc), chi_sq_st(nsrc) )                        !
   chi_sq_st = 0.0
   isrc = src                                  !
                                               !   calculate x = x0 + K*(y - H*x0)
   src = isrc + matmul( KALM , (dat - matmul( pdm , isrc )) )
                                               !
   chi_sq_st(:) = (src(:) - isrc(:))**2									   
 deallocate( isrc )                            !
!----------------------------------------------------
   write(*,*) 'calculate new csrc'             !
 allocate( KH(nsrc,nsrc) )                     !
                                               !
   KH = matmul( KALM , pdm )                   !
 deallocate( KALM )                            !   calculate (I - K*H)
   do i=1,nsrc                                 !
       KH(i,i) = 1.0 - KH(i,i)                 !
   enddo                                       !
!----------------------------------------------------
 allocate( icsrc(nsrc,nsrc) )                  !
                                               !
   icsrc = matmul( KH , csrc )                 !
 deallocate( KH )                              !   calculate B = (I - K*H)*B
                                               !
   csrc = icsrc                                !
 deallocate( icsrc )                           !

 do i=1,nsrc
     chi_sq_st(i) = chi_sq_st(i)/csrc(i,i)
 enddo
 
 end subroutine kalman_1
!==================================================================================
 subroutine kalman_2 ! if nsrc < ndat

 use global_data

 implicit none

! Variables

 integer i

 real*8, dimension( : , : ), allocatable:: icsrc, HOHB, KALM, HO, icdat
 real*8, dimension(   :   ), allocatable:: isrc

! icsrc( : , : )                         ! [ nsrc , nsrc ]
!  HOHB( : , : )                         ! [ nsrc , nsrc ]
!  KALM( : , : )                         ! [ nsrc , ndat ]
!    HO( : , : )                         ! [ ndat , nsrc ]
! icdat( : , : )                         ! [ ndat , ndat ]
!  isrc(   :   )                         ! [    nsrc     ]

! Body of kalman_2

 allocate( icdat(ndat,ndat) )
   icdat = 0.0
   do i=1,ndat
       icdat(i,i) = 1./cdat(i,i) ! R^-1
   enddo

 allocate( icsrc(nsrc,nsrc) )
   icsrc = csrc
 call mtxchi( icsrc , nsrc )     ! B^-1

 allocate( HO(nsrc,ndat) , HOHB(nsrc,nsrc) )
   HO = matmul( transpose(pdm) , icdat )
 deallocate( icdat )
   HOHB = matmul( HO , pdm ) + icsrc
 deallocate( icsrc )

 call mtxchi( HOHB , nsrc )

 allocate( KALM(nsrc,ndat) )
   KALM = matmul( HOHB , HO )
 deallocate( HO )
 
 allocate( isrc(nsrc), chi_sq_st(nsrc) )                        !
 chi_sq_st = 0.0
   isrc = src 
 
   write(*,*) 'calculate new src'
   src = isrc + matmul( KALM , (dat - matmul( pdm , isrc )) )
   
 do i=1,nsrc
     chi_sq_st(i) = ((src(i) - isrc(i))**2)/csrc(i,i)
 enddo   
 
   write(*,*) 'calculate new csrc' 
   csrc = HOHB

 deallocate( HOHB , isrc )

 end subroutine kalman_2
