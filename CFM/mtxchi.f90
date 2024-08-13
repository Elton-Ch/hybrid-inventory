subroutine mtxchi(A,n) ! from "Data Analysis Statistical and Computational
                       ! Methods for Scientists and Engineers", Siegmund Brandt
implicit none

integer i, j, k, n
real*8, allocatable:: U( : , : )
real*8 A( n , n ), s

! Step 1: Cholesky decomposition
allocate( U(n,n) )

U = 0.0
do i=1,n
    s = 0.0
    do j=i,n
        if(i.gt.1) then
            s = 0.0
            do k=1,i-1
                s = s + U(k,i)*U(k,j)
            enddo
        endif
        U(i,j) = A(i,j) - s
        if(i.eq.j) then
            U(i,j) = SQRT( ABS( U(i,j) ) )
        else
            U(i,j) = U(i,j)/U(i,i)
        endif
    enddo
enddo

do i=1,n
! Step 2: Forward Substitution
    do j=i,n
        if(j.eq.i) then
            A(n,j) = (1.0)/U(j,j)
        else
            A(n,j) = 0.0
            do k=i,j-1
                A(n,j) = A(n,j) - U(k,j)*A(n,k)
            enddo
            A(n,j) = A(n,j)/U(j,j)
        endif
    enddo
! Step 3: Back Substitution
    do j=n,i,-1
        if(j.eq.n) then
            A(i,j) = A(n,j)/U(j,j)
        else
            A(i,j) = A(n,j)
            do k=n,j+1,-1
                A(i,j) = A(i,j) - U(j,k)*A(i,k)
            enddo
                A(i,j) = A(i,j)/U(j,j)
        endif
    enddo
enddo
deallocate( U )
! Fill lower triangle symmetrically
if(n.gt.1) then
    do i=1,n
        do j=1,i-1
            A(i,j) = A(j,i)
        enddo
    enddo
endif

end subroutine mtxchi
