    module henondim
!---------------------------------------------------------------------------------
!  HENONDIM - Implement Algorithm B & D to determine the dimension of the basin of
!  infinity of the Henon map.
!
!
!  REQUIRED DEPENDENCIES
!    precision - to define KINDs for single and double floating-point.
!    And any other modules that you define.
!
!  REVISION HISTORY
!   10/28/2018 - First Implementation
!  PROGRAMMERS
!   Rahul Rathnakumar, rrathnak@asu.edu, 1215337552
!	Karan Maniar,	ksmaniar@asu.edu, 1213179589
!---------------------------------------------------------------------------------
	use precision
    implicit none
    real(DP), parameter:: LOCKOUT=100, MAXITER=100
    real(DP), parameter:: BOXMIN=-3.0, BOXXMAX=3.0  ! box limits
!   integer, parameter:: NGRID=10  ! number of points on a side
    contains
    subroutine henon_map(basin, m, x) 
!----------------------------------------------------------------------------------
!  Iterate the Henon map for one or more initial conditions up to
!  MAXITER times or until the orbit of one of the points is further than
!  LOCKOUT units away from the origin.  
!
!  Argument declarations
	integer , intent(in) :: m
	real (DP) , intent(in) :: x(m,2)
	integer, intent(inout) :: basin(m) !This vector maps 1 to those points which remain INSIDE the lockout limit and 0 to all other points
	integer :: k, i
	real (DP), parameter :: a=2.12, b=-0.3
	real (DP) :: f(m,2), t(m)
	t=0
	f=x
	!HENON ALGORITHM B STARTS
	!$OMP PARALLEL DO
	do k=1, MAXITER
		where(abs(f(:,1))<LOCKOUT .AND. abs(f(:,2))<LOCKOUT)
			t(:)=f(:,2)
			f(:,2)=f(:,1)
			f(:,1)=a-(f(:,1)**2)+(b*t(:))
		endwhere		
	enddo
	!$OMP END PARALLEL DO
	do i=1,m	
		if(abs(f(i,1))<LOCKOUT .AND. abs(f(i,2))<LOCKOUT) then
			basin(i)=1
		else
			basin(i)=0
		endif
	enddo
	return 
end subroutine henon_map
!-------------------------------------------------------------------------------
    subroutine basin_alg(eps, neps, m)
	real(DP) , intent(inout) :: eps, neps
	integer, intent(inout) :: m
!	LOCAL VARIABLES	
	integer, allocatable :: basin(:), basinp(:), basinm(:)
	real(DP), allocatable:: x(:,:), xeps(:,:)
	real(DP) :: ymin, ymax, xmin, xmax, spx, spy
	integer :: i, j ,k
	allocate(basin(m))
	allocate(basinp(m))
	allocate(basinm(m))
	allocate(x(m,2))
	allocate(xeps(m,2))
	ymin=-3
	ymax=3
	xmin=-3
	xmax=3
	spx=(xmax-xmin)/(m-1)
	spy=(ymax-ymin)/(m-1)
	x=0
	basin=0
	basinp=0
	basinm=0
	do j=1, m
		x(:,1)=xmin + (j-1)*spx !Create X coordinates
		do k=1, m !Create Y coordinates
			x(k,2)=ymin+(k-1)*spy
		enddo	
!	Calculate the basin of infinity for eps=0 first
		call henon_map(basin, m, x) !Returns basin
! 	BASIN_ALG - implements Algorithm D.
		xeps(:,1)=x(:,1) + eps
		xeps(:,2)=x(:,2)
		call henon_map(basinp, m, xeps) !Returns basinp
		xeps(:,1)=x(:,1) - eps
		xeps(:,2)=x(:,2)
		call henon_map(basinm, m, xeps) !Returns basinm
		do i=1,m
			if(basin(i)/=basinm(i) .OR. basin(i)/=basinp(i)) then
				neps=neps+1
			endif
		enddo
	enddo
    return
    end subroutine basin_alg
!-------------------------------------------------------------------------------
	end module henondim
