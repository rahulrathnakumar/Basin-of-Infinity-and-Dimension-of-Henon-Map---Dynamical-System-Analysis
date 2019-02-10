!  DIMENSION OF THE BASIN OF INFINTIY OF THE HENON MAP
!
!  SYNOPSIS
!  This code calculates the basin of infinity of the henon map.
!
!  DESCRIPTION
!	1. Program main contains the modules lsq and henondim which perform the least squares fit and the henon-map calcations respectively.
!	2. Main transfers the array of epsilons and NEpsilon, one by one, into basin_alg
!	3. Basin_alg first calls the henon_map subroutine for eps=0 and then calls henon_map for the 
!	eps+ and eps- values. Comparing these with the base calculation, the number epsilon uncertain points for this column is stored.
!	The next iteration consists of the (j+1)th column.
!	This logic makes the program parallelizable at the highest level, ie, at the main program.
!	5. The number of epsilon uncertain points for each case of epsilon is calculated by the end of the do loop in the main program. 
!	6. The next step is to find the power law fit for epsilon and N(epsilon) using LSQ and output the dimension
!	
!	KNOWN ISSUES
!	Implementation of OpenMP does not produce performance scaling for a 4096x4096 grid. 
!	
! 
!   REVISION HISTORY
!   11/02/2018 - Second Implementation
!   PROGRAMMERS
!   Rahul Rathnakumar, rrathnak@asu.edu, 1215337552
!	Karan Maniar,	ksmaniar@asu.edu, 1213179589
!	
program main
	use lsq
	use precision
	use henondim
	implicit none
	integer ::  m, i, j, ierr
	real(DP) :: ymin, ymax, xmin, xmax, spx, spy,param(2), eps(10), neps(10)
	character(1) :: model
	integer, external :: omp_get_thread_num
	real(DP) :: t1, t2, dt
	print*, "Enter the grid resolution:"
	read*, m
	Neps=0
	ierr=0
	print*, "Do you want to use a linear/power fit?"
	read*, model
	do j=1,10
		eps(j)=2**(-(11.0+j))
	enddo
	call cpu_time(t1)
	!$OMP PARALLEL DO 
	do j=1, 10
		print*, "Thread:" ,omp_get_thread_num()
		call basin_alg(eps(j), neps(j), m)
	enddo
	!$OMP END PARALLEL DO
	call cpu_time(t2)
	dt=t2-t1
	print*, "CPU time:",dt
	call linfit(model, 10, eps, neps, param, ierr)
	if(ierr==0) print*, "Dimension of the basin:" , 1+ param(2)	
end program main