   module lsq
!  LSQ - Interface module for LAPACK QR routines, and a special-purpose
!  subroutine for simple Y versus X linear regressions.
!
!  SYNOPSIS
!    use lsq
!
!  DESCRIPTION
!    Write a one-sentence description of LINFIT.
!
!  PUBLIC ROUTINES DEFINED
!    linfit - Y versus X linear regression.
!
!  DEPENDENCIES
!    precision - defines KINDs for single- and double-precision floating point.
!
!  REVISION HISTORY
!    10/28/18 - Second implementation.
!
!  PROGRAMMER
!    Your name, email@asu.edu
!
   use precision
   implicit none
!
!  Generic interface (for single- or double-precision) for LAPACK
!  QR factorizations SGELS and DGELS.
!
   interface xgels
      subroutine sgels(tr,n,k,nrhs,x,ldx,y,ldy,work,lwork,info)
         import
         character, intent(in):: tr
         integer, intent(in):: n, k, ldx, ldy, nrhs
         integer, intent(inout):: lwork
         real(SP), intent(inout):: x(ldx,k), y(ldy,nrhs), work(lwork)
         integer, intent(out):: info
      end subroutine sgels
!
      subroutine dgels(tr, n, k, nrhs, x, ldx, y, ldy, work, lwork, info)
         import
         character, intent(in):: tr
         integer, intent(in):: n, k, ldx, ldy, nrhs
         integer, intent(inout):: lwork
         real(DP), intent(inout):: x(ldx,k), y(ldy,nrhs), work(lwork)
         integer, intent(out):: info
      end subroutine dgels
   end interface
!
   contains
!----------
   subroutine linfit(model, n, x, y, param, ierr)
!  LINFIT - briefly describe the input and output arguments here.
!
   character(*), intent(in):: model
   integer, intent(in):: n
   real(DP), intent(in):: x(n)
   real(DP) ::y(n,1)
   real(DP), intent(out):: param(2)
   integer, intent(out):: ierr

   real(DP) :: xmatrix(n,2)
   real (DP) :: work(n)
   integer :: lwork, info
   real (DP) :: lx(n,2) , ly(n,1)   
   
   !integer :: lwork,k,ldx,ldy
  
   ! integer :: info
   ! real(DP) :: lx(n) , ly(n)
   ! integer :: nrhs
   ! real(DP) :: Ymatrix(n,1)
   ! nrhs = 1
   ! k=2
   lwork= max( 1, 2*n + max(2*n, 1) )
  ! real(DP) ::x1(ldx,k), Y1(ldy,nrhs), work(1,lwork)
!  your code here
select case (model(1:1))
case ("L", "l")
xmatrix(:,1) = 1
xmatrix(:,2) = x
call xgels("N", n, 2, 1, xmatrix,n, y, n, work, lwork, info)
!straight-line model
param(1) = y(1,1)
param(2) = y(2,1)

case ("P", "p")
ly=log(y)
lx(:,1) =1
lx(:,2)=log(x)
!power-law model
call xgels("N", n, 2, 1, lx, n, ly, n, work, lwork, info)
param(1) = exp(y(1,1))
param(2) = ly(2,1)
case default
ierr=-1
end select
return
!unknown model; set ierr appropriately

!
   end subroutine linfit
!-----------------------
   end module lsq
