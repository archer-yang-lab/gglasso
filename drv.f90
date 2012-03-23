!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine logdrv(nobs,nvars,x,y,r,vl)
      IMPLICIT NONE
      integer :: nobs
      integer :: nvars
      double precision :: y(nobs)
      double precision :: r(nobs)
      double precision :: x(nobs,nvars)
      double precision :: vl(nvars)                                                                                                                                                                                                                 
	  vl = - matmul(y/(1.0D0+exp(r)), x) / nobs
      END                                                 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine hubdrv(delta,nobs,nvars,x,y,r,vl)
      IMPLICIT NONE
      integer :: nobs
      integer :: nvars
      integer :: i
      double precision :: delta
      double precision :: dl(nobs)
      double precision :: y(nobs)
      double precision :: r(nobs)
      double precision :: x(nobs,nvars)
      double precision :: vl(nvars)                                                                                                                                                                                                                 
	  vl=0.0
	   DO i = 1, nobs
	       IF (r(i) > 1.0D0) THEN
	          dl (i) = 0.0D0
	       ELSEIF (r(i) <= (1-delta)) THEN
	          dl (i) = - 1.0D0
	       ELSE
	          dl (i) = (r(i)-1.0D0) / delta
	       ENDIF
       ENDDO
	       vl = matmul(dl*y, x) / nobs
       END                                                 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sqdrv(nobs,nvars,x,y,r,vl)
      IMPLICIT NONE
      integer :: nobs
      integer :: nvars
      integer :: i
      double precision :: dl(nobs)
      double precision :: y(nobs)
      double precision :: r(nobs)
      double precision :: x(nobs,nvars)
      double precision :: vl(nvars)                                                                                                                                                                                                                 
	  vl=0.0
	  DO i = 1, nobs
	       dl (i) = dim (1.0D0, r(i))
      ENDDO
	  vl = - 2.0D0 * matmul(dl*y, x) / nobs
      END                                                 







