MODULE helper_new


USE iso_Fortran_env, ONLY: wp => real64
IMPLICIT NONE


REAL(wp), parameter :: sigma = 2.0, sigmag = 2.0, gamma = 2.0, chi = 1.0
REAL(wp), parameter :: r = 1.01, theta = 0.0385, shi = 1.0, n = 10.0, m = 10.0
REAL(wp), parameter :: tolq = 1.0D-5, tolv = 1.0D-5
INTEGER, parameter :: maxiter = 1000, maximum = 100000, allow = 90001, zsiz = 201, bsiz = 100, tsiz = 100, k = 1,&
period = 40



CONTAINS

  SUBROUTINE penalty(shock,d0, d1, d2, pen)
  REAL(wp), INTENT(IN) :: shock, d0, d1, d2
  REAL(wp), INTENT(OUT) :: pen
  !pen = shock - d0*shock**d1
  pen = shock-max(0.0,d0+d1*shock+d2*shock**2.0)
  END SUBROUTINE penalty



  SUBROUTINE utility(c,l,g,x,pi,util)
  REAL(wp), INTENT(IN) :: c, l, g, x, pi
  REAL(wp), INTENT(OUT) :: util

  IF (c>0.0 .and. l>0.0 .and. g>0.0) THEN
    IF (k == 2) THEN
      util = (1.0/(1.0-sigma))*(c+shi*x-chi*(l**(1.0+gamma))/(1.0+gamma))**(1.0-sigma)
    ELSE
      util = ((1.0-pi)/(1.0-sigma))*((c-chi*(l**(1.0+gamma))/(1.0+gamma))**(1.0-sigma))+ &
      (pi/(1.0-sigmag))*g**(1.0-sigmag) + shi*x
    END IF
  ELSE
  util = -999.0
  END IF

  END SUBROUTINE utility

  SUBROUTINE mean(xx,x_bar)
    REAL(wp), dimension(allow),INTENT(IN) :: xx
    REAL(wp), INTENT(OUT) :: x_bar

    x_bar = SUM(xx, xx .NE. 1000.0)/(SIZE(xx)-COUNT(xx==1000.0))

  END SUBROUTINE mean

  SUBROUTINE mean1(xx,x_bar)
    REAL(wp), dimension(period),INTENT(IN) :: xx
    REAL(wp), INTENT(OUT) :: x_bar

    x_bar = SUM(xx)/SIZE(xx)

  END SUBROUTINE mean1


  SUBROUTINE variance(xx,x_var)
    REAL(wp), dimension(period),INTENT(IN) :: xx
    REAL(wp), INTENT(OUT) :: x_var
    REAL(wp) :: x_bar

    CALL mean1(xx,x_bar)
    x_var = SUM((xx-x_bar)**2.0, xx .NE. 1000.0)/(SIZE(xx)-COUNT(xx==1000.0))

  END SUBROUTINE variance

  SUBROUTINE cov(xx,yy,xy_cov)
    REAL(wp), dimension(period),INTENT(IN) :: xx, yy
    REAL(wp), INTENT(OUT) :: xy_cov
    REAL(wp) :: x_bar, y_bar

    CALL mean1(xx,x_bar)
    CALL mean1(yy,y_bar)
    xy_cov = SUM((xx-x_bar)*(yy-y_bar), xx .NE. 1000.0)/(SIZE(xx)-COUNT(xx==1000.0))

  END SUBROUTINE cov

  SUBROUTINE corr(xx,yy,cor)
    REAL(wp), dimension(period),INTENT(IN) :: xx, yy
    REAL(wp), INTENT(OUT) :: cor
    REAL(wp) :: xy_cov,x_var,y_var

    CALL variance(xx,x_var)
    CALL variance(yy,y_var)
    CALL cov(xx,yy,xy_cov)

    cor = xy_cov/(SQRT(x_var)*SQRT(y_var))

  END SUBROUTINE corr


END MODULE helper_new
