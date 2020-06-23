PROGRAM new_loss

  USE iso_Fortran_env, ONLY: wp => real64
  !USE omp_lib
  use helper_new
  use main_linear_g

  double precision ub(4), lb(4)
  integer*8 opt, q, maxeval
  double precision x(4), minf, stopval
  integer ires
  include 'nlopt.f'

  q = 4

  call nlo_create(opt, NLOPT_LN_COBYLA, q)
  call nlo_get_upper_bounds(ires, opt, ub)
  ub(1) = -0.6
  ub(2) = 0.9
  ub(3) = 0.6
  ub(4) = 0.99

  call nlo_set_upper_bounds(ires, opt, ub)
  call nlo_get_lower_bounds(ires, opt, lb)
  lb(1) = -0.8
  lb(2) = 0.6
  lb(3) = 0.3
  lb(4) = 0.9
  call nlo_set_lower_bounds(ires, opt, lb)

  call nlo_set_min_objective(ires, opt, loss_linear, 0)

  call nlo_set_xtol_rel(ires, opt, 100.0)

  maxeval = 1000

  call nlo_set_maxeval(ires, opt, maxeval)
  call nlo_get_maxeval(maxeval, opt)

  stopval = 10.0
  call nlo_set_stopval(ires, opt, stopval)
  call nlo_get_stopval(stopval, opt)


  x(1) = -0.64208
  x(2) = 0.70235
  x(3) = 0.50994
  x(4) = 0.94568

  call nlo_optimize(ires, opt, x, minf)
  if (ires.lt.0) then
    write(*,*) 'nlopt failed!'
  else
    write(*,*) 'd1 ', x(1)
    write(*,*) 'd2 ', x(2)
    write(*,*) 'pi ', x(3)
    write(*,*) 'beta ', x(4)
    write(*,*) 'min val = ', minf
   endif



call nlo_force_stop(ires, opt)


  call nlo_destroy(opt)




END PROGRAM new_loss
