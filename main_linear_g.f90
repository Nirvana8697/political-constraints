MODULE main_linear_g

  USE iso_Fortran_env, ONLY: wp => real64
  !USE omp_lib
  use helper_new



  !IMPLICIT NONE

CONTAINS

  SUBROUTINE loss_linear(loss_value, q, x, grad, need_gradient, f_data)


  INTEGER :: q, need_gradient
  REAL(wp), dimension(q), INTENT(IN) :: x
  REAL(wp), dimension(q) :: grad
  REAL(wp), INTENT(OUT) :: loss_value
  REAL(wp) :: m_spds, sd_spds, rho_tb_y, &
  avg_g_to_y, avg_b_to_y, d0, d1, d2, pi, beta, sd_c_y, data_sd_c_y
  REAL(wp) :: c, l, g, errq, errv, pen, util, Po, unif, data_m_spds,&
   data_sd_spds, data_avg_g_to_y, data_avg_b_to_y, data_rho_spds_y, data_sd_c, data_sd_y
  REAL(wp) :: rho_g_y, rho_c_y, rho_spds_y, rho_tax_y, sd_tax, sd_GDP, &
  sd_TB, sd_c, var_c, var_spds, var_tax, var_GDP, var_TB, m_tax, m_GDP, m_TB, m_c
  INTEGER :: zix, tix, bix, bpix, iter
  REAL(wp), dimension(tsiz) :: V_tax, t
  REAL(wp), dimension(bsiz) :: b, Vpos
  REAL(wp), ALLOCATABLE, dimension(:, :, :, :) :: B_rep
  REAL(wp), dimension(zsiz, tsiz) :: B_def
  REAL(wp), dimension(zsiz, zsiz) :: Zprob
  REAL(wp), dimension(zsiz) :: Z, sdef, gov_def
  REAL(wp), dimension(4) :: loss_array


  ! Policies
  INTEGER, dimension(zsiz) :: tax_def
  INTEGER, ALLOCATABLE, dimension(:, :, :) :: tax_rep
  REAL(wp), ALLOCATABLE, dimension(:, :, :) :: gov_rep
  !INTEGER, dimension(zsiz, bsiz, n) :: d
  INTEGER, dimension(zsiz, bsiz) :: b_pol, d
  REAL(wp), dimension(zsiz, bsiz) :: q1, Vc
  REAL(wp), dimension(zsiz, bsiz) :: Brep, tax, Bpol, q0, Tr
  REAL(wp), dimension(zsiz, bsiz) ::  Jc1, taxrep, srep
  REAL(wp), dimension(zsiz, bsiz) ::Prop_rep, Prop_def, Prop
  REAL(wp), dimension(zsiz, bsiz) :: J0, J1, gov, labor, GDP, cons
  REAL(wp), dimension(zsiz) :: Jd0, Jd1, Vd
  REAL(wp), dimension(zsiz) :: taxdef, Bdef
  REAL(wp), ALLOCATABLE, dimension(:) :: c_simul, g_simul, GDP_simul, spread_simul, &
  TB_simul, labor_simul, q_simul, g_to_y, b_to_y, tax_simul, Tr_simul
  INTEGER, ALLOCATABLE, dimension(:) :: b_simul_loc, d_simul
  INTEGER, ALLOCATABLE, dimension(:) :: z_simul_loc
  REAL(wp), ALLOCATABLE, dimension(:) :: spds_rel, g_y_rel, b_y_rel, sd_spds_rel, sd_c_rel, &
  sd_GDP_rel, sd_TB_rel, sd_tax_rel, rho_TB_y_rel, rho_g_y_rel, rho_c_y_rel, rho_spds_y_rel, rho_tax_y_rel
  INTEGER, ALLOCATABLE, dimension(:) :: d_simulate, b_simulate_loc
  REAL(wp), ALLOCATABLE, dimension(:) :: g_simulate, c_simulate, GDP_simulate, q_simulate, &
  TB_simulate, spread_simulate, tax_simulate, b_to_y_simulate, g_to_y_simulate, Tr_simulate
  ALLOCATE(B_rep(zsiz, bsiz, bsiz, tsiz), tax_rep(zsiz, bsiz, bsiz), gov_rep(zsiz, bsiz, bsiz))
  ALLOCATE(c_simul(maximum), g_simul(maximum), GDP_simul(maximum), spread_simul(maximum), &
  TB_simul(maximum), b_simul_loc(maximum), labor_simul(maximum), &
  q_simul(maximum), d_simul(maximum), b_to_y(maximum), g_to_y(maximum), tax_simul(maximum), Tr_simul(maximum))
  ALLOCATE(z_simul_loc(maximum))
  ALLOCATE(spds_rel(allow), g_y_rel(allow), b_y_rel(allow), sd_spds_rel(allow), sd_c_rel(allow),&
  sd_GDP_rel(allow), sd_TB_rel(allow), rho_TB_y_rel(allow), rho_g_y_rel(allow), rho_c_y_rel(allow), &
  rho_spds_y_rel(allow), rho_tax_y_rel(allow), sd_tax_rel(allow))
  ALLOCATE(d_simulate(allow), b_simulate_loc(allow), g_simulate(allow), c_simulate(allow), GDP_simulate(allow), &
  q_simulate(allow), TB_simulate(allow), spread_simulate(allow), tax_simulate(allow), b_to_y_simulate(allow), &
   g_to_y_simulate(allow), Tr_simulate(allow))


  OPEN(unit=100, file="Z.txt", status="old", action="read", position="rewind")
  READ (100, *) Z
  CLOSE(100)

  OPEN(unit=100, file="Zprob.txt", status="old", action="read", position="rewind")
  READ (100, *) Zprob
  CLOSE(100)

  OPEN(unit=100, file="b.txt", status="old", action="read", position="rewind")
  READ (100, *) b
  CLOSE(100)

  OPEN(unit=100, file="t.txt", status="old", action="read", position="rewind")
  READ (100, *) t
  CLOSE(100)

  OPEN(unit=100, file="z_simul_loc.txt", status="old", action="read", position="rewind")
  READ (100, *) z_simul_loc
  CLOSE(100)

  d1 = x(1)
  d2 = x(2)
  pi = x(3)
  beta = x(4)
  d0 = 0.0

  if (need_gradient.ne.0) then
     grad(1) = 0.0
     grad(2) = 0.0
     grad(3) = 0.0
   end if

  !CALL omp_set_num_threads(8)

  Zprob = transpose(Zprob)


  ! Solving for taxes in default

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(zix, pen, tix, c, l, g, V_tax, util)
  DO zix = 1, zsiz
    CALL penalty(Z(zix),d0, d1, d2, pen)
    DO tix = 1, tsiz
      c = ((1.0/chi)**(1.0/gamma))*((1.0-t(tix))*pen)**(1.0+1.0/gamma)
      l = ((1.0-t(tix))*pen/chi)**(1.0/gamma)
      g = (pi*m)**(1.0/sigmag)
      B_def(zix, tix) = n*t(tix)*pen*l - g
      CALL utility(c,l,g,B_def(zix, tix)/m, pi, util)
      V_tax(tix) = util
    END DO
      tax_def(zix) = MAXLOC(V_tax, 1)
      gov_def(zix) = g
      IF (B_def(zix, tax_def(zix)) < 0.0) THEN
        DO tix = 1, tsiz
          c = ((1.0/chi)**(1.0/gamma))*((1.0-t(tix))*pen)**(1.0+1.0/gamma)
          l = ((1.0-t(tix))*pen/chi)**(1.0/gamma)
          g = n*t(tix)*pen*l
          B_def(zix, tix) = 0.0
        CALL utility(c,l,g,B_def(zix, tix)/m, pi,util)
        V_tax(tix) = util
      END DO
      tax_def(zix) = MAXLOC(V_tax, 1)
      gov_def(zix) = n*t(tax_def(zix))*pen*((1.0-t(tax_def(zix)))*pen/chi)**(1.0/gamma)
    END IF
  END DO
!$OMP END PARALLEL DO

    errq = 1.0
    errv = 1.0
    J0 = 0.0
    Jd0 = 0.0

    q0 = 1.0/r

    iter = 0

    ! The main computation loop


    DO WHILE ((errq>tolq .or. errv>tolv) .and. (iter<maxiter))


      ! Solving for taxes in repayment

!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(zix, bix, bpix, tix, c, l, g, V_tax, util)
      DO zix = 1, zsiz
        DO bix = 1, bsiz
          DO bpix = 1, bsiz
            DO tix = 1, tsiz
             c = ((1.0/chi)**(1.0/gamma))*((1.0-t(tix))*Z(zix))**(1.0+1.0/gamma)
             l = ((1.0-t(tix))*Z(zix)/chi)**(1.0/gamma)
             g = (pi*m)**(1.0/sigmag)
             B_rep(zix, bix, bpix, tix) = n*t(tix)*Z(zix)*l - g +&
             b(bix) - b(bpix)*q0(zix, bpix)
             CALL utility(c, l, g, B_rep(zix, bix, bpix, tix)/m, pi, util)
             V_tax(tix) = util
           END DO
           tax_rep(zix, bix, bpix) = MAXLOC(V_tax,1)
           gov_rep(zix, bix, bpix) = g
            IF (B_rep(zix, bix, bpix, tax_rep(zix, bix, bpix)) < 0.0) THEN
              DO tix = 1, tsiz
                c = ((1.0/chi)**(1.0/gamma))*((1.0-t(tix))*Z(zix))**(1.0+1.0/gamma)
                l = ((1.0-t(tix))*Z(zix)/chi)**(1.0/gamma)
                g = n*t(tix)*Z(zix)*l +&
                b(bix) - b(bpix)*q0(zix, bpix)
                B_rep(zix, bix, bpix, tix) = 0.0
             CALL utility(c, l, g, B_rep(zix, bix, bpix, tix)/m, pi, util)
             V_tax(tix) = util
           END DO
             tax_rep(zix, bix, bpix) = MAXLOC(V_tax,1)
             gov_rep(zix, bix, bpix) = n*t(tax_rep(zix, bix, bpix))*Z(zix)*&
             ((1.0-t(tax_rep(zix, bix, bpix)))*Z(zix)/chi)**(1.0/gamma) +&
             b(bix) - b(bpix)*q0(zix, bpix)
           END IF
         END DO
       END DO
     END DO
!$OMP END PARALLEL DO

     ! Value function iteration for default values

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(zix, pen, c, l, g, util)
     DO zix = 1, zsiz
       CALL penalty(Z(zix),d0, d1, d2, pen)
       c = ((1.0/chi)**(1.0/gamma))*((1.0-t(tax_def(zix)))*pen)**(1.0+1.0/gamma)
       l = ((1.0-t(tax_def(zix)))*pen/chi)**(1.0/gamma)
       g = gov_def(zix)
       CALL utility(c, l, g, B_def(zix, tax_def(zix))/m, pi, util)
       Vd(zix) = util + beta*DOT_PRODUCT(Zprob(zix,:), theta*J0(:, bsiz)+(1.0-theta)*Jd0)
     END DO
!$OMP END PARALLEL DO

     ! Value function iteration for repayment values

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix, bpix, c, l, g, Vpos, util)
     DO zix = 1, zsiz
       DO bix = 1, bsiz
         DO bpix = 1, bsiz
           c = ((1.0/chi)**(1.0/gamma))*((1.0-t(tax_rep(zix, bix, bpix)))*Z(zix))**(1.0+1.0/gamma)
           l = ((1.0-t(tax_rep(zix, bix, bpix)))*Z(zix)/chi)**(1.0/gamma)
           g = gov_rep(zix, bix, bpix)
           CALL utility(c, l, g, B_rep(zix, bix, bpix, tax_rep(zix, bix, bpix))/m, pi, util)
           Vpos(bpix) = util + beta*DOT_PRODUCT(Zprob(zix,:), J0(:, bpix))
        END DO
         b_pol(zix, bix) = MAXLOC(Vpos, 1)
         Vc(zix, bix) = MAXVAL(Vpos)
       END DO
     END DO
     !$OMP END PARALLEL DO

  ! Default decision

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix)
  DO zix = 1, zsiz
  DO bix = 1, bsiz
    IF (Vc(zix, bix)>Vd(zix)) THEN
      d(zix, bix) = 0
    ELSE
      d(zix, bix) = 1
    END IF
  END DO
  END DO
  !$OMP END PARALLEL DO

  ! Constructing the continuation values
  ! In default

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(zix, pen, c, l, g, util)
  DO zix = 1, zsiz
  CALL penalty(Z(zix),d0, d1, d2, pen)
  c = ((1.0/chi)**(1.0/gamma))*((1.0-t(tax_def(zix)))*pen)**(1.0+1.0/gamma)
  l = ((1.0-t(tax_def(zix)))*pen/chi)**(1.0/gamma)
  g = gov_def(zix)
  CALL utility(c, l, g, B_def(zix, tax_def(zix))/n, pi, util)
  Jd1(zix) = util + beta*DOT_PRODUCT(Zprob(zix,:), theta*J0(:, bsiz)+(1.0-theta)*Jd0)
  END DO
  !$OMP END PARALLEL DO

  ! In repayment

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix, c, l, g, util)
  DO zix = 1, zsiz
   DO bix = 1, bsiz
    c = ((1.0/chi)**(1.0/gamma))*((1.0-t(tax_rep(zix, bix, b_pol(zix, bix))))*Z(zix))**(1.0+1.0/gamma)
    l = ((1.0-t(tax_rep(zix, bix, b_pol(zix, bix))))*Z(zix)/chi)**(1.0/gamma)
    g = gov_rep(zix, bix, b_pol(zix, bix))
    CALL utility(c, l, g, B_rep(zix, bix, b_pol(zix, bix), tax_rep(zix, bix, b_pol(zix, bix)))/n, pi, util)
    Jc1(zix, bix) = util + beta*DOT_PRODUCT(Zprob(zix,:), J0(:,b_pol(zix, bix)))
   END DO
  END DO
  !$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix)
  DO zix = 1, zsiz
   DO bix = 1, bsiz
    IF (d(zix, bix) == 0) THEN
      J1(zix, bix) = Jc1(zix, bix)
    ELSE
      J1(zix, bix) = Jd1(zix)
    END IF
   END DO
  END DO
!$OMP END PARALLEL DO

  ! The bond price schedule

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix)
  DO zix = 1, zsiz
   DO bix = 1, bsiz
    q1(zix, bix) = (1.0/r)*DOT_PRODUCT(Zprob(zix,:),(1.0-d(:,bix)))
   END DO
  END DO
!$OMP END PARALLEL DO

  errq = MAXVAL(abs(q0-q1))
  errv = MAX(MAXVAL(abs(J1-J0)),MAXVAL(abs(Jd1-Jd0)))

  WRITE(*,*) "iter ", iter, " errq = ",errq,  "J1-J0", MAXVAL(abs(J1-J0)), &
  "Jd1-Jd0",MAXVAL(abs(Jd1-Jd0))

  iter = iter+1

  J0 = J1
  q0 = 0.5*q1+0.5*q0
  Jd0 = Jd1

  END DO

  ! Print out policies for comparison with the other code

  ! Bond Policy

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix)
  DO zix = 1, zsiz
   DO bix = 1, bsiz
    IF (d(zix, bix) == 0) THEN
      Bpol(zix, bix) = b(b_pol(zix, bix))
    ELSE
      Bpol(zix, bix) = 0.0
    END IF
   END DO
  END DO
  !$OMP END PARALLEL DO

  ! Transfers in repayment

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix)
  DO zix = 1, zsiz
   DO bix = 1, bsiz
    Brep(zix, bix) = B_rep(zix, bix, b_pol(zix, bix), tax_rep(zix, bix, b_pol(zix, bix)))
   END DO
  END DO
  !$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(zix)
  DO zix = 1, zsiz
    Bdef(zix) = B_def(zix, tax_def(zix))
  END DO
!$OMP END PARALLEL DO

  ! Taxes in repayment

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix)
  DO zix = 1, zsiz
   DO bix = 1, bsiz
    taxrep(zix, bix) = t(tax_rep(zix, bix, b_pol(zix, bix)))
   END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(zix)
  DO zix = 1, zsiz
    taxdef(zix) = t(tax_def(zix))
  END DO
!$OMP END PARALLEL DO

  ! Adjust policies for default

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix)
  DO zix = 1, zsiz
   DO bix = 1, bsiz
    IF (d(zix, bix)==1) THEN
      tax(zix, bix) = taxdef(zix)
      Tr(zix, bix) = Bdef(zix)
      gov(zix, bix) = gov_def(zix)
    ELSE
      tax(zix, bix) = taxrep(zix, bix)
      Tr(zix, bix) = Brep(zix, bix)
      gov(zix, bix) = gov_rep(zix, bix, b_pol(zix, bix))
    END IF
   END DO
  END DO
!$OMP END PARALLEL DO

  ! Print out transfers in repayment and default

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix, pen, c, l, g, Po, util)
  DO zix = 1, zsiz
   DO bix = 1, bsiz
    CALL penalty(Z(zix),d0, d1, d2, pen)
    c = ((1.0/chi)**(1.0/gamma))*((1.0-t(tax_rep(zix, bix, b_pol(zix, bix))))*Z(zix))**(1.0+1.0/gamma)
    l = ((1.0-t(tax_rep(zix, bix, b_pol(zix, bix))))*Z(zix)/chi)**(1.0/gamma)
    g = gov_rep(zix, bix, b_pol(zix, bix))
    Po = 0.0
    CALL utility(c, l, g, Po, pi, util)
    srep(zix, bix) =1.0/shi*(J1(zix, bix)- &
    beta*DOT_PRODUCT(Zprob(zix,:), J0(:,b_pol(zix, bix)))-util)
    c = ((1.0/chi)**(1.0/gamma))*((1.0-t(tax_def(zix)))*pen)**(1.0+1.0/gamma)
    l = ((1.0-t(tax_def(zix)))*pen/chi)**(1.0/gamma)
    g = gov_def(zix)
    CALL utility(c, l, g, Po, pi, util)
    sdef(zix) = Jd0(zix)- &
    beta*DOT_PRODUCT(Zprob(zix,:), theta*J0(:,bsiz)+(1.0-theta)*Jd0)-util
   END DO
  END DO
  !$OMP END PARALLEL DO

  ! Proposer's transfers

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix)
  DO zix = 1, zsiz
   DO bix = 1, bsiz
    Prop_rep(zix, bix) = Brep(zix, bix)-(m-1.0)*srep(zix, bix)
    Prop_def(zix, bix) = Bdef(zix)-(m-1.0)*sdef(zix)
   END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix)
  DO zix = 1, zsiz
   DO bix = 1, bsiz
    IF (d(zix, bix) == 1) THEN
     Prop(zix, bix) = Prop_def(zix, bix)
    ELSE
     Prop(zix, bix) = Prop_rep(zix, bix)
    END IF
   END DO
  END DO
!$OMP END PARALLEL DO

! GDP, consumption and labor

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(zix, bix, pen)
DO zix = 1, zsiz
  DO bix = 1, bsiz
    IF (d(zix, bix) == 1) THEN
      CALL penalty(Z(zix),d0, d1, d2, pen)
    labor(zix, bix) = ((1.0-t(tax_def(zix)))*pen/chi)**(1.0/gamma)
    GDP(zix, bix) = pen*labor(zix, bix)*n
    cons(zix, bix) = (1.0-t(tax_def(zix)))*pen*labor(zix, bix)
  ELSE
    labor(zix, bix) = ((1.0-t(tax_rep(zix, bix, b_pol(zix, bix))))*Z(zix)/chi)**(1.0/gamma)
    GDP(zix, bix) = Z(zix)*labor(zix, bix)*n
    cons(zix, bix) = (1.0-t(tax_rep(zix, bix, b_pol(zix, bix))))*Z(zix)*labor(zix, bix)
  END IF
 END DO
END DO
!$OMP END PARALLEL DO

! Simulations

d_simul(1) = 1
b_simul_loc(1) = bsiz


DO iter = 2, maximum
  IF (d_simul(iter-1)==1) THEN
    CALL RANDOM_NUMBER(unif)
    IF (unif>=theta) THEN
      d_simul(iter) = 1
      b_simul_loc(iter) = bsiz
    ELSE
      d_simul(iter) = 0
      b_simul_loc(iter) = bsiz
    END IF
  ELSE
    b_simul_loc(iter) = b_pol(z_simul_loc(iter-1),b_simul_loc(iter-1))
    d_simul(iter) = d(z_simul_loc(iter),b_simul_loc(iter))
  END IF
END DO

! GDP, consumption, government spending  and labor simulations

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(iter)
DO iter = 1, maximum
  c_simul(iter) = cons(z_simul_loc(iter), b_simul_loc(iter))*n
  g_simul(iter) = gov(z_simul_loc(iter), b_simul_loc(iter))
  labor_simul(iter) = labor(z_simul_loc(iter), b_simul_loc(iter))*n
  GDP_simul(iter) = GDP(z_simul_loc(iter), b_simul_loc(iter))
  tax_simul(iter) = tax(z_simul_loc(iter), b_simul_loc(iter))
  Tr_simul(iter) = Tr(z_simul_loc(iter), b_simul_loc(iter))
END DO
!$OMP END PARALLEL DO

! Spread simulation

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(iter)
DO iter = 1, maximum-1
  IF (d_simul(iter) == 1) THEN
    q_simul(iter) = 100.0
  ELSE
    q_simul(iter) = q0(z_simul_loc(iter),b_simul_loc(iter+1))
  END IF
END DO
!$OMP END PARALLEL DO

! Price to spreads

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(iter)
DO iter = 1, maximum-1
  IF (d_simul(iter) == 1) THEN
    spread_simul(iter) = 100.0
  ELSE
    spread_simul(iter) = (((1.0/q_simul(iter))**4.0)-r) ! Annualized Percentage
  END IF
END DO
!$OMP END PARALLEL DO

! Simulation of Trade Balance to GDP

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(iter)
DO iter = 1, maximum-1
  IF (d_simul(iter) == 1) THEN
    TB_simul(iter) = 0.0
  ELSE
  TB_simul(iter) = (q_simul(iter)*b(b_simul_loc(iter+1))-b(b_simul_loc(iter)))/GDP_simul(iter)
 END IF
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(iter)
DO iter = 1, maximum
  b_to_y(iter) = -b(b_simul_loc(iter))/(4.0*GDP_simul(iter))
END DO
!$OMP END PARALLEL DO

! Government Spending to GDP

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(iter)
DO iter = 1, maximum
  g_to_y(iter) = g_simul(iter)/GDP_simul(iter)
END DO
!$OMP END PARALLEL DO

d_simulate = d_simul(5000:95000)
b_simulate_loc = b_simul_loc(5000:95000)
g_simulate = g_simul(5000:95000)
c_simulate = c_simul(5000:95000)
GDP_simulate = GDP_simul(5000:95000)
q_simulate = q_simul(5000:95000)
TB_simulate = TB_simul(5000:95000)
spread_simulate = spread_simul(5000:95000)
tax_simulate = tax_simul(5000:95000)
b_to_y_simulate = b_to_y(5000:95000)
g_to_y_simulate = g_to_y(5000:95000)
Tr_simulate = Tr_simul(5000:95000)

! Work with 30 quarters prior to a default

    spds_rel = 1000.0
    g_y_rel = 1000.0
    b_y_rel = 1000.0
    sd_spds_rel = 1000.0
    sd_c_rel = 1000.0
    sd_GDP_rel = 1000.0
    sd_tax_rel = 1000.0
    sd_TB_rel = 1000.0
    rho_tb_y_rel = 1000.0
    rho_g_y_rel = 1000.0
    rho_c_y_rel = 1000.0
    rho_spds_y_rel = 1000.0
    rho_tax_y_rel = 1000.0

    DO iter = allow, 1+period, -1
      IF (d_simulate(iter) == 1 .and. SUM(d_simulate(iter-period:iter-1)) == 0) THEN
        CALL mean1(spread_simulate(iter-period:iter-1), m_spds)
        spds_rel(iter) = m_spds
        CALL mean1(g_to_y_simulate(iter-period:iter-1), avg_g_to_y)
        g_y_rel(iter) = avg_g_to_y
        CALL mean1(b_to_y_simulate(iter-period:iter-1), avg_b_to_y)
        b_y_rel(iter) = avg_b_to_y
        CALL variance(spread_simulate(iter-period:iter-1),var_spds)
        sd_spds_rel(iter) = SQRT(var_spds)
        CALL variance(log(c_simulate(iter-period:iter-1)),var_c)
        sd_c_rel(iter) = SQRT(var_c)
        CALL variance(log(GDP_simulate(iter-period:iter-1)),var_GDP)
        sd_GDP_rel(iter) = SQRT(var_GDP)
        CALL variance(log(tax_simulate(iter-period:iter-1)),var_tax)
        sd_tax_rel(iter) = SQRT(var_tax)
        CALL variance(TB_simulate(iter-period:iter-1),var_TB)
        sd_TB_rel(iter) = SQRT(var_TB)
        CALL corr(TB_simulate(iter-period:iter-1),log(GDP_simulate(iter-period:iter-1)),rho_tb_y)
        rho_tb_y_rel(iter) = rho_TB_y
        CALL corr(log(g_simulate(iter-period:iter-1)),log(GDP_simulate(iter-period:iter-1)),rho_g_y)
        rho_g_y_rel(iter) = rho_g_y
        CALL corr(log(c_simulate(iter-period:iter-1)),log(GDP_simulate(iter-period:iter-1)),rho_c_y)
        rho_c_y_rel(iter) = rho_c_y
        CALL corr(spread_simulate(iter-period:iter-1),log(GDP_simulate(iter-period:iter-1)),rho_spds_y)
        rho_spds_y_rel(iter) = rho_spds_y
        CALL corr(log(tax_simulate(iter-period:iter-1)),log(GDP_simulate(iter-period:iter-1)),rho_tax_y)
        rho_tax_y_rel(iter) = rho_tax_y
      ELSE
        spds_rel(iter) = 1000.0
        g_y_rel(iter) = 1000.0
        b_y_rel(iter) = 1000.0
        sd_spds_rel(iter) = 1000.0
        sd_c_rel(iter) = 1000.0
        sd_GDP_rel(iter) = 1000.0
        sd_tax_rel(iter) = 1000.0
        sd_TB_rel(iter) = 1000.0
        rho_tb_y_rel(iter) = 1000.0
        rho_g_y_rel(iter) = 1000.0
        rho_c_y_rel(iter) = 1000.0
        rho_spds_y_rel(iter) = 1000.0
        rho_tax_y_rel(iter) = 1000.0
      END IF
    END DO

    CALL mean(spds_rel, m_spds)
    CALL mean(g_y_rel, avg_g_to_y)
    CALL mean(b_y_rel, avg_b_to_y)
    CALL mean(sd_spds_rel, sd_spds)
    CALL mean(sd_c_rel, sd_c)
    CALL mean(sd_tax_rel, sd_tax)
    CALL mean(sd_GDP_rel, sd_GDP)
    CALL mean(sd_TB_rel, sd_TB)
    CALL mean(rho_tb_y_rel, rho_tb_y)
    CALL mean(rho_g_y_rel, rho_g_y)
    CALL mean(rho_c_y_rel, rho_c_y)
    CALL mean(rho_spds_y_rel, rho_spds_y)
    CALL mean(rho_tax_y_rel, rho_tax_y)
    sd_c_y = sd_c/sd_GDP


  data_m_spds = 0.0815
  data_sd_spds = 0.0443
  data_avg_g_to_y = 0.17
  data_rho_tb_y  = -0.88 ! 0.62 if detrended
  data_avg_b_to_y = 0.053
  data_rho_spds_y = -0.79
  data_sd_c_y = 1.09
  data_sd_c = 0.083
  data_sd_y = 0.076


  loss_array = (/0.5*(m_spds-data_m_spds)**2.0,&
  0.4*(sd_spds-data_sd_spds)**2.0,0.05*(avg_g_to_y-data_avg_g_to_y)**2.0,&
  0.05*(data_avg_b_to_y-avg_b_to_y)**2.0/)

  loss_value = SUM(loss_array)
  WRITE(*,*) "this step beta=", x(4)
  WRITE(*,*) "this step d2=", x(2)
  WRITE(*,*) "this step d1=", x(1)
  WRITE(*,*) "this step pi=", x(3)
  WRITE(*,*) "this step minf=", loss_value

  WRITE(*,*) "sd_spds", sd_spds
  WRITE(*,*) "g/y", avg_g_to_y
  WRITE(*,*) "m_spds", m_spds
  WRITE(*,*) "sd_c_y", sd_c_y
  WRITE(*,*) "sd_c", sd_c
  WRITE(*,*) "sd_y", sd_GDP
  WRITE(*,*) "b/y", avg_b_to_y
  WRITE(*,*) "rho_spds_y", rho_spds_y
  WRITE(*,*) "rho_tb_y", rho_tb_y
  WRITE(*,*) "rho_c_y", rho_c_y
  WRITE(*,*) "rho_g_y", rho_g_y
  WRITE(*,*) "rho_tax_y", rho_tax_y
  WRITE(*,*) "sd_tax", sd_tax
  WRITE(*,*) "sd_TB_by_y", sd_TB/sd_GDP



  OPEN(unit=100, file="Bpol_linear.txt")
  WRITE (100, "(F20.10)") Bpol
  CLOSE(100)

  OPEN(unit=100, file="tax_linear.txt")
  WRITE (100, "(F20.10)") tax
  CLOSE(100)

  OPEN(unit=100, file="GDP_linear.txt")
  WRITE (100, "(F20.10)") GDP
  CLOSE(100)

  OPEN(unit=100, file="q_linear.txt")
  WRITE (100, "(F20.10)") q0
  CLOSE(100)

  OPEN(unit=100, file="cons_linear.txt")
  WRITE (100, "(F20.10)") cons
  CLOSE(100)

  OPEN(unit=100, file="labor_linear.txt")
  WRITE (100, "(F20.10)") labor
  CLOSE(100)

  OPEN(unit=100, file="Tr_linear.txt")
  WRITE (100, "(F20.10)") Tr
  CLOSE(100)

  OPEN(unit=100, file="g_linear.txt")
  WRITE (100, "(F20.10)") gov
  CLOSE(100)

  OPEN(unit=100, file="d_linear.txt")
  WRITE (100, "(I20)") d
  CLOSE(100)

  OPEN(unit=100, file="GDP_simul_linear.txt")
  WRITE (100, "(F20.10)") GDP_simulate
  CLOSE(100)

  OPEN(unit=100, file="Tr_simul_linear.txt")
  WRITE (100, "(F20.10)") Tr_simulate
  CLOSE(100)

  OPEN(unit=100, file="cons_simul_linear.txt")
  WRITE (100, "(F20.10)") c_simulate
  CLOSE(100)

  OPEN(unit=100, file="gov_simul_linear.txt")
  WRITE (100, "(F20.10)") g_simulate
  CLOSE(100)

  OPEN(unit=100, file="tax_simul_linear.txt")
  WRITE (100, "(F20.10)") tax_simulate
  CLOSE(100)

  OPEN(unit=100, file="TB_simul_linear.txt")
  WRITE (100, "(F20.10)") TB_simulate
  CLOSE(100)

  OPEN(unit=100, file="spread_simul_linear.txt")
  WRITE (100, "(F20.10)") spread_simulate
  CLOSE(100)

  OPEN(unit=100, file="d_simul_linear.txt")
  WRITE (100, "(I20)") d_simulate
  CLOSE(100)

  OPEN(unit=100, file="Vc_linear.txt")
  WRITE (100, "(F20.10)") Vc
  CLOSE(100)

  OPEN(unit=100, file="Vd_linear.txt")
  WRITE (100, "(F20.10)") Vd
  CLOSE(100)

  OPEN(unit=100, file="Bdef_linear.txt")
  WRITE (100, "(F20.10)") Bdef
  CLOSE(100)

  OPEN(unit=100, file="taxdef_linear.txt")
  WRITE (100, "(F20.10)") taxdef
  CLOSE(100)

  OPEN(unit=100, file="govdef_linear.txt")
  WRITE (100, "(F20.10)") gov_def
  CLOSE(100)




END SUBROUTINE loss_linear

END MODULE main_linear_g
