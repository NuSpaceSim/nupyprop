module constants

implicit none
integer, parameter :: dp=kind(0.d0)
real(dp), parameter :: pi = 3.1415927_dp
real(dp), parameter :: N_A = 6.0221409e+23_dp
real(dp), parameter :: rho_water = 1.02_dp ! g/cm^3
real(dp), parameter :: rho_rock = 2.65_dp ! g/cm^3
real(dp), parameter :: rho_iron = 7.87_dp ! g/cm^3
real(dp), parameter :: R_earth = 6371.0_dp ! radius of the earth in km
real(dp), parameter :: step_size = 4500.0_dp ! step size for continuous energy loss, in cm
real(dp), parameter :: E_nu(91) = (/1.00000000e+03, 1.25892541e+03, 1.58489319e+03, 1.99526231e+03,&
       & 2.51188643e+03, 3.16227766e+03, 3.98107171e+03, 5.01187234e+03,&
       & 6.30957344e+03, 7.94328235e+03, 1.00000000e+04, 1.25892541e+04,&
       & 1.58489319e+04, 1.99526231e+04, 2.51188643e+04, 3.16227766e+04,&
       & 3.98107171e+04, 5.01187234e+04, 6.30957344e+04, 7.94328235e+04,&
       & 1.00000000e+05, 1.25892541e+05, 1.58489319e+05, 1.99526231e+05,&
       & 2.51188643e+05, 3.16227766e+05, 3.98107171e+05, 5.01187234e+05,&
       & 6.30957344e+05, 7.94328235e+05, 1.00000000e+06, 1.25892541e+06,&
       & 1.58489319e+06, 1.99526231e+06, 2.51188643e+06, 3.16227766e+06,&
       & 3.98107171e+06, 5.01187234e+06, 6.30957344e+06, 7.94328235e+06,&
       & 1.00000000e+07, 1.25892541e+07, 1.58489319e+07, 1.99526231e+07,&
       & 2.51188643e+07, 3.16227766e+07, 3.98107171e+07, 5.01187234e+07,&
       & 6.30957344e+07, 7.94328235e+07, 1.00000000e+08, 1.25892541e+08,&
       & 1.58489319e+08, 1.99526231e+08, 2.51188643e+08, 3.16227766e+08,&
       & 3.98107171e+08, 5.01187234e+08, 6.30957344e+08, 7.94328235e+08,&
       & 1.00000000e+09, 1.25892541e+09, 1.58489319e+09, 1.99526231e+09,&
       & 2.51188643e+09, 3.16227766e+09, 3.98107171e+09, 5.01187234e+09,&
       & 6.30957344e+09, 7.94328235e+09, 1.00000000e+10, 1.25892541e+10,&
       & 1.58489319e+10, 1.99526231e+10, 2.51188643e+10, 3.16227766e+10,&
       & 3.98107171e+10, 5.01187234e+10, 6.30957344e+10, 7.94328235e+10,&
       & 1.00000000e+11, 1.25892541e+11, 1.58489319e+11, 1.99526231e+11,&
       & 2.51188643e+11, 3.16227766e+11, 3.98107171e+11, 5.01187234e+11,&
       & 6.30957344e+11, 7.94328235e+11, 1.00000000e+12/)
real(dp), parameter :: E_lep(121) = (/1.00000000e+00, 1.25892541e+00, 1.58489319e+00, 1.99526231e+00,&
       & 2.51188643e+00, 3.16227766e+00, 3.98107171e+00, 5.01187234e+00,&
       & 6.30957344e+00, 7.94328235e+00, 1.00000000e+01, 1.25892541e+01,&
       & 1.58489319e+01, 1.99526231e+01, 2.51188643e+01, 3.16227766e+01,&
       & 3.98107171e+01, 5.01187234e+01, 6.30957344e+01, 7.94328235e+01,&
       & 1.00000000e+02, 1.25892541e+02, 1.58489319e+02, 1.99526231e+02,&
       & 2.51188643e+02, 3.16227766e+02, 3.98107171e+02, 5.01187234e+02,&
       & 6.30957344e+02, 7.94328235e+02, 1.00000000e+03, 1.25892541e+03,&
       & 1.58489319e+03, 1.99526231e+03, 2.51188643e+03, 3.16227766e+03,&
       & 3.98107171e+03, 5.01187234e+03, 6.30957344e+03, 7.94328235e+03,&
       & 1.00000000e+04, 1.25892541e+04, 1.58489319e+04, 1.99526231e+04,&
       & 2.51188643e+04, 3.16227766e+04, 3.98107171e+04, 5.01187234e+04,&
       & 6.30957344e+04, 7.94328235e+04, 1.00000000e+05, 1.25892541e+05,&
       & 1.58489319e+05, 1.99526231e+05, 2.51188643e+05, 3.16227766e+05,&
       & 3.98107171e+05, 5.01187234e+05, 6.30957344e+05, 7.94328235e+05,&
       & 1.00000000e+06, 1.25892541e+06, 1.58489319e+06, 1.99526231e+06,&
       & 2.51188643e+06, 3.16227766e+06, 3.98107171e+06, 5.01187234e+06,&
       & 6.30957344e+06, 7.94328235e+06, 1.00000000e+07, 1.25892541e+07,&
       & 1.58489319e+07, 1.99526231e+07, 2.51188643e+07, 3.16227766e+07,&
       & 3.98107171e+07, 5.01187234e+07, 6.30957344e+07, 7.94328235e+07,&
       & 1.00000000e+08, 1.25892541e+08, 1.58489319e+08, 1.99526231e+08,&
       & 2.51188643e+08, 3.16227766e+08, 3.98107171e+08, 5.01187234e+08,&
       & 6.30957344e+08, 7.94328235e+08, 1.00000000e+09, 1.25892541e+09,&
       & 1.58489319e+09, 1.99526231e+09, 2.51188643e+09, 3.16227766e+09,&
       & 3.98107171e+09, 5.01187234e+09, 6.30957344e+09, 7.94328235e+09,&
       & 1.00000000e+10, 1.25892541e+10, 1.58489319e+10, 1.99526231e+10,&
       & 2.51188643e+10, 3.16227766e+10, 3.98107171e+10, 5.01187234e+10,&
       & 6.30957344e+10, 7.94328235e+10, 1.00000000e+11, 1.25892541e+11,&
       & 1.58489319e+11, 1.99526231e+11, 2.51188643e+11, 3.16227766e+11,&
       & 3.98107171e+11, 5.01187234e+11, 6.30957344e+11, 7.94328235e+11,&
       & 1.00000000e+12/)
real(dp), parameter, dimension(31) :: v2 = (/0.,-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1. ,&
       & -1.1, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.8, -1.9, -2. , -2.1,&
       & -2.2, -2.3, -2.4, -2.5, -2.6, -2.7, -2.8, -2.9, -3./) ! x (interpolating arr for find_y)
!real(dp), parameter, dimension(2) :: E_nu = (/1.,2./)

!real(dp), parameter, dimension(10) :: A=(/(j,j=1,10)/)

!real(dp) :: E_nu(:)
!real(dp):: E_lep(:)
!integer :: i, j
!
!do i = 0, 90
!    E_nu(i) = 10.**(3 + i*(12 - 3)/dble(91-1))
!end do
!
!do j = 0, 120
!    E_lep(i) = 10.**(0 + i*(12 - 0)/dble(121-1))
!end do

end module constants

module geometry

use constants
implicit none
contains

subroutine PREMdensity(Rin, idepth, edens)

    implicit none
    !!f2py threadsafe
    integer, intent(in) :: idepth
    real(dp), intent(in) :: Rin
    !f2py intent(out) :: edens
    real(dp) edens
    real(dp) :: x,y
    real(dp), dimension(10) :: Rlay

    Rlay = (/1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0, 6346.6, 6356.0, 6368.0, 6371.0/) ! # PREM layers based on R_earth
!    print *, E_nu
    Rlay(9) = 6368.0_dp + (3.0_dp-dble(idepth))
!    print *,'Rlay=',Rlay
    x = Rin
    y = x/R_earth

!    print *,'x=',x

    if (x<=Rlay(1)) then
        edens = 13.0885_dp-8.8381_dp*y**2
    else if (x<=Rlay(2)) then
        edens = 12.5815_dp-1.2638_dp*y-3.6426_dp*y**2-5.5281_dp*y**3
    else if (x<=Rlay(3)) then
        edens = 7.9565_dp-6.4761_dp*y+5.5283_dp*y**2-3.0807_dp*y**3
    else if (x<=Rlay(4)) then
        edens = 5.3197_dp-1.4836_dp*y
    else if (x<=Rlay(5)) then
        edens = 11.2494_dp-8.0298_dp*y
    else if (x<=Rlay(6)) then
        edens = 7.1089_dp-3.8045_dp*y
    else if (x<=Rlay(7)) then
        edens = 2.6910_dp+0.6924_dp*y
    else if (x<=Rlay(8)) then
        edens = 2.900_dp
    else if (x<=Rlay(9)) then
        edens = 2.600_dp
    else if (x<=Rlay(10)) then
        edens = 1.020_dp
    else if (x<=Rlay(10)*1.001_dp) then ! too close to call!
        edens = 1.020_dp
    else
        edens=0._dp
    end if
end subroutine

subroutine densityatx(x, beta_deg, idepth, r, rho_at_x)

    implicit none
    !!f2py threadsafe
    integer, intent(in) :: idepth
    real(dp), intent(in) :: x, beta_deg
    !f2py intent(out) :: r, rho_at_x
    real(dp) r, rho_at_x
    real(dp) :: tnadir, ell, r2

    tnadir = (90.0_dp-beta_deg)*(pi/180.0_dp)
    ell = R_earth*dcos(tnadir)*2
    r2 = x**2 - (ell*x) + R_earth**2

    if (beta_deg < 5.0_dp) then
        r = R_earth*(1.0_dp + 0.5_dp*(x**2-ell*x)/R_earth**2)
    else
        r = dsqrt(r2)
    end if

!    print *,'r=',r

    call PREMdensity(r, idepth, rho_at_x)
end subroutine

end module geometry

module interpolation
use constants
implicit none
contains

subroutine interp_linear_internal(x, y, xout, yout)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in)  :: x(2), y(2), xout
    !f2py intent(out) :: yout
    real(dp) yout
    real(dp) :: alph

    if ( xout .lt. x(1) .or. xout .gt. x(2) ) then
        write(*,*) "interp1: xout < x0 or xout > x1 !"
        write(*,*) "xout = ",xout
        write(*,*) "x0   = ",x(1)
        write(*,*) "x1   = ",x(2)
        stop
    end if

    alph = (xout - x(1)) / (x(2) - x(1))
    yout = y(1) + alph*(y(2) - y(1))

end subroutine interp_linear_internal

subroutine interp_linear_pt(x, y, xout, yout)
    ! Interpolate y from ordered x to ordered xout positions

    implicit none
    !!f2py threadsafe
    real(dp), intent(in), dimension(:) :: x, y
    real(dp), intent(in) :: xout
    !f2py intent(out) :: yout
    real(dp) yout
    integer :: j, n

    n = size(x)

    if (xout < x(1)) then
        yout = y(1)
    else if (xout > x(n)) then
        yout = y(n)
    else
        do j = 1, n
            if (x(j) >= xout) exit
        end do

        if (j == 1) then
            yout = y(1)
        else if (j == n+1) then
            yout = y(n)
        else
            call interp_linear_internal(x(j-1:j),y(j-1:j),xout,yout)
        end if
    end if

end subroutine interp_linear_pt

subroutine interpol(x_val, x, y, interp_val)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(in) :: x_val
    !f2py intent(out) :: interp_val
    real(dp) interp_val

    call interp_linear_pt(x, y, x_val, interp_val)
end subroutine interpol

subroutine cd2distd(xalong, cdalong, col_depth, out_val)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: xalong(:), cdalong(:)
    real(dp), intent(in) :: col_depth
    !f2py intent(out) :: out_val
    real(dp) out_val
!    print *, xalong(1)

    if (col_depth < minval(cdalong)) then
        out_val = (col_depth/cdalong(1))*xalong(1)
    else if (col_depth > maxval(cdalong)) then
        out_val = maxval(xalong)
    else
        call interpol(col_depth, cdalong, xalong, out_val)
    end if

end subroutine cd2distd

subroutine int_xc_nu(energy, nu_xc, sig_cc, sig_nc)

    implicit none
    !!f2py threadsafe
!    real(dp), intent(in) :: E_nu(:)
    real(dp), intent(in) :: nu_xc(:,:)
    real(dp), intent(in) :: energy
    !f2py intent(out) :: sig_cc, sig_nc
    real(dp) sig_cc, sig_nc

    call interpol(energy, E_nu, nu_xc(:,1), sig_cc)
    sig_cc = N_A * sig_cc
    call interpol(energy, E_nu, nu_xc(:,2), sig_nc)
    sig_nc = N_A * sig_nc

end subroutine int_xc_nu

subroutine int_xc_lep(energy, xc_arr, frac, frac_pn, sig_brem, sig_pair, sig_pn)

    implicit none
    !!f2py threadsafe
!    real(dp), intent(in) :: E_lep(:)
    real(dp), intent(in) :: xc_arr(:,:)
    real(dp), intent(in) :: energy, frac, frac_pn
    !f2py intent(out) :: sig_brem, sig_pair, sig_pn
    real(dp) sig_brem, sig_pair, sig_pn

    call interpol(energy, E_lep, xc_arr(:,1), sig_brem)
    call interpol(energy, E_lep, xc_arr(:,2), sig_pair)
    call interpol(energy, E_lep, xc_arr(:,3), sig_pn)

    sig_brem = frac * sig_brem
    sig_pair = frac * sig_pair
    sig_pn = frac_pn * sig_pn

end subroutine int_xc_lep

subroutine int_alpha(energy, alpha_sig, alpha)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in):: energy
    real(dp), intent(in) :: alpha_sig(:)
    !f2py intent(out) :: alpha
    real(dp) alpha

    call interpol(energy, E_lep, alpha_sig, alpha)

end subroutine int_alpha

subroutine int_beta(energy, beta_arr, frac, frac_pn, tot)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in):: energy, frac, frac_pn
    real(dp), intent(in) :: beta_arr(:,:)
    real(dp) :: brem, pair, pn
    !f2py intent(out) :: tot
    real(dp) tot

    call interpol(energy, E_lep, beta_arr(:,1), brem)
    call interpol(energy, E_lep, beta_arr(:,2), pair)
    call interpol(energy, E_lep, beta_arr(:,3), pn)

    tot = (frac*brem)+(frac*pair)+(frac_pn*pn)

end subroutine int_beta

subroutine searchsorted(array, search_value, binarysearch)
    ! Given an array and a value, returns the index of the element that
    ! is closest to, but less than, the given value.
    ! Uses a binary search algorithm.
    ! "delta" is the tolerance used to determine if two values are equal
    ! if ( abs(x1 - x2) <= delta) then
    !    assume x1 = x2
    ! endif

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: array(:), search_value
!    real(dp), intent(in), optional :: delta
    !f2py intent(out) :: binarysearch
    integer binarysearch
    integer :: length, left, middle, right
    real(dp) :: d

    length = size(array)
!    if (present(delta) .eqv. .true.) then
!        d = delta
!    else
!        d = 1e-9
!    endif

    d = 1e-9
    left = 1
    right = length
    do
        if (left > right) then
            exit
        endif
        middle = nint((left+right) / 2.0_dp)
        if ( abs(array(middle) - search_value) <= d) then
            binarySearch = middle
            return
        else if (array(middle) > search_value) then
            right = middle - 1
        else
            left = middle + 1
        end if
    end do
    binarysearch = right

end subroutine searchsorted

end module interpolation

module transport

use constants
use geometry
use interpolation
implicit none
contains

subroutine random_no(r)

    implicit none
    !!f2py threadsafe
    !f2py intent(out) :: r
    real(dp) r

    call random_number(r)

!    r = 0.5 ! for debugging only!

end subroutine random_no

subroutine idecay(energy, distance, m_le, c_tau, decay)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: energy, distance, m_le, c_tau
    !f2py intent(out) :: decay
    integer decay
    real(dp) :: gamma_val, prob_decay, dy

    gamma_val = energy/m_le
    prob_decay = 1.0_dp - dexp(-distance/(gamma_val*c_tau))
    call random_no(dy)

    if (dy < prob_decay) then
        decay = 0
    else
        decay = 1
    end if

end subroutine idecay

subroutine em_cont_part(E_init, alpha_val, beta_val, x, m_le, E_fin)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: E_init, alpha_val, beta_val, x, m_le
    !f2py intent(out) :: E_fin
    real(dp) E_fin

    if (beta_val * x < 1e-6_dp) then
        E_fin = E_init * (1._dp-beta_val*x) - alpha_val*x
    else
        E_fin = E_init * dexp(-beta_val*x) - alpha_val/beta_val*(1-dexp(-beta_val*x))
    end if

    if (E_fin < 0) then
        E_fin = m_le
    end if

end subroutine em_cont_part

subroutine int_length_nu(energy, nu_xc, fac_nu, x_int) ! interaction lengths for neutrinos; in cm
!     int_length_nu = 1/(N_A * sigma) = 1/((1/cm^3)*cm) = cm

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: energy, fac_nu
    real(dp), intent(in) :: nu_xc(:,:)
    !f2py intent(out) :: x_int
    real(dp) x_int
    real(dp) :: sig_cc, sig_nc

    call int_xc_nu(energy, nu_xc, sig_cc, sig_nc) ! initialize CC & NC xc interpolations
    x_int = 1/(((sig_cc + sig_nc))*fac_nu) ! check the * or / fac_nu!

end subroutine int_length_nu

subroutine int_length_lep(energy, xc_arr, rho, m_le, c_tau, frac, frac_pn, x_int) ! interaction lengths for leptons; in cm

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: energy, rho, m_le, c_tau, frac, frac_pn
    real(dp), intent(in) :: xc_arr(:,:)
    !f2py intent(out) :: x_int
    real(dp) x_int
    real (8) :: sig_cc, sig_nc, gamma_val, sig_brem, sig_pair, sig_pn

    sig_cc = 0 ! placeholder for CC lepton interactions
    sig_nc = 0 ! placeholder for NC lepton interactions
    gamma_val = energy/m_le
    call int_xc_lep(energy, xc_arr, frac, frac_pn, sig_brem, sig_pair, sig_pn) ! initialize brem, pair & pn xc interpolations
    x_int = 1/((sig_brem + sig_pair + sig_pn + (1/(gamma_val*c_tau*rho)) + sig_cc + sig_nc))

end subroutine int_length_lep

subroutine interaction_type_nu(energy, nu_xc, fac_nu, int_type)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: energy, fac_nu
    real(dp), intent(in) :: nu_xc(:,:)
    !f2py intent(out) :: int_type
    integer int_type
    real(dp) :: sig_cc, sig_nc, x, int_nu, tot_frac, cc_frac

    call int_xc_nu(energy, nu_xc, sig_cc, sig_nc)
    call int_length_nu(energy, nu_xc, fac_nu, int_nu)

    tot_frac = 1/int_nu
    cc_frac = sig_cc/tot_frac

    call random_no(x)

    if (x <= cc_frac) then
        int_type = 0 ! CC
    else
        int_type = 1 ! NC
    end if

end subroutine interaction_type_nu

subroutine interaction_type_lep(energy, xc_arr, rho, m_le, c_tau, frac, frac_pn, int_type)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: energy, rho, m_le, c_tau, frac, frac_pn
    real(dp), intent(in) :: xc_arr(:,:)
    !f2py intent(out) :: int_type
    integer int_type
    real(dp) :: sig_cc, sig_nc, sig_brem, sig_pair, sig_pn, gamma_val, int_lep
    real(dp) :: tot_frac, decay_frac, cc_frac, nc_frac, brem_frac, pair_frac, pn_frac, x

    sig_cc = 0 ! placeholder for CC lepton interactions
    sig_nc = 0 ! placeholder for NC lepton interactions
    call int_xc_lep(energy, xc_arr, frac, frac_pn, sig_brem, sig_pair, sig_pn)

    gamma_val = energy/m_le
    call int_length_lep(energy, xc_arr, rho, m_le, c_tau, frac, frac_pn, int_lep)

    tot_frac = 1/int_lep
    decay_frac = (1/(gamma_val*c_tau*rho))/tot_frac
    cc_frac = sig_cc/tot_frac ! placeholder for CC lepton interactions
    nc_frac = sig_nc/tot_frac ! placeholder for NC lepton interactions
    brem_frac = sig_brem/tot_frac
    pair_frac = sig_pair/tot_frac
    pn_frac = sig_pn/tot_frac

    call random_no(x)
!    x = 0.33 # for debugging only!

    if (x < decay_frac) then
        int_type = 2 ! decay
    else if (decay_frac < x .and. x < decay_frac+brem_frac) then
        int_type = 3 ! brem
    else if (decay_frac+brem_frac < x .and. x < decay_frac+brem_frac+pair_frac) then
        int_type = 4 ! pair
    else if (decay_frac+brem_frac+pair_frac < x .and. x < decay_frac+brem_frac+pair_frac+pn_frac) then
        int_type = 5 ! pn
    else
        int_type = 6 ! lep_nc => cc/nc (both need CDFs) - placeholder for now
    end if

end subroutine interaction_type_lep

subroutine find_y(energy, ixc_arr, ip, y)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: energy
    real(dp), intent(in) :: ixc_arr(:,:,:)
    integer, intent(in) :: ip ! type of interaction
    !f2py intent(out) :: y
    real(dp) y
    real(dp) :: dy, dlv, search_arr(31) ! interpolating_arr is f(x)
!    real(dp) :: loc
    integer :: ip_id, energy_index

    call random_no(dy)
!    dy = 0.33 ! for debugging only!

    if (ip == 0 .or. ip == 1) then ! basically, for neutrinos
        call searchsorted(E_nu, energy, energy_index)

        if (ip == 0) then ! CC
            ip_id = 1
        else
            ip_id = 2
        end if

    else
        call searchsorted(E_lep, energy, energy_index)

        if (ip == 3) then ! brem
            ip_id = 1
        else if (ip == 4) then ! pair
            ip_id = 2
        else if (ip == 5) then ! pn
            ip_id = 3
        else ! lep_nc?
            ip_id = 4 ! shouldn't happen for now
        end if
    end if

!    print *,"energy_index=",energy_index
!    search_arr = ixc_arr(ip_id,energy_index,1:)
    search_arr = ixc_arr(:,energy_index,ip_id)
!    print *,search_arr
!    y = 1
    call interpol(dy, search_arr, v2, dlv) ! interpolate in v2
!    print *,"dlv=",dlv
    y = 10.**dlv

    if (y > 1._dp) then
        y = 1.0_dp
    end if

!    loc = -10.*dlog10(y)
!    print *,'loc=',loc

end subroutine find_y

subroutine get_frac(rho, frac, frac_pn)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: rho
    !f2py intent(out) :: frac, frac_pn
    real(dp) frac, frac_pn
    real(dp) :: f_rock

    if (rho > rho_rock .and. rho < rho_iron) then
        f_rock = (rho_iron - rho)/(rho_iron - rho)
        frac = 1.97_dp - 0.97_dp * f_rock
        frac_pn = 0.91_dp + 0.09_dp * f_rock
    else
        frac = 1._dp
        frac_pn = 1._dp
    end if

end subroutine get_frac

subroutine propagate_nu(e_init, nu_xc, nu_ixc, depth_max, fac_nu, part_type, d_travel, e_fin)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: nu_xc(:,:)
    real(dp), intent(in) :: e_init, nu_ixc(:,:,:), depth_max, fac_nu
    !f2py intent(out) :: d_travel, e_fin, part_type
    real(dp) d_travel, e_fin
    real(dp) :: e_nu, int_len, x, x_f, x_0, y, col_depth_total, dy
    integer :: int_type
    integer part_type

    part_type = 0 ! starting off as a neutrino
    col_depth_total = depth_max*1e5_dp
    e_nu = e_init
    e_fin = e_init
    x_0 = 0.0_dp ! starting depth in cm
    d_travel = depth_max ! added this in case there is a problem and needed for breaking out when E<1e3

    do while (e_nu > 1e3_dp)
        call random_no(dy)
        call int_length_nu(e_nu, nu_xc, fac_nu, int_len)
        x = -int_len*dlog(dy)
        x_f = x_0 + x

        if (x_f > col_depth_total) then
            return
        end if

        if (x_f > col_depth_total) then
            print *, 'This should never happen'
        end if

        x_0 = x_f
        call interaction_type_nu(e_nu, nu_xc, fac_nu, int_type)

!        print *,'int_type=',int_type

        if (part_type == 0 .and. int_type == 1) then
            part_type = 0
        else if (part_type == 0 .and. int_type == 0) then
            part_type = 1
        end if

        call find_y(e_nu, nu_ixc, int_type, y)

        e_fin = e_nu*(1-y)

        if (part_type == 1) then
            d_travel = x_0*1e-5_dp
            return
            exit ! whoops! forgot this!
        else
            e_nu = e_fin
            d_travel = x_0*1e-5_dp
        end if

        if (e_nu <= 1e3_dp) then
            return
        end if
    end do

    if (e_nu<1e3_dp .or. part_type == 1) then
        return
    end if

end subroutine propagate_nu

subroutine propagate_lep_water(e_init, xc_water, lep_ixc, alpha_water, beta_water, d_in, lepton, prop_type, part_id, d_fin, e_fin)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in):: xc_water(:,:), beta_water(:,:)
    real(dp), intent(in) :: lep_ixc(:,:,:)
    real(dp), intent(in) :: alpha_water(:)
    real(dp), intent(in) :: e_init, d_in
    integer, intent(in) :: lepton, prop_type

    !f2py intent(out) :: d_fin, e_fin, part_id
    real(dp) d_fin, e_fin
    integer part_id

    real(dp) :: e_min, cd_left, e_lep, x_0, m_le, c_tau, dy, int_len, frac, frac_pn
    real(dp) :: x, x_f, d_rem, alpha, beta, e_int, e_avg, y, d0, delta_d, delta_x
    integer :: int_type
    integer :: cnt, j_max, i

    e_min = 1e3_dp ! minimum tau energy, in GeV

    part_id = 1 ! start with tau that's obviously not decayed

    cd_left = d_in*1e5_dp ! how much to go, in cm.w.e
    e_lep = e_init
    e_fin = e_init ! in case the first interaction is too far
    x_0 = 0._dp ! haven't gone anywhere yet
    cnt = 0

    frac = 1._dp
    frac_pn = 1._dp

    if (lepton == 1) then
        m_le = 1.77682_dp ! m_tau in GeV
        c_tau = 8.703e-3_dp ! c*lifetime, in cm, for taus (taken from PDB)
    else
        m_le = 0.10565837550000001_dp ! m_mu in GeV
        c_tau = 6.586384e4_dp ! c*lifetime, in cm, for muons (taken from PDB 2020)
    end if

    if (prop_type == 1) then

!        do while (e_lep > e_min)

        do i = 1,1000
            if (e_lep <= e_min) then
                exit
            end if
!            print *,'e_lep=',e_lep
            cnt = cnt + 1
!            print *,"cnt=",cnt
            call random_no(dy)
            call int_length_lep(e_lep, xc_water, rho_water, m_le, c_tau, frac, frac_pn, int_len)
!            print *,'int_len=', int_len
            x = -int_len*dlog(dy) ! t is rho*L and has units g/cm^2
            x_f = x_0 + x ! going 1D; how far have we traveled here
            d_fin = x_f/1e5_dp ! make sure it is not past the old number, in km.w.e

            if (x_f >= cd_left) then ! already past maximum depth but still a tau; go to 30
                d_rem = cd_left - x_0
                call int_alpha(e_lep, alpha_water, alpha) ! changed 12/9/2020
                call int_beta(e_lep, beta_water, frac, frac_pn, beta) ! changed 12/9/2020
                e_fin = e_lep - (e_lep*beta + alpha)*d_rem
                d_fin = d_in

                if (e_fin > e_init) then ! sanity check
                    e_fin = e_init
                end if

                if (e_fin <= e_min) then ! tau has decayed
                    d_fin = d_in ! just in case; added 12/9/2020
                    e_fin = e_min
                    part_id = 2
                end if
!                print *,'x_f=',x_f
                return

            end if

            x_0 = x_f ! update x_0 and keep going
            call int_alpha(e_lep, alpha_water, alpha)
            call int_beta(e_lep, beta_water, frac, frac_pn, beta)
!            if (cnt == 2) then
!                print *,'e_lep=',e_lep
!                print *,'x_0=',x_0
!                print *,'alpha=', alpha
!                print *,'beta=',beta
!            end if
            e_int = e_lep - (e_lep*beta + alpha)*x ! find some intermediate energy to get reasonable values of energy between initial and final energy, a la MUSIC

            if (e_int <= e_min) then
                e_int = e_min
            end if

            e_avg = 10._dp**((dlog10(e_lep)+dlog10(e_int))/2._dp) ! does this work?; changed 12/9/2020

            call int_alpha(e_avg, alpha_water, alpha)
            call int_beta(e_avg, beta_water, frac, frac_pn, beta)

            call em_cont_part(e_lep, alpha, beta, x, m_le, e_int) ! get the continuous energy loss

            if (e_int <= e_min) then ! below minimum energy; go to 20
                ! 20 continue
                d_fin = d_in ! changed 12/9/2020
                e_fin = e_min
                part_id = 2 ! don't count this
                return
            end if

            call interaction_type_lep(e_int, xc_water, rho_water, m_le, c_tau, frac, frac_pn, int_type)

            if (int_type == 2) then ! tau has decayed
                part_id = 0
                e_fin = e_int
                d_fin = x_f/1e5_dp
                ! go to 50
                ! 50 continue
                return ! basically what d_fin does this return?
            end if

            ! tau didn't decay. Now how much energy does it have after interaction?

            call find_y(e_int, lep_ixc, int_type, y)

!            if (cnt == 1) then
!                print *,'int_type=',int_type
!                print *,'e_int=',e_int
!                print *,'y=',y
!!                print *,'x_0=',x_0
!!                print *,'alpha=', alpha
!!                print *,'beta=',beta
!            end if
            ! outgoing tau energy is old e_lep*(1-y)
            e_lep = e_int*(1._dp-y) ! this is the energy for the next interaction
            e_fin = e_lep
            ! go to 10
        end do

        ! Outside the while loop, e_lep has to be < e_min
        if (e_lep <= e_min) then ! only continuous energy loss; go to 20
            d_fin = d_in
            e_fin = e_min
            part_id = 2 ! don't count this
            return
        end if

    else ! continuous energy loss
!        print *,"continuous"
        d0 = 0._dp
        delta_d = step_size ! for now, not adaptive, distance into decay, cm
!        delta_d = 5000._dp ! test
        j_max = int(cd_left/(delta_d*rho_water)) ! we will get close to exit.
!        print *,'j_max=',j_max

        do cnt = 1, j_max+1
            if (e_lep < e_min) then
                exit
            end if

            d0 = d0 + delta_d ! where we are in xalong - use in rock to find rho
            delta_x = delta_d * rho_water ! distance goes into decay
            x_f = x_0 + delta_x
            ! does the particle decay over this distance?
            call idecay(e_lep, delta_d, m_le, c_tau, part_id)

!            if (cnt==1) then
!                print *,'part_id=',part_id
!                print *,'e_lep=',e_lep
!                print *,'delta_d=',delta_d
!                print *,'m_le=',m_le
!                print *,'c_tau=',c_tau
!            end if

            if (part_id == 0) then ! we are all done
                e_fin = e_lep
                d_fin = d_in
                return
            else ! find the new energy; assume alpha and beta are total values, not cut values
                call int_alpha(e_lep, alpha_water, alpha) ! changed 12/9/2020
                call int_beta(e_lep, beta_water, frac, frac_pn, beta) ! changed 12/9/2020

                e_fin = e_lep - (e_lep*beta + alpha)*delta_x
                d_fin = x_f/1e+5_dp
                x_0 = x_f
                e_lep = e_fin

                if (cnt >= j_max) then
                    exit ! if so, break out of while loop
                end if

            end if

        end do

        if (cnt >= j_max) then
            delta_x = cd_left - x_f
            if (delta_x > 0._dp) then !last little energy loss
                e_fin = e_lep - (e_lep*beta + alpha)*delta_x
                d_fin = d_in ! take care of that last little delta-x
                return
            else
                if (e_fin <= e_min) then
                    e_fin = e_min
                    d_fin = d_in
                    part_id = 2
                else if (e_fin > e_init) then
                    e_fin = e_init ! sanity check
                return
                end if

            end if

        end if

        ! Outside the while e_lep has to be < e_min
        if (e_lep <= e_min) then ! only continuous energy loss; go to 20
            d_fin = d_in
            e_fin = e_min
            part_id = 2 ! don't count this
            return
        end if
    end if
end subroutine propagate_lep_water

subroutine propagate_lep_rock(angle, e_init, xc_rock, lep_ixc, alpha_rock, beta_rock, d_entry,&
& d_in, xalong, cdalong, idepth, lepton, prop_type, part_id, d_fin, e_fin)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: xc_rock(:,:), beta_rock(:,:)
    real(dp), intent(in) :: lep_ixc(:,:,:)
    real(dp), intent(in) :: alpha_rock(:)
    real(dp), intent(in) :: xalong(:), cdalong(:)
    real(dp), intent(in) :: angle, e_init, d_entry, d_in
    integer, intent(in) :: idepth
    integer, intent(in) :: lepton, prop_type
    !f2py intent(out) :: d_fin, e_fin, part_id
    real(dp) d_fin, e_fin
    integer part_id

    real(dp) :: e_min, col_depth, d_max, e_lep, x_0, m_le, c_tau, x_interp, r, frac, frac_pn
    real (dp) :: rho, dy, int_len, x, x_f, alpha, beta, e_int, e_avg, y, d_0, delta_x, d_rem, delta_d
    integer :: cnt, j_max
    integer :: int_type

    e_min = 1e3_dp ! minimum tau energy, in GeV
    part_id = 1 ! start with tau
    col_depth = d_entry*1e5_dp ! how far in
    d_max = d_in*1e5_dp ! how much to go, in cm.w.e
    e_lep = e_init
    x_0 = 0._dp
    cnt = 0
    ! cd_left = d_in*1e5

    if (lepton == 1) then
        m_le = 1.77682_dp ! m_tau in GeV
        c_tau = 8.703e-3_dp ! c*lifetime, in cm, for taus (taken from PDB)
    else
        m_le = 0.10565837550000001_dp ! m_mu in GeV
        c_tau = 6.586384e4_dp ! c*lifetime, in cm, for muons (taken from PDB 2020)
    end if

    if (prop_type == 1) then ! stochastic energy loss

        do while (e_lep > e_min)
            cnt = cnt + 1
            call cd2distd(xalong, cdalong, col_depth, x_interp) ! find how far we are along the chord for given beta
            call densityatx(x_interp, angle, idepth, r, rho) ! find the density at x

            call get_frac(rho, frac, frac_pn)

!            rho = rho_rock ! USE FOR TESTING P_SURV FOR ROCK ONLY!
!            frac = 1._dp
!            frac_pn = 1._dp

            call random_no(dy)
            call int_length_lep(e_lep, xc_rock, rho, m_le, c_tau, frac, frac_pn, int_len)

            x = -int_len * dlog(dy)
            col_depth = col_depth + x ! update along trajectory, from the start of the chord
            x_f = x_0 + x ! going 1D
            d_fin = x_f/1e5_dp ! make sure it is not past the old number, in km.w.e

            if (x_f > d_max) then ! already past max depth
                ! go to 30
                d_rem = d_max - x_0
                call int_alpha(e_lep, alpha_rock, alpha)
                call int_beta(e_lep, beta_rock, frac, frac_pn, beta)
                e_fin = e_lep - (e_lep*beta + alpha)*d_rem
                d_fin = d_max/1e5_dp
                if (e_fin > e_init) then
                    e_fin=e_init
                end if
                return

            end if

            x_0 = x_f ! update x_0 and keep going
            call int_alpha(e_lep, alpha_rock, alpha)
            call int_beta(e_lep, beta_rock, frac, frac_pn, beta)
            e_int = e_lep - (e_lep*beta + alpha)*x ! find some intermediate energy to get reasonable values of energy between initial and final energy, a la MUSIC

            if (e_int <= e_min) then
                e_int = e_min
            end if

            e_avg = 10._dp**((dlog10(e_lep)+dlog10(e_int))/2) ! does this work?

            call int_alpha(e_avg, alpha_rock, alpha)
            call int_beta(e_avg, beta_rock, frac, frac_pn, beta)

            call em_cont_part(e_lep, alpha, beta, x, m_le, e_int) ! get the continuous energy

            if (e_int <= e_min) then ! is it below minimum energy now?
!                file.write("e_int < 1e3, so considering this to be a decay" + "\n")
                e_fin = e_int
                ! go to 20
                d_fin = d_max/1e5_dp ! in km.w.e
                e_fin = e_min
                part_id = 0
                return
            end if

            call interaction_type_lep(e_int, xc_rock, rho, m_le, c_tau, frac, frac_pn, int_type)

            if (int_type == 2) then ! tau has decayed
!                file.write(str("decayed; no. of stochastic interactions before decay = %d" % stoch_int) + "\n")
                part_id = 0
                e_fin = e_int
                ! go to 50
                ! 50 continue
                return
            end if

            ! tau didn't decay. Now how much energy does it have after interaction?
!            stoch_int += 1
            call find_y(e_int, lep_ixc, int_type, y)
!            file.write(str("interaction type = %s, y = %f" % (int_type,y)) + "\n")

            ! outgoing tau energy is old e_lep*(1-y)
            e_lep = e_int*(1-y) ! this is the energy for the next interaction
            e_fin = e_lep

        end do

        ! Outside the while loop, e_lep has to be < e_min
        if (e_lep <= e_min) then ! only continuous energy loss; go to 20
!            file.write("e_int < 1e3, so considering this to be a decay" + "\n")
            d_fin = d_max/1e5_dp
            e_fin = e_min
            part_id = 0 ! don't count this; decayed or no_count??? should be decayed
            return
        end if

    else ! continuous energy loss

        d_0 = 0._dp
        delta_x = step_size ! for now, not adaptive, distance into decay, cm; works for taus
!        delta_d = 5000._dp ! test

        call cd2distd(xalong, cdalong, col_depth, x_interp) ! find how far we are along the chord for given beta
        call densityatx(x_interp, angle, idepth, r, rho) ! find the density at x
        call get_frac(rho, frac, frac_pn)

!        rho = rho_rock ! FOR TESTING P_SURV ONLY!!
!        frac = 1._dp
!        frac_pn = 1._dp

        j_max = dint(d_max/(rho*delta_x))

        if (j_max == 0) then
            e_fin = e_init
            d_fin = d_in
            return
        end if

        do cnt = 1, j_max+1

            if (e_lep < e_min) then
                exit
            end if
!            cnt += 1

            call cd2distd(xalong, cdalong, col_depth, x_interp) ! find how far we are along the chord for given beta
            call densityatx(x_interp, angle, idepth, r, rho) ! find the density at x
            call get_frac(rho, frac, frac_pn)

            delta_d = delta_x * rho

            x_0 = x_0 + delta_x
            d_0 = d_0 + delta_d

            ! does the particle decay over this distance?
            call idecay(e_lep, delta_x, m_le, c_tau, part_id)

            if (part_id == 0) then ! we are all done
                e_fin = e_lep
                d_fin = d_0/1e5_dp
                return
            else ! find the new energy; assume alpha and beta are total values, not cut values
                call int_alpha(e_lep, alpha_rock, alpha) ! changed 12/9/2020
                call int_beta(e_lep, beta_rock, frac, frac_pn, beta) ! changed 12/9/2020
                e_lep = e_lep - (e_lep*beta + alpha)*delta_d
                d_fin = d_0/1e5_dp ! updating the d_final
                e_fin = e_lep
                col_depth = d_entry*1e5_dp + d_0 ! in order to update rho

                if (cnt >= j_max) then
                    exit ! if so, break out of while loop
                end if

            end if

        end do

        if (cnt >= j_max) then
            if (delta_x>0._dp) then ! last little energy loss
                d_fin = d_in ! take care of that last little delta-x
                return
            else ! this never runs
                if (e_fin <= e_min) then
                    e_fin = e_min
                    d_fin = d_in
                    part_id = 2
                else if (e_fin > e_init) then
                    e_fin = e_init ! sanity check
                end if

                return

            end if

        end if

        if (e_lep <= e_min) then ! only continuous energy loss; go to 20
            d_fin = d_in
            e_fin = e_min
            part_id = 2 ! don't count this
            return
        else if (cnt >= j_max) then
            d_fin = d_in ! take care of that last little delta-x; naaaah!
            return
        end if

    end if

end subroutine propagate_lep_rock

subroutine tau_thru_layers(angle, depth, d_water, depth_traj, e_lep_in, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,&
& alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth, lepton, prop_type, part_type, d_fin, e_fin)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: angle, depth, d_water, xc_water(:,:), xc_rock(:,:), lep_ixc_water(:,:,:)
    real(dp), intent(in) :: lep_ixc_rock(:,:,:), alpha_water(:), alpha_rock(:), beta_water(:,:), beta_rock(:,:)
    real(dp), intent(in) :: xalong(:), cdalong(:), e_lep_in, depth_traj
    integer, intent(in) :: idepth, lepton, prop_type
!    real(dp), intent(inout) :: e_lep_in, depth_traj ! see if this creates problems
    !f2py intent(out) :: d_fin, e_fin, part_type
    real(dp) d_fin, e_fin
    integer part_type

    real(dp) :: col_depth, rho, x, r, d_in, d_f, e_lep_ch, depth_traj_ch

    depth_traj_ch = depth_traj ! so depth_traj doesn't change
    d_fin = depth_traj_ch
    col_depth = depth_traj_ch*1e5_dp ! g/cm^2
    e_lep_ch = e_lep_in ! so e_lep_in doesn't change
    e_fin = e_lep_ch
    part_type = 1 ! tau going in

    if (e_lep_ch < 1e3_dp) then ! just in case
        part_type = 0
        d_fin = depth
        return
    end if

    if (angle <= 1.5_dp .or. depth-depth_traj_ch < d_water) then
        rho = rho_water ! water
    else
        call cd2distd(xalong, cdalong, col_depth, x)
        call densityatx(x, angle, idepth, r, rho)

        if (rho <= 0._dp) then ! round off error happening here; went too far
            print *,"col_depth = ", col_depth
            print *,"x = ", x
            print *,'rho is 0'
        end if
        if (rho <= 1.5_dp .and. r < 6365.0_dp) then
            print *,'rho too small!'
        end if
    end if

    if (rho > 1.5_dp) then ! we aren't in water yet
!    if (rho > rho_water) then ! we aren't in water yet
        d_in = depth - depth_traj_ch - d_water ! propagate this far in rock

        call propagate_lep_rock(angle, e_lep_ch, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, depth_traj_ch,&
        & d_in, xalong, cdalong, idepth, lepton, prop_type, part_type, d_f, e_fin)

        if (part_type == 1 .and. idepth /= 0) then ! still a tau; added .and. clause on 3/18

            e_lep_ch = e_fin
            d_in = d_water
            depth_traj_ch = depth_traj_ch + d_f ! now propagate through final layer of water
            d_fin = depth_traj_ch

            call propagate_lep_water(e_lep_ch, xc_water, lep_ixc_water, alpha_water, beta_water, d_in,&
            & lepton, prop_type, part_type, d_f, e_fin)

        else ! neutrino; or just make it out
            return
        end if

    else
        d_in = depth - depth_traj_ch

        call propagate_lep_water(e_lep_ch, xc_water, lep_ixc_water, alpha_water, beta_water, d_in, lepton,&
        & prop_type, part_type, d_f, e_fin)

        depth_traj_ch = depth_traj_ch + d_f
        d_fin = depth_traj_ch ! needed to add this since return was depth_traj here
        return
    end if

    if (part_type == 0) then ! tau decayed
        return
    end if

    return

end subroutine tau_thru_layers

!function fnu(y)
!
!    implicit none
!    real(dp) :: fnu
!    real(dp), intent(in) :: y
!!    real(dp), intent(out) :: fnu_val
!
!    fnu = y/3 * (5 - 3* y**2 + y**3) - y/3 * (1 - 3 * y**2 + 2 * y**3)
!end function fnu

subroutine distnu(r, ithird, dist_val)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: r
    integer, intent(in) :: ithird
    !f2py intent(out) :: dist_val
    real(dp) dist_val
    real(dp) :: fnu, y, fm, ff, y0, y1

    fnu(y) = y/3._dp * (5._dp - 3._dp* y**2 + y**3) - y/3._dp * (1._dp - 3._dp * y**2 + 2._dp * y**3)

    if (ithird /= 1) then
        fm = 1._dp ! max value of distribution
        ff = r*fm
        y0 = 0._dp
        y1 = fm

        do while (abs(y1-y0) > 0.1e-2_dp)
            y = (y0+y1)/2._dp
            if (fnu(y) < ff) then
                y0 = y
            else
                y1 = y
            end if
        dist_val = (y0+y1)/2._dp
        end do

    else
        dist_val = 1._dp/3._dp
    end if

end subroutine distnu

subroutine regen(angle, e_lep, depth, d_water, d_lep, nu_xc, nu_ixc, ithird, xc_water, xc_rock,&
& ixc_water, ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth,&
& lepton, fac_nu, prop_type, part_type, d_exit, e_fin)

    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: angle, e_lep, depth, d_water, nu_xc(:,:), nu_ixc(:,:,:), d_lep
    real(dp), intent(in) :: xc_water(:,:), xc_rock(:,:), ixc_water(:,:,:), ixc_rock(:,:,:)
    real(dp), intent(in) :: alpha_water(:), alpha_rock(:), beta_water(:,:), beta_rock(:,:)
    real(dp), intent(in) :: xalong(:), cdalong(:), fac_nu
    integer, intent(in) :: idepth, lepton, prop_type, ithird
!    real(dp), intent(inout) :: d_lep ! see if this creates problems
    !f2py intent(out) :: d_exit, e_fin, part_type
    real(dp) d_exit, e_fin
    integer part_type

    real(dp) :: r, frac, e_nu, d_left, dtr, etau2, d_lep_ch
    integer :: int_part

    call random_no(r)
    call distnu(r, ithird, frac)
    e_nu = frac * e_lep

    d_lep_ch = d_lep
    d_left = depth-d_lep_ch ! this is how far the neutrino can go
    e_fin = e_nu ! in case we need to exit
    part_type = 3 ! in case we need to exit; change later to string (HLS = 2)

    if (d_left <= 0._dp) then ! past the point of interactions allowed
        d_exit = depth
        ! go to 60
        ! 60 continue
        return
    end if

    d_exit = d_lep_ch ! we are starting this far into the Earth with a neutrino
    int_part = 0 ! starting with a neutrino with energy e_nu; change later to string

    ! tnu follows NC to the end, or gives results if CC interactions
    call propagate_nu(e_nu, nu_xc, nu_ixc, d_left, fac_nu, int_part, dtr, etau2) ! does the neutrino interact?

    if (int_part /= 1) then ! neutrinos at the end
        d_exit = depth
        part_type = 0 ! (HLS = 0); changed 22/12/2020
        e_fin = etau2 ! final neutrino energy
        ! go to 60; all done
        ! 60 continue
        return
    end if

    ! otherwise we have a tau
    d_lep_ch = d_lep_ch + dtr
    d_left = d_left - dtr

    if (d_left <= 0._dp) then
        d_exit = depth
        e_fin = etau2
        part_type = 0 ! went too far to make a tau, so don't count (HLS = 0); changed 22/12/2020
        ! go to 60; no, still a neutrino
        ! 60 continue
        return
    end if

    ! we have a tau with room to travel for tauthrulayers

    call tau_thru_layers(angle, depth, d_water, d_lep_ch, etau2, xc_water, xc_rock, ixc_water, ixc_rock,&
    & alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth, lepton, prop_type,&
    & part_type, d_exit, e_fin)

    return

end subroutine regen

subroutine p_exit_w(xc_water, lep_ixc_water, alpha_water, beta_water)
    implicit none
    !!f2py threadsafe
    real(dp), intent(in) :: xc_water(:,:), lep_ixc_water(:,:,:), alpha_water(:), beta_water(:,:)
    real(dp) :: energy(104) = (/1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07,&
       & 1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07,&
       & 1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07, 1.e+07,&
       & 1.e+07, 1.e+07, 1.e+08, 1.e+08, 1.e+08, 1.e+08, 1.e+08, 1.e+08,&
       & 1.e+08, 1.e+08, 1.e+08, 1.e+08, 1.e+08, 1.e+08, 1.e+08, 1.e+08,&
       & 1.e+08, 1.e+08, 1.e+08, 1.e+08, 1.e+08, 1.e+08, 1.e+08, 1.e+08,&
       & 1.e+08, 1.e+08, 1.e+08, 1.e+08, 1.e+09, 1.e+09, 1.e+09, 1.e+09,&
       & 1.e+09, 1.e+09, 1.e+09, 1.e+09, 1.e+09, 1.e+09, 1.e+09, 1.e+09,&
       & 1.e+09, 1.e+09, 1.e+09, 1.e+09, 1.e+09, 1.e+09, 1.e+09, 1.e+09,&
       & 1.e+09, 1.e+09, 1.e+09, 1.e+09, 1.e+09, 1.e+09, 1.e+10, 1.e+10,&
       & 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10,&
       & 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10,&
       & 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10, 1.e+10/)

    real(dp) :: dm(104) = (/0.      ,  0.096604,  0.19321 ,  0.28981 ,  0.38642 ,  0.48302 ,&
        & 0.57962 ,  0.67623 ,  0.77283 ,  0.86943 ,  0.96604 ,  1.0626  ,&
        & 1.1592  ,  1.2558  ,  1.3525  ,  1.4491  ,  1.5457  ,  1.6423  ,&
        & 1.7389  ,  1.8355  ,  1.9321  ,  2.0287  ,  2.1253  ,  2.2219  ,&
        & 2.3185  ,  2.4151  ,  0.      ,  0.96604 ,  1.9321  ,  2.8981  ,&
        & 3.8642  ,  4.8302  ,  5.7962  ,  6.7623  ,  7.7283  ,  8.6943  ,&
        & 9.6604  , 10.626   , 11.592   , 12.558   , 13.525   , 14.491   ,&
       & 15.457   , 16.423   , 17.389   , 18.355   , 19.321   , 20.287   ,&
       & 21.253   , 22.219   , 23.185   , 24.151   ,  0.      ,  1.9623  ,&
        & 3.9245  ,  5.8868  ,  7.8491  ,  9.8113  , 11.774   , 13.736   ,&
       & 15.698   , 17.66    , 19.623   , 21.585   , 23.547   , 25.509   ,&
       & 27.472   , 29.434   , 31.396   , 33.358   , 35.321   , 37.283   ,&
       & 39.245   , 41.208   , 43.17    , 45.132   , 47.094   , 49.057   ,&
        & 0.      ,  3.0189  ,  6.0377  ,  9.0566  , 12.075   , 15.094   ,&
       & 18.113   , 21.132   , 24.151   , 27.17    , 30.189   , 33.208   ,&
       & 36.226   , 39.245   , 42.264   , 45.283   , 48.302   , 51.321   ,&
       & 54.34    , 57.358   , 60.377   , 63.396   , 66.415   , 69.434   ,&
       & 72.453   , 75.472/)
    real(dp) :: surv, e_fin, df, imax, ic, energy_single, dm_single
    integer :: j, k
    integer :: p_id

!    energy_single = 1e10
!    dm_single = 50
!    surv = 0.0
!    ic = 0._dp
!    call propagate_lep_water(energy_single, xc_water, lep_ixc_water, alpha_water,&
!    & beta_water, dm_single, 1, 1, p_id, df, e_fin)
!    print *,p_id,df,e_fin

    imax = 1e4_dp
    do j = 1, 104
        surv = 0.0_dp
        ic = 0._dp
        do k = 1, int(imax)
            call propagate_lep_water(energy(j), xc_water, lep_ixc_water, alpha_water, beta_water, dm(j), 1, 1, p_id, df, e_fin)
            if (p_id == 1 .and. e_fin > 50._dp) then
                ic = ic + 1._dp
            end if
        end do
        surv = ic/imax
        print *,energy(j), dm(j), surv, df
    end do

end subroutine p_exit_w

end module transport

module run

use constants
use geometry
use interpolation
use transport
implicit none
contains

subroutine single_stat(energy, angle, nu_xc, nu_ixc, depth, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,&
& alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu, prop_type,&
& u, no_regen_tot, regen_tot)

    implicit none
    !!f2py threadsafe

    real(dp), intent(in) :: energy, angle, nu_xc(:,:), nu_ixc(:,:,:), dwater, xc_water(:,:), xc_rock(:,:)
    real(dp), intent(in) :: lep_ixc_water(:,:,:), lep_ixc_rock(:,:,:), alpha_water(:), alpha_rock(:)
    real(dp), intent(in) :: depth, depthE, beta_water(:,:), beta_rock(:,:), xalong(:), cdalong(:), fac_nu
    integer, intent(in) :: ithird, idepth, lepton, prop_type
    integer(dp), intent(in) :: u

    !f2py intent(inout) :: no_regen_tot, regen_tot
    integer:: no_regen_tot, regen_tot

    real(dp) :: depth0, dtr, ef, etauin, dfinal, etauf, dleft, dtau2, ef2
    integer :: regen_cnt, ip, ipp, ipp3

    depth0 = 0.0_dp ! start with this each time

    ! 80 continue

    ! tnu goes until neutrino either goes to dtot, or converts to a tau
!        print *,'depth=',depth
    call propagate_nu(energy, nu_xc, nu_ixc, depth, fac_nu, ip, dtr, ef)
!        print *,"After propagate_nu..."
!        print *,ip, dtr, ef

    ! how far did the neutrino go? dtr is how far traveled

    depth0 = depth0 + dtr ! how far is the neutrino on trajectory?

    dleft = depth - depth0 ! how far is left for the neutrino to travel?

    if (ip == 0) then ! still a neutrino at the end of the road
        ! go to 10
        return ! break outside stat; continue is correct here
    end if

    ! continue here: we have a tau

    regen_cnt = 1 ! tau out after first interaction


    etauin = ef
    ! still need to propagate the tau, column depth to go


    call tau_thru_layers(angle, depth, dwater, depth0, etauin, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water,&
    & alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth, lepton, prop_type, ipp, dfinal, etauf)
!        print *,"After tau_thru_layers..."
!        print *,ipp, dfinal, etauf

    dleft = depth-dfinal

    if (ipp == 1 .and. dleft <= 0.0_dp) then ! a tau has emerged through column depth
        no_regen_tot = no_regen_tot + 1
        regen_tot = regen_tot + 1 ! update the regen tau array once
        !e_out.append(etauf)
        write(u, *) etauf
        ! go to 10; we are done with the loop
        return ! break outside stat; continue is correct here
    end if

    ! 11 continue; beginning of regeneration loop
    ! must be a neutrino. Is there still column depth to propagate?

    ipp3 = 99 ! dummy value
    do while (dfinal < depthE .and. ipp3 /= 1 .and. regen_cnt <= 6) ! tau has decayed before the end

        etauin = etauf ! regen finds neutrino energy


        call regen(angle, etauin, depth, dwater, dfinal, nu_xc, nu_ixc, ithird, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,&
        & alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth, lepton, fac_nu, prop_type, ipp3, dtau2, ef2)

!            print *,"After regen..."
!            print *,ipp3, dtau2, ef2

        regen_cnt = regen_cnt + 1

        if (ipp3 == 1) then ! then we are back to a tau at the end of the road
            regen_tot = regen_tot + 1
            ! e_out.append(ef2)
            write(u, *) ef2
            ! go to 10; we are done with the loop
            return ! need to check if this breaks out of stat loop or not. Yes??
        end if

        if (regen_cnt > 6) then ! 6 rounds of regeneration
            return ! only if regen > 6, break and go to run_stat for next iteration
        end if

        etauf = ef2
        dfinal = dtau2 ! go to 11

    end do

end subroutine single_stat


subroutine run_stat_single(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,&
& alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu, stats, prop_type,&
& no_regen_tot, regen_tot)

    implicit none
    !!f2py threadsafe

    real(dp), intent(in) :: energy, angle, nu_xc(:,:), nu_ixc(:,:,:), depthE, dwater, xc_water(:,:), xc_rock(:,:)
    real(dp), intent(in) :: lep_ixc_water(:,:,:), lep_ixc_rock(:,:,:), alpha_water(:), alpha_rock(:)
    real(dp), intent(in) :: beta_water(:,:), beta_rock(:,:), xalong(:), cdalong(:), fac_nu
    integer, intent(in) :: ithird, idepth, lepton, prop_type
    integer(kind=8), intent(in) :: stats

    !f2py intent(out) :: no_regen_tot, regen_tot
    integer :: no_regen_tot, regen_tot

    real(dp) :: depth, depth0, dtr, ef, etauin, dfinal, etauf, dleft, dtau2, ef2
    integer(kind=8):: regen_cnt, i
    integer :: ip, ipp, ipp3
!    integer :: regen_cnt, i
    integer(kind=8) :: u
    character(25) filename

!    write(6,*)'stat_val=',stats

    if (angle < 10._dp) then
        write(filename,'(a,es8.2,a,F4.2)') 'e_out_',energy,'_',angle ! filename is e_out_energy_angle
    else
        write(filename,'(a,es8.2,a,F5.2)') 'e_out_',energy,'_',angle ! filename is e_out_energy_angle
    end if

    open(newunit=u, file=trim(filename), status="replace")
!    call omp_set_nested(.false.)

    depth = depthE
    regen_cnt = 0
    no_regen_tot = 0
    regen_tot = 0
!    e_out = 0.0_dp ! initialize e_out
!$OMP PARALLEL DO
    do i = 1, stats
        call single_stat(energy, angle, nu_xc, nu_ixc, depth, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,&
        & alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu, prop_type, &
        & u, no_regen_tot, regen_tot)
    end do
!$OMP END PARALLEL DO
    close(u)
    return
end subroutine run_stat_single


!subroutine init_main(E_prop, angles, nu_xc, xc_water, xc_rock, alpha_water, alpha_rock, beta_water, beta_rock, nu_ixc, lep_ixc_water, lep_ixc_rock, idepth, lepton, fac_nu, stats, prop_type)
!
!    implicit none
!    !!f2py threadsafe
!
!    real(dp), intent(in) :: E_prop(:), angles(:), nu_xc(:,:), nu_ixc(:,:,:), xc_water(:,:), xc_rock(:,:)
!    real(dp), intent(in) :: lep_ixc_water(:,:,:), lep_ixc_rock(:,:,:), alpha_water(:), alpha_rock(:)
!    real(dp), intent(in) :: beta_water(:,:), beta_rock(:,:), xalong(:), cdalong(:), fac_nu
!    integer, intent(in) :: idepth, lepton, stats, prop_type
!
!    real(dp) ::
!    integer :: ithird
!
!    ithird = 0
!
!end subroutine init_main

end module run

!module data
!use constants
!!use h5fortran
!implicit none
!contains
!
!subroutine get_col_traj()
!
!    implicit none
!
!    integer :: idepth
!    character(16) :: filename
!    character(35) :: group
!
!    filename = 'lookup_tables.h5'
!    idepth = 4
!    write(group,'(a,i0,a)') 'Earth/traj_',idepth',/Column_Trajectories'
!    print *,filename
!    print *,group
!
!end subroutine get_col_traj
!
!end module data

module test
use constants
use interpolation
implicit none
contains

!subroutine multidim_test(ixc_arr, a)
!
!    implicit none
!    real(dp), intent(in) :: ixc_arr(:,:,:)
!    real(dp), intent(out):: a(121)
!
!    a =  xc_arr(:,3)
!end subroutine multidim_test

subroutine test_y(energy, ixc_arr, ip, search_arr)

    implicit none

    real(dp), intent(in) :: energy
    real(dp), intent(in) :: ixc_arr(:,:,:)
    integer, intent(in) :: ip ! type of interaction
    real(dp), intent(out) :: search_arr(31)
    integer :: ip_id, energy_index

    if (ip == 0 .or. ip == 1) then ! basically, for neutrinos
        call searchsorted(E_nu, energy, energy_index)

        if (ip == 0) then ! CC
            ip_id = 1
        else
            ip_id = 2
        end if

    else
        call searchsorted(E_lep, energy, energy_index)

        if (ip == 3) then ! brem
            ip_id = 1
        else if (ip == 4) then ! pair
            ip_id = 2
        else if (ip == 5) then ! pn
            ip_id = 3
        else ! lep_nc?
            ip_id = 4 ! shouldn't happen for now
        end if
    end if

!    print *,"energy_index=",energy_index
!    search_arr = ixc_arr(ip_id,energy_index,1:)
    search_arr = ixc_arr(:,energy_index,ip_id)

end subroutine test_y

subroutine test_write()
    implicit none

    integer(kind=16) :: u
    integer :: i

    open(newunit=u, file="log.txt", position="append", status="old")

    do i = 1,10
        write(u, *) i
    end do

    close(u)

end subroutine test_write

subroutine test_format()

    implicit none
    character(50) :: filename

!    write(filename,'(a,ES6.0E2,a,F3.1)') 'e_out_',energy,'_',angle ! filename is e_out_energy_angle
    write(filename,'(a,es8.2,a,F3.1)') 'e_out_',10.**10,'_',10 ! filename is e_out_energy_angle

    print *, filename
end subroutine test_format

end module test

module OTmod
  !$ use omp_lib
  implicit none

  public :: get_threads

contains

  function get_threads() result(nt)
    integer :: nt

    nt = 0
    !$ nt = omp_get_max_threads()
    print *, nt
  end function get_threads

end module OTmod