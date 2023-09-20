module constants
  
  implicit none
  integer, parameter :: dp=kind(0.d0)
  real(dp), parameter :: pi = 3.1415927_dp
  
  real(dp), parameter :: N_A = 6.0221409e+23_dp ! 1/mole
  real(dp), parameter :: rho_water = 1.02_dp ! g/cm^3
  real(dp), parameter :: rho_rock = 2.65_dp ! g/cm^3
  real(dp), parameter :: rho_iron = 7.87_dp ! g/cm^3
  real(dp), parameter :: R_earth = 6371.0_dp ! radius of the Earth, in km
  real(dp), parameter :: step_size = 4500.0_dp ! step size for continuous energy loss, in cm
  real(dp), parameter :: E_nu(91) = (/1.00000000e+03, 1.25892541e+03, 1.58489319e+03, 1.99526231e+03,& ! bins for neutrino energies
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
  real(dp), parameter :: E_lep(121) = (/1.00000000e+00, 1.25892541e+00, 1.58489319e+00, 1.99526231e+00,& ! bins for charged lepton energies
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
  ! real(dp), parameter, dimension(31) :: v2 = (/0.,-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1. ,&
  !        & -1.1, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.8, -1.9, -2. , -2.1,&
!        & -2.2, -2.3, -2.4, -2.5, -2.6, -2.7, -2.8, -2.9, -3./) ! x (interpolating arr for find_y)
  real(dp), parameter, dimension(31) :: yvals = (/1., 0.79432823, 0.63095734, 0.50118723, 0.39810717,0.31622777,&
       & 0.25118864, 0.19952623, 0.15848932, 0.12589254, 0.1, 0.07943282, 0.06309573, 0.05011872, 0.03981072,&
       & 0.03162278, 0.02511886, 0.01995262, 0.01584893, 0.01258925, 0.01, 0.00794328, 0.00630957, 0.00501187,&
       & 0.00398107,0.00316228, 0.00251189, 0.00199526, 0.00158489, 0.00125893, 0.001/) ! x (interpolating arr for find_y)
  
  ! arrays needed for interplation of polarization values
  real(dp), parameter, dimension(99) :: ypol = (/ .01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,&
       & 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 , 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27,&
       & 0.28, 0.29, 0.3 , 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4 , 0.41, 0.42, 0.43, 0.44,&
       & 0.45, 0.46, 0.47, 0.48, 0.49, 0.5 , 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6 , 0.61,&
       & 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7 , 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78,&
       & 0.79, 0.8 , 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9 , 0.91, 0.92, 0.93, 0.94, 0.95,&
       & 0.96, 0.97, 0.98, 0.99  /) 
  real(dp), parameter, dimension(99) :: Pcthp = (/9.99987641e-01,  9.99942103e-01,  9.99855593e-01,  9.99722311e-01,&
       & 9.99537246e-01,  9.99295721e-01,  9.98993123e-01,  9.98624698e-01,&
       & 9.98185400e-01,  9.97669756e-01,  9.97071768e-01,  9.96384828e-01,&
       & 9.95601646e-01,  9.94714200e-01,  9.93713699e-01,  9.92590548e-01,&
       & 9.91334341e-01,  9.89933842e-01,  9.88376993e-01,  9.86650920e-01,&
       & 9.84741945e-01,  9.82635613e-01,  9.80316721e-01,  9.77769353e-01,&
       & 9.74976927e-01,  9.71922248e-01,  9.68587574e-01,  9.64954681e-01,&
       & 9.61004952e-01,  9.56719464e-01,  9.52079095e-01,  9.47064635e-01,&
       & 9.41656911e-01,  9.35836925e-01,  9.29585993e-01,  9.22885909e-01,&
       & 9.15719096e-01,  9.08068786e-01,  8.99919190e-01,  8.91255675e-01,&
       & 8.82064945e-01,  8.72335218e-01,  8.62056401e-01,  8.51220253e-01,&
       & 8.39820542e-01,  8.27853191e-01,  8.15316401e-01,  8.02210765e-01,&
       & 7.88539349e-01,  7.74307762e-01,  7.59524188e-01,  7.44199402e-01,&
       & 7.28346753e-01,  7.11982116e-01,  6.95123820e-01,  6.77792550e-01,&
       & 6.60011214e-01,  6.41804799e-01,  6.23200191e-01,  6.04225984e-01,&
       & 5.84912275e-01,  5.65290436e-01,  5.45392883e-01,  5.25252844e-01,&
       & 5.04904116e-01,  4.84380825e-01,  4.63717198e-01,  4.42947334e-01,&
       & 4.22104990e-01,  4.01223381e-01,  3.80334990e-01,  3.59471402e-01,&
       & 3.38663150e-01,  3.17939580e-01,  2.97328741e-01,  2.76857286e-01,&
       & 2.56550401e-01,  2.36431743e-01,  2.16523401e-01,  1.96845879e-01,&
       & 1.77418081e-01,  1.58257324e-01,  1.39379360e-01,  1.20798406e-01,&
       & 1.02527192e-01,  8.45770141e-02,  6.69578032e-02,  4.96781973e-02,&
       & 3.27456283e-02,  1.61664186e-02, -5.41065269e-05, -1.59114852e-02,&
       & -3.14019518e-02, -4.65222127e-02, -6.12691149e-02, -7.56390959e-02,&
       & -8.96271198e-02, -1.03224155e-01, -1.16408787e-01 /)
  real(dp), parameter, dimension(99) :: P = (/0.99998943, 0.99995156, 0.99988105, 0.99977423, 0.99962816,&
       & 0.99944032, 0.99920841, 0.99893022, 0.9986035 , 0.99822587,&
       & 0.99779476, 0.99730737, 0.99676061, 0.9961511 , 0.99547513,&
       & 0.99472867, 0.99390739, 0.99300667, 0.9920216 , 0.99094703,&
       & 0.98977762, 0.98850782, 0.987132  , 0.98564442, 0.98403932,&
       & 0.982311  , 0.98045382, 0.97846234, 0.97633133, 0.9740559 ,&
       & 0.97163153, 0.9690542 , 0.96632043, 0.96342741, 0.96037307,&
       & 0.95715618, 0.95377643, 0.95023452, 0.94653223, 0.94267251,&
       & 0.93865957, 0.93449888, 0.93019725, 0.92576287, 0.92120527,&
       & 0.91653539, 0.91176549, 0.90690911, 0.90198105, 0.89699721,&
       & 0.89197452, 0.88693079, 0.88188454, 0.87685483, 0.87186109,&
       & 0.86692288, 0.86205967, 0.85729067, 0.85263456, 0.84810926,&
       & 0.84373177, 0.83951795, 0.83548232, 0.83163791, 0.82799614,&
       & 0.82456669, 0.82135745, 0.81837441, 0.81562174, 0.81310168,&
       & 0.81081468, 0.80875941, 0.80693288, 0.80533053, 0.80394637,&
       & 0.80277309, 0.80180226, 0.80102446, 0.80042943, 0.80000624,&
       & 0.79974345, 0.79962922, 0.7996515 , 0.79979813, 0.80005696,&
       & 0.80041596, 0.8008633 , 0.80138744, 0.80197721, 0.80262182,&
       & 0.80331098, 0.80403482, 0.80478404, 0.80554978, 0.80632371,&
       & 0.80709795, 0.80786499, 0.80861753, 0.80934773 /)

end module constants

module geometry
  
  use constants
  implicit none
contains
  
  subroutine PREMdensity(Rin, idepth, edens)
    !! Calculates the density at a radius inside the Earth according the PREM model.
    ! hep-ph/9512364v1 eq. 25
    
    implicit none
    
    integer, intent(in) :: idepth
    !! Depth of water layer, in km.
    real(dp), intent(in) :: Rin
    !! Point at which density is to be calculated.
    
    real(dp), intent(out) :: edens
    !! Density at Rin, in g/cm^3.

    real(dp) :: x,y
    real(dp), dimension(10) :: Rlay
    
    Rlay = (/1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0, 6346.6, 6356.0, 6368.0, 6371.0/) ! # PREM layers based on R_earth

    Rlay(9) = 6368.0_dp + (3.0_dp-dble(idepth)) ! sets the last layer as water layer
   
    x = Rin
    y = x/R_earth
    
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
  end subroutine PREMdensity

  subroutine densityatx(x, beta_deg, idepth, r, rho_at_x)
    !! Calculates the density at a distance x, for a given Earth emergence angle.
    ! 1905.13223v2 fig. 2 for chord length
    ! x is the distance at which we need to find the radial distance from the
    ! center of the Earth in order to find the density there.
    implicit none
    
    integer, intent(in) :: idepth
    !! Depth of water layer, in km.
    real(dp), intent(in) :: x
    !! Distance along the chord of the trajectory, in km.
    real(dp), intent(in) :: beta_deg
    !! Earth emergence angle, in degrees.
    
    real(dp), intent(out) :: r
    !! Radial distance from the center of the Earth to x, in km.
    real(dp), intent(out) :: rho_at_x
    !! Density at x, in g/cm^3
    
    real(dp) :: chord_length, r2
    
    chord_length = 2*R_earth*dsin(beta_deg*(pi/180.0_dp)) ! 2 R_E sin(beta)
    r2 = (chord_length-x)**2 + R_earth**2 - 2*R_earth*(chord_length-x)*dsin(beta_deg*(pi/180.0_dp)) ! just the law of cosines
    
    if (beta_deg < 5.0_dp) then
       r = R_earth*(1.0_dp + 0.5_dp*(x**2-chord_length*x)/R_earth**2)
    else
       r = dsqrt(r2)
    end if
    
    call PREMdensity(r, idepth, rho_at_x)
end subroutine

end module geometry

module interpolation
  use constants
  implicit none
contains
  
  subroutine interp_linear_internal(x, y, xout, yout)
    !! Interpolates between x(1), x(2), y(1) & y(2) for a given x value
    
    implicit none
    
    real(dp), intent(in) :: x(2)
    !! 2 element x array for interpolation
    real(dp), intent(in) :: y(2)
    !! 2 element y array for interpolation
    real(dp), intent(in) :: xout
    !! Point of interpolation
    
    real(dp), intent(out) :: yout
    !! Interpolated value
    
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
    !! Interpolates y from ordered x to ordered xout positions
    
    implicit none
    
    real(dp), intent(in), dimension(:) :: x
    !! x array for interpolation
    real(dp), intent(in), dimension(:) :: y
    !! y array for interpolation
    real(dp), intent(in) :: xout
    !! Point of interpolation
    
    real(dp), intent(out) :: yout
    !! Interpolated value
    
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
    !! Custom interpolation subroutine suited for all following subroutines
    
    implicit none
    
    real(dp), intent(in) :: x(:)
    !! x array for interpolation
    real(dp), intent(in) :: y(:)
    !! x array for interpolation
    real(dp), intent(in) :: x_val
    !! Point of interpolation
    
    real(dp), intent(out) :: interp_val
    !! Interpolated value
    
    call interp_linear_pt(x, y, x_val, interp_val)
  end subroutine interpol
  
  subroutine cd2distd(xalong, cdalong, col_depth, out_val)
    !! Interpolate between column depth & distance in water.
    
    implicit none
    
    real(dp), intent(in) :: xalong(:)
    !! Array of distance in water, in km.
    real(dp), intent(in) :: cdalong(:)
    !! Array of column depth at xalong, in g/cm^2.
    real(dp), intent(in) :: col_depth
    !! Column depth to interpolate at, in g/cm^2.
    
    real(dp), intent(out) :: out_val
    !! Interpolated distance in water, in km.
    
    if (col_depth < minval(cdalong)) then
       out_val = (col_depth/cdalong(1))*xalong(1)
    else if (col_depth > maxval(cdalong)) then
       out_val = maxval(xalong)
    else
       call interpol(col_depth, cdalong, xalong, out_val)
    end if
    
  end subroutine cd2distd
  
  subroutine get_rho_frac(rho, frac, frac_pn)
    !! Calculate the correction/scaling fraction for material density between rock & iron.
    
    implicit none

    real(dp), intent(in) :: rho
    !! Density of material, in g/cm^3.
    
    real(dp), intent(out) :: frac
    !! Scaling factor for density change between rock & iron (for Bremmstrahlung & pair production).
    real(dp), intent(out) :: frac_pn
    !! Scaling factor for density change between rock & iron (for photonuclear).
    
    real(dp) :: f_rock
    
    if (rho > rho_rock .and. rho < rho_iron) then
       f_rock = (rho_iron - rho)/(rho_iron - rho)
       frac = 1.97_dp - 0.97_dp * f_rock
       frac_pn = 0.91_dp + 0.09_dp * f_rock
    else ! for rho <= rho_water or rho>=iron (that shouldn't happen!) 
       frac = 1._dp
       frac_pn = 1._dp
     end if
     
   end subroutine get_rho_frac
   
   subroutine int_xc_nu(energy, nu_xc, fac_nu, sig_cc, sig_nc)
     !! Interpolate between neutrino energy & cross-section values.
     
     implicit none
     
     real(dp), intent(in) :: nu_xc(:,:)
     !! 2D array of neutrino cross-section values.
     real(dp), intent(in) :: energy
     !! Energy value to interpolate at, in GeV.
     real(dp), intent(in) :: fac_nu
     !! Rescaling factor for SM neutrino cross-sections.
     
     real(dp), intent(out) :: sig_cc, sig_nc
     !! Interpolated CC & NC cross-section values, in cm^2.
     
     call interpol(energy, E_nu, nu_xc(:,1), sig_cc)
     call interpol(energy, E_nu, nu_xc(:,2), sig_nc)

     sig_cc = fac_nu * sig_cc
     sig_nc = fac_nu * sig_nc
     
   end subroutine int_xc_nu
   
   subroutine int_xc_lep(energy, xc_arr, rho, sig_brem, sig_pair, sig_pn)
     !! Interpolate between lepton energy & cross-section values.
     
     implicit none
     
     real(dp), intent(in) :: xc_arr(:,:)
     !! 2D array of N_A/A*charged lepton cross-section values, in cm^2/g.
     real(dp), intent(in) :: energy
     !! Energy value to interpolate at, in GeV.
     real(dp), intent(in) :: rho
     !! Desnity of the material, in g/cm^3.
     
     real(dp), intent(out) :: sig_brem, sig_pair, sig_pn
     !! Interpolated cross-section values of charged lepton for energy losses, in cm^2/g.

     real(dp) :: frac, frac_pn
     
     call get_rho_frac(rho, frac, frac_pn)
     
     ! Note: The lookup tables already have N_A multiplied by lep_xc!
     call interpol(energy, E_lep, xc_arr(:,1), sig_brem)
     call interpol(energy, E_lep, xc_arr(:,2), sig_pair)
     call interpol(energy, E_lep, xc_arr(:,3), sig_pn)
     
     sig_brem = frac * sig_brem
     sig_pair = frac * sig_pair
     sig_pn = frac_pn * sig_pn
     
   end subroutine int_xc_lep
   
   subroutine int_alpha(energy, alpha_sig, alpha)
  !! Interpolate between charged lepton energy & ionization energy loss values.
     
     implicit none
     
     real(dp), intent(in):: energy
     !! Energy value to interpolate at, in GeV.
     real(dp), intent(in) :: alpha_sig(:)
    !! 1D array of ionization energy loss, in (GeV*cm^2)/g.
     
     real(dp), intent(out) :: alpha
     !! Interpolated ionization energy loss value, in (GeV*cm^2)/g.
     
     call interpol(energy, E_lep, alpha_sig, alpha)
     
   end subroutine int_alpha
   
   subroutine int_beta(energy, beta_arr, rho, tot)
     !! Interpolate between charged lepton energy & beta (energy loss parameter) values.

     implicit none
     
     real(dp), intent(in):: energy
     !! Energy value to interpolate at, in GeV.
     real(dp), intent(in) :: beta_arr(:,:)
     !! 2D array of beta values, in cm^2/g.
     real(dp), intent(in) :: rho
     !! Density of material, in g/cm^3.
     
     real(dp), intent(out) :: tot
     !! Interpolated (& summed) value of beta, in cm^2/g.
    
     real(dp) :: brem, pair, pn, frac, frac_pn
     
     call get_rho_frac(rho, frac, frac_pn) ! get the density scaling fractions
     
     call interpol(energy, E_lep, beta_arr(:,1), brem)
     call interpol(energy, E_lep, beta_arr(:,2), pair)
     call interpol(energy, E_lep, beta_arr(:,3), pn)

     tot = (frac*brem)+(frac*pair)+(frac_pn*pn)
     
   end subroutine int_beta

   subroutine searchsorted(array, search_value, binarysearch)
     !! Given an array and a value, returns the index of the element that is closest to, but less than, the given value.
     
     ! Uses a binary search algorithm.
     ! "delta" is the tolerance used to determine if two values are equal
     ! if ( abs(x1 - x2) <= delta) then
     !    assume x1 = x2
     ! endif
     
     implicit none
     
     real(dp), intent(in) :: array(:)
     !! Input array.
     real(dp), intent(in) :: search_value
     !! Value to search for.
     
     integer, intent(out) :: binarysearch
     !! Output result index.

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
     !! Generate random number in the range [0,1).

     implicit none
     
     real(dp), intent(out) :: r
     !! Random number.
     
     call random_number(r)
     
!    r = 0.5 ! for debugging only!
     
   end subroutine random_no
   
   subroutine random(r)
     
     implicit none
     
     real(dp), intent(out) :: r
     !! Random number.
     
     call random_no(r)
     
     if (r<0.5) then
        r = 1_dp
     else
        r = -1_dp
        
     end if
     
   end subroutine random
   
   subroutine idecay(energy, distance, m_le, c_tau, decay)
     !! Calculate decay probability of lepton.
     
     implicit none

     real(dp), intent(in) :: energy
     !! Charged lepton energy, in GeV.
     real(dp), intent(in) :: distance
     !! Distance of charged lepton travel, in cm.
     real(dp), intent(in) :: m_le
     !! Mass of charged lepton, in GeV.
     real(dp), intent(in) :: c_tau
     !! Decay length of charged lepton, in cm.
     
     integer, intent(out) :: decay
     !! Decay = 0 means the charged lepton decayed.
     
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
    !! Calculate the charged lepton electromagnetic energy loss (continuous part) a la MUSIC.

    implicit none
    
    real(dp), intent(in) :: E_init
    !! Initial charged lepton energy, in GeV.
    real(dp), intent(in) :: alpha_val
    !! Ionization energy loss value, in (GeV*cm^2)/g.
    real(dp), intent(in) :: beta_val
    !! Energy loss parameter (brem + pair + pn), in cm^2/g.
    real(dp), intent(in) :: x
    !! Distance (column depth) of charged lepton travel, in g/cm^2.
    real(dp), intent(in) :: m_le
    !! Mass of charged lepton, in GeV.
    
    
    real(dp), intent(out) :: E_fin
    !! Final charged lepton energy, in GeV.
    
    if (beta_val * x < 1e-6_dp) then
       E_fin = E_init * (1._dp-beta_val*x) - alpha_val*x
    else
       E_fin = E_init * dexp(-beta_val*x) - alpha_val/beta_val*(1._dp-dexp(-beta_val*x))
    end if
    
    if (E_fin < 0) then
       E_fin = m_le
    end if
    
  end subroutine em_cont_part

  subroutine int_depth_nu(energy, nu_xc, fac_nu, x_int)
    !! Calculate neutrino interaction depth.
    ! int_depth = M/(N_A*sigma_tot)
    
    implicit none
    
    real(dp), intent(in) :: energy
    !! Neutrino energy, in GeV.
    real(dp), intent(in) :: fac_nu
    !! Rescaling factor for SM neutrino cross-sections.
    real(dp), intent(in) :: nu_xc(:,:)
    !! 2D array containing neutrino CC & NC cross-section values, in cm^2.
    
    real(dp), intent(out) :: x_int
    !! Neutrino interaction depth, in cm^2/g.
    
    real(dp) :: sig_cc, sig_nc, sig_weak
    
    call int_xc_nu(energy, nu_xc, fac_nu, sig_cc, sig_nc) ! initialize CC & NC xc interpolations; moved fac_nu to here
    sig_weak = sig_cc + sig_nc ! weak interactions
    x_int = 1._dp/(N_A*sig_weak)
    
  end subroutine int_depth_nu

  subroutine int_depth_lep(energy, xc_arr, rho, m_le, c_tau, x_int)
    !! Calculate charged lepton interaction depth.
    ! int_depth = M/((N_A/A)*sigma_tot + 1/(gamma*c*tau*rho)); here we need rho to convert decay distance to decay depth

    implicit none
    
    real(dp), intent(in) :: energy
    !! Charged charged lepton energy, in GeV.
    real(dp), intent(in) :: rho
    !! Density of material, in g/cm^3.
    real(dp), intent(in) :: m_le
    !! Mass of charged lepton, in GeV.
    real(dp), intent(in) :: c_tau
    !! Decay length of charged lepton, in cm.
    real(dp), intent(in) :: xc_arr(:,:)
    !! 2D array containing N_A/A*charged lepton-nucleon cross-section values, in cm^2/g.
    
    real(dp), intent(out) :: x_int
    !! Charged lepton interaction length, in g/cm^2.
    
    real(dp) :: sig_cc, sig_nc, sig_brem, sig_pair, sig_pn, sig_em, sig_weak
    real(dp) :: decay_length, decay_depth_inv
    
    sig_cc = 0 ! placeholder for CC (charged) lepton interactions; can be read in from the lookup table in the future
    sig_nc = 0 ! placeholder for NC (charged) lepton interactions; can be read in from the lookup table in the future
    
    decay_length = (energy/m_le)*c_tau
    decay_depth_inv = 1._dp/(decay_length*rho)
    
    call int_xc_lep(energy, xc_arr, rho, sig_brem, sig_pair, sig_pn) ! initialize Brem, pair & pn xc interpolations
    
    sig_em = sig_brem + sig_pair + sig_pn ! EM interactions
    sig_weak = sig_cc + sig_nc ! weak interactions
    x_int = 1/(sig_em + sig_weak + decay_depth_inv)
    ! Note: N_A has already been multiplied by xc's in the lookup tables for taus and muons!
    
  end subroutine int_depth_lep
  
  subroutine interaction_type_nu(energy, nu_xc, fac_nu, int_type)
    !! Determine/calculate the type of neutrino-nucleon interaction.
    
    implicit none
    
    real(dp), intent(in) :: energy
    !! Neutrino energy, in GeV.
    real(dp), intent(in) :: fac_nu
    !! Rescaling factor for SM neutrino cross-sections.
    real(dp), intent(in) :: nu_xc(:,:)
    !! 2D array containing neutrino CC & NC cross-section values, in cm^2.
    
    integer, intent(out) :: int_type
    !! Type of neutrino interaction. 0=CC; 1=NC.
    
    real(dp) :: sig_cc, sig_nc, x, tot_frac, cc_frac
    
    call int_xc_nu(energy, nu_xc, fac_nu, sig_cc, sig_nc) ! initialize CC & NC xc interpolations
    
    tot_frac = sig_cc + sig_nc
    cc_frac = sig_cc/tot_frac
    
    call random_no(x)
    
    if (x <= cc_frac) then
       int_type = 0 ! CC
    else
       int_type = 1 ! NC
    end if
    
  end subroutine interaction_type_nu

  subroutine interaction_type_lep(energy, xc_arr, rho, m_le, c_tau, int_type)
    !! Determine/calculate the type of charged lepton-nucleon interaction.
    
    implicit none
    
    real(dp), intent(in) :: energy
    !! Charged lepton energy, in  GeV.
    real(dp), intent(in) :: rho
    !! Density of material, in g/cm^3.
    real(dp), intent(in) :: m_le
    !! Mass of charged lepton, in GeV.
    real(dp), intent(in) :: c_tau
    !! Decay length of charged lepton, in cm.
    real(dp), intent(in) :: xc_arr(:,:)
    !! 2D array containing N_A/A*charged lepton-nucleon cross-section values, in cm^2/g.
    
    integer, intent(out) :: int_type
    !! Type of lepton interaction. 2=decay; 3=Bremmstrahlung; 4=pair-production; 5=photonuclear; 6=CC/NC (placeholder).
    
    real(dp) :: sig_cc, sig_nc, sig_brem, sig_pair, sig_pn, decay_length, decay_depth_inv, int_lep
    real(dp) :: tot_frac, decay_frac, cc_frac, nc_frac, brem_frac, pair_frac, pn_frac, y
    real(dp) :: sig_em, sig_weak
    
    sig_cc = 0 ! placeholder for CC lepton interactions
    sig_nc = 0 ! placeholder for NC lepton interactions
    call int_xc_lep(energy, xc_arr, rho, sig_brem, sig_pair, sig_pn)
    
    decay_length = (energy/m_le)*c_tau
    decay_depth_inv = 1._dp/(decay_length*rho)
    ! call int_depth_lep(energy, xc_arr, rho, m_le, c_tau, int_lep)
    sig_em = sig_brem + sig_pair + sig_pn ! EM interactions
    sig_weak = sig_cc + sig_nc ! weak interactions
    int_lep = 1/(sig_em + sig_weak + decay_depth_inv)
    
    tot_frac = 1._dp/int_lep
    decay_frac = decay_depth_inv/tot_frac
    cc_frac = sig_cc/tot_frac ! placeholder for CC (charged) lepton interactions
    nc_frac = sig_nc/tot_frac ! placeholder for NC (charged) lepton interactions
    brem_frac = sig_brem/tot_frac
    pair_frac = sig_pair/tot_frac
    pn_frac = sig_pn/tot_frac

    call random_no(y)
    
    if (y <= decay_frac) then
       int_type = 2 ! decay
    else if (decay_frac < y .and. y <= decay_frac+brem_frac) then
       int_type = 3 ! Bremmstrahlung
    else if (decay_frac+brem_frac < y .and. y <= decay_frac+brem_frac+pair_frac) then
       int_type = 4 ! pair-production
    else if (decay_frac+brem_frac+pair_frac < y .and. y <= decay_frac+brem_frac+pair_frac+pn_frac) then
       int_type = 5 ! photonuclear
    else
       int_type = 6 ! lep_nc => cc/nc (both need CDFs) - placeholder for now
    end if
    
  end subroutine interaction_type_lep
  
  subroutine find_y(energy, ixc_arr, ip, y)  
    !! Stochastic determination of neutrino/lepton inelasticity.

    implicit none
    
    real(dp), intent(in) :: energy
    !! Neutrino or charged lepton energy, in GeV.
    real(dp), intent(in) :: ixc_arr(:,:,:)
    !! Neutrino or charged lepton integrated cross-section CDF values.
    integer, intent(in) :: ip                                            
    !! Type of neutrino-nucleon or lepton-nucleon interaction.
    
    real(dp), intent(out) :: y
    !! Inelasticity, y = (E_init-E_final)/E_initial.

    real(dp) :: dy, search_arr(31) ! interpolating_arr is f(x)
    integer :: ip_id, energy_index

    call random_no(dy)
    
    if (ip == 0 .or. ip == 1) then ! basically, for neutrinos
       call searchsorted(E_nu, energy, energy_index)
       
       if (ip == 0) then ! CC
          ip_id = 1
       else
          ip_id = 2
       end if
       
    else ! for charged leptons
       call searchsorted(E_lep, energy, energy_index)
       
       if (ip == 3) then ! Brem
          ip_id = 1
       else if (ip == 4) then ! pair
          ip_id = 2
       else if (ip == 5) then ! pn
          ip_id = 3
       else ! lep_nc?
          ip_id = 4 ! shouldn't happen, for now
       end if
    end if

    search_arr = ixc_arr(:,energy_index,ip_id)
    
    call interpol(dy, search_arr, yvals, y) ! interpolate in yvals directly
    ! dy is the randomly sampled cross-section CDF value (between 0 & 1)
    ! search_arr = cross-section CDF value array for energy_index
    ! yvals = array of min. y values from which the cross-section CDF is calculated (see models.py for calculation details) 
    ! y is the interpolated (yvals) value corresponding to the cross-section CDF value = dy; this y is responsible for stochastic energy losses

    if (y > 1._dp) then
       y = 1._dp
    end if
    
  end subroutine find_y
  
  subroutine polarization(y, pin, theta_in, pout, theta_out)
  
    implicit none
    
    real(dp), intent(in) :: pin, theta_in, y
    real(dp), intent(out) :: pout, theta_out
    
    real(dp) :: pzout, theta, p0, rs, cth
  
    if (y<0.01) then
       theta_out = theta_in
       pout = pin
       return
    else
       call interpol(y, ypol, Pcthp, pzout)  !!interpolating Pcthp value for y 
       call interpol(y, ypol, P, pout)   !!interpolating P value for y 
       
    end if

    cth = pzout/pout
    theta_out = acos(cth)
    p0 = pin*pout
    pout = p0
    
    call random(rs)
    
    theta = theta_in + rs*theta_out
    theta_out = theta
    
  end subroutine polarization

  subroutine propagate_nu(e_init, nu_xc, nu_ixc, depth_max, fac_nu, part_type, d_travel, e_fin)
    !! Propagates a neutrino inside the Earth.
    
    implicit none
    
    real(dp), intent(in) :: e_init
    !! Initial neutrino energy, in GeV.
    real(dp), intent(in) :: nu_xc(:,:)
    !! 2D array containing neutrino CC & NC cross-section values, in cm^2.
    real(dp), intent(in) :: nu_ixc(:,:,:)
    !! 3D array containing neutrino integrated cross-section CDF values.
    real(dp), intent(in) :: depth_max
    !! Maximum column depth for neutrino propagation, in kmwe.
    real(dp), intent(in) :: fac_nu
    !! Rescaling factor for SM neutrino cross-sections.

    integer, intent(out) :: part_type
    !! Type of outgoing particle. 0=neutrino; 1=charged lepton.
    real(dp), intent(out) :: d_travel
    !! Distance traveled until converted to charged lepton or total distance traveled by neutrino (if no conversion to charged lepton), in kmwe.
    real(dp), intent(out) :: e_fin
    !! Final neutrino energy, in GeV.

    real(dp) :: e_nu, int_depth, x, x_f, x_0, y, col_depth_total, r    
    integer :: int_type !int_type=0:CC, 1:NC
    
    part_type = 0 ! starting off as a neutrino
    col_depth_total = depth_max*1e5_dp ! kmwe to cmwe
    e_nu = e_init
    e_fin = e_init
    x_0 = 0.0_dp ! starting depth in g/cm^2
    d_travel = depth_max ! added this in case there is a problem and needed for breaking out when E<1e3

    do while (e_nu > 1e3_dp)
       call random_no(r)
       
       call int_depth_nu(e_nu, nu_xc, fac_nu, int_depth)
       x = -int_depth*dlog(r) ! prob of interaction=exp(-x/int_depth)
       x_f = x_0 + x ! x_f is tracking total column depth traveled
       
       if (x_f > col_depth_total) then ! total col depth exceeded
          return ! returns (depth_max, e_init, neutrino)
       end if
       
       x_0 = x_f
       call interaction_type_nu(e_nu, nu_xc, fac_nu, int_type) ! CC or NC interaction?
       
       if (part_type == 0 .and. int_type == 1) then
            part_type = 0
         else if (part_type == 0 .and. int_type == 0) then
            part_type = 1
         end if
         
        call find_y(e_nu, nu_ixc, int_type, y)
        
        e_fin = e_nu*(1._dp-y) ! energy transfer; y = (E_init-E_final)/E_initial
        
        if (part_type == 1) then ! converted to tau
           d_travel = x_0*1e-5_dp ! cmwe to kmwe
           return ! returns (d_travel, e_nu<=1e3, tau)
        end if
        
        ! still a neutrino; keep going
        e_nu = e_fin
        d_travel = x_0*1e-5_dp ! cmwe to kmwe

        if (e_nu <= 1e3_dp) then ! our E_nu lookup tables only go to 10^3 GeV
           return ! returns (d_travel, e_nu<=1e3, neutrino)
        end if
     end do
     
     
   end subroutine propagate_nu

   subroutine propagate_lep_water(e_init, xc_water, lep_ixc, alpha_water, beta_water, d_in, lepton,&
        & prop_type, cthi, Pi, part_id, d_fin, e_fin, pcthf)
    !! Propagates a charged lepton in water inside the Earth.
     
     implicit none
     
     real(dp), intent(in) :: e_init
     !! Initial energy of the charged lepton, in GeV.
     real(dp), intent(in) :: xc_water(:,:)
     !! 2D array containing N_A/A*charged lepton-nucleon cross-section values in water, in cm^2/g.
     real(dp), intent(in) :: lep_ixc(:,:,:)
     !! 3D array containing charged lepton integrated cross-section CDF values in water.
     real(dp), intent(in) :: alpha_water(:)
     !! 1D array containing ionization energy loss values in water, in (GeV*cm^2)/g.
     real(dp), intent(in) :: beta_water(:,:)
    !! 2D array of beta values in water, in cm^2/g.
     real(dp), intent(in) :: d_in
     !! Maximum distance for charged lepton to propagate in water, in kmwe.
     integer, intent(in) :: lepton
     !! Type of charged lepton. 1=tau; 2=muon.
     integer, intent(in) :: prop_type
     !! Type of energy loss propagation. 1=stochastic, 2=continuous.
     real(dp), intent(in) :: cthi
    !! costheta value obtained from tau EM interaction with rock
    real(dp), intent(in) :: Pi
    !! degree of polarization obtained from tau EM interaction with rock
    
     integer, intent(out) :: part_id
     !! Type of outgoing charged lepton. 0=decayed; 1=not decayed; 2=don't count.
     real(dp), intent(out) :: d_fin
     !! Distance traveled before charged lepton decays or total distance traveled by charged lepton, in kmwe.
     real(dp), intent(out) :: e_fin
     !! Final energy of the charged lepton, in GeV.
     real(dp), intent(out) :: pcthf
     !! Final polarization after EM interaction of the tau lepton
     

     real(dp) :: e_min, x_total, e_lep, x_0, m_le, c_tau, r, int_depth
     real(dp) :: x, x_f, x_step, alpha, beta, e_int, e_avg, y, delta_d, delta_x
     real(dp) :: Pin, theta_in, Pout, theta_out  !Varibales needed to calculate polarization for PN interaction
     integer :: int_type, j_max
     integer(dp) :: cnt
    
     e_min = 1e3_dp ! minimum tau energy, in GeV
     
     part_id = 1 ! start with tau that's obviously not decayed
     
     x_total = d_in*1e5_dp ! kmwe to cmwe
     e_lep = e_init
     e_fin = e_init ! in case the first interaction is too far
     x_0 = 0._dp ! haven't gone anywhere yet; initiate tracker
     
     if (lepton == 1) then
        m_le = 1.77682_dp ! m_tau in GeV
        c_tau = 8.703e-3_dp ! c*lifetime, in cm, for taus (taken from PDB)
     else
        m_le = 0.10565837550000001_dp ! m_mu in GeV
        c_tau = 6.586384e4_dp ! c*lifetime, in cm, for muons (taken from PDB 2020)
     end if
     
     if (prop_type == 1) then ! stochastic energy loss
        Pin = Pi
        theta_in = acos(cthi)
        
        do while (e_lep > e_min)
           
           if (e_lep <= e_min) then
              exit ! taken care of outside the do while loop
           end if
           
           call random_no(r)
           call int_depth_lep(e_lep, xc_water, rho_water, m_le, c_tau, int_depth)
           x = -int_depth*dlog(r) ! basically the step size
           ! prob of interaction=exp(-x/int_depth)
           
           x_f = x_0 + x ! how far have we traveled here
           d_fin = x_f/1e5_dp ! make sure it is not past the old number, in kmwe
           
           if (x_f >= x_total) then ! already past maximum depth but still a tau
              x_step = x_total - x_0 ! backtrack one step
              
              call int_alpha(e_lep, alpha_water, alpha) ! changed 12/9/2020
              call int_beta(e_lep, beta_water, rho_water, beta) ! changed 12/9/2020
              
              e_fin = e_lep - (e_lep*beta + alpha)*x_step
              d_fin = d_in
              
              if (e_fin <= e_min) then ! tau has decayed
                 d_fin = d_in ! just in case; added 12/9/2020
                 e_fin = e_min
                 part_id = 2 ! don't count
              end if
              pcthf = Pin*cos(theta_in)
              return ! d_fin = d_in
           end if
           
           x_0 = x_f ! update x_0 and keep going
           call int_alpha(e_lep, alpha_water, alpha)
           call int_beta(e_lep, beta_water, rho_water, beta)

           e_int = e_lep - (e_lep*beta + alpha)*x ! find some intermediate energy to get reasonable values of energy between initial and final energy, a la MUSIC
           
           if (e_int <= e_min) then
              e_int = e_min
           end if
           
           e_avg = 10._dp**((dlog10(e_lep)+dlog10(e_int))/2._dp) ! avg. of 10^7 & 10^8 is 10^(7.5)
           
           call int_alpha(e_avg, alpha_water, alpha)
           call int_beta(e_avg, beta_water, rho_water, beta)
            
           call em_cont_part(e_lep, alpha, beta, x, m_le, e_int) ! get the continuous energy loss for this average energy
           
           if (e_int <= e_min) then ! below minimum energy
              exit ! taken care of outside the do while loop
           end if
           
           call interaction_type_lep(e_int, xc_water, rho_water, m_le, c_tau, int_type)
            
           if (int_type == 2) then ! tau has decayed
              part_id = 0
              e_fin = e_int
              d_fin = x_f/1e5_dp
              pcthf = Pin*cos(theta_in)
              return ! basically what d_fin does this return?
           end if
           
           ! tau didn't decay. Now how much energy does it have after interaction?
           
            call find_y(e_int, lep_ixc, int_type, y) ! stochastic energy loss sampling

            ! outgoing tau energy is old e_lep*(1-y)
            e_lep = e_int*(1._dp-y) ! this is the energy for the next interaction
            e_fin = e_lep

            !! This routine will run for only PN interaction 
            if (int_type == 5) then          
             call polarization(y, Pin, theta_in, Pout, theta_out)
             Pin = Pout
             theta_in = theta_out
          end if
        end do
        
        ! Outside the while loop, e_lep has to be <= e_min
        if (e_lep <= e_min) then
           d_fin = d_in ! max. distance in water
           e_fin = e_min
           part_id = 2 ! don't count this
           !print *, "I have less energy than emin!"
           pcthf = Pin*cos(theta_in)
           return
        end if
        
     else ! continuous energy loss
        
        j_max = int(x_total/(step_size*rho_water)) ! we will get close to exit.
        
        do cnt = 1, j_max+1 ! takes care of the integer truncation issue
           if (e_lep < e_min) then
              exit ! taken care of outside the do while loop
           end if
           
            delta_x = step_size * rho_water ! distance goes into decay
            x_f = x_0 + delta_x
            ! does the particle decay over this distance?
            call idecay(e_lep, step_size, m_le, c_tau, part_id)
            
            if (part_id == 0) then ! we are all done
               e_fin = e_lep
               d_fin = d_in
               return
            else ! find the new energy; assume alpha and beta are total values, not cut values
               call int_alpha(e_lep, alpha_water, alpha) ! changed 12/9/2020
               call int_beta(e_lep, beta_water, rho_water, beta) ! changed 12/9/2020
               
               e_fin = e_lep - (e_lep*beta + alpha)*delta_x
               d_fin = x_f/1e5_dp
               x_0 = x_f
               e_lep = e_fin
               
               if (cnt >= j_max) then
                  exit ! if so, break out of while loop
               end if

            end if
            
         end do

         if (cnt >= j_max) then
            x_step = x_total - x_f ! backtrack a step
            if (x_step > 0._dp) then !last little energy loss
               e_fin = e_lep - (e_lep*beta + alpha)*x_step
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
       & d_in, xalong, cdalong, idepth, lepton, prop_type, part_id, d_fin, e_fin, cthf, Pf)
    !! Propagates a charged lepton in rock & iron inside the Earth.
    
    implicit none
    
    real(dp), intent(in) :: angle
    !! Earth emergence angle (beta), in degrees.
    real(dp), intent(in) :: e_init
    !! Initial energy of charged lepton, in GeV.
    real(dp), intent(in) :: xc_rock(:,:)
    !! 2D array containing N_A/A*charged lepton-nucleon cross-section values in rock, in cm^2/g.
    real(dp), intent(in) :: lep_ixc(:,:,:)
    !! 3D array containing charged lepton integrated cross-section CDF values in rock.
    real(dp), intent(in) :: alpha_rock(:)
    !! 1D array containing ionization energy loss values in rock, in (GeV*cm^2)/g.
    real(dp), intent(in) :: beta_rock(:,:)
    !! 2D array of beta values in rock, in cm^2/g.
    real(dp), intent(in) :: d_entry
    !! Column depth along the chord for a given Earth emergence angle, in kmwe.
    real(dp), intent(in) :: d_in
    !! How much distance in rock/iron the charged lepton is supposed to travel, in kmwe.
    real(dp), intent(in) :: xalong(:)
    !! 1D array containing distance in water, in km.
    real(dp), intent(in) :: cdalong(:)
    !! 1D array containing column depth at xalong, in g/cm^2.
    integer, intent(in) :: idepth
    !! Depth of water layer in km.
    integer, intent(in) :: lepton
    !! Type of charged lepton. 1=tau; 2=muon.
    integer, intent(in) :: prop_type
    !! Type of energy loss propagation. 1=stochastic, 2=continuous.
    
    integer, intent(out) :: part_id
    !! Type of outgoing charged lepton. 0=decayed; 1=not decayed; 2=don't count.
    real(dp), intent(out) :: d_fin
    !! Distance traveled before charged lepton decays or total distance traveled by charged lepton, in kmwe.
    real(dp), intent(out) :: e_fin
    !! Final charged lepton energy, in GeV.
    real(dp), intent(out) :: cthf
    !! Final costheta after EM interaction of the tau lepton
    real(dp), intent(out) :: Pf
    !! Final degree of polarization after EM interaction of the tau lepton
    
    real(dp) :: e_min, col_depth, d_max, e_lep, x_0, m_le, c_tau, x_interp, r
    real(dp) :: rho, dy, int_len, x, x_f, alpha, beta, e_int, e_avg, y, d_0, delta_x, d_rem, delta_d
    real(dp) :: Pin, theta_in, Pout, theta_out  !Varibales needed to calculate polarization for PN interaction
    integer :: cnt, j_max
    integer :: int_type
    
    e_min = 1e3_dp ! minimum tau energy, in GeV
    part_id = 1 ! start with tau
    col_depth = d_entry*1e5_dp ! how far in
    d_max = d_in*1e5_dp ! how much to go, in cm.w.e
    e_lep = e_init
    x_0 = 0._dp
    cnt = 0
    ! x_total = d_in*1e5
    !rho = rho_rock !g/cm^3  ! USE FOR TESTING P_SURV FOR ROCK ONLY!
  
    if (lepton == 1) then
       m_le = 1.77682_dp ! m_tau in GeV
       c_tau = 8.703e-3_dp ! c*lifetime, in cm, for taus (taken from PDB)
    else
       m_le = 0.10565837550000001_dp ! m_mu in GeV
       c_tau = 6.586384e4_dp ! c*lifetime, in cm, for muons (taken from PDB 2020)
    end if
    
    if (prop_type == 1) then ! stochastic energy loss
       Pin = 1._dp
       theta_in = 0._dp
       
       do while (e_lep > e_min)
          cnt = cnt + 1
          call cd2distd(xalong, cdalong, col_depth, x_interp) ! find how far we are along the chord for given beta
          call densityatx(x_interp, angle, idepth, r, rho) ! find rho at x
          
          call random_no(dy)
          call int_depth_lep(e_lep, xc_rock, rho, m_le, c_tau, int_len)
        
          x = -int_len * dlog(dy)
          col_depth = col_depth + x ! update along trajectory, from the start of the chord
          x_f = x_0 + x ! going 1D
          d_fin = x_f/1e5_dp ! make sure it is not past the old number, in km.w.e
          
          if (x_f > d_max) then ! already past max depth
             ! go to 30
             d_rem = d_max - x_0
             call int_alpha(e_lep, alpha_rock, alpha)
             call int_beta(e_lep, beta_rock, rho, beta)
             e_fin = e_lep - (e_lep*beta + alpha)*d_rem
             d_fin = d_max/1e5_dp !this might be the problem here!!!!!!!!!!
             if (e_fin > e_init) then
                e_fin=e_init  
             end if
             !print *,"Passed Earth already!"
             Pf = Pin
             cthf = cos(theta_in)
             return
             
          end if
        
          x_0 = x_f ! update x_0 and keep going
          call int_alpha(e_lep, alpha_rock, alpha)
          call int_beta(e_lep, beta_rock, rho, beta)
          e_int = e_lep - (e_lep*beta + alpha)*x ! find some intermediate energy to get reasonable values of energy between initial and final energy, a la MUSIC
          
          if (e_int <= e_min) then
             e_int = e_min
          end if
          
          e_avg = 10._dp**((dlog10(e_lep)+dlog10(e_int))/2) ! does this work?
          
          call int_alpha(e_avg, alpha_rock, alpha)
          call int_beta(e_avg, beta_rock, rho, beta)
          
          call em_cont_part(e_lep, alpha, beta, x, m_le, e_int) ! get the continuous energy
          
          if (e_int <= e_min) then ! is it below minimum energy now?
             e_fin = e_int
             ! go to 20
             d_fin = d_max/1e5_dp ! in km.w.e
             e_fin = e_min
             part_id = 0
             print *, "too small to handle"
             !pcthf = Pin*cos(theta_in)
             Pf = Pin
             cthf = cos(theta_in)
             return
          end if
        
          call interaction_type_lep(e_int, xc_rock, rho, m_le, c_tau, int_type)
        
          if (int_type == 2) then ! tau has decayed
             part_id = 0
             e_fin = e_int
             ! go to 50
             ! 50 continue
             Pf = Pin
             cthf = cos(theta_in)
             return
          end if
          
          ! tau didn't decay. Now how much energy does it have after interaction?
          call find_y(e_int, lep_ixc, int_type, y)
          
          ! outgoing tau energy is old e_lep*(1-y)
          e_lep = e_int*(1._dp-y) ! this is the energy for the next interaction
          e_fin = e_lep
          
          !! This routine will run for only PN interaction 
          if (int_type == 5) then          
             !print *, "Is PN"
             call polarization(y, Pin, theta_in, Pout, theta_out)
             Pin = Pout
             theta_in = theta_out
          end if
                
       end do
       
     ! Outside the while loop, e_lep has to be < e_min
       if (e_lep <= e_min) then ! only continuous energy loss; go to 20
          d_fin = d_max/1e5_dp
          e_fin = e_min
          part_id = 0 ! decayed or no_count??? should be decayed
          Pf = Pin
          cthf = cos(theta_in)
          return
       end if
       
    else ! continuous energy loss
       
       d_0 = 0._dp
       delta_x = step_size ! for now, not adaptive, distance into decay, cm; works for taus
       !        delta_d = 5000._dp ! test
       
       call cd2distd(xalong, cdalong, col_depth, x_interp) ! find how far we are along the chord for given beta
       call densityatx(x_interp, angle, idepth, r, rho) ! find rho at x
       
       !        rho = rho_rock ! FOR TESTING P_SURV ONLY!!
       
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
          
          call cd2distd(xalong, cdalong, col_depth, x_interp) ! find how far we are along the chord for given beta
          call densityatx(x_interp, angle, idepth, r, rho) ! find rho at x
          
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
             call int_beta(e_lep, beta_rock, rho, beta) ! changed 12/9/2020
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
       & alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth, lepton, prop_type, part_type, d_fin, e_fin, pcthf)
    
    implicit none
    
    real(dp), intent(in) :: angle
    !! Earth emergence angle (beta), in degrees.
    real(dp), intent(in) :: depth
    !! Total column depth of the chord, in kmwe.
    real(dp), intent(in) :: d_water
    !! Column depth of the final layer of water (or full distance in water if only water layer), in kmwe.
    real(dp), intent(inout) :: depth_traj
    !! Column depth along the chord for a given Earth emergence angle, in kmwe.
    real(dp), intent(in) :: e_lep_in
    !! Ingoing charged lepton energy, in GeV.
    real(dp), intent(in) :: xc_water(:,:)
    !! 2D array containing N_A/A*charged lepton-nucleon cross-section values in water, in cm^2/g.
    real(dp), intent(in) :: xc_rock(:,:)
    !! 2D array containing N_A/A*charged lepton-nucleon cross-section values in rock, in cm^2/g.
    real(dp), intent(in) :: lep_ixc_water(:,:,:)
    !! 3D array containing charged lepton integrated cross-section CDF values in water.
    real(dp), intent(in) :: lep_ixc_rock(:,:,:)
    !! 3D array containing charged lepton integrated cross-section CDF values in rock.
    real(dp), intent(in) :: alpha_water(:)
    !! 1D array containing ionization energy loss values in water, in (GeV*cm^2)/g.
    real(dp), intent(in) :: alpha_rock(:)
    !! 1D array containing ionization energy loss values in rock, in (GeV*cm^2)/g.
    real(dp), intent(in) :: beta_water(:,:)
    !! 2D array of beta values in water, in cm^2/g.
    real(dp), intent(in) :: beta_rock(:,:)
    !! 2D array of beta values in rock, in cm^2/g.
    real(dp), intent(in) :: xalong(:)
    !! 1D array containing distance in water, in km.
    real(dp), intent(in) :: cdalong(:)
    !! 1D array containing column depth at xalong, in g/cm^2.
    integer, intent(in) :: idepth
    !! Depth of water layer in km.
    integer, intent(in) :: lepton
    !! Type of charged lepton. 1=tau; 2=muon.
    integer, intent(in) :: prop_type
    !! Type of energy loss propagation. 1=stochastic, 2=continuous.
    
    integer, intent(out) :: part_type
    !! Type of outgoing charged lepton. 0=decayed, 1=not decayed.
    real(dp), intent(out) :: d_fin
    !! Distance traveled before charged lepton decays or total distance traveled by charged lepton, in kmwe.
    real(dp), intent(out) :: e_fin
    !! Outgoing charged lepton energy, in GeV.
    real(dp), intent(out) :: pcthf
    !! Final polarization after EM interaction of the tau lepton
    
    real(dp) :: col_depth, rho, x, r, d_in, d_f, e_lep
    real(dp) :: Pi, cthi, Pf, cthf
    
    d_fin = depth_traj
    col_depth = depth_traj*1e5_dp ! g/cm^2
    e_lep = e_lep_in ! so e_lep_in doesn't change
    e_fin = e_lep
    part_type = 1 ! tau going in
    
    if (e_lep < 1e3_dp) then ! just in case
       part_type = 0
       d_fin = depth
       return
    end if
    
    if (depth-depth_traj < d_water) then
       rho = rho_water ! water
    else
       call cd2distd(xalong, cdalong, col_depth, x)
       call densityatx(x, angle, idepth, r, rho) ! find rho at x
       
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
       d_in = depth - depth_traj - d_water ! propagate this far in rock
       call propagate_lep_rock(angle, e_lep, xc_rock, lep_ixc_rock, alpha_rock, beta_rock, depth_traj,&
            & d_in, xalong, cdalong, idepth, lepton, prop_type, part_type, d_f, e_fin, cthf, Pf)

       depth_traj = depth_traj + d_f
       
       if (part_type == 1 .and. idepth /= 0) then ! still a tau; added .and. clause on 3/18
          
          e_lep = e_fin
          d_in = d_water
          !depth_traj = depth_traj + d_f ! now propagate through final layer of water
          !d_fin = depth_traj
          Pi = Pf
          cthi = cthf
          
          call propagate_lep_water(e_lep, xc_water, lep_ixc_water, alpha_water, beta_water, d_in,&
               & lepton, prop_type, cthi, Pi, part_type, d_f, e_fin, pcthf)

          d_fin = depth_traj + d_f
          
       else !either tau decayed or idepth==0
          pcthf = Pf*cthf
          d_fin = depth_traj
          return
          
       end if

    else
       d_in = depth - depth_traj
       Pi = 1._dp
       cthi = cos(0._dp)

       call propagate_lep_water(e_lep, xc_water, lep_ixc_water, alpha_water, beta_water, d_in, lepton,&
            & prop_type, cthi, Pi, part_type, d_f, e_fin, pcthf)
       
       depth_traj = depth_traj + d_f
       d_fin = depth_traj ! needed to add this since return was depth_traj here
       return
    end if
    
    return
    
  end subroutine tau_thru_layers

  subroutine distnu(r, ithird, Pin, dist_val)
    !! Determines the neutrino energy from tau decay.
    ! The energy fraction is determined by tau energy CDF, approximated by leptonic decay channel.
    ! Approximated by left-handed leptonic decay channel.
    ! Approximate (good enough) for taus; exact for muons

    implicit none
    
    real(dp), intent(in) :: r
    !! Random number.
    integer, intent(in) :: ithird
    !! Choice for neutrino -> charged lepton energy fraction selection.
    real(dp), intent(in) :: Pin
    !! tau's polarization
    
    real(dp), intent(out) :: dist_val
    !! Energy fraction, y = E_nu_tau/E_tau.

    real(dp) :: fnu, y, fm, ff, y0, y1
    !real(dp) :: P  

    !P = -1._dp * Pin ! -1 = fully LH polarized; 0 = fully unpolarized tau
    
    fnu(y) = y/3._dp * (5._dp - 3._dp* y**2 + y**3) + (- 1._dp * Pin) * (y/3._dp * (1._dp - 3._dp * y**2 + 2._dp * y**3))
    
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
       
    else ! if ithird == 1; use 1/3 of energy of 3 body decay.
       dist_val = 1._dp/3._dp
    end if

  end subroutine distnu

  subroutine regen(angle, e_lep, depth, d_water, d_lep, nu_xc, nu_ixc, ithird, xc_water, xc_rock,&
       & ixc_water, ixc_rock, alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth,&
       & lepton, fac_nu, prop_type, Pin, part_type, d_exit, e_fin, Pout)
    !! Regeneration loop.   !!this should take a pin and also throw out pout. pin will be used by distnu and the pout it will throw will be used by regen again as input in the single_stat(). 
    
    implicit none
    
    real(dp), intent(in) :: angle
    !! Earth emergence angle (beta), in degrees.
    real(dp), intent(in) :: e_lep
    !! Incoming charged lepton energy, in GeV.
    real(dp), intent(in) :: depth
    !! Total column depth of the chord, in kmwe.
    real(dp), intent(in) :: d_water
    !! Column depth of the final layer of water (or full distance in water if only water layer), in kmwe.
    real(dp), intent(inout) :: d_lep
    !! Column depth along the chord for a given Earth emergence angle, in kmwe.
    real(dp), intent(in) :: nu_xc(:,:)
    !! 2D array containing neutrino CC & NC cross-section values, in cm^2.
    real(dp), intent(in) :: nu_ixc(:,:,:)
    !! 3D array containing neutrino integrated cross-section CDF values.
    integer, intent(in) :: ithird
    !! Choice for neutrino -> charged lepton energy fraction selection.
    real(dp), intent(in) :: xc_water(:,:)
    !! 2D array containing N_A/A*charged lepton-nucleon cross-section values in water, in cm^2/g.
    real(dp), intent(in) :: xc_rock(:,:)
    !! 2D array containing N_A/A*charged lepton-nucleon cross-section values in rock, in cm^2/g.
    real(dp), intent(in) :: ixc_water(:,:,:)
    !! 3D array containing charged lepton integrated cross-section CDF values in water.
    real(dp), intent(in) :: ixc_rock(:,:,:)
    !! 3D array containing charged lepton integrated cross-section CDF values in rock.
    real(dp), intent(in) :: alpha_water(:)
    !! 1D array containing ionization energy loss values in water, in (GeV*cm^2)/g.
    real(dp), intent(in) :: alpha_rock(:)
    !! 1D array containing ionization energy loss values in rock, in (GeV*cm^2)/g.
    real(dp), intent(in) :: beta_water(:,:)
    !! 2D array of beta values in water, in cm^2/g.
    real(dp), intent(in) :: beta_rock(:,:)
    !! 2D array of beta values in rock, in cm^2/g.
    real(dp), intent(in) :: xalong(:)
    !! 1D array containing distance in water, in km.
    real(dp), intent(in) :: cdalong(:)
    !! 1D array containing column depth at xalong, in g/cm^2.
    integer, intent(in) :: idepth
    !! Depth of water layer in km.
    integer, intent(in) :: lepton
    !! Type of charged lepton. 1=tau; 2=muon.
    real(dp), intent(in) :: fac_nu
    !! Rescaling factor for SM neutrino cross-sections.
    integer, intent(in) :: prop_type
    !! Type of energy loss propagation. 1=stochastic, 2=continuous.
    real(dp), intent(in) :: Pin
    !!tau's polarization input from tau passing thru layers
    
    integer, intent(out) :: part_type
    !! Type of outgoing particle. 0=neutrino; 3=exit.
    real(dp), intent(out) :: d_exit
    !! Distance traveled before charged lepton decays or total distance traveled by charged lepton, in kmwe.
    real(dp), intent(out) :: e_fin
    !! Final particle energy, in GeV.
    real(dp), intent(out) :: Pout
    !!tau's polarization output from regen function
    
    real(dp) :: r, frac, e_nu, d_left, dtr, etau2
    real(dp) :: Pi !intermediate polarization of tau
    integer :: int_part

    call random_no(r)
    call distnu(r, ithird, Pin, frac)  !!**that pola input goes here as input and cal. the neutrino energy **!!
    e_nu = frac * e_lep

    d_left = depth-d_lep ! this is how far the neutrino can go
    e_fin = e_nu ! in case we need to exit
    part_type = 3 ! in case we need to exit
    
    if (d_left <= 0._dp) then ! past the point of interactions allowed
       d_exit = depth
       Pout = Pin
       ! go to 60
       ! 60 continue
       return
    end if

    d_exit = d_lep ! we are starting this far into the Earth with a neutrino
    int_part = 0 ! starting with a neutrino with energy e_nu; change later to string
    
    ! tnu follows NC to the end, or gives results if CC interactions
    call propagate_nu(e_nu, nu_xc, nu_ixc, d_left, fac_nu, int_part, dtr, etau2) ! does the neutrino interact?
    
    if (int_part /= 1) then ! neutrinos at the end
       d_exit = depth
       part_type = 0 ! (HLS = 0); changed 22/12/2020
       e_fin = etau2 ! final neutrino energy
       ! go to 60; all done
       ! 60 continue
       Pout = -2._dp  !fake polarization value which will help us filter out the neutrino events in python
       return
    end if

    ! otherwise we have a tau
    d_lep = d_lep + dtr
    d_left = d_left - dtr
    
    if (d_left <= 0._dp) then
       d_exit = depth
       e_fin = etau2
       part_type = 0 ! went too far to make a tau, so don't count (HLS = 0); changed 22/12/2020
       ! go to 60; no, still a neutrino
       ! 60 continue
       Pout = Pin
       return
    end if

    ! we have a tau with room to travel for tauthrulayers

    call tau_thru_layers(angle, depth, d_water, d_lep, etau2, xc_water, xc_rock, ixc_water, ixc_rock,&
         & alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth, lepton, prop_type,&
         & part_type, d_exit, e_fin, Pi)

    Pout = Pi
    
    return

  end subroutine regen
  
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
       & u, w, no_regen_tot, regen_tot)
    !! Propagates a single ingoing neutrino event.
    
    implicit none
    
    real(dp), intent(in) :: energy
    !! Incoming neutrino energy, in GeV.
    real(dp), intent(in) :: angle
    !! Earth emergence angle (beta), in degrees.
    real(dp), intent(in) :: nu_xc(:,:)
    !! 2D array containing neutrino CC & NC cross-section values, in cm^2.
    real(dp), intent(in) :: nu_ixc(:,:,:)
    !! 3D array containing neutrino integrated cross-section CDF values.
    real(dp), intent(in) :: depth
    !! Maximum column depth for neutrino propagation, in kmwe.
    real(dp), intent(in) :: depthE
    !! Total column depth for a given Earth emergence angle, in kmwe.
    real(dp), intent(in) :: dwater
    !! Column depth along the chord for a given Earth emergence angle, in kmwe.
    real(dp), intent(in) :: xc_water(:,:)
    !! 2D array containing N_A/A*charged lepton-nucleon cross-section values in water, in cm^2/g.
    real(dp), intent(in) :: xc_rock(:,:)
    !! 2D array containing N_A/A*charged lepton-nucleon cross-section values in rock, in cm^2/g.
    real(dp), intent(in) :: lep_ixc_water(:,:,:)
    !! 3D array containing charged lepton integrated cross-section CDF values in water.
    real(dp), intent(in) :: lep_ixc_rock(:,:,:)
    !! 3D array containing charged lepton integrated cross-section CDF values in rock.
    real(dp), intent(in) :: alpha_water(:)
    !! 1D array containing ionization energy loss values in water, in (GeV*cm^2)/g.
    real(dp), intent(in) :: alpha_rock(:)
    !! 1D array containing ionization energy loss values in rock, in (GeV*cm^2)/g.
    real(dp), intent(in) :: beta_water(:,:)
    !! 2D array of beta values in water, in cm^2/g.
    real(dp), intent(in) :: beta_rock(:,:)
    !! 2D array of beta values in rock, in cm^2/g.
    real(dp), intent(in) :: xalong(:)
    !! 1D array containing distance in water, in km.
    real(dp), intent(in) :: cdalong(:)
    !! 1D array containing column depth at xalong, in g/cm^2.
    integer, intent(in) :: ithird
    !! Choice for neutrino -> charged lepton energy fraction selection.
    integer, intent(in) :: idepth
    !! Depth of water layer in km.
    integer, intent(in) :: lepton
    !! Type of charged lepton. 1=tau; 2=muon.
    real(dp), intent(in) :: fac_nu
    !! Rescaling factor for SM cross-sections.
    integer, intent(in) :: prop_type
    !! Type of energy loss propagation. 1=stochastic, 2=continuous.
    integer(dp), intent(in) :: u
    !! Filename for Etau_out character size.
    integer(dp), intent(in) :: w
    !! Filename for Polarization_out character size.

    integer(kind=8), intent(inout) :: no_regen_tot
    !! No. of outgoing leptons without regeneration.
    integer(kind=8), intent(inout) :: regen_tot
    !! No. of outgoing leptons with regeneration.
    
    real(dp) :: depth0, dtr, ef, etauin, dfinal, etauf, dleft, dtau2, ef2
    integer :: ip, ipp, ipp3
    integer(kind=8) :: regen_cnt
    real(dp) :: Pi, Pint, Pout
    CHARACTER(LEN=10) :: Format, Format1
    Format = "(F5.2)"
    Format1 = "(F8.5)"

    depth0 = 0.0_dp ! start with this each time
    
    ! 80 continue
    
    ! tnu goes until neutrino either goes to dtot, or converts to a tau
    !        print *,'depth=',depth
    call propagate_nu(energy, nu_xc, nu_ixc, depth, fac_nu, ip, dtr, ef)
    
    ! how far did the neutrino go? dtr is how far traveled
    
    depth0 = depth0 + dtr ! how far is the neutrino on trajectory?

    dleft = depth - depth0 ! how far is left for the neutrino to travel?
    
    if (ip == 0) then ! still a neutrino at the end of the road
       !print *, "First go is a neutrino"
       ! go to 10
       return ! break outside stat; continue is correct here
    end if

    ! continue here: we have a tau
    
    regen_cnt = 1 ! tau out after first interaction
    
    etauin = ef
    ! still need to propagate the tau, column depth to go
    
    call tau_thru_layers(angle, depth, dwater, depth0, etauin, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock, alpha_water,&
         & alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth, lepton, prop_type, ipp, dfinal, etauf, Pi)

    dleft = depth-dfinal
    
    if (ipp == 1 .and. dleft <= 0.0_dp) then ! a tau has emerged through column depth
       Pout = Pi
       no_regen_tot = no_regen_tot + 1
       regen_tot = regen_tot + 1 ! update the regen tau array once
       write(u, Format) dlog10(etauf)
       write(w, Format1) (Pout)
       ! go to 10; we are done with the loop
       return ! break outside stat; continue is correct here
    end if 

    ! 11 continue; beginning of regeneration loop
    ! must be a neutrino. Is there still column depth to propagate?
     
    ipp3 = 99 ! dummy value
    do while (dfinal < depthE .and. ipp3 /= 1 .and. regen_cnt <= 6) ! tau has decayed before the end
       
       etauin = etauf ! regen finds neutrino energy
       
       call regen(angle, etauin, depth, dwater, dfinal, nu_xc, nu_ixc, ithird, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,&
            & alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, idepth, lepton,&
            & fac_nu, prop_type, Pi, ipp3, dtau2, ef2, Pint)  

       !print *, "After regeneration, Pint = ", Pint
       !print *, "ipp3 = ",ipp3, " ipp3=1: Tau; ipp3=0: neutrino"
       
       regen_cnt = regen_cnt + 1
       
       if (ipp3 == 1) then ! then we are back to a tau at the end of the road
          regen_tot = regen_tot + 1
          Pout = Pint
          write(u, Format) dlog10(ef2)
          write(w, Format1) (Pout)
          ! go to 10; we are done with the loop
          return ! need to check if this breaks out of stat loop or not. Yes??
       end if

       if (regen_cnt > 6) then ! 6 rounds of regeneration
          !print *, "regeneration exceeded 6"
          return ! only if regen > 6, break and go to run_stat for next iteration
       end if
          
       etauf = ef2
       dfinal = dtau2 ! go to 11
       Pi = Pint
       
    end do
    
  end subroutine single_stat
  
  subroutine run_stat_single(energy, angle, nu_xc, nu_ixc, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,&
       & alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu, stats, prop_type,&
       & no_regen_tot, regen_tot)
    !! Run a loop for all ingoing neutrinos.
    
    implicit none
    
    real(dp), intent(in) :: energy
    !! Incoming neutrino energy, in GeV.
    real(dp), intent(in) :: angle
    !! Earth emergence angle (beta), in degrees.
    real(dp), intent(in) :: nu_xc(:,:)
    !! 2D array containing neutrino CC & NC cross-section values, in cm^2.
    real(dp), intent(in) :: nu_ixc(:,:,:)
    !! 3D array containing neutrino integrated cross-section CDF values.
    real(dp), intent(in) :: depthE
    !! Total column depth for a given Earth emergence angle, in kmwe.
    real(dp), intent(in) :: dwater
    !! Column depth along the chord for a given Earth emergence angle, in kmwe.
    real(dp), intent(in) :: xc_water(:,:)
    !! 2D array containing N_A/A*charged lepton-nucleon cross-section values in water, in cm^2/g.
    real(dp), intent(in) :: xc_rock(:,:)
    !! 2D array containing N_A/A*charged lepton-nucleon cross-section values in rock, in cm^2/g.
    real(dp), intent(in) :: lep_ixc_water(:,:,:)
    !! 3D array containing charged lepton integrated cross-section CDF values in water.
    real(dp), intent(in) :: lep_ixc_rock(:,:,:)
    !! 3D array containing lepton integrated cross-section CDF values in rock.
    real(dp), intent(in) :: alpha_water(:)
    !! 1D array containing ionization energy loss values in water, in (GeV*cm^2)/g.
    real(dp), intent(in) :: alpha_rock(:)
    !! 1D array containing ionization energy loss values in rock, in (GeV*cm^2)/g.
    real(dp), intent(in) :: beta_water(:,:)
    !! 2D array of beta values in water, in cm^2/g.
    real(dp), intent(in) :: beta_rock(:,:)
    !! 2D array of beta values in rock, in cm^2/g.
    real(dp), intent(in) :: xalong(:)
    !! 1D array containing distance in water, in km.
    real(dp), intent(in) :: cdalong(:)
    !! 1D array containing column depth at xalong, in g/cm^2.
    integer, intent(in) :: ithird
    !! Choice for neutrino -> charged lepton energy fraction selection.
    integer, intent(in) :: idepth
    !! Depth of water layer in km.
    integer, intent(in) :: lepton
    !! Type of charged lepton. 1=tau; 2=muon.
    real(dp), intent(in) :: fac_nu
    !! Rescaling factor for SM neutrino cross-sections.
    integer(kind=8), intent(in) :: stats
    !! Statistics or no. of ingoing neutrinos.
    integer, intent(in) :: prop_type
    !! Type of energy loss propagation. 1=stochastic, 2=continuous.

    integer(kind=8), intent(out) :: no_regen_tot
    !! No. of outgoing charged leptons without regeneration.
    integer(kind=8), intent(out) :: regen_tot
    !! No. of outgoing charged leptons with regeneration.
   
    real(dp) :: depth
    integer(kind=8):: i
    integer(kind=8) :: u, w
    character(25) Efilename, Pfilename
    
    write(Efilename,'(a,F0.2,a,F4.1,a)') 'eout_',dlog10(energy),'_',angle,'.dat' ! filename is eout_energy_angle.dat
    open(newunit=u, file=trim(Efilename), status="replace")

    write(Pfilename,'(a,F0.2,a,F4.1,a)') 'Pout_',dlog10(energy),'_',angle,'.dat' ! filename is Pout_energy_angle.dat
    open(newunit=w, file=trim(Pfilename), status="replace")
    
    depth = depthE
    no_regen_tot = 0
    regen_tot = 0
    
    !$OMP PARALLEL DO
    do i = 1, stats
       call single_stat(energy, angle, nu_xc, nu_ixc, depth, depthE, dwater, xc_water, xc_rock, lep_ixc_water, lep_ixc_rock,&
            & alpha_water, alpha_rock, beta_water, beta_rock, xalong, cdalong, ithird, idepth, lepton, fac_nu, prop_type, &
            & u,w, no_regen_tot, regen_tot)
       
    end do
    !$OMP END PARALLEL DO
    close(u)
    close(w)
    return
  end subroutine run_stat_single

end module run

  
