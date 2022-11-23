module coare36_parameters

  use params, only: &
       Rgas => rgas, & ! gas constant air ( J/kg/K )
       cpa => cp, & ! spec heat air at p ( J/kg/K )
       pi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1  

  public

  integer, parameter :: kind_phys = KIND(1.)

  real (kind=kind_phys), parameter :: &
    !%***********  set constants **********************************************
    zref = 10., &
    Beta = 1.2, &
    von  = 0.4, &
    fdg  = 1.00, & !% Turbulent Prandtl number
    T2K  = 273.15, &
    sigma_r = 5.67e-8, & ! Stephan-Boltzmann constant in W/m2/K4
    !%***********  air constants **********************************************
    eps = 0.622, &
    !%***********  cool skin constants  ***************************************
    !%%% includes salinity dependent thermal expansion coeff for water
    !%%%%%%%%%%%%%%%%%%%
    bets = 7.5e-4, & !% salintity expansion coeff; assumes beta independent of
                     !  temperature
    !%%%%  see "Computing the seater expansion coefficients directly from the
    !%%%%  1980 equation of state".  J. Lillibridge, J.Atmos.Oceanic.Tech, 1980.
    cpw  = 4000., &
    rhow = 1022., &
    visw = 1.e-6, &
    tcw  = 0.6, &
    dT_skin_min = 0.001, & ! minimum value for cool-skin temperature
                           ! depression change
                           ! use it for abs(dT_skin-dT_skin_pre) criterion
                           ! for especaping loop

    !bloss(TODO): Should these be linked to the radiation scheme and perhaps
    !   allow a dirnal cycle in sw albedo??
    lw_emiss = 0.97, &  ! longwave emissivity
    sw_onema = 0.945    ! 1-[shortwave albedo]

end module coare36_parameters


