subroutine nudging()
	
use vars
use params
use microphysics, only: micro_field, index_water_vapor, nmicro_fields, &
     flag_nudging, mknudge, mk0, mkobs
implicit none

real coef, coef1
integer i,j,k,n
	
!bloss: nudging that adjusts to inversion height
real :: nudge_ramp(nzm), nudge_test(nzm)
real :: pii

!bloss: option to apply strong nudging for a limited time period
real :: itau_transient_z(nzm), time_ramp

real :: itauz_t(nzm), itauz_q(nzm)
real :: tmp_zinv_obs
real :: dist_zinv

! tendency for CGILS qfloor
!   This is a strictly positive nudging that prevents the above-inversion
!   air from becoming unrealistically dry (mostly due to horizontal advective
!   drying that is applied above the inversion).
real dqdt_qfloor

real, external :: get_inversion_height

call t_startf ('nudging')

tnudge = 0.
qnudge = 0.
unudge = 0.
vnudge = 0.
mknudge = 0.

if(donudging_transient) then
  ! initialize transient nudging strength (i.e., 1/tau_transient_z(k)) to zero
  itau_transient_z(:) = 0.

  if((day.gt.transient_nudging_start).AND.(day.lt.transient_nudging_end)) then
    ! switch on nudging of moisture/temperature
    donudging_tq = .true. ! switch

    ! compute profile of transient nudging
    tmp_zinv_obs = get_inversion_height(nzm,z,pres,tg0+gamaz,tg0,qg0) ! inversion height in meters

    call define_transient_nudging(nzm,z,tmp_zinv_obs, & 
                                tau_transient_nudging, & ! nudging timescale below inversion (sec)
                                600., & ! nudging timescale above inversion (sec)
                                itau_transient_z)

    !smooth start/finish to transient nudging
    time_ramp = 0.5*(1.-cos(pi*MAX(0., MIN(1., (day-transient_nudging_start)/transient_nudging_ramp ) ) ) ) &
         *0.5*(1.-cos(pi*MAX(0., MIN(1., (transient_nudging_end-day)/transient_nudging_ramp ) ) ) ) 

    itau_transient_z(:) = itau_transient_z(:)*time_ramp

    if(icycle.eq.1.AND.mod(nstep,10).eq.1) then
      if(masterproc) write(*,*) 'Transient nudging on with itau_BL(day) = ', itau_transient_z(1), &
           ' itau_aloft(day) = ', itau_transient_z(nzm)
    end if

  end if

end if

coef = 1./tauls

if(donudging_uv) then
  if(nudge_to_sounding_winds) then
    do k=1,nzm
      if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
        unudge(k)=unudge(k) - (u0(k)-usounding0(k))*coef
        vnudge(k)=vnudge(k) - (v0(k)-vsounding0(k))*coef
        do j=1,ny
          do i=1,nx
             dudt(i,j,k,na)=dudt(i,j,k,na)-(u0(k)-usounding0(k))*coef
             dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vsounding0(k))*coef
          end do
        end do
      end if
    end do
  else
    do k=1,nzm
      if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
        unudge(k)=unudge(k) - (u0(k)-ug0(k))*coef
        vnudge(k)=vnudge(k) - (v0(k)-vg0(k))*coef
        do j=1,ny
          do i=1,nx
             dudt(i,j,k,na)=dudt(i,j,k,na)-(u0(k)-ug0(k))*coef
             dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vg0(k))*coef
          end do
        end do
      end if
    end do
  end if
endif

if(donudging_tq.or.donudging_t.or.donudging_q) then

  ! set up nudging heights automatically by tracking inversion height
  if(dovariable_tauz) then
    nudging_t_z1 = get_inversion_height(nzm,z,pres,t0,tabs0,q0) &
         + variable_tauz_offset_above_inversion
    nudging_t_z1 = MAX(nudging_t_z1, variable_tauz_minimum_height)
    nudging_t_zramp = variable_tauz_thickness_of_onset

    nudging_q_z1 = nudging_t_z1
    nudging_q_zramp = nudging_t_zramp
  end if

  ! vertically-varying nudging amplitude for temperature
  call define_nudge_ramp(nzm,z,nudging_t_z1,nudging_t_z2, &
                             nudging_t_zramp,nudge_ramp)

  !bloss: Use itauz_t(1:nzm) to keep track of vertically-varying inverse nudging timescale
  coef = 1./tautqls
  itauz_t(:) = coef*nudge_ramp(:)

  ! vertically-varying nudging amplitude for moisture
  call define_nudge_ramp(nzm,z,nudging_q_z1,nudging_q_z2, &
                               nudging_q_zramp,nudge_ramp)

  !bloss: Use itauz_q(1:nzm) to keep track of vertically-varying inverse nudging timescale
  coef = 1./tautqls
  itauz_q(:) = coef*nudge_ramp(:)


  if(donudging_transient) then
    ! If transient nudging is enabled, choose strongest nudging at each level.
    do k = 1,nzm
      itauz_t(k) = MAX(itauz_t(k), itau_transient_z(k))
      itauz_q(k) = MAX(itauz_q(k), itau_transient_z(k))
    end do
  end if

end if


if(donudging_tq.or.donudging_t) then
    do k=1,nzm
!bloss: Apply nudging for all levels -- itauz_t(k) will be zero where no nudging is applied.
        tnudge(k)=tnudge(k) -(t0(k)-tg0(k)-gamaz(k))*itauz_t(k)

        do j=1,ny
          do i=1,nx
            t(i,j,k)=t(i,j,k)-(t0(k)-tg0(k)-gamaz(k))*dtn*itauz_t(k)
          end do
        end do
    end do
endif

if(donudging_tq.or.donudging_q) then
    do k=1,nzm
!bloss: Apply nudging for all levels -- itauz_t(k) will be zero where no nudging is applied.
      qnudge(k)=qnudge(k) -(q0(k)-qg0(k))*itauz_q(k)
      do j=1,ny
        do i=1,nx
          micro_field(i,j,k,index_water_vapor)=micro_field(i,j,k,index_water_vapor)-(q0(k)-qg0(k))*dtn*itauz_q(k)
        end do
      end do
    end do

    do n = 1,nmicro_fields
      if((n.ne.index_water_vapor).AND.(flag_nudging(n).eq.1)) then
        !bloss: Nudge isotopic water vapor in free troposphere
        do k=1,nzm
          mknudge(k,n) = mknudge(k,n) - (mk0(k,n)-mkobs(k,n))*itauz_q(k)*float(nx*ny)
          do j=1,ny
            do i=1,nx
              micro_field(i,j,k,n) = &
                   micro_field(i,j,k,n) - (mk0(k,n)-mkobs(k,n))*dtn*itauz_q(k)
            end do ! do i = 1,nx
          end do ! do j = 1,ny
        end do ! do k = 1,nzm
      endif ! if(flag_nudging(n).eq.1).AND. not index_water_vapor
    end do ! do n = 1,nmicro_fields

endif

!bloss: option for CGILS cases to avoid over-drying of sounding above BL
if(doenforce_cgils_qfloor) then
  do k = 1,nzm
    if((z(k).lt.ztop_qfloor).AND.(q0(k).lt.qfloor)) then
      dqdt_qfloor = (qfloor-q0(k))/tau_qfloor
      qnudge(k) = qnudge(k) + dqdt_qfloor
      do j=1,ny
        do i=1,nx
          micro_field(i,j,k,index_water_vapor) &
               = micro_field(i,j,k,index_water_vapor) + dtn*dqdt_qfloor
        end do
      end do
    end if
  end do
end if

mknudge(:,index_water_vapor) = qnudge(:)*float(nx*ny)

call t_stopf('nudging')

end subroutine nudging

subroutine define_nudge_ramp(N,xi,xi1,xi2,xi_ramp,nudge_ramp)
  implicit none
 
  integer, intent(in) :: N
  real, intent(in) :: xi(N), xi1, xi2, xi_ramp
  real, intent(out) :: nudge_ramp(N)

  integer :: k
  real :: pii

  ! define function that is 0 where there is no nudging 
  !   and increases smoothly to 1 in regions of full nudging.
  nudge_ramp(:) = 0.

  pii = acos(-1.)
  do k=1,N
    if(xi(k).ge.xi1.and.xi(k).le.xi2) then
      !nudging will be applied
      nudge_ramp(k) = 1.

      if(xi(k).lt.xi1+xi_ramp) then
        ! gradual onset of nudging between z1 and z1+zramp
        nudge_ramp(k) = 0.5*(1-cos(pii*(xi(k)-xi1) / xi_ramp ) )
      elseif(xi(k).gt.xi2 - xi_ramp) then
        ! gradual falloff of nudging between z2 and z2-zramp
        nudge_ramp(k) = 0.5*(1-cos(pii*(xi2 - xi(k)) / xi_ramp ) )
      end if
    end if
  end do

end subroutine define_nudge_ramp
  
subroutine define_transient_nudging(N,xi,xi0,tau_bl,tau_aloft,itau_z)
  implicit none
 
  integer, intent(in) :: N
  real, intent(in) :: xi(N), xi0, tau_bl, tau_aloft
  real, intent(out) :: itau_z(N)

  integer :: k
  real :: pii, nudge_ramp(N), dist_xi0

  pii = acos(-1.)

  nudge_ramp(:) = 1.  
  do k = 1,N
    dist_xi0 = ABS(xi(k)-xi0) 
    IF (dist_xi0.lt.50.) then
      ! no nudging within 50m of xi0
      nudge_ramp(k) = 0.
    elseif (dist_xi0.lt.100.) then
      ! ramp up to full nudging between 50 and 100m from xi0
      nudge_ramp(k) = 0.5*(1. - cos(pii*(dist_xi0-50.)/50.))
    end IF
    if(xi(k).lt.xi0) then
      ! specified nudging timescale below xi0 (i.e., within boundary layer)
      itau_z(k) = nudge_ramp(k)/tau_bl
    else
      ! different, shorter nudging time above inversion
      itau_z(k) = nudge_ramp(k)/tau_aloft
    end if
  end do

end subroutine define_transient_nudging
  
