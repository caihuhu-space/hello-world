!
      subroutine newton(tn,o2,o1,n2,no,barm,cp,cool_implicit,
     |  cool_explicit,lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate implicit and explicit cooling terms.
!
      use params_module,only: dz
      use cons_module,only: rmass_co2,p0,expz,expzmid_inv,boltz,
     |  avto,avo,rmassinv_o2,rmassinv_o1,rmassinv_n2,rmassinv_no,
     |  gask
      use addfld_module,only: addfld
      use diags_module,only: mkdiag_CO2COOL,mkdiag_NOCOOL
      use init_module,only: sfeps,iday,iyear,secs,istep
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o2,      ! molecular oxygen (mmr)
     |  o1,      ! atomic oxygen (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  no,      ! nitric oxide (mmr)
     |  cp,      ! specific heat
     |  barm     ! mean molecular weight
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  cool_implicit, cool_explicit ! output
!
! Local:
      real :: co2u
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  xmco2,
     |  nco2,aco2,bco2,  ! n(co2), a(co2), b(co2) (s10,s11,s12)
     |  ano,             ! a(no) (s11)
     |  co2_cool,        ! co2 contribution to cooling
     |  no_cool          ! no contribution to cooling
      integer :: k,i
      real :: uttime,yfrac
!
! set CO2 concentration
! calculate the fractional year value for the model time
      uttime= secs/3600.
      if (mod(iyear,4) == 0) then
         yfrac=iyear+(iday-1+(uttime/24.))/366.
      else
         yfrac=iyear+(iday-1+(uttime/24.))/365.
      endif
      if (yfrac > 1953.0) then 
         co2u=248.0+30.0*exp(0.014*(yfrac-1900.))
      else
         co2u=280.0+0.30*(yfrac-1850.0)
      endif
      if (co2u < 280.0) co2u=280.0
      co2u=co2u*1.e-6
!
! Integral(m(co2)*dz/mbar) (s10).
! Low boundary: 
      do i=lon0,lon1
        xmco2(1,i) = .25*rmass_co2*(1./barm(1,i)+1./barm(2,i))*dz
!
! Integrate upwards:
        do k=lev0+1,lev1-1
          xmco2(k,i) = xmco2(k-1,i)+rmass_co2/barm(k,i)*dz
        enddo ! k=lev0+1,lev1-1
!
! n(co2) (s10):
        do k=lev0,lev1-1
          nco2(k,i) = co2u*p0*expz(1)*expzmid_inv/(boltz*tn(k,i))*
     |      exp(-xmco2(k,i))
!
! a(co2) (s11):
          if (tn(k,i) >= 200.) then
            aco2(k,i) = 2.5e-15*(1.+0.03*(tn(k,i)-200.))
          else
            aco2(k,i) = 2.5e-15
          endif
!
! b(co2) (s12):
          if (tn(k,i) < 260.) then
             bco2(k,i) = 1.56e-12
          elseif (tn(k,i) >= 260..and.tn(k,i) <= 300.) then
             bco2(k,i) = (2.6-0.004*tn(k,i))*1.0e-12
          else
             bco2(k,i) = 1.4e-12
          endif
!
! Co2 cooling (s12):
          co2_cool(k,i) = 2.65e-13*nco2(k,i)*exp(-960./tn(k,i))*
     |      avo*((o2(k,i)*rmassinv_o2 + n2(k,i)*rmassinv_n2)*
     |      aco2(k,i)+o1(k,i)*rmassinv_o1*bco2(k,i))
!
! a(no) (s11):
          ano(k,i) = 
     |      .5*(barm(k,i)+barm(k+1,i))*avo*p0*expz(k)/(gask*tn(k,i))*
     |      (4.2e-11*o1(k,i)*rmassinv_o1+
     |       2.4e-14*o2(k,i)*rmassinv_o2)
!
! No cooling (s11):
          no_cool(k,i) = 4.956e-12*(avo*no(k,i)*rmassinv_no)*
     |      (ano(k,i)/(ano(k,i)+13.3))*exp(-2700./tn(k,i))
!
! Total cooling terms.
          cool_implicit(k,i) = (2700.*no_cool(k,i)+960.*co2_cool(k,i))/
     |      (.5*(cp(k,i)+cp(k+1,i))*tn(k,i)**2)
          cool_explicit(k,i) = ((1.-2700./tn(k,i))*no_cool(k,i)+
     |      (1.-960./tn(k,i))*co2_cool(k,i))*expz(k)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
!     call addfld('XMCO2',' ',' ',xmco2(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('CO2_COOL',' ',' ',co2_cool(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('NO_COOL' ,' ',' ',no_cool(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('COOL_IMP' ,' ',' ',cool_implicit(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('COOL_EXP' ,' ',' ',cool_explicit(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Save as diagnostics on secondary histories (define top level as k-1):
      co2_cool(lev1,:) = co2_cool(lev1-1,:)
      no_cool (lev1,:) = no_cool (lev1-1,:)
      call mkdiag_CO2COOL('CO2_COOL',co2_cool,lev0,lev1,lon0,lon1,lat)
      call mkdiag_NOCOOL ('NO_COOL' ,no_cool ,lev0,lev1,lon0,lon1,lat)

      end subroutine newton
!-----------------------------------------------------------------------
#include <defs.h>
      subroutine newto3p(tn,o1,cp,cool_implicit,cool_explicit,
     |  lev0,lev1,lon0,lon1,lat)
!
! Add O(3P) cooling to implicit and explicit cooling terms.
! This is called from dynamics after newton.
!
      use params_module,only: dz,nlevp1,zibot
      use cons_module,only: avo,rmassinv_o1,expz
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o1,      ! atomic oxygen (mmr)
     |  cp       ! specific heat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(inout) ::
     |  cool_implicit,
     |  cool_explicit  
!
! Local:
      integer :: k,i
!
! Change O(3P) cooling factor from 0.5 to 1.0 of the Bates expression.
      real,parameter :: 
     |  an(3) = (/0.835E-18,0.6,0.2/),
     |  bn(3) = (/228.,228.,325./)
      real :: xfac(lev0:lev1)
      real,dimension(lev0:lev1,lon0:lon1) :: fac1,fac2
      integer::ih
!
!#if (NLEV==32)
      if (dz==0.5) then
        xfac(1:29) = (/
     |  1.0000e-02, 1.0000e-02, 1.0000e-02, 5.0000e-02, 1.0000e-01,
     |  2.0000e-01, 4.0000e-01, 5.5000e-01, 7.0000e-01, 7.5000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01/)
! extending the xfac when ZITOP > 7 (ZITOP=7,dz=0.5,nlev=28) --Yihui Cai, Oct. 14, 2020
        xfac(29:nlevp1) = 8.0000e-01 
!#elif (NLEV==64)
      elseif (dz==0.25) then
        xfac(1:57) = (/
     |  1.0000e-02, 1.0000e-02, 1.0000e-02, 1.0000e-02, 1.0000e-02,
     |  3.0000e-02, 5.0000e-02, 7.5000e-02, 1.0000e-01, 1.5000e-01,
     |  2.0000e-01, 3.0000e-01, 4.0000e-01, 4.7500e-01, 5.5000e-01,
     |  6.2500e-01, 7.0000e-01, 7.2500e-01, 7.5000e-01, 7.7500e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01/)
! extending the xfac when ZITOP > 7 (ZITOP=7,dz=0.25,nlev=56) --Yihui Cai, Oct. 14, 2020
        xfac(57:nlevp1) = 8.0000e-01 
!#elif (NLEV==128)
      elseif (dz==0.125) then
        xfac(1:113) = (/
     |  1.0000e-02, 1.0000e-02, 1.0000e-02, 1.0000e-02, 1.0000e-02,
     |  1.0000e-02, 1.0000e-02, 1.0000e-02, 1.0000e-02, 2.0000e-02,
     |  3.0000e-02, 4.0000e-02, 5.0000e-02, 6.2500e-02, 7.5000e-02,
     |  8.7500e-02, 1.0000e-01, 1.2500e-01, 1.5000e-01, 1.7500e-01,
     |  2.0000e-01, 2.5000e-01, 3.0000e-01, 3.5000e-01, 4.0000e-01,
     |  4.3750e-01, 4.7500e-01, 5.1250e-01, 5.5000e-01, 5.8750e-01,
     |  6.2500e-01, 6.6250e-01, 7.0000e-01, 7.1250e-01, 7.2500e-01,
     |  7.3750e-01, 7.5000e-01, 7.6250e-01, 7.7500e-01, 7.8750e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01/)
! extending the xfac when ZITOP > 7 (ZITOP=7,dz=0.125,nlev=112) --Yihui Cai, Oct. 14, 2020
        xfac(113:nlevp1) = 8.0000e-01 
!#elif (NLEV==256)
      elseif (dz==0.0625) then
        xfac(1:225) = (/
     |  1.0000e-02, 1.0000e-02, 1.0000e-02, 1.0000e-02, 1.0000e-02,
     |  1.0000e-02, 1.0000e-02, 1.0000e-02, 1.0000e-02, 1.0000e-02,
     |  1.0000e-02, 1.0000e-02, 1.0000e-02, 1.0000e-02, 1.0000e-02,
     |  1.0000e-02, 1.0000e-02, 1.5000e-02, 2.0000e-02, 2.5000e-02,
     |  3.0000e-02, 3.5000e-02, 4.0000e-02, 4.5000e-02, 5.0000e-02,
     |  5.6250e-02, 6.2500e-02, 6.8750e-02, 7.5000e-02, 8.1250e-02,
     |  8.7500e-02, 9.3750e-02, 1.0000e-01, 1.1250e-01, 1.2500e-01,
     |  1.3750e-01, 1.5000e-01, 1.6250e-01, 1.7500e-01, 1.8750e-01,
     |  2.0000e-01, 2.2500e-01, 2.5000e-01, 2.7500e-01, 3.0000e-01,
     |  3.2500e-01, 3.5000e-01, 3.7500e-01, 4.0000e-01, 4.1875e-01,
     |  4.3750e-01, 4.5625e-01, 4.7500e-01, 4.9375e-01, 5.1250e-01,
     |  5.3125e-01, 5.5000e-01, 5.6875e-01, 5.8750e-01, 6.0625e-01,
     |  6.2500e-01, 6.4375e-01, 6.6250e-01, 6.8125e-01, 7.0000e-01,
     |  7.0625e-01, 7.1250e-01, 7.1875e-01, 7.2500e-01, 7.3125e-01,
     |  7.3750e-01, 7.4375e-01, 7.5000e-01, 7.5625e-01, 7.6250e-01,
     |  7.6875e-01, 7.7500e-01, 7.8125e-01, 7.8750e-01, 7.9375e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01,
     |  8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01, 8.0000e-01/)
! extending the xfac when ZITOP > 7 (ZITOP=7,dz=0.0625,nlev=224) --Yihui Cai, Oct. 14, 2020
        xfac(225:nlevp1) = 8.0000e-01 
      else
!#else
        write(6,"('>>> WARNING newto3p: unsupported dz=',f10.3)") dz
        call shutdown('newto3p')
      endif
!#endif
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!       
! fac1 = AN(1)*XFAC*N0*PSI2/M2*EXP(-BN(1)/T) (s11)
          fac1(k,i) = an(1)*.5*(xfac(k)+xfac(k+1))*avo*rmassinv_o1
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          fac1(k,i) = fac1(k,i)*o1(k,i)*exp(-bn(1)/tn(k,i)) ! s11
!
! fac2 = 1.+AN(2)*EXP(-BN(2)/T)+AN(3)*EXP(-BN(3)/T) (s12)

          fac2(k,i) = 1.+an(2)*exp(-bn(2)/tn(k,i))+
     |                   an(3)*exp(-bn(3)/tn(k,i))
          fac1(k,i) = fac1(k,i)/fac2(k,i)
!
! fac2 = DF/DT=S11/S12/T**2*(BN(1)+(BN(1)-BN(3))*AN(3)*EXP(-BN(3)/T)) (s12)
          fac2(k,i) = fac1(k,i)/fac2(k,i)/tn(k,i)**2*
     |      (bn(1)+(bn(1)-bn(3))*an(3)*exp(-bn(3)/tn(k,i)))
!
! Add implicit contributions to cool_implicit:
          cool_implicit(k,i) = cool_implicit(k,i)+fac2(k,i)/
     |      (.5*(cp(k,i)+cp(k+1,i)))
!
! Add explicit contributions to cool_explicit:
          fac1(k,i) = fac1(k,i)-tn(k,i)*fac2(k,i)
          cool_explicit(k,i) = cool_explicit(k,i)+fac1(k,i)*expz(k)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
!     call addfld('COOLIMP' ,' ',' ',cool_implicit(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('COOLEXP' ,' ',' ',cool_explicit(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      end subroutine newto3p
