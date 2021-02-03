!
      subroutine pdrag(tn,te,ti,ne,op,o1,o2,n2,he,barm,diffs,Fe,Fn,
     | lev0,lev1,lon0,lon1,lat0,lat1)
! para ion-drag
! Jiuhou Lei
      use cons_module,only: evergs,gask,grav,dipmin,rmass_op,boltz
      use magfield_module,only: csdec,sndec,dipmag
      use oplus_module,only: rrk
      use addfld_module,only: addfld

#ifdef MPI
!      use mpi_module,only: mp_bndlons_f3d,! mp_bndlats_f3d, 
      use mpi_module,only: mp_bndlons_f3d, mp_periodic_f3d
#endif
      implicit none
!
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
! Input fields (full 3d task subdomain):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
! Lei
     |  diffs,       ! (D/(H*DZ)*2.*TP+M*G/R)*N(O+) (s7,s8,s9)
     |  tn, te, ti, ! neutral, electron, and ion temperatures (deg K)
     |  ne,         ! electron density
     |  op,         ! O+ ion
     |  o2, o1,     ! o2, o mass mixing ratios
     |  he,         ! he mass mixing ratio
     |  barm       ! mean molecular mass
!
! N2 (mmr) is input, but is also updated from new o2,o,he.
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) :: n2 
!
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) :: 
     |  Fe,      ! east-component of field line ion-drag;
     |  Fn       ! north-component of field line ion-drag;
! local
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  tpj,        ! (Te+Ti)=2*Tp
     |  Da,Wd,Vni     
! local       
      integer :: jj, k, lat, i
      real :: sini,cosi,sind,cosd
!
       do lat=lat0,lat1
          do i=lon0,lon1
            do k=lev0,lev1-1
               tpj(k,i,lat  ) = 0.5*(tn(k,i,lat  )+ti(k,i,lat  ))
            enddo
          enddo
        enddo
!
! diffs=d(2*n*tp)/(h*dz)+mgn/R
!
! calculte Vni=Vin*Ni*mi/rho 
      do lat=lat0,lat1
         call rrkvni(
     |   tn (:,lon0:lon1,lat),barm(:,lon0:lon1,lat),
     |   o2 (:,lon0:lon1,lat),o1(:,lon0:lon1,lat),
     |   tpj(:,lon0:lon1,lat),
     |   vni(:,lon0:lon1,lat),
     |   lon0,lon1,lev0,lev1,lat)
!
         do i=lon0,lon1
            Vni(:,i,lat)=Vni(:,i,lat)*ne(:,i,lat)
         enddo
!
! Da=(k/(mi*Vin);
       call rrk(
     |   tn (:,lon0:lon1,lat),barm(:,lon0:lon1,lat),
     |   o2 (:,lon0:lon1,lat),o1(:,lon0:lon1,lat),
     |   n2(:,lon0:lon1,lat),he(:,lon0:lon1,lat),
     |   tpj(:,lon0:lon1,lat),
     |   Da(:,lon0:lon1,lat),
     |   lon0,lon1,lev0,lev1,lat)
! Lei
       call addfld('Da',' ',' ',Da(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
       call addfld('DIFFS',' ',' ',diffs(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
       call addfld('VNI',' ',' ',Vni(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

      enddo
!
! compute parr ion-drag  including pressure gradient and gravite contributions
      do lat=lat0,lat1
          do i=lon0,lon1
!
             sind=sndec(i,lat)
	     cosd=csdec(i,lat)
	     sini=sin(dipmag(i,lat))
	     cosi=cos(dipmag(i,lat))
! Wd=diffs*Da/Ne; southward positive
             Wd(:,i,lat)=-diffs(:,i,lat)*Da(:,i,lat)/ne(:,i,lat)*sini
! drag=-Vni*Wd,download/northward posive and then have postive component in northward;
! drag forcing: eastward
	     Fe(:,i,lat)=-Vni(:,i,lat)*Wd(:,i,lat)*sind*cosi
! northward                       
	     Fn(:,i,lat)=-Vni(:,i,lat)*Wd(:,i,lat)*cosd*cosi
!             
             Fe(lev1,i,lat)=1.5* Fe(lev1-1,i,lat)-0.5* Fe(lev1-2,i,lat)
             Fn(lev1,i,lat)=1.5* Fn(lev1-1,i,lat)-0.5* Fn(lev1-2,i,lat)
!                        
	    enddo
! Lei
       call addfld('WD',' ',' ',Wd(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)

       call addfld('FE',' ',' ',Fe(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
       call addfld('FN',' ',' ',Fn(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
       enddo
!
      end subroutine pdrag

!-----------------------------------------------------------------------
      subroutine rrkvni(t,rms,ps1,ps2,tp,ans,lon0,lon1,lev0,lev1,lat)
!
! Returns Vni/ne in ans.
! Vin=Vin*Ni*mi/rho
! Jiuhou lei
!
      use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,boltz,
     |  p0,expz
      use input_module,only: colfac
!
! Args:
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat
      real,dimension(lev0:lev1,lon0:lon1),intent(in) ::
     |  t,rms,ps1,ps2,tp
      real,dimension(lev0:lev1,lon0:lon1),intent(out) :: ans
!
! Local:
      integer :: k,i
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! rho*Vni=sum(Vin*Ni*mi);
! and then simpified as Vni=Vo+-n*Ne*mo+/rho
! ans=mo+*Vo+-n (without Ne)
! atom mass is not used because it is cancaled by that in rho;
          ans(k,i) = 16*3.53E-11*(p0*expz(k)*.5*(rms(k,i)+
     |      rms(k+1,i))*(ps2(k,i)*rmassinv_o1*sqrt(tp(k,i))*(1.-0.064*
     |      alog10(tp(k,i)))**2*colfac+18.6*(1.-ps1(k,i)-ps2(k,i))*
     |      rmassinv_n2+18.1*ps1(k,i)*rmassinv_o2))/(boltz*t(k,i))
!
! Vni=Vin*Ni*mi/rho, where rho=mbar*P/kT
          ans(k,i) = ans(k,i)*(boltz*t(k,i))/(p0*expz(k)*.5*(rms(k,i)+
     |      rms(k+1,i)))
!
        enddo ! k=lev0,lev1
!
        ans(lev1,i)=1.5*ans(lev1-1,i)-0.5*ans(lev1-2,i)
!
      enddo ! i=lon0,lon1
      end subroutine rrkvni
