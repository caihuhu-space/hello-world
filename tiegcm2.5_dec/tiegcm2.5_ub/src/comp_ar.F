!
      module ar_module
!
! Advance argon by one time step.
!
! Boundary conditions, production and loss for argon are defined
! by comp_ar, and referenced by minor_ar. Comp_ar is called
! from a latitude loop in dynamics. After comp_ar, dynamics calls
! minor_ar, which passes this module data to sub minor. Sub
! minor contains 3d mpi calls and its own latitude loops. 
!
      use params_module,only: nlevp1
      use addfld_module,only: addfld
      implicit none
!
! Boundary conditions and production and loss terms are allocated
! on subdomains by sub alloc_ar (called from allocdata.F).
!
      real,allocatable,dimension(:,:)   :: ar_ubc ! upper boundary (i,j)
      real,allocatable,dimension(:,:,:) :: ar_lbc ! lower boundary (i,3,j)
      real,allocatable,dimension(:,:,:) ::
     |  ar_prod,     ! production of argon (k,i,j)
     |  ar_loss      ! loss of argon       (k,i,j)
      real :: ar_glbm(nlevp1)
!
      contains
!-----------------------------------------------------------------------
      subroutine comp_ar(ar,lev0,lev1,lon0,lon1,lat)
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  ar           ! argon from previous step
!
! Local:
      integer :: i
!
! ar_lbc, ar_ubc, ar_glbm, ar_prod, ar_loss are module data above.
!
      do i=lon0,lon1 
        ar_lbc(i,1,lat) = 0.
        ar_lbc(i,2,lat) = 1.
        ar_lbc(i,3,lat) = -sqrt(ar_glbm(1)*ar_glbm(2))
        ar_ubc(i,lat)   = 0.
      enddo
!     write(6,"('comp_ar: lat=',i3,' ar_lbc(3)=',/,(6e12.4))") 
!    |  lat,ar_lbc(:,3,lat)
!
! Production and loss of argon are zero:
      ar_prod = 0.
      ar_loss = 0.

      end subroutine comp_ar
!-----------------------------------------------------------------------
      subroutine minor_ar(tn,o2,o1,n2,he,ar,ar_nm,ar_out,arnm_out,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
      use cons_module,only: rmass_ar
      implicit none
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  tn,     ! neutral temperature (deg K)
     |  o2,     ! molecular oxygen (mmr)
     |  o1,     ! oxygen family (mmr)
     |  n2,     ! molecular nitrogen (mmr)
     |  he,     ! helium (mmr)
     |  ar,     ! argon (mmr)
     |  ar_nm   ! argon at time n-1
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  ar_out,    ! ar output
     |  arnm_out   ! ar output at time n-1
!
! Local:
      real,parameter :: phi_ar(3) = (/1.042,1.509,1.176/) ! from tgcm24
      real,parameter :: alfa_ar   = 0.17    ! thermal diffusion coefficient (from tgcm24)
      real,parameter :: xyar      = 1.e-10  ! from tgcm24
!
      call minor(tn,o2,o1,n2,he,ar,ar_nm,ar_out,arnm_out,ar_loss,
     |  ar_prod,ar_lbc,ar_ubc,rmass_ar,phi_ar,alfa_ar,lev0,lev1,
     |  lon0,lon1,lat0,lat1,0,'AR')
!
      end subroutine minor_ar
!-----------------------------------------------------------------------
#include <defs.h>
      subroutine alloc_ar(lon0,lon1,lat0,lat1)
!
! Allocate subdomains (without ghost cells) to module data for boundary
! conditions and production and loss terms. This is called once per run
! from sub allocdata (allocdata.F).
!
      use params_module,only: nlev, dz,zibot,nlevp1
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1
!
! Local:
      integer :: istat,ih
      real::ind
!
! Allocate subdomains to boundary conditions:
      allocate(ar_ubc(lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_ar: error allocating',
     |  ' ar_ubc: stat=',i3)") istat
      allocate(ar_lbc(lon0:lon1,3,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_ar: error allocating',
     |  ' ar_lbc: stat=',i3)") istat
      ar_ubc = 0. ; ar_lbc = 0.
!
! Allocate subdomains to production and loss:
      allocate(ar_prod(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_ar: error allocating',
     |  ' ar_prod: stat=',i3)") istat
      allocate(ar_loss(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_ar: error allocating',
     |  ' ar_loss: stat=',i3)") istat
      ar_prod = 0. ; ar_loss = 0.
!
      write(6,"('alloc_ar: allocated module data')")
!
!#if (NLEV==32)
      if (dz==0.5) then
        ar_glbm(1:29) = (/
     |  1.1640e-02, 1.0740e-02, 9.4430e-03, 7.9050e-03, 6.3420e-03,
     |  4.9270e-03, 3.7370e-03, 2.7760e-03, 2.0170e-03, 1.4300e-03,
     |  9.8550e-04, 6.5750e-04, 4.2290e-04, 2.6130e-04, 1.5480e-04,
     |  8.7970e-05, 4.8060e-05, 2.5350e-05, 1.2980e-05, 6.4880e-06,
     |  3.1840e-06, 1.5410e-06, 7.3810e-07, 3.5110e-07, 1.6620e-07,
     |  7.8400e-08, 3.6890e-08, 1.7330e-08, 7.4360e-09/)
! extending the ar_glbm when ZITOP > 7 (ZITOP=7,dz=0.5,nlev=28) --Yihui Cai, Oct. 14, 2020
        do ih=29,nlevp1
           ind = zibot+(ih-1)*dz
           ar_glbm(ih) = 10**(-0.0255*ind**3+0.4354*ind**2
     |                   -3.1191*ind+1.1204)
        enddo
        ar_glbm(29) = 7.4360e-09
!#elif (NLEV==64)
      elseif (dz==0.25) then
        ar_glbm(1:57) = (/
     |  1.1640e-02, 1.1190e-02, 1.0740e-02, 1.0091e-02, 9.4430e-03,
     |  8.6740e-03, 7.9050e-03, 7.1235e-03, 6.3420e-03, 5.6345e-03,
     |  4.9270e-03, 4.3320e-03, 3.7370e-03, 3.2565e-03, 2.7760e-03,
     |  2.3965e-03, 2.0170e-03, 1.7235e-03, 1.4300e-03, 1.2078e-03,
     |  9.8550e-04, 8.2150e-04, 6.5750e-04, 5.4020e-04, 4.2290e-04,
     |  3.4210e-04, 2.6130e-04, 2.0805e-04, 1.5480e-04, 1.2139e-04,
     |  8.7970e-05, 6.8015e-05, 4.8060e-05, 3.6705e-05, 2.5350e-05,
     |  1.9165e-05, 1.2980e-05, 9.7340e-06, 6.4880e-06, 4.8360e-06,
     |  3.1840e-06, 2.3625e-06, 1.5410e-06, 1.1395e-06, 7.3810e-07,
     |  5.4460e-07, 3.5110e-07, 2.5865e-07, 1.6620e-07, 1.2230e-07,
     |  7.8400e-08, 5.7645e-08, 3.6890e-08, 2.7110e-08, 1.7330e-08,
     |  1.2383e-08, 7.4360e-09/)
! extending the zr_glbm when ZITOP > 7 (ZITOP=7,dz=0.25,nlev=56) --Yihui Cai, Oct. 14, 2020
        do ih=57,nlevp1
           ind = zibot+(ih-1)*dz
           ar_glbm(ih) = 10**(-0.0255*ind**3+0.4354*ind**2
     |                   -3.1191*ind+1.1204)
        enddo
        ar_glbm(57) = 7.4360e-09
!#elif (NLEV==128)
      elseif (dz==0.125) then
        ar_glbm(1:113) = (/
     |  1.1640e-02, 1.1415e-02, 1.1190e-02, 1.0965e-02, 1.0740e-02,
     |  1.0416e-02, 1.0091e-02, 9.7672e-03, 9.4430e-03, 9.0585e-03,
     |  8.6740e-03, 8.2895e-03, 7.9050e-03, 7.5143e-03, 7.1235e-03,
     |  6.7328e-03, 6.3420e-03, 5.9883e-03, 5.6345e-03, 5.2807e-03,
     |  4.9270e-03, 4.6295e-03, 4.3320e-03, 4.0345e-03, 3.7370e-03,
     |  3.4967e-03, 3.2565e-03, 3.0162e-03, 2.7760e-03, 2.5863e-03,
     |  2.3965e-03, 2.2067e-03, 2.0170e-03, 1.8703e-03, 1.7235e-03,
     |  1.5768e-03, 1.4300e-03, 1.3189e-03, 1.2078e-03, 1.0966e-03,
     |  9.8550e-04, 9.0350e-04, 8.2150e-04, 7.3950e-04, 6.5750e-04,
     |  5.9885e-04, 5.4020e-04, 4.8155e-04, 4.2290e-04, 3.8250e-04,
     |  3.4210e-04, 3.0170e-04, 2.6130e-04, 2.3467e-04, 2.0805e-04,
     |  1.8143e-04, 1.5480e-04, 1.3809e-04, 1.2139e-04, 1.0468e-04,
     |  8.7970e-05, 7.7993e-05, 6.8015e-05, 5.8037e-05, 4.8060e-05,
     |  4.2382e-05, 3.6705e-05, 3.1028e-05, 2.5350e-05, 2.2257e-05,
     |  1.9165e-05, 1.6072e-05, 1.2980e-05, 1.1357e-05, 9.7340e-06,
     |  8.1110e-06, 6.4880e-06, 5.6620e-06, 4.8360e-06, 4.0100e-06,
     |  3.1840e-06, 2.7733e-06, 2.3625e-06, 1.9517e-06, 1.5410e-06,
     |  1.3403e-06, 1.1395e-06, 9.3882e-07, 7.3810e-07, 6.4135e-07,
     |  5.4460e-07, 4.4785e-07, 3.5110e-07, 3.0487e-07, 2.5865e-07,
     |  2.1242e-07, 1.6620e-07, 1.4425e-07, 1.2230e-07, 1.0035e-07,
     |  7.8400e-08, 6.8022e-08, 5.7645e-08, 4.7268e-08, 3.6890e-08,
     |  3.2000e-08, 2.7110e-08, 2.2220e-08, 1.7330e-08, 1.4857e-08,
     |  1.2383e-08, 9.9095e-09, 7.4360e-09/)
! extending the ar_glbm when ZITOP > 7 (ZITOP=7,dz=0.125,nlev=112) --Yihui Cai, Oct. 14, 2020
        do ih=113,nlevp1
           ind = zibot+(ih-1)*dz
           ar_glbm(ih) = 10**(-0.0255*ind**3+0.4354*ind**2
     |                   -3.1191*ind+1.1204)
        enddo
        ar_glbm(113) = 7.4360e-09
!#elif (NLEV==256)
      elseif (dz==0.0625) then
        ar_glbm(1:225) = (/
     |  1.1640e-02, 1.1527e-02, 1.1415e-02, 1.1303e-02, 1.1190e-02,
     |  1.1077e-02, 1.0965e-02, 1.0852e-02, 1.0740e-02, 1.0578e-02,
     |  1.0416e-02, 1.0254e-02, 1.0091e-02, 9.9294e-03, 9.7672e-03,
     |  9.6051e-03, 9.4430e-03, 9.2508e-03, 9.0585e-03, 8.8663e-03,
     |  8.6740e-03, 8.4817e-03, 8.2895e-03, 8.0973e-03, 7.9050e-03,
     |  7.7096e-03, 7.5143e-03, 7.3189e-03, 7.1235e-03, 6.9281e-03,
     |  6.7328e-03, 6.5374e-03, 6.3420e-03, 6.1651e-03, 5.9883e-03,
     |  5.8114e-03, 5.6345e-03, 5.4576e-03, 5.2807e-03, 5.1039e-03,
     |  4.9270e-03, 4.7782e-03, 4.6295e-03, 4.4808e-03, 4.3320e-03,
     |  4.1832e-03, 4.0345e-03, 3.8857e-03, 3.7370e-03, 3.6169e-03,
     |  3.4967e-03, 3.3766e-03, 3.2565e-03, 3.1364e-03, 3.0162e-03,
     |  2.8961e-03, 2.7760e-03, 2.6811e-03, 2.5863e-03, 2.4914e-03,
     |  2.3965e-03, 2.3016e-03, 2.2067e-03, 2.1119e-03, 2.0170e-03,
     |  1.9436e-03, 1.8703e-03, 1.7969e-03, 1.7235e-03, 1.6501e-03,
     |  1.5768e-03, 1.5034e-03, 1.4300e-03, 1.3744e-03, 1.3189e-03,
     |  1.2633e-03, 1.2078e-03, 1.1522e-03, 1.0966e-03, 1.0411e-03,
     |  9.8550e-04, 9.4450e-04, 9.0350e-04, 8.6250e-04, 8.2150e-04,
     |  7.8050e-04, 7.3950e-04, 6.9850e-04, 6.5750e-04, 6.2818e-04,
     |  5.9885e-04, 5.6952e-04, 5.4020e-04, 5.1088e-04, 4.8155e-04,
     |  4.5223e-04, 4.2290e-04, 4.0270e-04, 3.8250e-04, 3.6230e-04,
     |  3.4210e-04, 3.2190e-04, 3.0170e-04, 2.8150e-04, 2.6130e-04,
     |  2.4799e-04, 2.3467e-04, 2.2136e-04, 2.0805e-04, 1.9474e-04,
     |  1.8143e-04, 1.6811e-04, 1.5480e-04, 1.4645e-04, 1.3809e-04,
     |  1.2974e-04, 1.2139e-04, 1.1303e-04, 1.0468e-04, 9.6324e-05,
     |  8.7970e-05, 8.2981e-05, 7.7993e-05, 7.3004e-05, 6.8015e-05,
     |  6.3026e-05, 5.8037e-05, 5.3049e-05, 4.8060e-05, 4.5221e-05,
     |  4.2382e-05, 3.9544e-05, 3.6705e-05, 3.3866e-05, 3.1028e-05,
     |  2.8189e-05, 2.5350e-05, 2.3804e-05, 2.2257e-05, 2.0711e-05,
     |  1.9165e-05, 1.7619e-05, 1.6072e-05, 1.4526e-05, 1.2980e-05,
     |  1.2168e-05, 1.1357e-05, 1.0545e-05, 9.7340e-06, 8.9225e-06,
     |  8.1110e-06, 7.2995e-06, 6.4880e-06, 6.0750e-06, 5.6620e-06,
     |  5.2490e-06, 4.8360e-06, 4.4230e-06, 4.0100e-06, 3.5970e-06,
     |  3.1840e-06, 2.9786e-06, 2.7733e-06, 2.5679e-06, 2.3625e-06,
     |  2.1571e-06, 1.9517e-06, 1.7464e-06, 1.5410e-06, 1.4406e-06,
     |  1.3403e-06, 1.2399e-06, 1.1395e-06, 1.0392e-06, 9.3882e-07,
     |  8.3846e-07, 7.3810e-07, 6.8972e-07, 6.4135e-07, 5.9298e-07,
     |  5.4460e-07, 4.9623e-07, 4.4785e-07, 3.9947e-07, 3.5110e-07,
     |  3.2799e-07, 3.0487e-07, 2.8176e-07, 2.5865e-07, 2.3554e-07,
     |  2.1242e-07, 1.8931e-07, 1.6620e-07, 1.5522e-07, 1.4425e-07,
     |  1.3328e-07, 1.2230e-07, 1.1132e-07, 1.0035e-07, 8.9375e-08,
     |  7.8400e-08, 7.3211e-08, 6.8022e-08, 6.2834e-08, 5.7645e-08,
     |  5.2456e-08, 4.7268e-08, 4.2079e-08, 3.6890e-08, 3.4445e-08,
     |  3.2000e-08, 2.9555e-08, 2.7110e-08, 2.4665e-08, 2.2220e-08,
     |  1.9775e-08, 1.7330e-08, 1.6093e-08, 1.4857e-08, 1.3620e-08,
     |  1.2383e-08, 1.1146e-08, 9.9095e-09, 8.6728e-09, 7.4360e-09/)
! extending the ar_glbm when ZITOP > 7 (ZITOP=7,dz=0.0625,nlev=224) --Yihui Cai, Oct. 14, 2020
        do ih=225,nlevp1
           ind = zibot+(ih-1)*dz
           ar_glbm(ih) = 10**(-0.0255*ind**3+0.4354*ind**2
     |                   -3.1191*ind+1.1204)
        enddo
        ar_glbm(225) = 7.4360e-09
      else
        write(6,"('>>> comp_ar: unsupported NLEV=',i4,' dz=',f10.4)")
     |    nlev,dz
        call shutdown('comp_ar dz')
      endif
!#endif
!
      end subroutine alloc_ar
!-----------------------------------------------------------------------
      end module ar_module
