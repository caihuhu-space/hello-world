!
#include <defs.h>
      module ringfilter_module
!
! Coded by Tong Dang, May, 2017
! This program is to build a new ring filter for the future high-resolution 
! TIE-GCM. The filter procedures are as follows:
! 1. Divide the longitude sessions into several chunks
! 2. Calculate the average values of each chunk
! 3. Do the restruction based on the mean value of each chunk
!
      use params_module, only: nlon, nlonp4, nlat
      use cons_module, only: dlamda
      implicit none
      real, parameter :: rnlon = real(nlon)
!      real, parameter :: dlamda = 2. * 3.14159 / rnlon
#if (NLAT==36 && NLON==72)
      integer, parameter :: nlat_filter = 6 ! 5 deg 6 in total (untested)
      integer, parameter, dimension(nlat_filter) :: chunk_array =
     | (/9,9,18,18,36,36/)
#elif (NLAT==72 && NLON==144)
      integer, parameter :: nlat_filter = 10 ! 2.5 deg 10 in total (untested)
      integer, parameter, dimension(nlat_filter) :: chunk_array =
     | (/9,9,18,18,36,36,72,72,72,72/)
#elif (NLAT==144 && NLON==288)
      integer, parameter :: nlat_filter = 18 ! 1.25 deg 18 in total
      integer, parameter, dimension(nlat_filter) :: chunk_array =
     | (/9,9,18,18,36,36,36,36,72,72,144,144,144,144,144,144,144,144/)
#elif (NLAT==288 && NLON==576)
      integer, parameter :: nlat_filter = 40 ! 0.625 deg 40 in total
      integer, parameter, dimension(nlat_filter) :: chunk_array =
     | (/9,9,9,9,18,18,18,18,36,36,36,36,36,36,36,36,
     | 72,72,72,72,72,72,72,72,144,144,144,144,144,144,144,144,
     | 288,288,288,288,288,288,288,288/)
#else
      UNKNOWN NLAT,NLON
#endif
!
      contains
!-----------------------------------------------------------------------
      subroutine ringfilter(f, lev0, lev1, lat)
!
! Ringfilter for the second order of FFT
! keep the first order of fourier series, and filter higher orders
! Args:
      integer, intent(in) :: lev0, lev1, lat
      real, intent(inout) :: f(nlonp4, lev0:lev1)
!
! Local:,
      integer :: k, i, n, nchnk, ichnk, im2, im1, ip1, ip2, m
      real :: a0, a1, b1, theta, fm2, fm1, ff, fp1, fp2,
     |  fL, fR, a, b, c
      real, dimension(nlon) :: w, wm, fx
      real, dimension(nlat) :: average
!
      if (lat<=nlat_filter .or. lat>=nlat-nlat_filter+1) then

        do k = lev0, lev1

! Load field data into w
! Fourier expansion: f(x) = a0+a1*cos(x)+b1*sin(x)
          a1 = 0.
          b1 = 0.
          do i = 1, nlon
            w(i) = f(i+2, k)
            theta = dlamda * i
            a1 = a1 + w(i)*cos(theta)
            b1 = b1 + w(i)*sin(theta)
          enddo
          a1 = a1 * 2. / rnlon
          b1 = b1 * 2. / rnlon
          a0 = sum(w) / rnlon

! Chunk numbers in this latitude
          if (lat <= nlat_filter) nchnk = chunk_array(lat)
          if (lat >= nlat-nlat_filter+1) nchnk = chunk_array(nlat-lat+1)

! w(i) = wm(i)+fx(i), then filter fx(i)
          do i = 1, nlon
            theta = dlamda * i
            wm(i) = a0 + a1*cos(theta) + b1*sin(theta)
!           fx(i) = w(i) - wm(i)
            fx(i) = w(i)
          enddo

! Start the ring average filtering

! Grid points in each chunk
          n = nlon / nchnk

! Calculate the average value in each chunk
          do ichnk = 1, nchnk
            average(ichnk) = sum(fx((ichnk-1)*n+1: ichnk*n)) / real(n)
          enddo

! Then do the linear interpolation between each fL, fR
          do ichnk = 1, nchnk

! Calculate fL, fR
            im2 = ichnk - 2
            if (im2 <= 0) im2 = im2 + nchnk
            im1 = ichnk - 1
            if (im1 <= 0) im1 = im1 + nchnk
            ip1 = ichnk + 1
            if (ip1 >= nchnk+1) ip1 = ip1 - nchnk
            ip2 = ichnk + 2
            if (ip2 >= nchnk+1) ip2 = ip2 - nchnk

            fm2 = average(im2)
            fm1 = average(im1)
            ff = average(ichnk)
            fp1 = average(ip1)
            fp2 = average(ip2)

            fL = (-fm2 + 7.*fm1 + 7.*ff - fp1) / 12.
            fR = (-fm1 + 7.*ff + 7.*fp1 - fp2) / 12.

            a = fL + fR - 2.*ff
            b = 3.*ff - fR - 2.*fL
            c = fL

! Calculate the filtered data at m + (ichnk-1)*n
            do m = 1, n
              fx(m + (ichnk-1)*n) =
!    |          fL*(2.*n-2.*m+1.)/(2.*n) + fR*(2.*m-1.)/(2.*n)
     |          a*real(3*m**2-3*m+1)/real(n**2)
     |          + b*real(2*m-1)/real(n) + c
            enddo

          enddo ! ichnk = 1, nchnk

! Save filtered field:
          do i = 1, nlon
!           f(i+2, k) = fx(i) + wm(i)
            f(i+2, k) = fx(i)
          enddo ! i = 1, nlon

        enddo ! k = lev0, lev1

      endif ! lat
      end subroutine ringfilter
!-----------------------------------------------------------------------
      subroutine ringfilter_uv(f, lev0, lev1, lat)
!
! Ringfilter for the fourth order of FFT
! keep first three orders of fourier series, and filter higher orders
! Args:
      integer, intent(in) :: lev0, lev1, lat
      real, intent(inout) :: f(nlonp4, lev0:lev1)
!
! Local:,
      integer :: k, i, n, nchnk, ichnk, im2, im1, ip1, ip2, m
      real :: a0, a1, b1, a2, b2, a3, b3, theta, fm2, fm1, ff, fp1, fp2,
     |  fL, fR, a, b, c
      real, dimension(nlon) :: w, wm, fx
      real, dimension(nlat) :: average
!
      if (lat<=nlat_filter .or. lat>=nlat-nlat_filter+1) then

        do k = lev0, lev1

! Load field data into w
! Fourier expansion: f(x) = a0+a1*cos(x)+b1*sin(x)
!                          +a2*cos(2x)+b2*sin(2x)+a3*cos(3x)+b3*sin(3x)
          a1 = 0.
          b1 = 0.
          a2 = 0.
          b2 = 0.
          a3 = 0.
          b3 = 0.
          do i = 1, nlon
            w(i) = f(i+2, k)
            theta = dlamda * i
            a1 = a1 + w(i)*cos(theta)
            b1 = b1 + w(i)*sin(theta)
            a2 = a2 + w(i)*cos(2*theta)
            b2 = b2 + w(i)*sin(2*theta)
            a3 = a3 + w(i)*cos(3*theta)
            b3 = b3 + w(i)*sin(3*theta)
          enddo
          a1 = a1 * 2. / rnlon
          b1 = b1 * 2. / rnlon
          a2 = a2 * 2. / rnlon
          b2 = b2 * 2. / rnlon
          a3 = a3 * 2. / rnlon
          b3 = b3 * 2. / rnlon
          a0 = sum(w) / rnlon

! Chunk numbers in this latitude
          if (lat <= nlat_filter) nchnk = chunk_array(lat)
          if (lat >= nlat-nlat_filter+1) nchnk = chunk_array(nlat-lat+1)

! w(i) = wm(i)+fx(i), then filter fx(i)
          do i = 1, nlon
            theta = dlamda * i
            wm(i) = a0 + a1*cos(theta) + b1*sin(theta)
     |        + a2*cos(2*theta) + b2*sin(2*theta)
     |        + a3*cos(3*theta) + b3*sin(3*theta)
            fx(i) = w(i) - wm(i)
          enddo

! Start the ring average filtering

! Grid points in each chunk
          n = nlon / nchnk

! Calculate the average value in each chunk
          do ichnk = 1, nchnk
            average(ichnk) = sum(fx((ichnk-1)*n+1: ichnk*n)) / real(n)
          enddo

! Then do the linear interpolation between each fL, fR
          do ichnk = 1, nchnk-1

! Calculate fL, fR
            im2 = ichnk - 2
            if (im2 <= 0) im2 = im2 + nchnk
            im1 = ichnk - 1
            if (im1 <= 0) im1 = im1 + nchnk
            ip1 = ichnk + 1
            if (ip1 >= nchnk+1) ip1 = ip1 - nchnk
            ip2 = ichnk + 2
            if (ip2 >= nchnk+1) ip2 = ip2 - nchnk

            fm2 = average(im2)
            fm1 = average(im1)
            ff = average(ichnk)
            fp1 = average(ip1)
            fp2 = average(ip2)

            fL = (-fm2 + 7.*fm1 + 7.*ff - fp1) / 12.
            fR = (-fm1 + 7.*ff + 7.*fp1 - fp2) / 12.

            a = fL + fR - 2.*ff
            b = 3.*ff - fR - 2.*fL
            c = fL

! Calculate the filtered data at m + (ichnk-1)*n
            do m = 1, n
              fx(m + (ichnk-1)*n) =
!    |          fL*(2.*n-2.*m+1.)/(2.*n) + fR*(2.*m-1.)/(2.*n)
     |          a*real(3*m**2-3*m+1)/real(n**2)
     |          + b*real(2*m-1)/real(n) + c
            enddo

          enddo ! ichnk = 1, nchnk

! Save filtered field:
          do i = 1, nlon
            f(i+2, k) = fx(i) + wm(i)
          enddo ! i = 1, nlon

        enddo ! k = lev0, lev1

      endif ! lat
      end subroutine ringfilter_uv
!-----------------------------------------------------------------------
      subroutine linearfilter(f, lev0, lev1, lat)
!
! Ringfilter for the second order of FFT
! keep the first order of fourier series, and filter higher orders
! Args:
      integer,intent(in) :: lev0, lev1, lat
      real,intent(inout) :: f(nlonp4, lev0:lev1)
!
! Local:
      integer :: k, i, n, nchnk, ichnk, m
      real :: a0, a1, b1, theta, fac
      real, dimension(nlon) :: w, wm, fx
      real, dimension(nlat) :: average
!
      if (lat<=nlat_filter .or. lat>=nlat-nlat_filter+1) then

        do k = lev0, lev1

! Load field data into w
! Fourier expansion: f(x) = a0+a1*cos(x)+b1*sin(x)
          a1 = 0.
          b1 = 0.
          do i = 1, nlon
            w(i) = f(i+2, k)
            theta = dlamda * i
            a1 = a1 + w(i)*cos(theta)
            b1 = b1 + w(i)*sin(theta)
          enddo
          a1 = a1 * 2. / rnlon
          b1 = b1 * 2. / rnlon
          a0 = sum(w) / rnlon

! Chunk numbers in this latitude
          if (lat <= nlat_filter) nchnk = chunk_array(lat)
          if (lat >= nlat-nlat_filter+1) nchnk = chunk_array(nlat-lat+1)

! w(i) = wm(i)+fx(i), then filter fx(i)
          do i = 1, nlon
            theta = dlamda * i
            wm(i) = a0 + a1*cos(theta) + b1*sin(theta)
            fx(i) = w(i) - wm(i)
          enddo

! Start the ring average filtering

! Grid points in each chunk
          n = nlon / nchnk

! Calculate the average value in each chunk
          do ichnk = 1, nchnk
            average(ichnk) = sum(fx((ichnk-1)*n+1: ichnk*n)) / real(n)
          enddo

! Then do the linear interpolation between each average point

          fac = 1. / real(2*n)
! Calculate first n/2 and last n/2  
          do m = 1, n/2
            fx(nlon - n/2 + m) = fac *
     |        (average(nchnk)*real(2*n-2*m+1) + average(1)*real(2*m-1))
            fx(m) = fac *
     |        (average(nchnk)*real(n-2*m+1) + average(1)*real(n+2*m-1))
          enddo

! Calculate the filtered data at the 1st grid of this chunk
          do m = 1, n
            fx(n/2 + m) = fac *
     |        (average(1)*real(2*n-2*m+1) + average(2)*real(2*m-1))
          enddo

          do ichnk = 2, nchnk-1

! Calculate the filtered data at the mth grid of this chunk
            do m = 1, n
              fx(ichnk*n - n/2 + m) = fac *
     |          (average(ichnk)*real(2*n-2*m+1) +
     |          average(ichnk+1)*real(2*m-1))
            enddo

          enddo ! ichnk = 2, nchnk-1

! Save filtered field:
          do i = 1, nlon
            f(i+2, k) = fx(i) + wm(i)
          enddo ! i = 1, nlon

        enddo ! k = lev0, lev1

      endif ! lat
      end subroutine linearfilter
!-----------------------------------------------------------------------
      subroutine linearfilter_uv(f, lev0, lev1, lat)
!
! Ringfilter for the fourth order of FFT
! keep first three orders of fourier series, and filter higher orders
! Args:
      integer, intent(in) :: lev0, lev1, lat
      real, intent(inout) :: f(nlonp4, lev0:lev1)
!
! Local:,
      integer :: k, i, n, nchnk, ichnk, m
      real :: a0, a1, b1, a2, b2, a3, b3, theta, fac
      real, dimension(nlon) :: w, wm, fx
      real, dimension(nlat) :: average
!
      if (lat<=nlat_filter .or. lat>=nlat-nlat_filter+1) then

        do k = lev0, lev1

! Load field data into w
! Fourier expansion: f(x) = a0+a1*cos(x)+b1*sin(x)
!                          +a2*cos(2x)+b2*sin(2x)+a3*cos(3x)+b3*sin(3x)
          a1 = 0.
          b1 = 0.
          a2 = 0.
          b2 = 0.
          a3 = 0.
          b3 = 0.
          do i = 1, nlon
            w(i) = f(i+2, k)
            theta = dlamda * i
            a1 = a1 + w(i)*cos(theta)
            b1 = b1 + w(i)*sin(theta)
            a2 = a2 + w(i)*cos(2*theta)
            b2 = b2 + w(i)*sin(2*theta)
            a3 = a3 + w(i)*cos(3*theta)
            b3 = b3 + w(i)*sin(3*theta)
          enddo
          a1 = a1 * 2. / rnlon
          b1 = b1 * 2. / rnlon
          a2 = a2 * 2. / rnlon
          b2 = b2 * 2. / rnlon
          a3 = a3 * 2. / rnlon
          b3 = b3 * 2. / rnlon
          a0 = sum(w) / rnlon

! Chunk numbers in this latitude
          if (lat <= nlat_filter) nchnk = chunk_array(lat)
          if (lat >= nlat-nlat_filter+1) nchnk = chunk_array(nlat-lat+1)

! w(i) = wm(i)+fx(i), then filter fx(i)
          do i = 1, nlon
            theta = dlamda * i
            wm(i) = a0 + a1*cos(theta) + b1*sin(theta)
     |        + a2*cos(2*theta) + b2*sin(2*theta)
     |        + a3*cos(3*theta) + b3*sin(3*theta)
            fx(i) = w(i) - wm(i)
          enddo

! Start the ring average filtering

! Grid points in each chunk
          n = nlon / nchnk

! Calculate the average value in each chunk
          do ichnk = 1, nchnk
            average(ichnk) = sum(fx((ichnk-1)*n+1: ichnk*n)) / real(n)
          enddo

! Then do the linear interpolation between each fL, fR

          fac = 1. / real(2*n)
! Calculate first n/2 and last n/2  
          do m = 1, n/2
            fx(nlon - n/2 + m) = fac *
     |        (average(nchnk)*real(2*n-2*m+1) + average(1)*real(2*m-1))
            fx(m) = fac *
     |        (average(nchnk)*real(n-2*m+1) + average(1)*real(n+2*m-1))
          enddo

! Calculate the filtered data at the 1st grid of this chunk
          do m = 1, n
            fx(n/2 + m) = fac *
     |        (average(1)*real(2*n-2*m+1) + average(2)*real(2*m-1))
          enddo

          do ichnk = 2, nchnk-1

! Calculate the filtered data at the mth grid of this chunk
            do m = 1, n
              fx(ichnk*n - n/2 + m) = fac *
     |          (average(ichnk)*real(2*n-2*m+1) +
     |          average(ichnk+1)*real(2*m-1))
            enddo

          enddo ! ichnk = 2, nchnk-1

! Save filtered field:
          do i = 1, nlon
            f(i+2, k) = fx(i) + wm(i)
          enddo ! i = 1, nlon

        enddo ! k = lev0, lev1

      endif ! lat
      end subroutine linearfilter_uv
!-----------------------------------------------------------------------
      end module ringfilter_module
