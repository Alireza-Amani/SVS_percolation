!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------

subroutine inisoili_svs(ni, trnch)
   use tdpack
   use sfcbus_mod
   use svs_configs
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>

   integer ni, trnch

   !@Author  Maria Abrahamowicz, Stephane Belair , Vincent Fortin (20xx)
   !@Object  Compute soil properties for given soil texture. Compute these properties on
   !         native levels of database providing soil texture, and then map properties
   !         unto SVS levels.
   !@Arguments
   !             - Input -
   ! NI          longueur d'une tranche horizontale

   integer :: i, k, kk, jj
   REAL b, usb, fb, crit1_wfcint, crit2_wfcint, ts

   ! "geo" variables are on the levels of the geophysical soil texture datbase
   REAL, dimension(ni,nl_stp) :: wsat_geo, wwilt_geo, wfc_geo, b_geo, psisat_geo, &
           ksat_geo, wfcint_geo, fb_geo, quartz_geo,rhosoil_geo,conddry_geo,condsld_geo , wunfrz_geo
   real, pointer, dimension(:) :: zcgsat, zgrkef, zdraindens, zslop

   ! variables on the levels of SVS
   real, pointer, dimension(:,:) :: zbcoef, zclay, zfbcof, zksat, zpsisat, zsand, zwfc, zwfcint, zwsat, zwwilt, &
                                         zconddry, zcondsld , zquartz, zrhosoil,zwunfrz

   ! Alireza
   INTEGER :: layer_i, iter1

   ! print info for each read param call
   LOGICAL :: read_param_message


#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

   MKPTR1D(zcgsat, cgsat)
   MKPTR1D(zdraindens, draindens)
   MKPTR1D(zgrkef, grkef)
   MKPTR1D(zslop, slop)

   MKPTR2D(zbcoef, bcoef)
   MKPTR2D(zclay, clay)
   MKPTR2D(zfbcof, fbcof)
   MKPTR2D(zksat, ksat)
   MKPTR2D(zpsisat , psisat)
   MKPTR2D(zsand, sand)
   MKPTR2D(zwfc, wfc)
   MKPTR2D(zwfcint, wfcint)
   MKPTR2D(zwsat, wsat)
   MKPTR2D(zwunfrz, wunfrz)
   MKPTR2D(zwwilt , wwilt)
   MKPTR2D(zconddry , conddry)
   MKPTR2D(zcondsld , condsld)
   MKPTR2D(zrhosoil , rhosoil)
   MKPTR2D(zquartz , quartz)


   ! store read parameters temporarily
   ! this helps with dealing with different number of layers
   ALLOCATE(temp_var(nl_svs))

   ! calculate soil parameters on native GEO layers, and then map them unto model layers.
   ! calculate weights to be used in phybusinit.... because here... we are
   ! re-doing the calculation for each row of domain...
   ! but the weights are the same !

   ! Compute averaged values for sand and clay for each layer when SOILGRID is used
   if ( soiltext == "SOILGRIDS" ) then
      do i=1,ni
         do k=1,(nl_stp-1)
            zsand(i,k) = (zsand(i,k)+zsand(i,k+1))/2.
            zclay(i,k) = (zclay(i,k)+zclay(i,k+1))/2.
         enddo
      enddo
   endif

   !     Computer soil properties for GEO layers
   do i=1,ni
      ! WRITE(*, "(/, 'inside inisoili_svs.f90', /)")
      do k=1,nl_stp
         wsat_geo  (i,k)  =  -0.00126   * zsand(i,k) + 0.489
         wwilt_geo (i,k)  =  37.1342e-3 * sqrt(max(1.,zclay(i,k)))
         wfc_geo   (i,k)  =  89.0467e-3 * max(1.,zclay(i,k))**0.3496
         psisat_geo(i,k)  =  0.01 * ( 10.0**(-0.0131 * zsand(i,k) + 1.88) )
         ksat_geo  (i,k)  =  ( 10.0**(0.0153 * zsand(i,k) - 0.884) ) * 7.0556E-6

         b                 =  0.137 * zclay(i,k)  + 3.501
         b_geo     (i,k)  =  b
         usb               =  1./b
         fb                =  b**usb/(b-1.) * ((3.*b+2.)**(1.-usb)-(2.*b+2.)**(1.-usb))
         fb_geo(i,k)      =  fb
         ! Compute water content at field capacity along sloping aquifer based on Soulis et al. 2012
         ! Ensure that wc at fc stays between wilting point and saturation

         crit1_wfcint   = 2.*zdraindens(i)*psisat_geo(i,k)*(wsat_geo(i,k)/wwilt_geo(i,k)*fb)**b
         crit2_wfcint   = 2.*zdraindens(i)*psisat_geo(i,k)*fb**b

         if (abs(zslop(i)).gt.crit1_wfcint) then
            wfcint_geo(i,k) = wwilt_geo(i,k)
         elseif (abs(zslop(i)).lt.crit2_wfcint) then
            wfcint_geo(i,k) = wsat_geo(i,k)
         elseif (zslop(i).ne.0.0) then
            wfcint_geo(i,k) = wsat_geo(i,k) * fb * &
                 ( psisat_geo(i,k)/ABS(zslop(i)) *2. * zdraindens(i) )**usb
         else
            wfcint_geo(i,k) = wfc_geo(i,k)
         endif

       ! Compute soil thermal properties for soil freezing

!       Quartz content (ref : NL95 & PL98)):
        quartz_geo(i,k)  = 0.038 + 0.0095*zsand(i,k)

!       Soil dry density (PL98):
        rhosoil_geo(i,k) = (1.0-wsat_geo(i,k))*2700.

!       Soil solid conductivity:
        if (quartz_geo(i,k).gt.0.20) then
           condsld_geo(i,k) = (7.7**quartz_geo(i,k)) *  &
                           (2.0**(1.0-quartz_geo(i,k)))
        endif
        if (quartz_geo(i,k).le.0.20) then
           condsld_geo(i,k) = (7.7**quartz_geo(i,k)) *  &
                           (3.0**(1.0-quartz_geo(i,k)))
        endif

!       Soil dry conductivity:
!       Calculate the value of the expression:
!       exp(4.31 * 10**-3 * clay_percentage + 0.821 * soil_dry_density - 0.543 * porosity - 2.372)
!        ref: `Modelling dry soil thermal conductivity`, He et al. 2019

        conddry_geo(i,k) = exp(&
            (4.31 * 10**-3 * zclay(i,k)) + &
            (0.821 * (rhosoil_geo(i,k) * 0.001)) - &
            (0.543 * wsat_geo(i,k)) - &
            2.372 &
        )



      !   conddry_geo(i,k) = (0.135*rhosoil_geo(i,k) + 64.7) / &
      !                   (2700. - 0.947*rhosoil_geo(i,k))

!       Unfrozen residual water content obtained from Niu and Yang (2006)
!       Average value between -10 and -2 deg C
        wunfrz_geo(i,k) = 0.
        do jj = 0,4
            ts =263.15+ 2.*jj
            wunfrz_geo(i,k)     =  wunfrz_geo(i,k)+ wsat_geo(i,k)*(CHLF*(ts-273.15)/(ts*(-1.0*psisat_geo(i,k))*9.81))**(-1.0*usb)
        enddo
        wunfrz_geo(i,k) =  wunfrz_geo(i,k)/5.

      enddo
   enddo

   ! Alireza
   !  save geo parameters in a file before they are recomputed for the soil profile
   IF (save_geoparams_to_file) THEN
   !  create and open file
      OPEN(UNIT=geoparam_file_unit, FILE="default_geoparams.csv", STATUS='REPLACE', ACTION='WRITE', FORM='FORMATTED')

      ! write headers
      WRITE(geoparam_file_unit, *) "Unit,Layer,Sand,Clay," // &
         "wsat_geo,wwilt_geo,wfc_geo, " // &
         "b_geo, " // &
         "fb_geo, wfcint_geo, quartz_geo, " // &
         "rhosoil_geo, condsld_geo, conddry_geo, " // &
         "wunfrz_geo,psisat_geo,ksat_geo"

      ! write data
      DO i=1, ni
      DO k=1, nl_stp
         WRITE(geoparam_file_unit, '(2(I5, ", "), 12(F10.4, ", "), F10.4, ", ", ES10.2, ", ", ES10.2)') &
            i, k, &
            zsand(i,k), zclay(i,k), &
            wsat_geo(i,k), wwilt_geo(i,k), wfc_geo(i,k), &
            b_geo(i,k), &
            fb_geo(i,k), wfcint_geo(i,k), quartz_geo(i,k), &
            rhosoil_geo(i,k), condsld_geo(i,k), conddry_geo(i,k), &
            wunfrz_geo(i,k), psisat_geo(i,k), ksat_geo(i,k)
      END DO
      END DO

      CLOSE(UNIT=geoparam_file_unit)

      ! set flag to false so that the file is not overwritten again
      save_geoparams_to_file = .FALSE.
   END IF







   ! Ali_mod: read parameters from parameters file
   IF (print_inisoili) THEN

      read_param_message = .TRUE.

      WRITE(*, "(/, A, /)") "Printing information from `inisoili_svs.F90` -------------------------"

   ELSE

      read_param_message = .FALSE.

   END IF

   temp_scalar_var = read_param("read_user_parameters", 1, .FALSE.)
   IF (temp_scalar_var(1) /= 666) THEN
      read_parameters = read_param("read_user_parameters", 1, .FALSE.)
   END IF

   IF (read_parameters(1) .EQ. 1) THEN

      !     Compute soil properties for GEO layers
      do i=1,ni

            temp_var = read_param("wsat", NL_SVS, .FALSE.)
            IF (temp_var(1) /= 666) THEN
               wsat_geo(i, :) = read_param("wsat", NL_SVS, read_param_message)
               DO layer_i = 1, nl_stp
                 !       Soil dry conductivity:
                  conddry_geo(i,layer_i) = exp(&
                     (4.31 * 10**-3 * zclay(i,layer_i)) + &
                     (0.821 * (rhosoil_geo(i,layer_i) * 0.001)) - &
                     (0.543 * wsat_geo(i,layer_i)) - &
                     2.372 &
               )
               END DO
            END IF

            temp_var = read_param("wwilt", NL_SVS, .FALSE.)
            IF (temp_var(1) /= 666) THEN
               wwilt_geo(i, :) = read_param("wwilt", NL_SVS, read_param_message)
            END IF

            temp_var = read_param("wfc", NL_SVS, .FALSE.)
            IF (temp_var(1) /= 666) THEN
               wfc_geo(i, :) = read_param("wfc", NL_SVS, read_param_message)
            END IF

            temp_var = read_param("psisat", NL_SVS, .FALSE.)
            IF (temp_var(1) /= 666) THEN
               psisat_geo(i, :) = read_param("psisat", NL_SVS, read_param_message)
            END IF

            temp_var = read_param("ksat", NL_SVS, .FALSE.)
            IF (temp_var(1) /= 666) THEN
               ksat_geo(i, :) = read_param("ksat", NL_SVS, read_param_message)
            END IF

            temp_var = read_param("bcoef", NL_SVS, .FALSE.)
            IF (temp_var(1) /= 666) THEN
               b_geo(i, :) = read_param("bcoef", NL_SVS, read_param_message)

               DO layer_i = 1, nl_stp
                  b = b_geo(i, layer_i)
                  usb = 1./b
                  fb =  b**usb/(b-1.) * ((3.*b+2.)**(1.-usb)-(2.*b+2.)**(1.-usb))
                  fb_geo(i, layer_i) =  fb
               END DO

            END IF

            DO layer_i = 1, nl_stp
               ! Compute water content at field capacity along sloping aquifer based on Soulis et al. 2012
               ! Ensure that wc at fc stays between wilting point and saturation
               crit1_wfcint = 2.*zdraindens(i)*psisat_geo(i,layer_i)* &
                              (wsat_geo(i,layer_i)/wwilt_geo(i,layer_i)*fb_geo(i, layer_i))**b_geo(i, layer_i)

               crit2_wfcint = 2.*zdraindens(i)*psisat_geo(i,layer_i)*fb_geo(i, layer_i)**b_geo(i, layer_i)

                if (abs(zslop(i)).gt.crit1_wfcint) then
                  wfcint_geo(i,layer_i) = wwilt_geo(i,layer_i)
               elseif (abs(zslop(i)).lt.crit2_wfcint) then
                  wfcint_geo(i,layer_i) = wsat_geo(i,layer_i)
               elseif (zslop(i).ne.0.0) then
                  wfcint_geo(i,layer_i) = wsat_geo(i,layer_i) * fb * &
                     ( psisat_geo(i,layer_i)/ABS(zslop(i)) *2. * zdraindens(i) )**(1./b_geo(i, layer_i))
               else
                  wfcint_geo(i,layer_i) = wfc_geo(i,layer_i)
               endif

            END DO

         ! Compute soil thermal properties for soil freezing

            temp_var = read_param("quartz", NL_SVS, .FALSE.)
            IF (temp_var(1) /= 666) THEN
               quartz_geo(i, :) = read_param("quartz", NL_SVS, read_param_message)

               DO layer_i = 1, nl_stp
                  ! I dont think I will over-write these values, let them be for now.
                  ! maybe late for sensitivity An.
                  !       Soil solid conductivity:
                  if (quartz_geo(i,layer_i).gt.0.20) then
                     condsld_geo(i,layer_i) = (7.7**quartz_geo(i,layer_i)) *  &
                                       (2.0**(1.0-quartz_geo(i,layer_i)))
                  endif
                  if (quartz_geo(i,layer_i).le.0.20) then
                     condsld_geo(i,layer_i) = (7.7**quartz_geo(i,layer_i)) *  &
                                       (3.0**(1.0-quartz_geo(i,layer_i)))
                  endif
               END DO

            END IF

            temp_var = read_param("rhosoil", NL_SVS, .FALSE.)
            IF (temp_var(1) /= 666) THEN
               rhosoil_geo(i, :) = read_param("rhosoil", NL_SVS, read_param_message)

               DO layer_i = 1, nl_stp
                 !       Soil dry conductivity:
                  ! conddry_geo(i,layer_i) = (0.135*rhosoil_geo(i,layer_i) + 64.7) / &
                  !                   (2700. - 0.947*rhosoil_geo(i,layer_i))
                  conddry_geo(i,layer_i) = exp(&
                     (4.31 * 10**-3 * zclay(i,layer_i)) + &
                     (0.821 * (rhosoil_geo(i,layer_i) * 0.001)) - &
                     (0.543 * wsat_geo(i,layer_i)) - &
                     2.372 &
               )
               END DO

            END IF

            ! read cond dry and solid
            temp_var = read_param("conddry", NL_SVS, .FALSE.)
            IF (temp_var(1) /= 666) THEN
               conddry_geo(i, :) = read_param("conddry", NL_SVS, read_param_message)
            END IF

            temp_var = read_param("condsld", NL_SVS, .FALSE.)
            IF (temp_var(1) /= 666) THEN
               condsld_geo(i, :) = read_param("condsld", NL_SVS, read_param_message)
            END IF


            temp_var = read_param("wunfrz", NL_SVS, .FALSE.)
            IF (temp_var(1) /= 666) THEN

               wunfrz_geo(i, :) = read_param("wunfrz", NL_SVS, read_param_message)

            ELSE

               DO layer_i = 1, nl_stp
                  !       Unfrozen residual water content obtained from Niu and Yang (2006)
                     !       Average value between -10 and -2 deg C
                  wunfrz_geo(i,layer_i) = 0.
                  do jj = 0,4
                        ts =263.15+ 2.*jj
                        wunfrz_geo(i,layer_i)     =  wunfrz_geo(i,layer_i)+ wsat_geo(i,layer_i)* &
                           (CHLF*(ts-273.15)/(ts*(-1.0*psisat_geo(i,layer_i))*9.81))**(-1.0*(1./b_geo(i, layer_i)))
                  enddo
                  wunfrz_geo(i,layer_i) =  wunfrz_geo(i,layer_i)/5.

               END DO

            END IF ! to read wunfrz


      enddo
   END IF




   ! "Map" GEO soil properties unto model soil layers
   Do i = 1 , ni
      Do k = 1, nl_svs
         do kk = 1 , nl_stp

            zwsat  (i,k)  = zwsat  (i,k) + wsat_geo  (i,kk)  * weights( k , kk)
            zwwilt (i,k)  = zwwilt (i,k) + wwilt_geo (i,kk)  * weights( k , kk)

            zwfc   (i,k)  = zwfc   (i,k) + wfc_geo   (i,kk)  * weights( k , kk)
            zbcoef (i,k)  = zbcoef (i,k) + b_geo     (i,kk)  * weights( k , kk)
            zfbcof (i,k)  = zfbcof (i,k) + fb_geo    (i,kk)  * weights( k , kk)
            zpsisat(i,k)  = zpsisat(i,k) + psisat_geo(i,kk)  * weights( k , kk)
            zksat  (i,k)  = zksat  (i,k) + ksat_geo  (i,kk)  * weights( k , kk)
            zwfcint(i,k)  = zwfcint(i,k) + wfcint_geo(i,kk)  * weights( k , kk)
            zconddry  (i,k)  = zconddry  (i,k) + conddry_geo  (i,kk)  * weights( k , kk)
            zcondsld  (i,k)  = zcondsld  (i,k) + condsld_geo  (i,kk)  * weights( k , kk)
            zquartz   (i,k)  = zquartz   (i,k) + quartz_geo   (i,kk)  * weights( k , kk)
            zrhosoil  (i,k)  = zrhosoil  (i,k) + rhosoil_geo  (i,kk)  * weights( k , kk)
            zwunfrz   (i,k)  = zwunfrz   (i,k) + wunfrz_geo   (i,kk)  * weights( k , kk)

         enddo
      enddo
      ! compute thermal coeff.
      ! for 1st model layer only --- here simply use 1st GEO soil texture !!! Do not map !
      zcgsat (i)  = ( -1.557e-2 * zsand(i,1) &
           -  1.441e-2 * zclay(i,1) + 4.7021 ) * 1.E-6

      ! Compute effective parameter for watdrain
      zgrkef(i)   = 2.* zdraindens(i) * zslop(i)

   enddo


   ! regardless of `read_user_parameters` values, read the following:

   temp_scalar_var = read_param("save_minutes_csv", 1, .FALSE.)
   IF (temp_scalar_var(1) /= 666) THEN
      save_minutes_csv = read_param("save_minutes_csv", 1, read_param_message)
   END IF

   temp_scalar_var = read_param("user_wfcdp", 1, .FALSE.)
   IF (temp_scalar_var(1) /= 666) THEN
      user_wfcdp = read_param("user_wfcdp", 1, read_param_message)
   END IF

   temp_scalar_var = read_param("OPT_LIQWAT", 1, .FALSE.)
   IF (temp_scalar_var(1) /= 666) THEN
      user_OPT_LIQWAT = read_param("OPT_LIQWAT", 1, read_param_message)
   END IF

   temp_scalar_var = read_param("OPT_FRAC", 1, .FALSE.)
   IF (temp_scalar_var(1) /= 666) THEN
      user_OPT_FRAC = read_param("OPT_FRAC", 1, read_param_message)
   END IF

   temp_scalar_var = read_param("OPT_SNOW", 1, .FALSE.)
   IF (temp_scalar_var(1) /= 666) THEN
      user_OPT_SNOW = read_param("OPT_SNOW", 1, read_param_message)
   END IF

   temp_scalar_var = read_param("water_ponding", 1, .FALSE.)
   IF (temp_scalar_var(1) /= 666) THEN
      user_water_ponding = read_param("water_ponding", 1, read_param_message)
      IF (user_water_ponding(1) .EQ. 1) THEN
         lwater_ponding_svs1 = .TRUE.
      ELSE IF (user_water_ponding(1) .EQ. 0) THEN
         lwater_ponding_svs1 = .FALSE.
      END IF
   END IF


   temp_scalar_var = read_param("maxpnd", 1, .FALSE.)
   IF (temp_scalar_var(1) /= 666) THEN
      user_maxpnd = read_param("maxpnd", 1, read_param_message)
   END IF

   temp_scalar_var = read_param("KFICE", 1, .FALSE.)
   IF (temp_scalar_var(1) /= 666) THEN
      user_KFICE = read_param("KFICE", 1, read_param_message)
   END IF

   temp_scalar_var = read_param("WAT_REDIS", 1, .FALSE.)
   IF (temp_scalar_var(1) /= 666) THEN
      user_WAT_REDIS = read_param("WAT_REDIS", 1, read_param_message)
   END IF


   IF (print_inisoili) THEN

      ! Ali_mod: printing info
      WRITE(*, "(A, I0)") "NL_SVS: ", NL_SVS
      WRITE(*, "(A, I0)") "NL_STP: ", nl_stp
      WRITE(*, "(A, I0, /)") "KHYD: ", KHYD


      ! main info for each layer
      WRITE(*, "('layer  #', 4x, 'Sand(%)', 4x, 'Clay(%)', 4x, 'Wsat(%)', 4x," // &
               "'Wfc(%)', 4x, 'Wwilt(%)', 4x, 'Ksat(m/s)', 4x, 'Psi_sat(kPa)')")

      DO iter1 = 1, NL_SVS

         WRITE(*, "('layer ', I2, 4x, F5.2, 6x, F5.2, 6x, F5.2, 6x, F5.2," // &
                  "5x, F5.2, 5x, ES10.2, 3x, ES10.2)") &

            iter1, zsand(ni,iter1), zclay(ni, iter1),  &
            zwsat(ni, iter1) * 100, zwfc(ni, iter1) * 100, zwwilt(ni, iter1) * 100, &
            zksat(ni, iter1), zpsisat(ni, iter1)

      END DO

      ! other info for each layer
      WRITE(*, "(/, 'layer  #', 4x, 'bcoef', 4x, 'rhosoil', 4x, 'quartz(%)', 4x," // &
               "'conddry', 4x, 'condsld', 4x, 'wunfrz(%)')")

      DO iter1 = 1, NL_SVS

         WRITE(*, "('layer ', I2, 3x, F5.2, 5x, F7.2, 4x, F5.2, 7x, F5.2," // &
                  "6x, F5.2, 7x, F5.2)") &

            iter1, zbcoef(ni,iter1), zrhosoil(ni, iter1),  &
            zquartz(ni, iter1) * 100, zconddry(ni, iter1), zcondsld(ni, iter1), &
            zwunfrz(ni, iter1) * 100

      END DO

      WRITE(*, "(/, A, /)") "-----------------------------------------------------------------------"

   END IF

   ! no more printing of info for next calls of this sub
   print_inisoili = .FALSE.

   !  save modified geoparams to file
   IF (save_modified_geoparams_to_file) THEN
   !  create and open file
      OPEN(UNIT=modified_geoparam_file_unit, FILE="modified_geoparams.csv", STATUS='REPLACE', ACTION='WRITE', FORM='FORMATTED')

      ! write headers
      WRITE(modified_geoparam_file_unit, *) "Unit,Layer,Sand,Clay," // &
         "wsat_geo,wwilt_geo,wfc_geo, " // &
         "b_geo, " // &
         "fb_geo, wfcint_geo, quartz_geo, " // &
         "rhosoil_geo, condsld_geo, conddry_geo, " // &
         "wunfrz_geo,psisat_geo,ksat_geo"

      ! write data
      DO i=1, ni
      DO k=1, nl_stp
         WRITE(modified_geoparam_file_unit, '(2(I5, ", "), 12(F10.4, ", "), F10.4, ", ", ES10.2, ", ", ES10.2)') &
            i, k, &
            zsand(i,k), zclay(i,k), &
            wsat_geo(i,k), wwilt_geo(i,k), wfc_geo(i,k), &
            b_geo(i,k), &
            fb_geo(i,k), wfcint_geo(i,k), quartz_geo(i,k), &
            rhosoil_geo(i,k), condsld_geo(i,k), conddry_geo(i,k), &
            wunfrz_geo(i,k), psisat_geo(i,k), ksat_geo(i,k)
      END DO
      END DO

      CLOSE(UNIT=modified_geoparam_file_unit)


      ! set flag to false so that the file is not overwritten again
      save_modified_geoparams_to_file = .FALSE.
   END IF


   return
 end subroutine inisoili_svs
