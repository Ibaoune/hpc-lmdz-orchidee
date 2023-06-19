        
        module m_simu_airs

        USE print_control_mod, ONLY: prt_level,lunout
          
        implicit none

        REAL, PARAMETER :: tau_thresh = 0.05 ! seuil nuages detectables
        REAL, PARAMETER :: p_thresh = 445.   ! seuil nuages hauts
        REAL, PARAMETER :: em_min = 0.2      ! seuils nuages semi-transp
        REAL, PARAMETER :: em_max = 0.85
        REAL, PARAMETER :: undef = 999.

        contains

        REAL function search_tropopause(P,T,alt,N) result(P_tropo)
! this function searches for the tropopause pressure in [hPa].
! The search is based on ideology described in
! Reichler et al., Determining the tropopause height from gridded data,
! GRL, 30(20) 2042, doi:10.1029/2003GL018240, 2003

        INTEGER N,i,i_lev,first_point,exit_flag,i_dir
        REAL P(N),T(N),alt(N),slope(N)
        REAL P_min, P_max, slope_limit,slope_2km, &
     & delta_alt_limit,tmp,delta_alt
        PARAMETER(P_min=75.0, P_max=470.0)   ! hPa
        PARAMETER(slope_limit=0.002)         ! 2 K/km converted to K/m
        PARAMETER(delta_alt_limit=2000.0)    ! 2000 meters

        do i=1,N-1
        slope(i)=-(T(i+1)-T(i))/(alt(i+1)-alt(i))
        end do
        slope(N)=slope(N-1)

        if (P(1).gt.P(N)) then
        i_dir= 1
        i=1
        else
        i_dir=-1
        i=N
        end if

        first_point=0
        exit_flag=0
        do while(exit_flag.eq.0)
        if (P(i).ge.P_min.and.P(i).le.P_max) then
        if (first_point.gt.0) then
        delta_alt=alt(i)-alt(first_point)
        if (delta_alt.ge.delta_alt_limit) then
        slope_2km=(T(first_point)-T(i))/delta_alt
        if (slope_2km.lt.slope_limit) then
        exit_flag=1
        else
        first_point=0
        end if
        end if
        end if
        if (first_point.eq.0.and.slope(i).lt.slope_limit) first_point=i
        end if
        i=i+i_dir
        if (i.le.1.or.i.ge.N) exit_flag=1
        end do

        if (first_point.le.0) P_tropo=65.4321
        if (first_point.eq.1) P_tropo=64.5432
        if (first_point.gt.1) then
        tmp=(slope_limit-slope(first_point))/(slope(first_point+1)- &
     & slope(first_point))*(P(first_point+1)-P(first_point))
        P_tropo=P(first_point)+tmp
        ! print*, 'P_tropo= ', tmp, P(first_point), P_tropo
        end if

! Ajout Marine
        if (P_tropo .lt. 60. .or. P_tropo .gt. 470.) then
        P_tropo = 999.
        endif

        return
        end function search_tropopause



        subroutine cloud_structure(len_cs, rneb_cs, temp_cs, &
     & emis_cs, iwco_cs, &
     & pres_cs, dz_cs, rhodz_cs, rad_cs, &
     & cc_tot_cs, cc_hc_cs, cc_hist_cs, &
     & cc_Cb_cs, cc_ThCi_cs, cc_Anv_cs, &
     & pcld_hc_cs, tcld_hc_cs, &
     & em_hc_cs, iwp_hc_cs, deltaz_hc_cs, &
     & pcld_Cb_cs, tcld_Cb_cs, em_Cb_cs, &
     & pcld_ThCi_cs, tcld_ThCi_cs, em_ThCi_cs, &
     & pcld_Anv_cs, tcld_Anv_cs, em_Anv_cs, &
     & em_hist_cs, iwp_hist_cs, deltaz_hist_cs, rad_hist_cs)


     
        INTEGER :: i, n, nss

        INTEGER, intent(in) :: len_cs
        REAL, DIMENSION(:), intent(in) :: rneb_cs, temp_cs
        REAL, DIMENSION(:), intent(in) :: emis_cs, iwco_cs, rad_cs
        REAL, DIMENSION(:), intent(in) :: pres_cs, dz_cs, rhodz_cs

        REAL, intent(out) :: cc_tot_cs, cc_hc_cs, cc_hist_cs, &
     & cc_Cb_cs, cc_ThCi_cs, cc_Anv_cs, &
     & pcld_hc_cs, tcld_hc_cs, em_hc_cs, iwp_hc_cs, &
     & pcld_Cb_cs, tcld_Cb_cs, em_Cb_cs, &
     & pcld_ThCi_cs, tcld_ThCi_cs, em_ThCi_cs, &
     & pcld_Anv_cs, tcld_Anv_cs, em_Anv_cs, &
     & em_hist_cs, iwp_hist_cs, &
     & deltaz_hc_cs, deltaz_hist_cs, rad_hist_cs

        REAL, DIMENSION(len_cs) :: rneb_ord
        REAL :: rneb_min

        REAL, DIMENSION(:), allocatable :: s, s_hc, s_hist, rneb_max
        REAL, DIMENSION(:), allocatable :: sCb, sThCi, sAnv
        REAL, DIMENSION(:), allocatable :: iwp_ss, pcld_ss, tcld_ss,&
     & emis_ss
        REAL, DIMENSION(:), allocatable :: deltaz_ss, rad_ss

        CHARACTER (len = 50)      :: modname = 'simu_airs.cloud_structure'
        CHARACTER (len = 160)     :: abort_message
        

        write(lunout,*) 'dans cloud_structure'

        call ordonne(len_cs, rneb_cs, rneb_ord)
       

! Definition des sous_sections

        rneb_min = rneb_ord(1)
        nss = 1

        do i = 2, size(rneb_cs)
        if (rneb_ord(i) .gt. rneb_min) then
        nss = nss+1
        rneb_min = rneb_ord(i)
        endif
        enddo

        allocate (s(nss))
        allocate (s_hc(nss))
        allocate (s_hist(nss))
        allocate (rneb_max(nss))
        allocate (emis_ss(nss))
        allocate (pcld_ss(nss))
        allocate (tcld_ss(nss))
        allocate (iwp_ss(nss))
        allocate (deltaz_ss(nss))
        allocate (rad_ss(nss))
        allocate (sCb(nss))
        allocate (sThCi(nss))
        allocate (sAnv(nss))

        rneb_min = rneb_ord(1)
        n = 1
        s(1) = rneb_ord(1)
        s_hc(1) = rneb_ord(1)
        s_hist(1) = rneb_ord(1)
        sCb(1) = rneb_ord(1)
        sThCi(1) = rneb_ord(1)
        sAnv(1) = rneb_ord(1)

        rneb_max(1) = rneb_ord(1)

        do i = 2, size(rneb_cs)
        if (rneb_ord(i) .gt. rneb_min) then
        n = n+1
        s(n) = rneb_ord(i)-rneb_min
        s_hc(n) = rneb_ord(i)-rneb_min
        s_hist(n) = rneb_ord(i)-rneb_min
        sCb(n) = rneb_ord(i)-rneb_min
        sThCi(n) = rneb_ord(i)-rneb_min
        sAnv(n) = rneb_ord(i)-rneb_min

        rneb_max(n) = rneb_ord(i)
        rneb_min = rneb_ord(i)
        endif
        enddo

! Appel de sous_section

        do i = 1, nss
         call sous_section(len_cs, rneb_cs, temp_cs, &
     &  emis_cs, iwco_cs, &
     &  pres_cs, dz_cs, rhodz_cs, rad_cs, rneb_ord, &
     &  rneb_max(i),s(i),s_hc(i),s_hist(i), &
     &  sCb(i), sThCi(i), sAnv(i), &
     &  emis_ss(i), &
     &  pcld_ss(i), tcld_ss(i), iwp_ss(i), deltaz_ss(i), rad_ss(i))
        enddo

! Caracteristiques de la structure nuageuse

        cc_tot_cs = 0.
        cc_hc_cs = 0.
        cc_hist_cs = 0.

        cc_Cb_cs = 0.
        cc_ThCi_cs = 0.
        cc_Anv_cs = 0.

        em_hc_cs = 0.
        iwp_hc_cs = 0.
        deltaz_hc_cs = 0.

        em_hist_cs = 0.
        iwp_hist_cs = 0.
        deltaz_hist_cs = 0.
        rad_hist_cs = 0.

        pcld_hc_cs = 0.
        tcld_hc_cs = 0.

        pcld_Cb_cs = 0.
        tcld_Cb_cs = 0.
        em_Cb_cs = 0.

        pcld_ThCi_cs = 0.
        tcld_ThCi_cs = 0.
        em_ThCi_cs = 0.

        pcld_Anv_cs = 0.
        tcld_Anv_cs = 0.
        em_Anv_cs = 0.

         do i = 1, nss

        cc_tot_cs = cc_tot_cs + s(i)
        cc_hc_cs = cc_hc_cs + s_hc(i)
        cc_hist_cs = cc_hist_cs + s_hist(i)

        cc_Cb_cs = cc_Cb_cs + sCb(i)
        cc_ThCi_cs = cc_ThCi_cs + sThCi(i)
        cc_Anv_cs = cc_Anv_cs + sAnv(i)

        iwp_hc_cs = iwp_hc_cs + s_hc(i)*iwp_ss(i)
        em_hc_cs = em_hc_cs + s_hc(i)*emis_ss(i)
        deltaz_hc_cs = deltaz_hc_cs + s_hc(i)*deltaz_ss(i)

        iwp_hist_cs = iwp_hist_cs + s_hist(i)*iwp_ss(i)
        em_hist_cs = em_hist_cs + s_hist(i)*emis_ss(i)
        deltaz_hist_cs = deltaz_hist_cs + s_hist(i)*deltaz_ss(i)
        rad_hist_cs = rad_hist_cs + s_hist(i)*rad_ss(i)

        pcld_hc_cs = pcld_hc_cs + s_hc(i)*pcld_ss(i)
        tcld_hc_cs = tcld_hc_cs + s_hc(i)*tcld_ss(i)

        pcld_Cb_cs = pcld_Cb_cs + sCb(i)*pcld_ss(i)
        tcld_Cb_cs = tcld_Cb_cs + sCb(i)*tcld_ss(i)
        em_Cb_cs = em_Cb_cs + sCb(i)*emis_ss(i)

        pcld_ThCi_cs = pcld_ThCi_cs + sThCi(i)*pcld_ss(i)
        tcld_ThCi_cs = tcld_ThCi_cs + sThCi(i)*tcld_ss(i)
        em_ThCi_cs = em_ThCi_cs + sThCi(i)*emis_ss(i)

        pcld_Anv_cs = pcld_Anv_cs + sAnv(i)*pcld_ss(i)
        tcld_Anv_cs = tcld_Anv_cs + sAnv(i)*tcld_ss(i)
        em_Anv_cs = em_Anv_cs + sAnv(i)*emis_ss(i)

        enddo

        deallocate(s)
        deallocate (s_hc)
        deallocate (s_hist)
        deallocate (rneb_max)
        deallocate (emis_ss)
        deallocate (pcld_ss)
        deallocate (tcld_ss)
        deallocate (iwp_ss)
        deallocate (deltaz_ss)
        deallocate (rad_ss)
        deallocate (sCb)
        deallocate (sThCi)
        deallocate (sAnv)

       call normal_undef(pcld_hc_cs,cc_hc_cs)
       call normal_undef(tcld_hc_cs,cc_hc_cs)
       call normal_undef(iwp_hc_cs,cc_hc_cs)
       call normal_undef(em_hc_cs,cc_hc_cs)
       call normal_undef(deltaz_hc_cs,cc_hc_cs)
       
       call normal_undef(iwp_hist_cs,cc_hist_cs)
       call normal_undef(em_hist_cs,cc_hist_cs)
       call normal_undef(deltaz_hist_cs,cc_hist_cs)
       call normal_undef(rad_hist_cs,cc_hist_cs)

       call normal_undef(pcld_Cb_cs,cc_Cb_cs)
       call normal_undef(tcld_Cb_cs,cc_Cb_cs)
       call normal_undef(em_Cb_cs,cc_Cb_cs)

       call normal_undef(pcld_ThCi_cs,cc_ThCi_cs)
       call normal_undef(tcld_ThCi_cs,cc_ThCi_cs)
       call normal_undef(em_ThCi_cs,cc_ThCi_cs)
    
       call normal_undef(pcld_Anv_cs,cc_Anv_cs)
       call normal_undef(tcld_Anv_cs,cc_Anv_cs)
       call normal_undef(em_Anv_cs,cc_Anv_cs)


! Tests

        if (cc_tot_cs .gt. maxval(rneb_cs) .and. &
     & abs(cc_tot_cs-maxval(rneb_cs)) .gt. 1.e-4 )  then
          WRITE(abort_message,*) 'cc_tot_cs > max rneb_cs', cc_tot_cs, maxval(rneb_cs)
          CALL abort_physic(modname,abort_message,1)
        endif

        if (iwp_hc_cs .lt. 0.) then
          abort_message= 'cloud_structure: iwp_hc_cs < 0'
          CALL abort_physic(modname,abort_message,1)
        endif
 
        end subroutine cloud_structure


        subroutine normal_undef(num, den)

        REAL, intent(in) :: den
        REAL, intent(inout) :: num

        if (den .ne. 0) then
        num = num/den
        else
        num = undef
        endif

        end subroutine normal_undef


        subroutine normal2_undef(res,num,den)

        REAL, intent(in) :: den
        REAL, intent(in) :: num
        REAL, intent(out) :: res

        if (den .ne. 0.) then
        res = num/den
        else
        res = undef
        endif

        end subroutine normal2_undef


        subroutine sous_section(len_cs, rneb_cs, temp_cs, &
     & emis_cs, iwco_cs, &
     & pres_cs, dz_cs, rhodz_cs, rad_cs, rneb_ord, &
     & rnebmax, stot, shc, shist, &
     & sCb, sThCi, sAnv, &
     & emis, pcld, tcld, iwp, deltaz, rad)

        INTEGER, intent(in) :: len_cs
        REAL, DIMENSION(len_cs), intent(in) :: rneb_cs, temp_cs
        REAL, DIMENSION(len_cs), intent(in) :: emis_cs, iwco_cs, &
     & rneb_ord
        REAL, DIMENSION(len_cs), intent(in) :: pres_cs, dz_cs, rad_cs
        REAL, DIMENSION(len_cs), intent(in) :: rhodz_cs
        REAL, DIMENSION(len_cs) :: tau_cs, w
        REAL, intent(in) :: rnebmax
        REAL, intent(inout) :: stot, shc, shist
        REAL, intent(inout) :: sCb, sThCi, sAnv
        REAL, intent(out) :: emis, pcld, tcld, iwp, deltaz, rad

        INTEGER :: i, ideb, ibeg, iend, nuage, visible
        REAL :: som, som_tau, som_iwc, som_dz, som_rad, tau

        CHARACTER (len = 50)      :: modname = 'simu_airs.sous_section'
        CHARACTER (len = 160)     :: abort_message


! Ponderation: 1 pour les nuages, 0 pour les trous

        do i = 1, len_cs
        if (rneb_cs(i) .ge. rnebmax) then
        w(i) = 1.
        else
        w(i) = 0.
        endif
        enddo

! Calcul des epaisseurs optiques a partir des emissivites

        som = 0.
        do i = 1, len_cs
        if (emis_cs(i) .eq. 1.) then
        tau_cs(i) = 10.
        else
        tau_cs(i) = -log(1.-emis_cs(i))
        endif
        som = som+tau_cs(i)
        enddo


        ideb = 1
        nuage = 0
        visible = 0


! Boucle sur les nuages
        do while (ideb .ne. 0 .and. ideb .le. len_cs)   


! Definition d'un nuage de la sous-section

        call topbot(ideb, w, ibeg, iend)
        ideb = iend+1

        if (ibeg .gt. 0) then

        nuage = nuage + 1

! On determine les caracteristiques du nuage
! (ep. optique, ice water path, pression, temperature)

        call caract(ibeg, iend, temp_cs, tau_cs, iwco_cs, &
     & pres_cs, dz_cs, rhodz_cs, rad_cs, pcld, tcld, &
     & som_tau, som_iwc, som_dz, som_rad)

! On masque le nuage s'il n'est pas detectable

        call masque(ibeg, iend, som_tau, visible, w)

        endif

! Fin boucle nuages
        enddo


! Analyse du nuage detecte

        call topbot(1, w, ibeg, iend)

        if (ibeg .gt. 0) then

        call caract(ibeg, iend, temp_cs, tau_cs, iwco_cs, &
     & pres_cs, dz_cs, rhodz_cs, rad_cs, pcld, tcld, &
     & som_tau, som_iwc, som_dz, som_rad)

        tau = som_tau
        emis = 1. - exp(-tau)
        iwp = som_iwc
        deltaz = som_dz
        rad = som_rad

        if (pcld .gt. p_thresh) then

        shc = 0.
        shist = 0.
        sCb = 0.
        sThCi = 0.
        sAnv = 0.

        else

        if (emis .lt. em_min .or. emis .gt. em_max  &
     & .or. tcld .gt. 230.) then
        shist = 0.
        endif

        if (emis .lt. 0.98) then
        sCb = 0.
        endif

        if (emis .gt. 0.5 .or. emis .lt. 0.1) then
        sThCi = 0.
        endif

        if (emis .le. 0.5 .or. emis .ge. 0.98) then
        sAnv = 0.
        endif

        endif

        else

        tau = 0.
        emis = 0.
        iwp = 0.
        deltaz = 0.
        pcld = 0.
        tcld = 0.
        stot = 0.
        shc = 0.
        shist = 0.
        rad = 0.
        sCb = 0.
        sThCi = 0.
        sAnv = 0.

        endif


! Tests

        if (iwp .lt. 0.) then
          WRITE(abort_message,*) 'ideb iwp =', ideb, iwp
          CALL abort_physic(modname,abort_message,1)
        endif

        if (deltaz .lt. 0.) then
          WRITE(abort_message,*)'ideb deltaz =', ideb, deltaz
          CALL abort_physic(modname,abort_message,1)
        endif

        if (emis .lt. 0.048 .and. emis .ne. 0.) then
          WRITE(abort_message,*) 'ideb emis =', ideb, emis
          CALL abort_physic(modname,abort_message,1)
        endif

        end subroutine sous_section


        subroutine masque (ibeg, iend, som_tau, &
     & visible, w)

        INTEGER, intent(in) :: ibeg, iend
        REAL, intent(in) :: som_tau

        INTEGER, intent(inout) :: visible
        REAL, DIMENSION(:), intent(inout) :: w

        INTEGER :: i



! Masque

! Cas ou il n'y a pas de nuage visible au-dessus

        if (visible .eq. 0) then

        if (som_tau .lt. tau_thresh) then
        do i = ibeg, iend
        w(i) = 0.
        enddo
        else
        visible = 1
        endif

! Cas ou il y a un nuage visible au-dessus

        else

        do i = ibeg, iend
        w(i) = 0.
        enddo

        endif


        end subroutine masque


         subroutine caract (ibeg, iend, temp_cs, tau_cs, iwco_cs, &
     & pres_cs, dz_cs, rhodz_cs, rad_cs, pcld, tcld, &
     & som_tau, som_iwc, som_dz, som_rad)

        INTEGER, intent(in) :: ibeg, iend
        REAL, DIMENSION(:), intent(in) :: tau_cs, iwco_cs, temp_cs
        REAL, DIMENSION(:), intent(in) :: pres_cs, dz_cs, rad_cs
        REAL, DIMENSION(:), intent(in) :: rhodz_cs
        REAL, intent(out) :: som_tau, som_iwc, som_dz, som_rad
        REAL , intent(out) :: pcld, tcld

        INTEGER :: i, ibase, imid

        CHARACTER (len = 50)      :: modname = 'simu_airs.caract'
        CHARACTER (len = 160)     :: abort_message

! Somme des epaisseurs optiques et des contenus en glace sur le nuage

        som_tau = 0.
        som_iwc = 0.
        som_dz = 0.
        som_rad = 0.
        ibase = -100

        do i = ibeg, iend

        som_tau = som_tau + tau_cs(i)

        som_dz = som_dz + dz_cs(i)
        som_iwc = som_iwc + iwco_cs(i)*1000*rhodz_cs(i)  ! en g/m2
        som_rad = som_rad + rad_cs(i)*dz_cs(i)

        if (som_tau .gt. 3. .and. ibase .eq. -100) then
        ibase = i-1
        endif

        enddo

        if (som_dz .ne. 0.) then
          som_rad = som_rad/som_dz
        else
          write(*,*) 'som_dez = 0 STOP'
          write(*,*) 'ibeg, iend =', ibeg, iend
          do i = ibeg, iend
             write(*,*) dz_cs(i), rhodz_cs(i)
          enddo
          abort_message='see above'
          CALL abort_physic(modname,abort_message,1)
        endif

! Determination de Pcld et Tcld

       if (ibase .lt. ibeg) then
       ibase = ibeg
       endif

       imid = (ibeg+ibase)/2

       pcld = pres_cs(imid)/100.        ! pcld en hPa
       tcld = temp_cs(imid)


       end subroutine caract
 
        subroutine topbot(ideb,w,ibeg,iend)

        INTEGER, intent(in) :: ideb
        REAL, DIMENSION(:), intent(in) :: w
        INTEGER, intent(out) :: ibeg, iend

        INTEGER :: i, itest

        itest = 0
        ibeg = 0
        iend = 0

        do i = ideb, size(w)

        if (w(i) .eq. 1. .and. itest .eq. 0) then
        ibeg = i
        itest = 1
        endif

        enddo


        i = ibeg
        do while (w(i) .eq. 1. .and. i .le. size(w))
        i = i+1
        enddo
        iend = i-1


        end subroutine topbot

        subroutine ordonne(len_cs, rneb_cs, rneb_ord)

        INTEGER, intent(in) :: len_cs
        REAL, DIMENSION(:), intent(in) :: rneb_cs
        REAL, DIMENSION(:), intent(out) :: rneb_ord

        INTEGER :: i, j, ind_min

        REAL, DIMENSION(len_cs) :: rneb
        REAL :: rneb_max


        do i = 1, size(rneb_cs)
        rneb(i) = rneb_cs(i)
        enddo

        do j = 1, size(rneb_cs)

        rneb_max = 100.

        do i = 1, size(rneb_cs)
        if (rneb(i) .lt. rneb_max) then
        rneb_max = rneb(i)
        ind_min = i
        endif
        enddo

        rneb_ord(j) = rneb_max
        rneb(ind_min) = 100.

        enddo
        
        end subroutine ordonne

 
        subroutine sim_mesh(rneb_1D, temp_1D, emis_1D, &
     & iwcon_1D, rad_1D, &
     & pres, dz, &
     & rhodz_1D, cc_tot_mesh, cc_hc_mesh, cc_hist_mesh, pcld_hc_mesh,&
     & tcld_hc_mesh, &
     & em_hc_mesh, iwp_hc_mesh, deltaz_hc_mesh, &
     & cc_Cb_mesh, cc_ThCi_mesh, cc_Anv_mesh, &
     & pcld_Cb_mesh, tcld_Cb_mesh, em_Cb_mesh, &
     & pcld_ThCi_mesh, tcld_ThCi_mesh, em_ThCi_mesh, &
     & pcld_Anv_mesh, tcld_Anv_mesh, em_Anv_mesh, &
     & em_hist_mesh, iwp_hist_mesh, deltaz_hist_mesh, rad_hist_mesh)

       USE dimphy

       REAL, DIMENSION(klev), intent(in) :: rneb_1D, temp_1D, emis_1D, &
     & iwcon_1D, rad_1D
        REAL, DIMENSION(klev), intent(in) :: pres, dz, rhodz_1D
        REAL, intent(out) :: cc_tot_mesh, cc_hc_mesh, cc_hist_mesh
        REAL, intent(out) :: cc_Cb_mesh, cc_ThCi_mesh, cc_Anv_mesh
        REAL, intent(out) :: em_hc_mesh, pcld_hc_mesh, tcld_hc_mesh, &
     & iwp_hc_mesh

        REAL, intent(out) :: pcld_Cb_mesh, tcld_Cb_mesh, em_Cb_mesh
        REAL, intent(out) :: pcld_ThCi_mesh, tcld_ThCi_mesh, &
     & em_ThCi_mesh
        REAL, intent(out) :: pcld_Anv_mesh, tcld_Anv_mesh, em_Anv_mesh

        REAL, intent(out) :: em_hist_mesh, iwp_hist_mesh, rad_hist_mesh
        REAL, intent(out) :: deltaz_hc_mesh, deltaz_hist_mesh

        REAL, DIMENSION(:), allocatable :: rneb_cs, temp_cs, emis_cs, &
     & iwco_cs
        REAL, DIMENSION(:), allocatable :: pres_cs, dz_cs, rad_cs, &
     & rhodz_cs

        INTEGER :: i,j,l
        INTEGER :: ltop, itop, ibot, num_cs, N_cs, len_cs, ics

        REAL :: som_emi_hc,som_pcl_hc,som_tcl_hc,som_iwp_hc,som_hc,&
     & som_hist
        REAL :: som_emi_hist, som_iwp_hist, som_deltaz_hc, &
     & som_deltaz_hist
        REAL :: som_rad_hist
        REAL :: som_Cb, som_ThCi, som_Anv
        REAL :: som_emi_Cb, som_tcld_Cb, som_pcld_Cb
        REAL :: som_emi_Anv, som_tcld_Anv, som_pcld_Anv
        REAL :: som_emi_ThCi, som_tcld_ThCi, som_pcld_ThCi
        REAL :: tsom_tot, tsom_hc, tsom_hist
        REAL :: prod, prod_hh

        REAL :: cc_tot_cs, cc_hc_cs, cc_hist_cs
        REAL :: cc_Cb_cs, cc_ThCi_cs, cc_Anv_cs
        REAL :: pcld_hc_cs, tcld_hc_cs
        REAL :: em_hc_cs, iwp_hc_cs, deltaz_hc_cs
        REAL :: pcld_Cb_cs, tcld_Cb_cs, em_Cb_cs
        REAL :: pcld_ThCi_cs, tcld_ThCi_cs, em_ThCi_cs
        REAL :: pcld_Anv_cs, tcld_Anv_cs, em_Anv_cs
        REAL :: em_hist_cs, iwp_hist_cs, deltaz_hist_cs, rad_hist_cs

        REAL, DIMENSION(klev) :: test_tot, test_hc, test_hist
        REAL, DIMENSION(klev) :: test_pcld, test_tcld, test_em, test_iwp

        CHARACTER (len = 50)      :: modname = 'simu_airs.sim_mesh'
        CHARACTER (len = 160)     :: abort_message
        

        do j = 1, klev
          WRITE(lunout,*) 'simu_airs, j, rneb_1D =', rneb_1D(j)
        enddo

! Definition des structures nuageuses, de la plus haute a la plus basse

        num_cs = 0
        ltop = klev-1

        prod = 1.

        som_emi_hc = 0.
        som_emi_hist = 0.
        som_pcl_hc = 0.
        som_tcl_hc = 0.
        som_iwp_hc = 0.
        som_iwp_hist = 0.
        som_deltaz_hc = 0.
        som_deltaz_hist = 0.
        som_rad_hist = 0.
        som_hc = 0.
        som_hist = 0.

        som_Cb = 0.
        som_ThCi = 0.
        som_Anv = 0.

        som_emi_Cb = 0.
        som_pcld_Cb = 0.
        som_tcld_Cb = 0.

        som_emi_ThCi = 0.
        som_pcld_ThCi = 0.
        som_tcld_ThCi = 0.

        som_emi_Anv = 0.
        som_pcld_Anv = 0.
        som_tcld_Anv = 0.

        tsom_tot = 0.
        tsom_hc = 0.
        tsom_hist = 0.

        do while (ltop .ge. 1)   ! Boucle sur toute la colonne

        itop = 0

        do l = ltop,1,-1

        if (itop .eq. 0 .and. rneb_1D(l) .gt. 0.001 ) then
        itop = l
        endif

        enddo

        ibot = itop

        do while (rneb_1D(ibot) .gt. 0.001 .and. ibot .ge. 1)
        ibot = ibot -1
        enddo


        ibot = ibot+1

        if (itop .gt. 0) then    ! itop > 0

        num_cs = num_cs +1
        len_cs = itop-ibot+1

! Allocation et definition des variables de la structure nuageuse
! le premier indice denote ce qui est le plus haut

        allocate (rneb_cs(len_cs))
        allocate (temp_cs(len_cs))
        allocate (emis_cs(len_cs))
        allocate (iwco_cs(len_cs))
        allocate (pres_cs(len_cs))
        allocate (dz_cs(len_cs))
        allocate (rad_cs(len_cs))
        allocate (rhodz_cs(len_cs))

        ics = 0

        do i = itop, ibot, -1
        ics = ics + 1
        rneb_cs(ics) = rneb_1D(i)
        temp_cs(ics) = temp_1D(i)
        emis_cs(ics) = emis_1D(i)
        iwco_cs(ics) = iwcon_1D(i)
        rad_cs(ics) = rad_1D(i)
        pres_cs(ics) = pres(i)
        dz_cs(ics) = dz(i)
        rhodz_cs(ics) = rhodz_1D(i)
        enddo

! Appel du sous_programme cloud_structure

        call cloud_structure(len_cs,rneb_cs,temp_cs,emis_cs,iwco_cs,&
     & pres_cs, dz_cs, rhodz_cs, rad_cs, &
     & cc_tot_cs, cc_hc_cs, cc_hist_cs, &
     & cc_Cb_cs, cc_ThCi_cs, cc_Anv_cs, &
     & pcld_hc_cs, tcld_hc_cs, &
     & em_hc_cs, iwp_hc_cs, deltaz_hc_cs, &
     & pcld_Cb_cs, tcld_Cb_cs, em_Cb_cs, &
     & pcld_ThCi_cs, tcld_ThCi_cs, em_ThCi_cs, &
     & pcld_Anv_cs, tcld_Anv_cs, em_Anv_cs, &
     & em_hist_cs, iwp_hist_cs, deltaz_hist_cs, rad_hist_cs)


        deallocate (rneb_cs)
        deallocate (temp_cs)
        deallocate (emis_cs)
        deallocate (iwco_cs)
        deallocate (pres_cs)
        deallocate (dz_cs)
        deallocate (rad_cs)
        deallocate (rhodz_cs)


! Pour la couverture nuageuse sur la maille

        prod_hh = prod

        prod = prod*(1.-cc_tot_cs)

! Pour les autres variables definies sur la maille

        som_emi_hc = som_emi_hc + em_hc_cs*cc_hc_cs*prod_hh
        som_iwp_hc = som_iwp_hc + iwp_hc_cs*cc_hc_cs*prod_hh
        som_deltaz_hc = som_deltaz_hc + deltaz_hc_cs*cc_hc_cs*prod_hh

        som_emi_Cb = som_emi_Cb + em_Cb_cs*cc_Cb_cs*prod_hh
        som_tcld_Cb = som_tcld_Cb + tcld_Cb_cs*cc_Cb_cs*prod_hh
        som_pcld_Cb = som_pcld_Cb + pcld_Cb_cs*cc_Cb_cs*prod_hh

        som_emi_ThCi = som_emi_ThCi + em_ThCi_cs*cc_ThCi_cs*prod_hh
        som_tcld_ThCi = som_tcld_ThCi + tcld_ThCi_cs*cc_ThCi_cs*prod_hh
        som_pcld_ThCi = som_pcld_ThCi + pcld_ThCi_cs*cc_ThCi_cs*prod_hh

        som_emi_Anv = som_emi_Anv + em_Anv_cs*cc_Anv_cs*prod_hh
        som_tcld_Anv = som_tcld_Anv + tcld_Anv_cs*cc_Anv_cs*prod_hh
        som_pcld_Anv = som_pcld_Anv + pcld_Anv_cs*cc_Anv_cs*prod_hh

        som_emi_hist = som_emi_hist + em_hist_cs*cc_hist_cs*prod_hh
        som_iwp_hist = som_iwp_hist + iwp_hist_cs*cc_hist_cs*prod_hh
        som_deltaz_hist = som_deltaz_hist + &
     & deltaz_hist_cs*cc_hist_cs*prod_hh
        som_rad_hist = som_rad_hist + rad_hist_cs*cc_hist_cs*prod_hh

        som_pcl_hc = som_pcl_hc + pcld_hc_cs*cc_hc_cs*prod_hh
        som_tcl_hc = som_tcl_hc + tcld_hc_cs*cc_hc_cs*prod_hh

        som_hc = som_hc + cc_hc_cs*prod_hh
        som_hist = som_hist + cc_hist_cs*prod_hh

        som_Cb = som_Cb + cc_Cb_cs*prod_hh
        som_ThCi = som_ThCi + cc_ThCi_cs*prod_hh
        som_Anv = som_Anv + cc_Anv_cs*prod_hh


! Pour test
        
        call test_bornes('cc_tot_cs     ',cc_tot_cs,1.,0.)
        call test_bornes('cc_hc_cs      ',cc_hc_cs,1.,0.)
        call test_bornes('cc_hist_cs    ',cc_hist_cs,1.,0.)
        call test_bornes('pcld_hc_cs    ',pcld_hc_cs,1200.,0.)
        call test_bornes('tcld_hc_cs    ',tcld_hc_cs,1000.,100.)
        call test_bornes('em_hc_cs      ',em_hc_cs,1000.,0.048)

        test_tot(num_cs) = cc_tot_cs
        test_hc(num_cs) = cc_hc_cs
        test_hist(num_cs) = cc_hist_cs
        test_pcld(num_cs) = pcld_hc_cs
        test_tcld(num_cs) = tcld_hc_cs
        test_em(num_cs) = em_hc_cs
        test_iwp(num_cs) = iwp_hc_cs

        tsom_tot = tsom_tot + cc_tot_cs
        tsom_hc = tsom_hc + cc_hc_cs
        tsom_hist = tsom_hist + cc_hist_cs


        endif                   ! itop > 0

        ltop = ibot -1

        enddo                   ! fin de la boucle sur la colonne

        N_CS = num_cs


! Determination des variables de sortie

        if (N_CS .gt. 0) then   ! if N_CS>0

        cc_tot_mesh = 1. - prod

        cc_hc_mesh = som_hc
        cc_hist_mesh = som_hist

        cc_Cb_mesh = som_Cb
        cc_ThCi_mesh = som_ThCi
        cc_Anv_mesh = som_Anv

        call normal2_undef(pcld_hc_mesh,som_pcl_hc, &
     & cc_hc_mesh)
        call normal2_undef(tcld_hc_mesh,som_tcl_hc, &
     & cc_hc_mesh)
        call normal2_undef(em_hc_mesh,som_emi_hc, &
     & cc_hc_mesh)
        call normal2_undef(iwp_hc_mesh,som_iwp_hc, &
     & cc_hc_mesh)
        call normal2_undef(deltaz_hc_mesh,som_deltaz_hc, &
     & cc_hc_mesh)

        call normal2_undef(em_Cb_mesh,som_emi_Cb, &
     & cc_Cb_mesh)
        call normal2_undef(tcld_Cb_mesh,som_tcld_Cb, &
     & cc_Cb_mesh)
        call normal2_undef(pcld_Cb_mesh,som_pcld_Cb, &
     & cc_Cb_mesh)

        call normal2_undef(em_ThCi_mesh,som_emi_ThCi, &
     & cc_ThCi_mesh)
        call normal2_undef(tcld_ThCi_mesh,som_tcld_ThCi, &
     & cc_ThCi_mesh)
        call normal2_undef(pcld_ThCi_mesh,som_pcld_ThCi, &
     & cc_ThCi_mesh)

       call normal2_undef(em_Anv_mesh,som_emi_Anv, &
     & cc_Anv_mesh)
        call normal2_undef(tcld_Anv_mesh,som_tcld_Anv, &
     & cc_Anv_mesh)
        call normal2_undef(pcld_Anv_mesh,som_pcld_Anv, &
     & cc_Anv_mesh)


        call normal2_undef(em_hist_mesh,som_emi_hist, &
     & cc_hist_mesh)
        call normal2_undef(iwp_hist_mesh,som_iwp_hist, &
     & cc_hist_mesh)
        call normal2_undef(deltaz_hist_mesh,som_deltaz_hist, &
     & cc_hist_mesh)
        call normal2_undef(rad_hist_mesh,som_rad_hist, &
     & cc_hist_mesh)


! Tests 

        ! Tests

       if (cc_tot_mesh .gt. tsom_tot .and. &
     & abs(cc_tot_mesh-tsom_tot) .gt. 1.e-4) then
           WRITE(abort_message,*)'cc_tot_mesh > tsom_tot', cc_tot_mesh, tsom_tot
           CALL abort_physic(modname,abort_message,1)
        endif

        if (cc_tot_mesh .lt. maxval(test_tot(1:N_CS)) .and. &
     & abs(cc_tot_mesh-maxval(test_tot(1:N_CS))) .gt. 1.e-4) then
           WRITE(abort_message,*) 'cc_tot_mesh < max', cc_tot_mesh, maxval(test_tot(1:N_CS))
           CALL abort_physic(modname,abort_message,1)
        endif

        if (cc_hc_mesh .gt. tsom_hc .and. &
     & abs(cc_hc_mesh-tsom_hc) .gt. 1.e-4) then
           WRITE(abort_message,*) 'cc_hc_mesh > tsom_hc', cc_hc_mesh, tsom_hc
           CALL abort_physic(modname,abort_message,1)
        endif

        if (cc_hc_mesh .lt. maxval(test_hc(1:N_CS)) .and. &
     & abs(cc_hc_mesh-maxval(test_hc(1:N_CS))) .gt. 1.e-4) then
           WRITE(abort_message,*) 'cc_hc_mesh < max', cc_hc_mesh, maxval(test_hc(1:N_CS))
           CALL abort_physic(modname,abort_message,1)
        endif

        if (cc_hist_mesh .gt. tsom_hist .and. &
     & abs(cc_hist_mesh-tsom_hist) .gt. 1.e-4) then
           WRITE(abort_message,*) 'cc_hist_mesh > tsom_hist', cc_hist_mesh, tsom_hist
           CALL abort_physic(modname,abort_message,1)
        endif

        if (cc_hist_mesh .lt. 0.) then
           WRITE(abort_message,*) 'cc_hist_mesh < 0', cc_hist_mesh
           CALL abort_physic(modname,abort_message,1)
        endif

        if ((pcld_hc_mesh .gt. maxval(test_pcld(1:N_CS)) .or. &
     & pcld_hc_mesh .lt. minval(test_pcld(1:N_CS))) .and. &
     & abs(pcld_hc_mesh-maxval(test_pcld(1:N_CS))) .gt. 1. .and. &
     & maxval(test_pcld(1:N_CS)) .ne. 999. &
     & .and. minval(test_pcld(1:N_CS)) .ne. 999.) then
           WRITE(abort_message,*) 'pcld_hc_mesh est faux', pcld_hc_mesh, maxval(test_pcld(1:N_CS)), &
     & minval(test_pcld(1:N_CS))
           CALL abort_physic(modname,abort_message,1)
        endif

       if ((tcld_hc_mesh .gt. maxval(test_tcld(1:N_CS)) .or. &
     & tcld_hc_mesh .lt. minval(test_tcld(1:N_CS))) .and. &
     & abs(tcld_hc_mesh-maxval(test_tcld(1:N_CS))) .gt. 0.1 .and. &
     & maxval(test_tcld(1:N_CS)) .ne. 999. &
     & .and. minval(test_tcld(1:N_CS)) .ne. 999.) then
           WRITE(abort_message,*) 'tcld_hc_mesh est faux', tcld_hc_mesh, maxval(test_tcld(1:N_CS)), &
                & minval(test_tcld(1:N_CS))
           CALL abort_physic(modname,abort_message,1)
        endif

        if ((em_hc_mesh .gt. maxval(test_em(1:N_CS)) .or. &
     & em_hc_mesh .lt. minval(test_em(1:N_CS))) .and. &
     & abs(em_hc_mesh-maxval(test_em(1:N_CS))) .gt. 1.e-4 .and. &
     & minval(test_em(1:N_CS)) .ne. 999. .and. &
     & maxval(test_em(1:N_CS)) .ne. 999. ) then
           WRITE(abort_message,*) 'em_hc_mesh est faux', em_hc_mesh, maxval(test_em(1:N_CS)), &
     & minval(test_em(1:N_CS))
           CALL abort_physic(modname,abort_message,1)
        endif

        else               ! if N_CS>0

        cc_tot_mesh = 0.
        cc_hc_mesh = 0.
        cc_hist_mesh = 0.

        cc_Cb_mesh = 0.
        cc_ThCi_mesh = 0.
        cc_Anv_mesh = 0.

        iwp_hc_mesh = undef 
        deltaz_hc_mesh = undef 
        em_hc_mesh = undef 
        iwp_hist_mesh = undef 
        deltaz_hist_mesh = undef 
        rad_hist_mesh = undef 
        em_hist_mesh = undef 
        pcld_hc_mesh = undef 
        tcld_hc_mesh = undef 

        pcld_Cb_mesh = undef 
        tcld_Cb_mesh = undef 
        em_Cb_mesh = undef 

        pcld_ThCi_mesh = undef 
        tcld_ThCi_mesh = undef 
        em_ThCi_mesh = undef 

        pcld_Anv_mesh = undef 
        tcld_Anv_mesh = undef 
        em_Anv_mesh = undef 


        endif                  ! if N_CS>0

        end subroutine sim_mesh

        subroutine test_bornes(sx,x,bsup,binf)

        REAL, intent(in) :: x, bsup, binf
        character*14, intent(in) :: sx
        CHARACTER (len = 50)      :: modname = 'simu_airs.test_bornes'
        CHARACTER (len = 160)     :: abort_message

        if (x .gt. bsup .or. x .lt. binf) then
          WRITE(abort_message,*) sx, 'est faux', sx, x
          CALL abort_physic(modname,abort_message,1)
        endif
 
        end subroutine test_bornes

        end module m_simu_airs


        subroutine simu_airs &
     & (itap, rneb_airs, temp_airs, cldemi_airs, iwcon0_airs, rad_airs, &
     & geop_airs, pplay_airs, paprs_airs, &
     & map_prop_hc,map_prop_hist,&
     & map_emis_hc,map_iwp_hc,map_deltaz_hc,map_pcld_hc,map_tcld_hc,&
     & map_emis_Cb,map_pcld_Cb,map_tcld_Cb,&
     & map_emis_ThCi,map_pcld_ThCi,map_tcld_ThCi,&
     & map_emis_Anv,map_pcld_Anv,map_tcld_Anv,&
     & map_emis_hist,map_iwp_hist,map_deltaz_hist,map_rad_hist,&
     & map_ntot,map_hc,map_hist,&
     & map_Cb,map_ThCi,map_Anv,alt_tropo )

        USE dimphy
        USE m_simu_airs

        IMPLICIT NONE

        include "YOMCST.h"

        INTEGER,intent(in) :: itap

        REAL, DIMENSION(klon,klev), intent(in) :: &
     & rneb_airs, temp_airs, cldemi_airs, iwcon0_airs, &
     & rad_airs, geop_airs, pplay_airs, paprs_airs

       REAL, DIMENSION(klon,klev) :: &
     & rhodz_airs, rho_airs, iwcon_airs

        REAL, DIMENSION(klon),intent(out) :: alt_tropo

        REAL, DIMENSION(klev) :: rneb_1D, temp_1D, &
     & emis_1D, rad_1D, pres_1D, alt_1D, &
     & rhodz_1D, dz_1D, iwcon_1D

        INTEGER :: i, j

        REAL :: cc_tot_mesh, cc_hc_mesh, cc_hist_mesh
        REAL :: cc_Cb_mesh, cc_ThCi_mesh, cc_Anv_mesh
        REAL :: pcld_hc_mesh, tcld_hc_mesh, em_hc_mesh, iwp_hc_mesh
        REAL :: em_hist_mesh, iwp_hist_mesh
        REAL :: deltaz_hc_mesh, deltaz_hist_mesh, rad_hist_mesh
        REAL :: pcld_Cb_mesh, tcld_Cb_mesh, em_Cb_mesh
        REAL :: pcld_ThCi_mesh, tcld_ThCi_mesh, em_ThCi_mesh
        REAL :: pcld_Anv_mesh, tcld_Anv_mesh, em_Anv_mesh

        REAL, DIMENSION(klon),intent(out) :: map_prop_hc, map_prop_hist
        REAL, DIMENSION(klon),intent(out) :: map_emis_hc, map_iwp_hc
        REAL, DIMENSION(klon),intent(out) :: map_deltaz_hc, map_pcld_hc
        REAL, DIMENSION(klon),intent(out) :: map_tcld_hc
        REAL, DIMENSION(klon),intent(out) :: map_emis_Cb,map_pcld_Cb,map_tcld_Cb 
        REAL, DIMENSION(klon),intent(out) :: &
     & map_emis_ThCi,map_pcld_ThCi,map_tcld_ThCi
        REAL, DIMENSION(klon),intent(out) :: &
     & map_emis_Anv,map_pcld_Anv,map_tcld_Anv
        REAL, DIMENSION(klon),intent(out) :: &
     & map_emis_hist,map_iwp_hist,map_deltaz_hist,&
     & map_rad_hist
        REAL, DIMENSION(klon),intent(out) :: map_ntot,map_hc,map_hist
        REAL, DIMENSION(klon),intent(out) :: map_Cb,map_ThCi,map_Anv
 
 
        write(*,*) 'simu_airs'
        write(*,*) 'itap, klon, klev', itap, klon, klev
        write(*,*) 'RG, RD =', RG, RD


! Definition des variables 1D

        do i = 1, klon
        do j = 1, klev-1
        rhodz_airs(i,j) = &
     & (paprs_airs(i,j)-paprs_airs(i,j+1))/RG
        enddo
        rhodz_airs(i,klev) = 0.
        enddo

        do i = 1, klon
        do j = 1,klev
        rho_airs(i,j) = &
     & pplay_airs(i,j)/(temp_airs(i,j)*RD)

        if (rneb_airs(i,j) .gt. 0.001) then
        iwcon_airs(i,j) = iwcon0_airs(i,j)/rneb_airs(i,j)
        else
        iwcon_airs(i,j) = 0.
        endif
 
        enddo
        enddo

!=============================================================================

        do i = 1, klon  ! boucle sur les points de grille

!=============================================================================
        
        do j = 1,klev

        rneb_1D(j) = rneb_airs(i,j)
        temp_1D(j) = temp_airs(i,j)
        emis_1D(j) = cldemi_airs(i,j)
        iwcon_1D(j) = iwcon_airs(i,j)
        rad_1D(j) = rad_airs(i,j)
        pres_1D(j) = pplay_airs(i,j)
        alt_1D(j) = geop_airs(i,j)/RG
        rhodz_1D(j) = rhodz_airs(i,j)
        dz_1D(j) = rhodz_airs(i,j)/rho_airs(i,j)

        enddo

        alt_tropo(i) = &
     & search_tropopause(pres_1D/100.,temp_1D,alt_1D,klev)


! Appel du ss-programme sim_mesh

!        if (itap .eq. 1 ) then

        call sim_mesh(rneb_1D, temp_1D, emis_1D, iwcon_1D, rad_1D, &
     & pres_1D, dz_1D, rhodz_1D, &
     & cc_tot_mesh, cc_hc_mesh, cc_hist_mesh, &
     & pcld_hc_mesh, tcld_hc_mesh, em_hc_mesh, iwp_hc_mesh, &
     & deltaz_hc_mesh,&
     & cc_Cb_mesh, cc_ThCi_mesh, cc_Anv_mesh, &
     & pcld_Cb_mesh, tcld_Cb_mesh, em_Cb_mesh, &
     & pcld_ThCi_mesh, tcld_ThCi_mesh, em_ThCi_mesh, &
     & pcld_Anv_mesh, tcld_Anv_mesh, em_Anv_mesh, &
     & em_hist_mesh, iwp_hist_mesh, deltaz_hist_mesh, rad_hist_mesh)

         write(*,*) '===================================='
         write(*,*) 'itap, i:', itap, i 
         write(*,*) 'cc_tot, cc_hc, cc_hist, pcld_hc, tcld_hc, em_hc, &
     & iwp_hc, em_hist, iwp_hist ='
         write(*,*) cc_tot_mesh, cc_hc_mesh, cc_hist_mesh
         write(*,*) pcld_hc_mesh, tcld_hc_mesh, em_hc_mesh, iwp_hc_mesh
         write(*,*)  em_hist_mesh, iwp_hist_mesh

!        endif

! Definition des variables a ecrire dans le fichier de sortie

        call normal2_undef(map_prop_hc(i),cc_hc_mesh, &
     & cc_tot_mesh)
        call normal2_undef(map_prop_hist(i),cc_hist_mesh, &
     & cc_tot_mesh)

       map_emis_hc(i) = em_hc_mesh
       map_iwp_hc(i) = iwp_hc_mesh
       map_deltaz_hc(i) = deltaz_hc_mesh
       map_pcld_hc(i) = pcld_hc_mesh
       map_tcld_hc(i) = tcld_hc_mesh

       map_emis_Cb(i) = em_Cb_mesh
       map_pcld_Cb(i) = pcld_Cb_mesh
       map_tcld_Cb(i) = tcld_Cb_mesh

       map_emis_ThCi(i) = em_ThCi_mesh
       map_pcld_ThCi(i) = pcld_ThCi_mesh
       map_tcld_ThCi(i) = tcld_ThCi_mesh

       map_emis_Anv(i) = em_Anv_mesh
       map_pcld_Anv(i) = pcld_Anv_mesh
       map_tcld_Anv(i) = tcld_Anv_mesh

       map_emis_hist(i) = em_hist_mesh
       map_iwp_hist(i) = iwp_hist_mesh
       map_deltaz_hist(i) = deltaz_hist_mesh
       map_rad_hist(i) = rad_hist_mesh

       map_ntot(i) = cc_tot_mesh
       map_hc(i) = cc_hc_mesh
       map_hist(i) = cc_hist_mesh

       map_Cb(i) = cc_Cb_mesh
       map_ThCi(i) = cc_ThCi_mesh
       map_Anv(i) = cc_Anv_mesh


        enddo         ! fin boucle sur les points de grille

        

        end subroutine simu_airs

