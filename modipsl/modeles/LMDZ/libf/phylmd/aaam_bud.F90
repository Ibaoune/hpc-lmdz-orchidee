
! $Id: aaam_bud.F90 2350 2015-08-25 11:40:19Z emillour $

SUBROUTINE aaam_bud(iam, nlon, nlev, rjour, rsec, rea, rg, ome, plat, plon, &
    phis, dragu, liftu, phyu, dragv, liftv, phyv, p, u, v, aam, torsfc)

  USE dimphy
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, klon_glo
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): F.Lott (LMD/CNRS) date: 20031020
  ! Object: Compute different terms of the axial AAAM Budget.
  ! No outputs, every AAM quantities are written on the IAM
  ! File.

  ! Modif : I.Musat (LMD/CNRS) date : 20041020
  ! Outputs : axial components of wind AAM "aam" and total surface torque
  ! "torsfc",
  ! but no write in the iam file.

  ! WARNING: Only valid for regular rectangular grids.
  ! REMARK: CALL DANS PHYSIQ AFTER lift_noro:
  ! CALL aaam_bud (27,klon,klev,rjourvrai,gmtime,
  ! C               ra,rg,romega,
  ! C               rlat,rlon,pphis,
  ! C               zustrdr,zustrli,zustrph,
  ! C               zvstrdr,zvstrli,zvstrph,
  ! C               paprs,u,v)

  ! ======================================================================
  ! Explicit Arguments:
  ! ==================
  ! iam-----input-I-File number where AAMs and torques are written
  ! It is a formatted file that has been opened
  ! in physiq.F
  ! nlon----input-I-Total number of horizontal points that get into physics
  ! nlev----input-I-Number of vertical levels
  ! rjour        -R-Jour compte depuis le debut de la simu (run.def)
  ! rsec         -R-Seconde de la journee
  ! rea          -R-Earth radius
  ! rg           -R-gravity constant
  ! ome          -R-Earth rotation rate
  ! plat ---input-R-Latitude en degres
  ! plon ---input-R-Longitude en degres
  ! phis ---input-R-Geopotential at the ground
  ! dragu---input-R-orodrag stress (zonal)
  ! liftu---input-R-orolift stress (zonal)
  ! phyu----input-R-Stress total de la physique (zonal)
  ! dragv---input-R-orodrag stress (Meridional)
  ! liftv---input-R-orolift stress (Meridional)
  ! phyv----input-R-Stress total de la physique (Meridional)
  ! p-------input-R-Pressure (Pa) at model half levels
  ! u-------input-R-Horizontal wind (m/s)
  ! v-------input-R-Meridional wind (m/s)
  ! aam-----output-R-Axial Wind AAM (=raam(3))
  ! torsfc--output-R-Total surface torque (=tmou(3)+tsso(3)+tbls(3))

  ! Implicit Arguments:
  ! ===================

  ! nbp_lon--common-I: Number of longitude intervals
  ! (nbp_lat-1)--common-I: Number of latitude intervals
  ! klon-common-I: Number of points seen by the physics
  ! nbp_lon*(nbp_lat-2)+2 for instance
  ! klev-common-I: Number of vertical layers
  ! ======================================================================
  ! Local Variables:
  ! ================
  ! dlat-----R: Latitude increment (Radians)
  ! dlon-----R: Longitude increment (Radians)
  ! raam  ---R: Wind AAM (3 Components, 1 & 2 Equatoriales; 3 Axiale)
  ! oaam  ---R: Mass AAM (3 Components, 1 & 2 Equatoriales; 3 Axiale)
  ! tmou-----R: Resolved Mountain torque (3 components)
  ! tsso-----R: Parameterised Moutain drag torque (3 components)
  ! tbls-----R: Parameterised Boundary layer torque (3 components)

  ! LOCAL ARRAY:
  ! ===========
  ! zs    ---R: Topographic height
  ! ps    ---R: Surface Pressure
  ! ub    ---R: Barotropic wind zonal
  ! vb    ---R: Barotropic wind meridional
  ! zlat  ---R: Latitude in radians
  ! zlon  ---R: Longitude in radians
  ! ======================================================================

  ! ARGUMENTS

  INTEGER iam, nlon, nlev
  REAL, INTENT (IN) :: rjour, rsec, rea, rg, ome
  REAL plat(nlon), plon(nlon), phis(nlon)
  REAL dragu(nlon), liftu(nlon), phyu(nlon)
  REAL dragv(nlon), liftv(nlon), phyv(nlon)
  REAL p(nlon, nlev+1), u(nlon, nlev), v(nlon, nlev)

  ! Variables locales:

  INTEGER i, j, k, l
  REAL xpi, hadley, hadday
  REAL dlat, dlon
  REAL raam(3), oaam(3), tmou(3), tsso(3), tbls(3)
  INTEGER iax
  ! IM ajout aam, torsfc
  ! aam = composante axiale du Wind AAM raam
  ! torsfc = composante axiale de (tmou+tsso+tbls)
  REAL aam, torsfc

  REAL zs(801, 401), ps(801, 401)
  REAL ub(801, 401), vb(801, 401)
  REAL ssou(801, 401), ssov(801, 401)
  REAL blsu(801, 401), blsv(801, 401)
  REAL zlon(801), zlat(401)

  CHARACTER (LEN=20) :: modname = 'aaam_bud'
  CHARACTER (LEN=80) :: abort_message



  ! PUT AAM QUANTITIES AT ZERO:

  IF (nbp_lon+1>801 .OR. nbp_lat>401) THEN
    abort_message = 'Pb de dimension dans aaam_bud'
    CALL abort_physic(modname, abort_message, 1)
  END IF

  xpi = acos(-1.)
  hadley = 1.E18
  hadday = 1.E18*24.*3600.
  IF(klon_glo.EQ.1) THEN
    dlat = xpi
  ELSE
    dlat = xpi/real(nbp_lat-1)
  ENDIF
  dlon = 2.*xpi/real(nbp_lon)

  DO iax = 1, 3
    oaam(iax) = 0.
    raam(iax) = 0.
    tmou(iax) = 0.
    tsso(iax) = 0.
    tbls(iax) = 0.
  END DO

  ! MOUNTAIN HEIGHT, PRESSURE AND BAROTROPIC WIND:

  ! North pole values (j=1):

  l = 1

  ub(1, 1) = 0.
  vb(1, 1) = 0.
  DO k = 1, nlev
    ub(1, 1) = ub(1, 1) + u(l, k)*(p(l,k)-p(l,k+1))/rg
    vb(1, 1) = vb(1, 1) + v(l, k)*(p(l,k)-p(l,k+1))/rg
  END DO

  zlat(1) = plat(l)*xpi/180.

  DO i = 1, nbp_lon + 1

    zs(i, 1) = phis(l)/rg
    ps(i, 1) = p(l, 1)
    ub(i, 1) = ub(1, 1)
    vb(i, 1) = vb(1, 1)
    ssou(i, 1) = dragu(l) + liftu(l)
    ssov(i, 1) = dragv(l) + liftv(l)
    blsu(i, 1) = phyu(l) - dragu(l) - liftu(l)
    blsv(i, 1) = phyv(l) - dragv(l) - liftv(l)

  END DO


  DO j = 2, nbp_lat-1

    ! Values at Greenwich (Periodicity)

    zs(nbp_lon+1, j) = phis(l+1)/rg
    ps(nbp_lon+1, j) = p(l+1, 1)
    ssou(nbp_lon+1, j) = dragu(l+1) + liftu(l+1)
    ssov(nbp_lon+1, j) = dragv(l+1) + liftv(l+1)
    blsu(nbp_lon+1, j) = phyu(l+1) - dragu(l+1) - liftu(l+1)
    blsv(nbp_lon+1, j) = phyv(l+1) - dragv(l+1) - liftv(l+1)
    zlon(nbp_lon+1) = -plon(l+1)*xpi/180.
    zlat(j) = plat(l+1)*xpi/180.

    ub(nbp_lon+1, j) = 0.
    vb(nbp_lon+1, j) = 0.
    DO k = 1, nlev
      ub(nbp_lon+1, j) = ub(nbp_lon+1, j) + u(l+1, k)*(p(l+1,k)-p(l+1,k+1))/rg
      vb(nbp_lon+1, j) = vb(nbp_lon+1, j) + v(l+1, k)*(p(l+1,k)-p(l+1,k+1))/rg
    END DO


    DO i = 1, nbp_lon

      l = l + 1
      zs(i, j) = phis(l)/rg
      ps(i, j) = p(l, 1)
      ssou(i, j) = dragu(l) + liftu(l)
      ssov(i, j) = dragv(l) + liftv(l)
      blsu(i, j) = phyu(l) - dragu(l) - liftu(l)
      blsv(i, j) = phyv(l) - dragv(l) - liftv(l)
      zlon(i) = plon(l)*xpi/180.

      ub(i, j) = 0.
      vb(i, j) = 0.
      DO k = 1, nlev
        ub(i, j) = ub(i, j) + u(l, k)*(p(l,k)-p(l,k+1))/rg
        vb(i, j) = vb(i, j) + v(l, k)*(p(l,k)-p(l,k+1))/rg
      END DO

    END DO

  END DO


  ! South Pole

  IF (nbp_lat-1>1) THEN
    l = l + 1
    ub(1, nbp_lat) = 0.
    vb(1, nbp_lat) = 0.
    DO k = 1, nlev
      ub(1, nbp_lat) = ub(1, nbp_lat) + u(l, k)*(p(l,k)-p(l,k+1))/rg
      vb(1, nbp_lat) = vb(1, nbp_lat) + v(l, k)*(p(l,k)-p(l,k+1))/rg
    END DO
    zlat(nbp_lat) = plat(l)*xpi/180.

    DO i = 1, nbp_lon + 1
      zs(i, nbp_lat) = phis(l)/rg
      ps(i, nbp_lat) = p(l, 1)
      ssou(i, nbp_lat) = dragu(l) + liftu(l)
      ssov(i, nbp_lat) = dragv(l) + liftv(l)
      blsu(i, nbp_lat) = phyu(l) - dragu(l) - liftu(l)
      blsv(i, nbp_lat) = phyv(l) - dragv(l) - liftv(l)
      ub(i, nbp_lat) = ub(1, nbp_lat)
      vb(i, nbp_lat) = vb(1, nbp_lat)
    END DO
  END IF


  ! MOMENT ANGULAIRE

  DO j = 1, nbp_lat-1
    DO i = 1, nbp_lon

      raam(1) = raam(1) - rea**3*dlon*dlat*0.5*(cos(zlon(i))*sin(zlat(j))*cos &
        (zlat(j))*ub(i,j)+cos(zlon(i))*sin(zlat(j+1))*cos(zlat(j+ &
        1))*ub(i,j+1)) + rea**3*dlon*dlat*0.5*(sin(zlon(i))*cos(zlat(j))*vb(i &
        ,j)+sin(zlon(i))*cos(zlat(j+1))*vb(i,j+1))

      oaam(1) = oaam(1) - ome*rea**4*dlon*dlat/rg*0.5*(cos(zlon(i))*cos(zlat( &
        j))**2*sin(zlat(j))*ps(i,j)+cos(zlon(i))*cos(zlat(j+ &
        1))**2*sin(zlat(j+1))*ps(i,j+1))

      raam(2) = raam(2) - rea**3*dlon*dlat*0.5*(sin(zlon(i))*sin(zlat(j))*cos &
        (zlat(j))*ub(i,j)+sin(zlon(i))*sin(zlat(j+1))*cos(zlat(j+ &
        1))*ub(i,j+1)) - rea**3*dlon*dlat*0.5*(cos(zlon(i))*cos(zlat(j))*vb(i &
        ,j)+cos(zlon(i))*cos(zlat(j+1))*vb(i,j+1))

      oaam(2) = oaam(2) - ome*rea**4*dlon*dlat/rg*0.5*(sin(zlon(i))*cos(zlat( &
        j))**2*sin(zlat(j))*ps(i,j)+sin(zlon(i))*cos(zlat(j+ &
        1))**2*sin(zlat(j+1))*ps(i,j+1))

      raam(3) = raam(3) + rea**3*dlon*dlat*0.5*(cos(zlat(j))**2*ub(i,j)+cos( &
        zlat(j+1))**2*ub(i,j+1))

      oaam(3) = oaam(3) + ome*rea**4*dlon*dlat/rg*0.5*(cos(zlat(j))**3*ps(i,j &
        )+cos(zlat(j+1))**3*ps(i,j+1))

    END DO
  END DO


  ! COUPLE DES MONTAGNES:


  DO j = 1, nbp_lat-1
    DO i = 1, nbp_lon
      tmou(1) = tmou(1) - rea**2*dlon*0.5*sin(zlon(i))*(zs(i,j)-zs(i,j+1))*( &
        cos(zlat(j+1))*ps(i,j+1)+cos(zlat(j))*ps(i,j))
      tmou(2) = tmou(2) + rea**2*dlon*0.5*cos(zlon(i))*(zs(i,j)-zs(i,j+1))*( &
        cos(zlat(j+1))*ps(i,j+1)+cos(zlat(j))*ps(i,j))
    END DO
  END DO

  DO j = 2, nbp_lat-1
    DO i = 1, nbp_lon
      tmou(1) = tmou(1) + rea**2*dlat*0.5*sin(zlat(j))*(zs(i+1,j)-zs(i,j))*( &
        cos(zlon(i+1))*ps(i+1,j)+cos(zlon(i))*ps(i,j))
      tmou(2) = tmou(2) + rea**2*dlat*0.5*sin(zlat(j))*(zs(i+1,j)-zs(i,j))*( &
        sin(zlon(i+1))*ps(i+1,j)+sin(zlon(i))*ps(i,j))
      tmou(3) = tmou(3) - rea**2*dlat*0.5*cos(zlat(j))*(zs(i+1,j)-zs(i,j))*( &
        ps(i+1,j)+ps(i,j))
    END DO
  END DO


  ! COUPLES DES DIFFERENTES FRICTION AU SOL:

  l = 1
  DO j = 2, nbp_lat-1
    DO i = 1, nbp_lon
      l = l + 1
      tsso(1) = tsso(1) - rea**3*cos(zlat(j))*dlon*dlat*ssou(i, j)*sin(zlat(j &
        ))*cos(zlon(i)) + rea**3*cos(zlat(j))*dlon*dlat*ssov(i, j)*sin(zlon(i &
        ))

      tsso(2) = tsso(2) - rea**3*cos(zlat(j))*dlon*dlat*ssou(i, j)*sin(zlat(j &
        ))*sin(zlon(i)) - rea**3*cos(zlat(j))*dlon*dlat*ssov(i, j)*cos(zlon(i &
        ))

      tsso(3) = tsso(3) + rea**3*cos(zlat(j))*dlon*dlat*ssou(i, j)*cos(zlat(j &
        ))

      tbls(1) = tbls(1) - rea**3*cos(zlat(j))*dlon*dlat*blsu(i, j)*sin(zlat(j &
        ))*cos(zlon(i)) + rea**3*cos(zlat(j))*dlon*dlat*blsv(i, j)*sin(zlon(i &
        ))

      tbls(2) = tbls(2) - rea**3*cos(zlat(j))*dlon*dlat*blsu(i, j)*sin(zlat(j &
        ))*sin(zlon(i)) - rea**3*cos(zlat(j))*dlon*dlat*blsv(i, j)*cos(zlon(i &
        ))

      tbls(3) = tbls(3) + rea**3*cos(zlat(j))*dlon*dlat*blsu(i, j)*cos(zlat(j &
        ))

    END DO
  END DO


  ! write(*,*) 'AAM',rsec,
  ! write(*,*) 'AAM',rjour+rsec/86400.,
  ! c      raam(3)/hadday,oaam(3)/hadday,
  ! c      tmou(3)/hadley,tsso(3)/hadley,tbls(3)/hadley

  ! write(iam,100)rjour+rsec/86400.,
  ! c      raam(1)/hadday,oaam(1)/hadday,
  ! c      tmou(1)/hadley,tsso(1)/hadley,tbls(1)/hadley,
  ! c      raam(2)/hadday,oaam(2)/hadday,
  ! c      tmou(2)/hadley,tsso(2)/hadley,tbls(2)/hadley,
  ! c      raam(3)/hadday,oaam(3)/hadday,
  ! c      tmou(3)/hadley,tsso(3)/hadley,tbls(3)/hadley
100 FORMAT (F12.5, 15(1X,F12.5))

  ! write(iam+1,*)((zs(i,j),i=1,nbp_lon),j=1,nbp_lat)
  ! write(iam+1,*)((ps(i,j),i=1,nbp_lon),j=1,nbp_lat)
  ! write(iam+1,*)((ub(i,j),i=1,nbp_lon),j=1,nbp_lat)
  ! write(iam+1,*)((vb(i,j),i=1,nbp_lon),j=1,nbp_lat)
  ! write(iam+1,*)((ssou(i,j),i=1,nbp_lon),j=1,nbp_lat)
  ! write(iam+1,*)((ssov(i,j),i=1,nbp_lon),j=1,nbp_lat)
  ! write(iam+1,*)((blsu(i,j),i=1,nbp_lon),j=1,nbp_lat)
  ! write(iam+1,*)((blsv(i,j),i=1,nbp_lon),j=1,nbp_lat)

  aam = raam(3)
  torsfc = tmou(3) + tsso(3) + tbls(3)

  RETURN
END SUBROUTINE aaam_bud
