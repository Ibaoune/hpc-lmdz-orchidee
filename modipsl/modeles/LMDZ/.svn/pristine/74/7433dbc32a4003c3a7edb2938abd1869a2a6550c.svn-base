! $Id$
module regr_pr_o3_m

  implicit none

contains

  subroutine regr_pr_o3(p3d, o3_mob_regr)

    ! "regr_pr_o3" stands for "regrid pressure ozone".
    ! This procedure reads Mobidic ozone mole fraction from
    ! "coefoz_LMDZ.nc" at the initial day of the run and regrids it in
    ! pressure.
    ! Ozone mole fraction from "coefoz_LMDZ.nc" at the initial day is
    ! a 2D latitude -- pressure variable.
    ! The target horizontal LMDZ grid is the "scalar" grid: "rlonv", "rlatu".
    ! The target vertical LMDZ grid is the grid of layer boundaries.
    ! We assume that the input variable is already on the LMDZ "rlatu"
    ! latitude grid.
    ! The input variable does not depend on longitude, but the
    ! pressure at LMDZ layers does.
    ! Therefore, the values on the LMDZ grid do depend on longitude.
    ! Regridding is by averaging, assuming a step function.
    ! We assume that, in the input file, the pressure levels are in
    ! hPa and strictly increasing.

    use netcdf95, only: nf95_open, nf95_close, nf95_inq_varid, handle_err
    use netcdf, only:  nf90_nowrite, nf90_get_var
    use assert_m, only: assert
    use regr_conserv_m, only: regr_conserv
    use press_coefoz_m, only: press_in_edg
    use time_phylmdz_mod, only: day_ref
    use mod_grid_phy_lmdz, only: nbp_lon, nbp_lat, nbp_lev

    REAL, intent(in):: p3d(:, :, :) ! pressure at layer interfaces, in Pa
    ! ("p3d(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)",
    ! for interface "l")

    real, intent(out):: o3_mob_regr(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! (ozone mole fraction from Mobidic adapted to the LMDZ grid)
    ! ("o3_mob_regr(i, j, l)" is at longitude "rlonv(i)", latitude
    ! "rlatu(j)" and pressure level "pls(i, j, l)")

    ! Variables local to the procedure:

    integer ncid, varid, ncerr ! for NetCDF
    integer i, j

    real r_mob(nbp_lat, size(press_in_edg) - 1)
    ! (ozone mole fraction from Mobidic at day "day_ref")
    ! (r_mob(j, k) is at latitude "rlatu(j)", in pressure interval
    ! "[press_in_edg(k), press_in_edg(k+1)]".)

    !------------------------------------------------------------

    print *, "Call sequence information: regr_pr_o3"
    call assert(shape(o3_mob_regr) == (/nbp_lon + 1, nbp_lat, nbp_lev/), &
         "regr_pr_o3 o3_mob_regr")
    call assert(shape(p3d) == (/nbp_lon + 1, nbp_lat, nbp_lev + 1/), "regr_pr_o3 p3d")

    call nf95_open("coefoz_LMDZ.nc", nf90_nowrite, ncid)

    call nf95_inq_varid(ncid, "r_Mob", varid)
    ! Get data at the right day from the input file:
    ncerr = nf90_get_var(ncid, varid, r_mob, start=(/1, 1, day_ref/))
    call handle_err("nf90_get_var r_Mob", ncerr)
    ! Latitudes are in ascending order in the input file while
    ! "rlatu" is in descending order so we need to invert order:
    r_mob = r_mob(nbp_lat:1:-1, :)

    call nf95_close(ncid)

    ! Regrid in pressure by averaging a step function of pressure:

    ! Poles:
    do j = 1, nbp_lat, nbp_lat-1
       call regr_conserv(1, r_mob(j, :), press_in_edg, &
            p3d(1, j, nbp_lev + 1:1:-1), o3_mob_regr(1, j, nbp_lev:1:-1))
       ! (invert order of indices because "p3d" is in descending order)
    end do

    ! Other latitudes:
    do j = 2, nbp_lat-1
       do i = 1, nbp_lon
          call regr_conserv(1, r_mob(j, :), press_in_edg, &
               p3d(i, j, nbp_lev + 1:1:-1), o3_mob_regr(i, j, nbp_lev:1:-1))
          ! (invert order of indices because "p3d" is in descending order)
       end do
    end do

    ! Duplicate pole values on all longitudes:
    o3_mob_regr(2:, 1, :) = spread(o3_mob_regr(1, 1, :), dim=1, ncopies=nbp_lon)
    o3_mob_regr(2:, nbp_lat, :) &
         = spread(o3_mob_regr(1, nbp_lat, :), dim=1, ncopies=nbp_lon)

    ! Duplicate first longitude to last longitude:
    o3_mob_regr(nbp_lon + 1, 2:nbp_lat-1, :) = o3_mob_regr(1, 2:nbp_lat-1, :)

  end subroutine regr_pr_o3

end module regr_pr_o3_m
