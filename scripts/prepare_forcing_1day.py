#!/usr/bin/env python3
"""
prepare_forcing_1day.py
=======================
Extract and format atmospheric forcing for a 1-day ORCHIDEE Offline test
(1979-02-01) from historical coupled LMDZ-ORCHIDEE simulation outputs.

Source dataset:
  /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist/ATM/Output/HF/sim35km-HIST_19790201_19790228_HF_histhf.nc

Target output:
  simulations/offline_1day/forcing_19790201.nc
"""

import sys
import os
import argparse
import numpy as np
import netCDF4 as nc

def prepare_forcing(input_path, output_path, start_idx=0, num_steps=8, force=False):
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"Input atmospheric file not found: {input_path}")

    if os.path.isfile(output_path) and not force:
        file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
        print(f"Le fichier de forçage existe déjà : {output_path} ({file_size_mb:.2f} MB)")
        print("Protection active : aucun remplacement par défaut. Utilisez l'option --force pour écraser et régénérer.")
        return

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    print(f"Opening source file: {input_path}")
    with nc.Dataset(input_path, 'r') as src:
        # Check required variables
        req_vars = ['lon', 'lat', 'time_counter', 't2m', 'q2m', 'psol', 
                    'wind10m', 'u10m', 'v10m', 'precip', 'snow', 
                    'SWdnSFC', 'LWdnSFC', 'contfracOR']
        missing = [v for v in req_vars if v not in src.variables]
        if missing:
            raise KeyError(f"Missing required variables in source file: {missing}")

        n_times = len(src.dimensions['time_counter'])
        end_idx = start_idx + num_steps
        if end_idx > n_times:
            raise IndexError(f"Requested {num_steps} steps from index {start_idx}, but file has {n_times} steps.")

        print(f"Extracting steps {start_idx} to {end_idx - 1} ({num_steps} steps of 3h = 24h)...")

        # Read coordinates
        lons = src.variables['lon'][:]
        lats = src.variables['lat'][:]
        time_var = src.variables['time_counter']
        times = time_var[start_idx:end_idx]
        time_units = time_var.units
        time_calendar = getattr(time_var, 'calendar', 'gregorian')

        n_lon = len(lons)
        n_lat = len(lats)
        print(f"Grid dimensions: lon={n_lon}, lat={n_lat}, time={num_steps}")

        print(f"Creating output forcing file: {output_path}")
        with nc.Dataset(output_path, 'w', format='NETCDF4') as dst:
            # Set global attributes
            dst.title = "ORCHIDEE Offline 1-Day Atmospheric Forcing (1979-02-01)"
            dst.source = f"Extracted from {os.path.basename(input_path)}"
            dst.institution = "UM6P / IWRI - Climate Modeling Team"
            dst.conventions = "CF-1.6"
            dst.history = "Prepared for ORCHIDEE standalone test"

            # Create dimensions
            dst.createDimension('lon', n_lon)
            dst.createDimension('lat', n_lat)
            dst.createDimension('time_counter', num_steps)

            # Coordinate variables
            v_lon = dst.createVariable('lon', 'f4', ('lon',))
            v_lon[:] = lons
            v_lon.units = "degrees_east"
            v_lon.long_name = "Longitude"
            v_lon.standard_name = "longitude"
            v_lon.axis = "X"

            v_lat = dst.createVariable('lat', 'f4', ('lat',))
            v_lat[:] = lats
            v_lat.units = "degrees_north"
            v_lat.long_name = "Latitude"
            v_lat.standard_name = "latitude"
            v_lat.axis = "Y"

            # 2D coordinates nav_lon, nav_lat for full IOIPSL/NEMO/ORCHIDEE compatibility
            v_nav_lon = dst.createVariable('nav_lon', 'f4', ('lat', 'lon'), zlib=True, complevel=4)
            v_nav_lat = dst.createVariable('nav_lat', 'f4', ('lat', 'lon'), zlib=True, complevel=4)
            lon_2d, lat_2d = np.meshgrid(lons, lats)
            v_nav_lon[:] = lon_2d
            v_nav_lon.units = "degrees_east"
            v_nav_lon.long_name = "Longitude"

            v_nav_lat[:] = lat_2d
            v_nav_lat.units = "degrees_north"
            v_nav_lat.long_name = "Latitude"

            v_time = dst.createVariable('time_counter', 'f8', ('time_counter',))
            v_time[:] = times
            v_time.units = time_units
            v_time.calendar = time_calendar
            v_time.long_name = "Time axis"
            v_time.standard_name = "time"
            v_time.axis = "T"

            # Extract meteorological variables
            print("Processing Tair (2m air temperature)...")
            v_tair = dst.createVariable('Tair', 'f4', ('time_counter', 'lat', 'lon'), zlib=True, complevel=4)
            v_tair[:] = src.variables['t2m'][start_idx:end_idx, :, :]
            v_tair.units = "K"
            v_tair.long_name = "Near-surface air temperature"
            v_tair.standard_name = "air_temperature"

            print("Processing Qair (2m specific humidity)...")
            v_qair = dst.createVariable('Qair', 'f4', ('time_counter', 'lat', 'lon'), zlib=True, complevel=4)
            v_qair[:] = src.variables['q2m'][start_idx:end_idx, :, :]
            v_qair.units = "kg/kg"
            v_qair.long_name = "Near-surface specific humidity"
            v_qair.standard_name = "specific_humidity"

            print("Processing PSurf (surface pressure)...")
            v_ps = dst.createVariable('PSurf', 'f4', ('time_counter', 'lat', 'lon'), zlib=True, complevel=4)
            v_ps[:] = src.variables['psol'][start_idx:end_idx, :, :]
            v_ps.units = "Pa"
            v_ps.long_name = "Surface atmospheric pressure"
            v_ps.standard_name = "surface_air_pressure"

            print("Processing Wind, Wind_E, Wind_N (10m wind speed & components)...")
            v_wind = dst.createVariable('Wind', 'f4', ('time_counter', 'lat', 'lon'), zlib=True, complevel=4)
            v_wind[:] = src.variables['wind10m'][start_idx:end_idx, :, :]
            v_wind.units = "m/s"
            v_wind.long_name = "10m wind speed"
            v_wind.standard_name = "wind_speed"

            v_u = dst.createVariable('Wind_E', 'f4', ('time_counter', 'lat', 'lon'), zlib=True, complevel=4)
            v_u[:] = src.variables['u10m'][start_idx:end_idx, :, :]
            v_u.units = "m/s"
            v_u.long_name = "Eastward 10m wind component"
            v_u.standard_name = "eastward_wind"

            v_v = dst.createVariable('Wind_N', 'f4', ('time_counter', 'lat', 'lon'), zlib=True, complevel=4)
            v_v[:] = src.variables['v10m'][start_idx:end_idx, :, :]
            v_v.units = "m/s"
            v_v.long_name = "Northward 10m wind component"
            v_v.standard_name = "northward_wind"

            print("Processing Rainf & Snowf (liquid and solid precipitation fluxes)...")
            raw_precip = src.variables['precip'][start_idx:end_idx, :, :]
            raw_snow = src.variables['snow'][start_idx:end_idx, :, :]
            raw_rain = np.maximum(0.0, raw_precip - raw_snow)

            v_rain = dst.createVariable('Rainf', 'f4', ('time_counter', 'lat', 'lon'), zlib=True, complevel=4)
            v_rain[:] = raw_rain
            v_rain.units = "kg/m^2/s"
            v_rain.long_name = "Rainfall flux (liquid precipitation rate)"
            v_rain.standard_name = "rainfall_flux"

            v_snow = dst.createVariable('Snowf', 'f4', ('time_counter', 'lat', 'lon'), zlib=True, complevel=4)
            v_snow[:] = raw_snow
            v_snow.units = "kg/m^2/s"
            v_snow.long_name = "Snowfall flux (solid precipitation rate)"
            v_snow.standard_name = "snowfall_flux"

            print("Processing SWdown & LWdown (surface radiative fluxes)...")
            v_sw = dst.createVariable('SWdown', 'f4', ('time_counter', 'lat', 'lon'), zlib=True, complevel=4)
            v_sw[:] = src.variables['SWdnSFC'][start_idx:end_idx, :, :]
            v_sw.units = "W/m^2"
            v_sw.long_name = "Surface downwelling solar radiation"
            v_sw.standard_name = "surface_downwelling_shortwave_flux_in_air"

            v_lw = dst.createVariable('LWdown', 'f4', ('time_counter', 'lat', 'lon'), zlib=True, complevel=4)
            v_lw[:] = src.variables['LWdnSFC'][start_idx:end_idx, :, :]
            v_lw.units = "W/m^2"
            v_lw.long_name = "Surface downwelling longwave radiation"
            v_lw.standard_name = "surface_downwelling_longwave_flux_in_air"

            print("Processing contfrac (continental land fraction)...")
            cf_data = src.variables['contfracOR'][start_idx:end_idx, :, :]
            v_cf = dst.createVariable('contfrac', 'f4', ('time_counter', 'lat', 'lon'), zlib=True, complevel=4)
            v_cf[:] = cf_data
            v_cf.units = "1"
            v_cf.long_name = "Continental land fraction"

    file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"Forcing file successfully prepared: {output_path} ({file_size_mb:.2f} MB)")

if __name__ == '__main__':
    default_input = "/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist/ATM/Output/HF/sim35km-HIST_19790201_19790228_HF_histhf.nc"
    default_output = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "simulations", "offline_1day", "forcing_19790201.nc"))

    parser = argparse.ArgumentParser(description="Prepare 1-day ORCHIDEE Offline forcing from LMDZ historical outputs.")
    parser.add_argument("--input", "-i", default=default_input, help="Path to histhf.nc input file")
    parser.add_argument("--output", "-o", default=default_output, help="Path to output forcing NetCDF file")
    parser.add_argument("--start", "-s", type=int, default=0, help="Starting time step index (default: 0 for 1979-02-01)")
    parser.add_argument("--steps", "-n", type=int, default=8, help="Number of 3-hour steps to extract (default: 8 = 24h)")
    parser.add_argument("--force", "-f", action="store_true", help="Force overwrite if target file already exists")

    args = parser.parse_args()
    prepare_forcing(args.input, args.output, args.start, args.steps, args.force)
