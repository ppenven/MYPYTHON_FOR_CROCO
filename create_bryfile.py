import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset  as netcdf
import time
import os
import sys
sys.path.insert(0,'/home/penven/Teletravail/MYPYTHON_FOR_ROMS')
import croco_vgrid as vgrd
#
#function create_bryfile(bryname,grdname,title,obc,...
#                        theta_s,theta_b,hc,N,...
#                        time,cycle,clobber,vtransform)
################################################################
#
# function create_bryfile(bryname,grdname,title,obc...
#                          theta_s,theta_b,hc,N,...
#                          time,cycle,clobber)
#
#   This function create the header of a Netcdf climatology 
#   file.
#
#   Input:
#
#   bryname      Netcdf climatology file name (character string).
#   grdname      Netcdf grid file name (character string).
#   obc          open boundaries flag (1=open , [S E N W]).
#   theta_s      S-coordinate surface control parameter.(Real)
#   theta_b      S-coordinate bottom control parameter.(Real)
#   hc           Width (m) of surface or bottom boundary layer 
#                where higher vertical resolution is required 
#                during stretching.(Real)
#   N            Number of vertical levels.(Integer)
#   time         time.(vector)
#   cycle        Length (days) for cycling the climatology.(Real)
#   clobber      Switch to allow or not writing over an existing
#                file.(character string)
# 
#  Further Information:  
#  http://www.croco-ocean.org
#  
#  This file is part of CROCOTOOLS
#
#  CROCOTOOLS is free software you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation either version 2 of the License,
#  or (at your option) any later version.
#
#  CROCOTOOLS is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
#  Copyright (c) 2001-2006 by Pierrick Penven 
#  e-mail:Pierrick.Penven@ird.fr  
#
################################################################
#
#
crocofiles_dir='/home/penven/Teletravail/SWAG/CROCO_FILES/'
bryname=crocofiles_dir+'swag4_bry.nc'
grdname=crocofiles_dir+'swag4_grd.nc'
title='Boundary file with GLORYS'
obc = [1,1,1,1] # open boundaries (1=open , [S E N W])
N = 32
theta_s    =  7.
theta_b    =  2.
hc         = 200.
vtransform =  2
brytime=[0.,1.,2.]
cycle=0.
#
#
################################################################
#
#
print(' ')
print(' Creating the file : '+bryname)
print(' ')
print(' VTRANSFORM = '+str(vtransform))
#
#  Read the grid file and check the topography
#
nc=netcdf(grdname, 'r')
h=nc.variables['h'][:]
maskr=nc.variables['mask_rho'][:]
[Lp,Mp]=np.shape(h)
nc.close()
#
#
#
hmin=np.min(h[np.where(maskr==1)])
if vtransform == 1:
  if hc > hmin:
    print('Error: hc ('+hc+' m) > hmin ('+hmin+' m)')


L=Lp-1
M=Mp-1
Np=N+1

Nt=0

#
#  Create the boundary file
#
type = 'BOUNDARY file'
history = 'CROCO' 

if os.path.exists(bryname):
    os.remove(bryname)

print('Create: ' +bryname)
nc=netcdf(bryname,'w',format='NETCDF4')

#
# set global attributes
#

nc.type='CROCO boundary file'
nc.history = 'Created '+str(time.ctime(time.time()))
nc.title = title
nc.clim_file = bryname
nc.grd_file = grdname

#
#  Create dimensions
#
nc_dim_xi_u=nc.createDimension('xi_u',L)
nc_dim_xi_v=nc.createDimension('xi_v',Lp)
nc_dim_xi_rho=nc.createDimension('xi_rho',Lp)
nc_dim_eta_u=nc.createDimension('eta_u',Mp)
nc_dim_eta_v=nc.createDimension('eta_v',M)
nc_dim_eta_rho=nc.createDimension('eta_rho',Mp)
nc_dim_s_rho=nc.createDimension('s_rho',N)
nc_dim_s_w=nc.createDimension('s_w',Np)
nc_dim_tracer=nc.createDimension('tracer',2)
nc_dim_bry_time=nc.createDimension('bry_time',Nt)
nc_dim_tclm_time=nc.createDimension('tclm_time',Nt)
nc_dim_temp_time=nc.createDimension('temp_time',Nt)
nc_dim_sclm_time=nc.createDimension('sclm_time',Nt)
nc_dim_salt_time=nc.createDimension('salt_time',Nt)
nc_dim_uclm_time=nc.createDimension('uclm_time',Nt)
nc_dim_vclm_time=nc.createDimension('vclm_time',Nt)
nc_dim_v2d_time=nc.createDimension('v2d_time',Nt)
nc_dim_v3d_time=nc.createDimension('v3d_time',Nt)
nc_dim_ssh_time=nc.createDimension('ssh_time',Nt)
nc_dim_zeta_time=nc.createDimension('zeta_time',Nt)
nc_dim_one=nc.createDimension('one',1)
#
#  Create variables and attributes
#
nc_spherical=nc.createVariable('spherical','S1', ('one',))
nc_spherical.long_name = 'grid type logical switch'
nc_spherical.flag_values = 'T, F'
nc_spherical.flag_meanings = 'spherical Cartesian'
#
nc_Vtransform=nc.createVariable('Vtransform','i4', ('one',))
nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
#
nc_Vstretching=nc.createVariable('Vstretching','i4', ('one',))
nc_Vstretching.long_name = 'vertical terrain-following stretching function'
#
nc_tstart=nc.createVariable('tstart',np.float64, ('one',))
nc_tstart.long_name = 'start processing day'
nc_tstart.units = 'day'
#
nc_tend=nc.createVariable('tend',np.float64, ('one',))
nc_tend.long_name = 'end processing day'
nc_tend.units = 'day'
#
nc_theta_s=nc.createVariable('theta_s',np.float64, ('one',))
nc_theta_s.long_name = 'S-coordinate surface control parameter'
nc_theta_s.units = 'nondimensional'
#
nc_theta_b=nc.createVariable('theta_b',np.float64, ('one',))
nc_theta_b.long_name = 'S-coordinate bottom control parameter'
nc_theta_b.units = 'nondimensional'
#
nc_Tcline=nc.createVariable('Tcline',np.float64, ('one',))
nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
nc_Tcline.units = 'meter'
#
nc_hc=nc.createVariable('hc',np.float64, ('one',))
nc_hc.long_name = 'S-coordinate parameter, critical depth'
nc_hc.units = 'meter'
#
nc_sc_r=nc.createVariable('sc_r',np.float64, ('s_rho',))
nc_sc_r.long_name = 'S-coordinate at RHO-points'
nc_sc_r.valid_min = -1.
nc_sc_r.valid_max = 0.
nc_sc_r.positive = 'up'

if vtransform == 1:
    nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
elif vtransform == 2:
    nc_sc_r.standard_name = 'ocean_s_coordinate_g2'     

nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
#
nc_sc_w=nc.createVariable('sc_w',np.float64, ('s_w',))
nc_sc_w.long_name = 'S-coordinate at W-points'
nc_sc_w.valid_min = -1. 
nc_sc_w.valid_max = 0. 
nc_sc_w.positive = 'up'

if vtransform == 1:
    nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
elif vtransform == 2:
    nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
#
nc_Cs_r=nc.createVariable('Cs_r',np.float64, ('s_rho',))
nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
nc_Cs_r.units = 'nondimensional'
nc_Cs_r.valid_min = -1
nc_Cs_r.valid_max = 0
#
nc_Cs_w=nc.createVariable('Cs_w',np.float64, ('s_w',))
nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
nc_Cs_w.units = 'nondimensional'
nc_Cs_w.valid_min = -1
nc_Cs_w.valid_max = 0
#
nc_bry_time=nc.createVariable('bry_time',np.float64, ('bry_time',))
nc_bry_time.long_name = 'time for boundary climatology'
nc_bry_time.units = 'day'
nc_bry_time.calendar = 'XXX days in every year'
nc_bry_time.cycle_length = cycle
#
nc_tclm_time=nc.createVariable('tclm_time',np.float64, ('tclm_time',))
nc_tclm_time.long_name = 'time for temperature climatology'
nc_tclm_time.units = 'day'
nc_tclm_time.calendar = 'XXX days in every year'
nc_tclm_time.cycle_length = cycle
#
nc_temp_time=nc.createVariable('temp_time',np.float64, ('temp_time',))
nc_temp_time.long_name = 'time for temperature climatology'
nc_temp_time.units = 'day'
nc_temp_time.calendar = 'XXX days in every year'
nc_temp_time.cycle_length = cycle
#
nc_sclm_time=nc.createVariable('sclm_time',np.float64, ('sclm_time',))
nc_sclm_time.long_name = 'time for salinity climatology'
nc_sclm_time.units = 'day'
nc_sclm_time.calendar = 'XXX days in every year'
nc_sclm_time.cycle_length = cycle
#
nc_salt_time=nc.createVariable('salt_time',np.float64, ('salt_time',))
nc_salt_time.long_name = 'time for salinity climatology'
nc_salt_time.units = 'day'
nc_salt_time.calendar = 'XXX days in every year'
nc_salt_time.cycle_length = cycle
#
nc_uclm_time=nc.createVariable('uclm_time',np.float64, ('uclm_time',))
nc_uclm_time.long_name = 'time climatological u'
nc_uclm_time.units = 'day'
nc_uclm_time.calendar = 'XXX days in every year'
nc_uclm_time.cycle_length = cycle
#
nc_vclm_time=nc.createVariable('vclm_time',np.float64, ('vclm_time',))
nc_vclm_time.long_name = 'time climatological v'
nc_vclm_time.units = 'day'
nc_vclm_time.calendar = 'XXX days in every year'
nc_vclm_time.cycle_length = cycle
#
nc_v2d_time=nc.createVariable('v2d_time',np.float64, ('v2d_time',))
nc_v2d_time.long_name = 'time for 2D velocity climatology'
nc_v2d_time.units = 'day'
nc_v2d_time.calendar = 'XXX days in every year'
nc_v2d_time.cycle_length = cycle
#
nc_v3d_time=nc.createVariable('v3d_time',np.float64, ('v3d_time',))
nc_v3d_time.long_name = 'time for 3D velocity climatology'
nc_v3d_time.units = 'day'
nc_v3d_time.calendar = 'XXX days in every year'
nc_v3d_time.cycle_length = cycle
#
nc_ssh_time=nc.createVariable('ssh_time',np.float64, ('ssh_time',))
nc_ssh_time.long_name = 'time for sea surface height'
nc_ssh_time.units = 'day'
nc_ssh_time.calendar = 'XXX days in every year'
nc_ssh_time.cycle_length = cycle
#
nc_zeta_time=nc.createVariable('zeta_time',np.float64, ('zeta_time',))
nc_zeta_time.long_name = 'time for sea surface height'
nc_zeta_time.units = 'day'
nc_zeta_time.calendar = 'XXX days in every year'
nc_zeta_time.cycle_length = cycle
#
if obc[0]==1:
#
#   Southern boundary
#
  nc_temp_south=nc.createVariable('temp_south',np.float64, ('temp_time','s_rho','xi_rho',))
  nc_temp_south.long_name = 'southern boundary potential temperature'
  nc_temp_south.units = 'Celsius'
  nc_temp_south.coordinates = 'lon_rho s_rho temp_time'
#
  nc_salt_south=nc.createVariable('salt_south',np.float64, ('salt_time','s_rho','xi_rho',))
  nc_salt_south.long_name = 'southern boundary salinity'
  nc_salt_south.units = 'PSU'
  nc_salt_south.coordinates = 'lon_rho s_rho salt_time'
#
  nc_u_south=nc.createVariable('u_south',np.float64, ('v3d_time','s_rho','xi_u',))
  nc_u_south.long_name = 'southern boundary u-momentum component'
  nc_u_south.units = 'meter second-1'
  nc_u_south.coordinates = 'lon_u s_rho u_time'
#
  nc_v_south=nc.createVariable('v_south',np.float64, ('v3d_time','s_rho','xi_rho',))
  nc_v_south.long_name = 'southern boundary v-momentum component'
  nc_v_south.units = 'meter second-1'
  nc_v_south.coordinates = 'lon_v s_rho vclm_time'
#
  nc_ubar_south=nc.createVariable('ubar_south',np.float64, ('v2d_time','xi_u',))
  nc_ubar_south.long_name = 'southern boundary vertically integrated u-momentum component'
  nc_ubar_south.units = 'meter second-1'
  nc_ubar_south.coordinates = 'lon_u uclm_time'
#
  nc_vbar_south=nc.createVariable('vbar_south',np.float64, ('v2d_time','xi_rho',))
  nc_vbar_south.long_name = 'southern boundary vertically integrated v-momentum component'
  nc_vbar_south.units = 'meter second-1'
  nc_vbar_south.coordinates = 'lon_v vclm_time'
#
  nc_zeta_south=nc.createVariable('zeta_south',np.float64, ('zeta_time','xi_rho',))
  nc_zeta_south.long_name = 'southern boundary sea surface height'
  nc_zeta_south.units = 'meter'
  nc_zeta_south.coordinates = 'lon_rho zeta_time'
#

if obc[1]==1:
#
#   Eastern boundary
#
  nc_temp_east=nc.createVariable('temp_east',np.float64, ('temp_time','s_rho','eta_rho',))
  nc_temp_east.long_name = 'eastern boundary potential temperature'
  nc_temp_east.units = 'Celsius'
  nc_temp_east.coordinates = 'lat_rho s_rho temp_time'
#
  nc_salt_east=nc.createVariable('salt_east',np.float64, ('salt_time','s_rho','eta_rho',))
  nc_salt_east.long_name = 'eastern boundary salinity'
  nc_salt_east.units = 'PSU'
  nc_salt_east.coordinates = 'lat_rho s_rho salt_time'
#
  nc_u_east=nc.createVariable('u_east',np.float64, ('v3d_time','s_rho','eta_rho',))
  nc_u_east.long_name = 'eastern boundary u-momentum component'
  nc_u_east.units = 'meter second-1'
  nc_u_east.coordinates = 'lat_u s_rho u_time'
#
  nc_v_east=nc.createVariable('v_east',np.float64, ('v3d_time','s_rho','eta_v',))
  nc_v_east.long_name = 'eastern boundary v-momentum component'
  nc_v_east.units = 'meter second-1'
  nc_v_east.coordinates = 'lat_v s_rho vclm_time'
#
  nc_ubar_east=nc.createVariable('ubar_east',np.float64, ('v2d_time','eta_rho',))
  nc_ubar_east.long_name = 'eastern boundary vertically integrated u-momentum component'
  nc_ubar_east.units = 'meter second-1'
  nc_ubar_east.coordinates = 'lat_u uclm_time'
#
  nc_vbar_east=nc.createVariable('vbar_east',np.float64, ('v2d_time','eta_v',))
  nc_vbar_east.long_name = 'eastern boundary vertically integrated v-momentum component'
  nc_vbar_east.units = 'meter second-1'
  nc_vbar_east.coordinates = 'lat_v vclm_time'
#
  nc_zeta_east=nc.createVariable('zeta_east',np.float64, ('zeta_time','eta_rho',))
  nc_zeta_east.long_name = 'eastern boundary sea surface height'
  nc_zeta_east.units = 'meter'
  nc_zeta_east.coordinates = 'lat_rho zeta_time'
#

if obc[2]==1:
#
#   Northern boundary
#
  nc_temp_north=nc.createVariable('temp_north',np.float64, ('temp_time','s_rho','xi_rho',))
  nc_temp_north.long_name = 'northern boundary potential temperature'
  nc_temp_north.units = 'Celsius'
  nc_temp_north.coordinates = 'lon_rho s_rho temp_time'
#
  nc_salt_north=nc.createVariable('salt_north',np.float64, ('salt_time','s_rho','xi_rho',))
  nc_salt_north.long_name = 'northern boundary salinity'
  nc_salt_north.units = 'PSU'
  nc_salt_north.coordinates = 'lon_rho s_rho salt_time'
#
  nc_u_north=nc.createVariable('u_north',np.float64, ('v3d_time','s_rho','xi_u',))
  nc_u_north.long_name = 'northern boundary u-momentum component'
  nc_u_north.units = 'meter second-1'
  nc_u_north.coordinates = 'lon_u s_rho u_time'
#
  nc_v_north=nc.createVariable('v_nort',np.float64, ('v3d_time','s_rho','xi_rho',))
  nc_v_north.long_name = 'northern boundary v-momentum component'
  nc_v_north.units = 'meter second-1'
  nc_v_north.coordinates = 'lon_v s_rho vclm_time'
#
  nc_ubar_north=nc.createVariable('ubar_north',np.float64, ('v2d_time','xi_u',))
  nc_ubar_north.long_name = 'northern boundary vertically integrated u-momentum component'
  nc_ubar_north.units = 'meter second-1'
  nc_ubar_north.coordinates = 'lon_u uclm_time'
#
  nc_vbar_north=nc.createVariable('vbar_north',np.float64, ('v2d_time','xi_rho',))
  nc_vbar_north.long_name = 'northern boundary vertically integrated v-momentum component'
  nc_vbar_north.units = 'meter second-1'
  nc_vbar_north.coordinates = 'lon_v vclm_time'

  nc_zeta_north=nc.createVariable('zeta_north',np.float64, ('zeta_time','xi_rho',))
  nc_zeta_north.long_name = 'northern boundary sea surface height'
  nc_zeta_north.units = 'meter'
  nc_zeta_north.coordinates = 'lon_rho zeta_time'
#

if obc[3]==1:
#
#   Western boundary
#
  nc_temp_west=nc.createVariable('temp_west',np.float64, ('temp_time','s_rho','eta_rho',))
  nc_temp_west.long_name = 'western boundary potential temperature'
  nc_temp_west.units = 'Celsius'
  nc_temp_west.coordinates = 'lat_rho s_rho temp_time'
#
  nc_salt_west=nc.createVariable('salt_west',np.float64, ('salt_time','s_rho','eta_rho',))
  nc_salt_west.long_name = 'western boundary salinity'
  nc_salt_west.units = 'PSU'
  nc_salt_west.coordinates = 'lat_rho s_rho salt_time'
#
  nc_u_west=nc.createVariable('u_west',np.float64, ('v3d_time','s_rho','eta_rho',))
  nc_u_west.long_name = 'western boundary u-momentum component'
  nc_u_west.units = 'meter second-1'
  nc_u_west.coordinates = 'lat_u s_rho u_time'
#
  nc_v_west=nc.createVariable('v_west',np.float64, ('v3d_time','s_rho','eta_v',))
  nc_v_west.long_name = 'western boundary v-momentum component'
  nc_v_west.units = 'meter second-1'
  nc_v_west.coordinates = 'lat_v s_rho vclm_time'
#
  nc_ubar_west=nc.createVariable('ubar_west',np.float64, ('v2d_time','eta_rho',))
  nc_ubar_west.long_name = 'western boundary vertically integrated u-momentum component'
  nc_ubar_west.units = 'meter second-1'
  nc_ubar_west.coordinates = 'lat_u uclm_time'
#
  nc_vbar_west=nc.createVariable('vbar_west',np.float64, ('v2d_time','eta_v',))
  nc_vbar_west.long_name = 'western boundary vertically integrated v-momentum component'
  nc_vbar_west.units = 'meter second-1'
  nc_vbar_west.coordinates = 'lat_v vclm_time'
#
  nc_zeta_west=nc.createVariable('zeta_west',np.float64, ('zeta_time','eta_rho',))
  nc_zeta_west.long_name = 'western boundary sea surface height'
  nc_zeta_west.units = 'meter'
  nc_zeta_west.coordinates = 'lat_rho zeta_time'
#

#
# Compute S coordinates
#
(sc_r,Cs_r,sc_w,Cs_w) = vgrd.scoordinate(theta_s,theta_b,N,hc,vtransform)
print('vtransform = '+str(vtransform))
#
# Write variables
#
nc_spherical[:]='T'
nc_Vtransform[:]=vtransform
nc_Vstretching[:]=1
nc_tstart[:] =   np.min(brytime) 
nc_tend[:] =     np.max(brytime) 
nc_theta_s[:] =  theta_s 
nc_theta_b[:] =  theta_b 
nc_Tcline[:] =  hc 
nc_hc[:] =  hc 
nc_sc_r[:] = sc_r
nc_sc_w[:] = sc_w
nc_Cs_r[:] = Cs_r  
nc_Cs_w[:] = Cs_w
nc_tclm_time[:] =  brytime 
nc_temp_time[:] =  brytime 
nc_sclm_time[:] =  brytime 
nc_salt_time[:] =  brytime 
nc_uclm_time[:] =  brytime 
nc_vclm_time[:] =  brytime 
nc_v2d_time[:]  =  brytime 
nc_v3d_time[:]  =  brytime 
nc_ssh_time[:]  =  brytime
nc_zeta_time[:] =  brytime
nc_bry_time[:]  =  brytime 

if obc[0]==1:
  nc_u_south[:] =  0. 
  nc_v_south[:] =  0. 
  nc_ubar_south[:] =  0. 
  nc_vbar_south[:] =  0. 
  nc_zeta_south[:] =  0. 
  nc_temp_south[:] =  0. 
  nc_salt_south[:] =  0.

if obc[1]==1:
  nc_u_east[:] =  0. 
  nc_v_east[:] =  0. 
  nc_ubar_east[:] =  0. 
  nc_vbar_east[:] =  0. 
  nc_zeta_east[:] =  0. 
  nc_temp_east[:] =  0. 
  nc_salt_east[:] =  0.

if obc[2]==1:
  nc_u_north[:] =  0. 
  nc_v_north[:] =  0. 
  nc_ubar_north[:] =  0. 
  nc_vbar_north[:] =  0. 
  nc_zeta_north[:] =  0. 
  nc_temp_north[:] =  0. 
  nc_salt_north[:] =  0.

if obc[3]==1:
  nc_u_west[:] =  0. 
  nc_v_west[:] =  0. 
  nc_ubar_west[:] =  0. 
  nc_vbar_west[:] =  0. 
  nc_zeta_west[:] =  0. 
  nc_temp_west[:] =  0. 
  nc_salt_west[:] =  0.
  
nc.close()



