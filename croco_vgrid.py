#
########################################################################
#
#
# A few python functions to deal with CROCO vertical grids :
#
# * CSF = get_csf(sc,theta_s,theta_b)
# Get CS-Curves for the new s-ccordinate system
#
# * (sc_r,Cs_r,sc_w,Cs_w) = scoordinate(theta_s,theta_b,N,hc,vtransform)
# Define S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
#
# * z=zlevs(h,zeta,theta_s,theta_b,hc,N,type,vtransform)
# Compute the depth of rho or w points for CROCO
#
# * vnew=vinterp(var,z,depth)
# Interpolate a 3D variable on a horizontal level of constant depth
#
# * (V,h0)=vintegr(var,zw,zr,z01,z02)
# Vertically integrate a CROCO variable (var) from a constant depth 
# z01 (ex z01=-4000 m) to a constant depth z02 (ex z02=-2000m).
# If z01 = NaN : perform the integration from the bottom to z02.
# If z02 = NaN : perform the integration from z01 to the surface.
# If they are both NaNs perform the integration from the bottom
# to the surface.
#
# * (LON,LAT,X,Z,VAR)=get_section(lonsec,latsec,var,lon,lat,
#                       rmask,h,zeta,theta_s,theta_b,hc,N,vtransform)
#  Extract a vertical slice in any direction (or along a curve).
#
# * (LON,LAT,X,Z,Un,Ut)=get_UV_section(lonsec,latsec,u,v, lon,lat,
#                            rmask,h,zeta,theta_s,theta_b,hc,N,vtransform)
#  Extract a vertical slice of cross-section and along section velocities
#  in any direction (or along a curve)
#
#
########################################################################
#
import sys
import numpy                as np
from scipy.interpolate import griddata
import tpx_tools as tpx
import croco_vgrid as vgrd
#
### FUNCTION GET_CSF ###################################################
#
# function CSF = get_csf(sc,theta_s,theta_b)
#
# Get CS-Curves for the new s-ccordinate system
#
                                          # NOTE: Mathematical 
def get_csf(sc,theta_s,theta_b):            # limits of CSF,csrf for
                                            # theta_s, theta_b --> 0
                                            # match that under "else"
                                            # logical branches.
  if (theta_s>0.0):
      csrf=(1.0 - np.cosh(theta_s * sc)) / (np.cosh(theta_s) - 1.0)
  else:
      csrf=-sc**2


  if (theta_b>0.0):
    CSF=(np.exp(theta_b * csrf) -1.0) / (1.0 - np.exp(-theta_b))  
  else:
    CSF=csrf

  return CSF

### END FUNCTION GET_CSF ###################################################


### FUNCTION SCOORDINATE ###################################################
#
# function [sc_r,Cs_r,sc_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform)
#
# Define S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
#
########################################################################
def scoordinate(theta_s,theta_b,N,hc,vtransform):
#
# Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
#
#  sc_r = np.zeros(N)
#  Cs_r = np.zeros(N)
#  sc_w = np.zeros(N+1)
#  Cs_w = np.zeros(N+1)

#
  if vtransform == 2:

    print('NEW_S_COORD')

    ds=1./N
 
    sc_r= ds*(np.arange(1,N+1)-N-0.5)
    Cs_r=get_csf(sc_r,theta_s,theta_b)
#
#    sc_w[0]   = -1.
#    sc_w[N]   =  0.
#    Cs_w[0]   = -1.
#    Cs_w[N]   =  0.

    sc_w = ds*(np.arange(0,N+1)-N)
    Cs_w=get_csf(sc_w, theta_s,theta_b);

  else:

    print('OLD_S_COORD')

    cff1 = 1.0 / np.sinh(theta_s)
    cff2 = 0.5 / np.tanh(0.5*theta_s)

    sc_w = (np.arange(0,N+1) - N)/N
    Cs_w = (1. - theta_b) * cff1 * np.sinh(theta_s * sc_w) + \
         theta_b * (cff2 * np.tanh(theta_s * (sc_w + 0.5)) - 0.5)
  
    sc_r= (np.arange(1,N+1)-N-0.5)/N
    Cs_r = (1. - theta_b) * cff1 * np.sinh(theta_s * sc_r) + \
           theta_b * (cff2 * np.tanh(theta_s * (sc_r + 0.5)) - 0.5)

  return sc_r,Cs_r,sc_w,Cs_w
### END FUNCTION SCOORDINATE ###################################################



### FUNCTION ZLEVS #####################################################
#
# function z=zlevs(h,zeta,theta_s,theta_b,hc,N,type,vtransform)
#
#  this function compute the depth of rho or w points for CROCO
#
#  On Input:
#
#    type    'r': rho point 'w': w point 
#    vtransform  1=> old v transform (Song, 1994); 
#                2=> new v transform (Shcheptekin, 2006)
#  On Output:
#
#    z       Depths (m) of RHO- or W-points (3D matrix).
#
########################################################################
def zlevs(h,zeta,theta_s,theta_b,hc,N,type,vtransform):

#
# Test the number of dimension for h
#
  
  Ndim=np.size(np.shape(h))
  if Ndim==2:  
    M, L = np.squeeze(np.shape(h))
  elif Ndim==1:
    L = np.squeeze(np.shape(h))
  else:
    print('zlevs: error - incorrect dimension for h')
    return
    
  hmin = h.min()
  hmax = h.max()

#
# Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
#

  ds = 1.0 / float(N)      

  if type=='w':
  
    sc = ds * (np.arange(0,N+1) - N)
    Nmax=N+1
    
  elif type=='r':
  
    sc = ds * (np.arange(1,N+1) - N - 0.5)
    Nmax=N
    
  else:
  
    print('Problem with type = ',type)
    sys.exit()

  if vtransform==1:  

    print('OLD_S_COORD')

    cff1 = 1.0 / np.sinh(theta_s)
    cff2 = 0.5 / np.tanh(0.5*theta_s)
    Cs = (1. - theta_b) * cff1 * np.sinh(theta_s * sc) + \
           theta_b * (cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5)


  elif vtransform==2:  

    print('NEW_S_COORD')

    Cs = get_csf(sc,theta_s,theta_b) 
         
  else:

    print('Problem with vtransform = ',vtransform)
    sys.exit()

#
#  Set vertical grid
#

  if Ndim==2:  
    z = np.zeros([Nmax, M, L])
  elif Ndim==1:
    z = np.zeros([Nmax, L])

  for k in np.arange(0,Nmax):

    if vtransform==1: 

      cff = hc*(sc[k]-Cs[k])
      cff1 = Cs[k]
      z0 = cff + cff1*h
      hinv=1./h  
      if Ndim==2:  
        z[k,:,:]=z0+zeta*(1.+z0*hinv)
      elif Ndim==1:
        z[k,:]=z0+zeta*(1.+z0*hinv)

    elif vtransform==2: 

      cff = hc*sc[k]
      cff1 = Cs[k]
      z0 = cff + cff1*h 
      hinv = 1./(h+hc)
      if Ndim==2:  
       z[k,:,:]=z0*h*hinv+zeta*(1.+z0*hinv)
      elif Ndim==1:
       z[k,:]=z0*h*hinv+zeta*(1.+z0*hinv)

  return z

### END FUNCTION ZLEVS #################################################



### FUNCTION VINTERP ###################################################
#
# function  vnew = vinterp(var,z,depth)
#
# This function interpolate a 3D variable on a horizontal level of constant
# depth
#
# On Input:  
# 
#    var     Variable to process (3D matrix).
#    z       Depths (m) of RHO- or W-points (3D matrix).
#    depth   Slice depth (scalar; meters, negative).
# 
# On Output: 
#
#    vnew    Horizontal slice (2D matrix). 
#
#  Further Information:  
#  http://www.croco-ocean.org
#
######################################################################
def  vinterp(var,z,depth):

  [N,M,L]=np.shape(z)
  i1=np.arange(0,L)
  j1=np.arange(0,M)

#
# Find the grid positions of the nearest vertical levels
#

  a=np.int_(z<depth)
  levs=np.sum(a,axis=0)
  levs[np.where(levs==N)] = N-1

  mask=0.*levs + 1.
  mask[np.where(levs==0)]=np.nan

  [i2,j2]=np.meshgrid(i1,j1)
  
  pos_up  = L*M*levs + L*j2 + i2
  pos_down= L*M*(levs-1) + L*j2 + i2
  
  za=np.reshape(z,N*M*L)
  z_up=za[pos_up]
  z_down=za[pos_down]

  va=np.reshape(var,N*M*L)
  v_up=va[pos_up]
  v_down=va[pos_down]

#
# Do the vertical interpolation
#

  vnew=mask * ( (v_up-v_down)*depth + v_down*z_up - v_up*z_down ) / (z_up-z_down)
  
  return vnew

### END FUNCTION VINTERP ###############################################



### FUNCTION VINTEGR ###################################################

def vintegr(var,zw,zr,z01,z02):

#function [V,h0]=vintegr2(var,zw,zr,z01,z02)
#
# Vertically integrate a CROCO variable (var) from a constant depth 
# z01 (ex z01=-4000 m) to a constant depth z02 (ex z02=-2000m).
#
# If z01 = NaN : perform the integration from the bottom to z02.
# If z02 = NaN : perform the integration from z01 to the surface.
# If they are both NaNs perform the integration from the bottom
# to the surface.
#
# Input :
#
# var : CROCO variable at RHO-points (3D matrix)
# zw  : Depth of the W-points (3D matrix)
# zr  : Depth of the RHO-points (3D matrix)
# z01 : lower limit of integration (scalar)
# z02 : upper limit of integration (scalar)
#
# Output :
#
# V   : intgrated value (2D matrix)
# h0  : layer thickness (2D matrix)
#
# Pierrick Penven 2005
#
  if z02 <= z01:
    print('vintegr2:  z02 <= z01')
    sys.exit()

  [Np,M,L]=np.shape(zw)
  N=Np-1
  mask=np.zeros([M,L]) + 1.
  
  i1=np.arange(0,L)
  j1=np.arange(0,M)
  [i2,j2]=np.meshgrid(i1,j1)

  za=np.reshape(zw,Np*M*L)
  vara=np.reshape(var,N*M*L)

#
# Get the matrix of the variables above of z01 and below z02
# 

  if (np.isfinite(z01) & np.isfinite(z02)):
    isgood=np.int_(zw[1:,:,:]>z01) * np.int_(zw[1:,:,:]<z02)

  elif np.isfinite(z01):
    isgood=np.int_(zw[1:,:,:]>z01)

  elif np.isfinite(z02):
    isgood=np.int_(zw[1:,:,:]<z02)

  else:
    isgood=np.int_(var==var)

#
  if np.isfinite(z01):
#
# Find the bottom limit of the corresponding grid cells
#
    a=np.int_(zw<z01)
    levs=np.sum(a,axis=0)-1
    mask=np.zeros([M,L]) + 1.
    mask[np.where(levs<0)]=np.nan 
    mask[np.where( (levs<0) | (levs==N) )]=np.nan 
    levs[np.where(levs==N)]=1
    
    pos = L*M*levs + L*j2 + i2
    z1=mask*za[pos]

#
# Compute the vertical size of the partial step to add at the bottom
#

    dzbot=z1-z01
    dzbot[np.where(np.isnan(mask))]=0.
#
# Get the value of the variable in the partial step to add at the bottom
#

    Vbot=vara[pos]

  else:

    dzbot=0 
    Vbot=0

  if np.isfinite(z02):
#
# Find the top positions
#
    a=np.int_(zw<z02)
    levs=np.sum(a,axis=0)-1
    mask=np.zeros([M,L]) + 1.
    mask[np.where( (levs<0) | (levs==N) )]=np.nan 
    levs[np.where(levs==N)]=1
    
    pos = L*M*levs + L*j2 + i2
    z2=mask*za[pos]

#
# Compute the vertical size of the partial step to add at the top
#

    dztop=z02-z2
    dztop[np.where(np.isnan(mask))]=0.

#
# Get the value of the variable in the partial step to add at the bottom
#

    Vtop=vara[pos]

  else:

    dztop=0
    Vtop=0

#
# Perform the vertical integration


  dz=zw[1:,:,:]-zw[:-1,:,:]
  V=np.sum(dz*isgood*var,axis=0) + dzbot*Vbot + dztop*Vtop

#
# Get the depth
#

  h0=np.sum(dz*isgood,axis=0) + dzbot + dztop
 
  V[np.where(h0==0)]=0
  h0[np.where(h0==0)]=0

  
  return V,h0

### END FUNCTION VINTEGR ###############################################


### FUNCTION GET_SECTION ###################################################

def get_section(lonsec,latsec,var,lon,lat,rmask,h,zeta,theta_s,theta_b,hc,N,vtransform):
#
#  Extract a vertical slice in any direction (or along a curve)
#  from a ROMS netcdf file.
#
# 
# On Input:
# 
#    lonsec      Longitudes of the points of the section. 
#    latsec      Latitudes of the points of the section.
#    var         CROCO 3D variable at rho point
#    lon         Longitudes at rho points
#    lat         Latitudes at rho points
#    rmask
#    h
#    zeta
#    theta_s
#    theta_b
#    hc
#    N
#    vtransform
#
#
# On Output:LON,LAT,X,Z,VAR
#
#    LON         Longitudes of the section (2D matrix).
#    LAT         Latitudes of the section (2D matrix).
#    X           Slice X-distances (km) from the first point (2D matrix).
#    Z           Slice Z-positions (matrix). 
#    VAR         Slice of the variable (matrix).
#
#

  Npts=np.size(lonsec)

#
# Get a first subgrid limited by the size of the section (extended by twice the resolution dl)
#


  (dlony,dlonx)=np.gradient(lon)
  (dlaty,dlatx)=np.gradient(lat)
  dl=2*np.max((np.max(np.abs(dlonx)),np.max(np.abs(dlony)),np.max(np.abs(dlatx)),np.max(np.abs(dlaty))))
  minlon=np.min(lonsec)-dl
  minlat=np.min(latsec)-dl
  maxlon=np.max(lonsec)+dl
  maxlat=np.max(latsec)+dl

  sub=( (lon>minlon) & (lon<maxlon) & (lat>minlat) & (lat<maxlat) )
  (jndx,indx)=np.where(sub)
  imin=np.min(indx)
  imax=np.max(indx)
  jmin=np.min(jndx)
  jmax=np.max(jndx)


  lon=lon[jmin:jmax,imin:imax]
  lat=lat[jmin:jmax,imin:imax]
  zeta=zeta[jmin:jmax,imin:imax]
  h=h[jmin:jmax,imin:imax]
  rmask=rmask[jmin:jmax,imin:imax]
  var=var[:,jmin:jmax,imin:imax]

  masksec=np.nan+0*lonsec
  zetasec=np.nan+0*lonsec
  HSEC=np.nan+0*lonsec
  VAR=np.nan+np.zeros((N,Npts))

#
# Do an interpolation for each point of the section
#

  for i in np.arange(0,Npts):

#
# Get a subgrid around each point
#

    sub=( (lon>lonsec[i]-dl) & (lon<lonsec[i]+dl) & \
          (lat>latsec[i]-dl) & (lat<latsec[i]+dl) )
    (jndx,indx)=np.where(sub)
    imin=np.min(indx)
    imax=np.max(indx)
    jmin=np.min(jndx)
    jmax=np.max(jndx)

    londata=lon[jmin:jmax,imin:imax].ravel()
    latdata=lat[jmin:jmax,imin:imax].ravel()
    maskdata=rmask[jmin:jmax,imin:imax].ravel()
    zetadata=zeta[jmin:jmax,imin:imax].ravel()
    hdata=h[jmin:jmax,imin:imax].ravel()

#
# Get the mask as nearest point
#

    masksec[i]=griddata((londata,latdata),maskdata, (lonsec[i], latsec[i]), method='nearest')

#
# If there is enough data, do the interpolation
#
 
    isgood=np.where(np.isfinite(zetadata))

    if np.size(isgood)>4:
      zetasec[i]=griddata((londata[isgood],latdata[isgood]),zetadata[isgood], (lonsec[i], latsec[i]), method='cubic')
      HSEC[i]=griddata((londata[isgood],latdata[isgood]),hdata[isgood], (lonsec[i], latsec[i]), method='cubic')
  
      for k in np.arange(0,N):
        vdata=var[k,jmin:jmax,imin:imax].ravel()
        VAR[k,i]=griddata((londata[isgood],latdata[isgood]),vdata[isgood], (lonsec[i], latsec[i]), method='cubic')

#
# Get the vertical positions
#
 
  zetasec[np.where(np.isnan(zetasec))]=0.
  HSEC[np.where(np.isnan(HSEC))]=hc
  Z=vgrd.zlevs(HSEC,zetasec,theta_s,theta_b,hc,N,'r',vtransform)

#
# Get the horizontal positions
#
 

# Get the horizontal positions
#
  dx=0.*lonsec
  for i in np.arange(1,Npts):
    dx[i]=1e-3*tpx.spheric_dist(latsec[i-1],latsec[i],lonsec[i-1],lonsec[i])

  dist=np.cumsum(dx)
  X=np.matlib.repmat(dist,N,1)
  LON=np.matlib.repmat(lonsec,N,1)
  LAT=np.matlib.repmat(latsec,N,1)

  return LON,LAT,X,Z,VAR




### FUNCTION GET_UV_SECTION ###################################################

def get_UV_section(lonsec,latsec,u,v,lon,lat,rmask,h,zeta,theta_s,theta_b,hc,N,vtransform):

#
#  Extract a vertical slice of cross-section and along section velocities
#  in any direction (or along a curve) from a ROMS netcdf file.
#
# 
# On Input:
# 
#    lonsec      Longitudes of the points of the section. 
#    latsec      Latitudes of the points of the section.
#    u           CROCO eastward velocities at rho point
#    v           CROCO northward velocities at rho point
#    lon         Longitudes at rho points
#    lat         Latitudes at rho points
#    rmask
#    h
#    zeta
#    theta_s
#    theta_b
#    hc
#    N
#    vtransform
#
#
# On Output: LON,LAT,X,Z,Un,Ut
#
#    LON         Longitudes of the section (2D matrix).
#    LAT         Latitudes of the section (2D matrix).
#    X           Slice X-distances (km) from the first point (2D matrix).
#    Z           Slice Z-positions (2D matrix). 
#    Un          Slice cross-section velocities (2D matrix). 
#    Ut          Slice along-section velocities (2D matrix). 
#
#

  Npts=np.size(lonsec)
    
  dlon=np.gradient(lonsec)
  dlat=np.gradient(latsec)
  Rearth=6367442.76
  deg2rad=np.pi/180
  dx=Rearth*deg2rad*dlon*np.cos(deg2rad*latsec)
  dy=Rearth*deg2rad*dlat
    
  dl=np.sqrt(dx*dx +dy*dy)
  
  tx=dx/dl
  ty=dy/dl
  nx=-ty
  ny=tx

#
# Get a first subgrid limited by the size of the section (extended by twice the resolution dl)
#


  (dlony,dlonx)=np.gradient(lon)
  (dlaty,dlatx)=np.gradient(lat)
  dl=2*np.max((np.max(np.abs(dlonx)),np.max(np.abs(dlony)),np.max(np.abs(dlatx)),np.max(np.abs(dlaty))))
  minlon=np.min(lonsec)-dl
  minlat=np.min(latsec)-dl
  maxlon=np.max(lonsec)+dl
  maxlat=np.max(latsec)+dl

  sub=( (lon>minlon) & (lon<maxlon) & (lat>minlat) & (lat<maxlat) )
  (jndx,indx)=np.where(sub)
  imin=np.min(indx)
  imax=np.max(indx)
  jmin=np.min(jndx)
  jmax=np.max(jndx)


  lon=lon[jmin:jmax,imin:imax]
  lat=lat[jmin:jmax,imin:imax]
  zeta=zeta[jmin:jmax,imin:imax]
  h=h[jmin:jmax,imin:imax]
  rmask=rmask[jmin:jmax,imin:imax]
  u=u[:,jmin:jmax,imin:imax]
  v=v[:,jmin:jmax,imin:imax]

  masksec=np.nan+0*lonsec
  zetasec=np.nan+0*lonsec
  HSEC=np.nan+0*lonsec
  Un=np.nan+np.zeros((N,Npts))
  Ut=np.nan+np.zeros((N,Npts))

#
# Do an interpolation for each point of the section
#

  for i in np.arange(0,Npts):

#
# Get a subgrid around each point
#

    sub=( (lon>lonsec[i]-dl) & (lon<lonsec[i]+dl) & \
          (lat>latsec[i]-dl) & (lat<latsec[i]+dl) )
    (jndx,indx)=np.where(sub)
    imin=np.min(indx)
    imax=np.max(indx)
    jmin=np.min(jndx)
    jmax=np.max(jndx)

    londata=lon[jmin:jmax,imin:imax].ravel()
    latdata=lat[jmin:jmax,imin:imax].ravel()
    maskdata=rmask[jmin:jmax,imin:imax].ravel()
    zetadata=zeta[jmin:jmax,imin:imax].ravel()
    hdata=h[jmin:jmax,imin:imax].ravel()

#
# Get the mask as nearest point
#

    masksec[i]=griddata((londata,latdata),maskdata, (lonsec[i], latsec[i]), method='nearest')

#
# If there is enough data, do the interpolation
#
 
    isgood=np.where(np.isfinite(zetadata))

    if np.size(isgood)>4:
      zetasec[i]=griddata((londata[isgood],latdata[isgood]),zetadata[isgood], (lonsec[i], latsec[i]), method='cubic')
      HSEC[i]=griddata((londata[isgood],latdata[isgood]),hdata[isgood], (lonsec[i], latsec[i]), method='cubic')
  
      for k in np.arange(0,N):
        vdata=u[k,jmin:jmax,imin:imax].ravel()
        U=griddata((londata[isgood],latdata[isgood]),vdata[isgood], (lonsec[i], latsec[i]), method='cubic')
        vdata=v[k,jmin:jmax,imin:imax].ravel()
        V=griddata((londata[isgood],latdata[isgood]),vdata[isgood], (lonsec[i], latsec[i]), method='cubic')
        Un[k,i]=U*nx[i]+V*ny[i]
        Ut[k,i]=U*tx[i]+V*ty[i]

#
# Get the vertical positions
#
 
  zetasec[np.where(np.isnan(zetasec))]=0.
  HSEC[np.where(np.isnan(HSEC))]=hc
  Z=vgrd.zlevs(HSEC,zetasec,theta_s,theta_b,hc,N,'r',vtransform)

#
# Get the horizontal positions
#
  dx=0.*lonsec
  for i in np.arange(1,Npts):
    dx[i]=1e-3*tpx.spheric_dist(latsec[i-1],latsec[i],lonsec[i-1],lonsec[i])

  dist=np.cumsum(dx)
  X=np.matlib.repmat(dist,N,1)
  LON=np.matlib.repmat(lonsec,N,1)
  LAT=np.matlib.repmat(latsec,N,1)

  return LON,LAT,X,Z,Un,Ut
