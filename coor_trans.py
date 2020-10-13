# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        GPS_file
# Purpose:     This library compute the geographic spherical coordinates from the geocentric cartesian coodinates
#              and the horizontal coodinates (azimut and elevation from station to satellite)
#
# Author:      AngelMartin
#
# Created:     12/08/2014
# Final Modification:
#              03/04/2020  
# Copyright:   (c) AngelMartin 2014
# References: 
#               Hofmann-Wellenhof et al. (2001): GPS theory and Practice, fifth edition
#-------------------------------------------------------------------------------
import numpy

###WGS84 elipsoid constants
a_wgs84 = 6378137.0
b_wgs84 = 6356752.3142
f_wgs84 = 1/298.257223563
e_wgs84 = numpy.sqrt((a_wgs84*a_wgs84-b_wgs84*b_wgs84)/(a_wgs84*a_wgs84))

def geod(x_approx):
    #x_approx is a vector containing the ECEF geocentric coordinates of the station in meters
    X_ECEF = x_approx[0]
    Y_ECEF = x_approx[1]
    Z_ECEF = x_approx[2]

    itera_lat = 10

    lon_ECEF = numpy.arctan(Y_ECEF/X_ECEF)
    p_ECEF = numpy.sqrt(X_ECEF*X_ECEF+Y_ECEF*Y_ECEF)
    lat_ECEF = numpy.arctan((Z_ECEF/p_ECEF)*((1-e_wgs84*e_wgs84)**(-1)))

    for i_ECEF in range (itera_lat):
        bajo = numpy.sqrt((a_wgs84*a_wgs84*(numpy.cos(lat_ECEF))**2)+(b_wgs84*b_wgs84*(numpy.sin(lat_ECEF))**2))
        N_ECEF = (a_wgs84*a_wgs84)/bajo
        h_ECEF = (p_ECEF/numpy.cos(lat_ECEF))-N_ECEF
        lat_ECEF = numpy.arctan((Z_ECEF/p_ECEF)*(1-(e_wgs84*e_wgs84)*(N_ECEF/(N_ECEF+h_ECEF)))**(-1))

    return(lat_ECEF,lon_ECEF,h_ECEF)

def horiz(lat_ECEF,lon_ECEF,h_ECEF,x_aprox2,sat_x_rot,sat_y_rot,sat_z_rot):
    #lat_ECEF is the geocentric latitude in radians of the station
    #lon_ECEF is the geocentric longitude in radians of the station
    #h_ECEF is the ellipsoidal height in meters of the station
    #x_approx is a vector containing the ECEF geocentric coordinates of the station in meters
    #sat_x_rot, sat_y_rot and sat_z_rot are the geocentric coordinates (meters) of the satellite in the recepcion time
    dx_ECEF = sat_x_rot-x_aprox2[0]
    dy_ECEF = sat_y_rot-x_aprox2[1]
    dz_ECEF = sat_z_rot-x_aprox2[2]
    
    a11=-numpy.sin(lat_ECEF)*numpy.cos(lon_ECEF)
    a12=-numpy.sin(lat_ECEF)*numpy.sin(lon_ECEF)
    a13=numpy.cos(lat_ECEF)
    a21=-numpy.sin(lon_ECEF)
    a22=numpy.cos(lon_ECEF)
    a23=0
    a31=numpy.cos(lat_ECEF)*numpy.cos(lon_ECEF)
    a32=numpy.cos(lat_ECEF)*numpy.sin(lon_ECEF)
    a33=numpy.sin(lat_ECEF)

    dN=a11*dx_ECEF+a12*dy_ECEF+a13*dz_ECEF
    dE=a21*dx_ECEF+a22*dy_ECEF+a23*dz_ECEF
    dh=a31*dx_ECEF+a32*dy_ECEF+a33*dz_ECEF

    dis_hor=numpy.sqrt(dN*dN+dE*dE)

    if dis_hor<0.0001:
        el_sat=numpy.pi/2
    else:
        el_sat=numpy.arctan2(dh,dis_hor)

    azi_sal=numpy.arctan2(dE,dN) 

    if azi_sal<0:
        azi_sal=azi_sal+2*numpy.pi

    return(el_sat,azi_sal) #satellite elevaion and azimut in radians

