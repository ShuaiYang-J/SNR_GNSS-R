# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        GPS_file
# Purpose:     This library compute the julian date, second in the GPS week and the check
#              parameter to account for beggingng and end of week crossovers
#
# Author:      AngelMartin
#
# Created:     15/07/2014
# Final Modification:
#              03/04/2020  
# Copyright:   (c) AngelMartin 2014
# References: 
#               Leick,A. (2004): GPS satellite surveying
#               https://gssc.esa.int/navipedia/index.php/Main_Page
#               Hofmann-Wellenhof et al. (2001): GPS theory and Practice, fifth edition
#-------------------------------------------------------------------------------

def jd(year,month,day,hour):
    if month<=2:
        year=year-1
        month=month+12
    hour_jd=hour/24
    jd=(int(365.25*year)+int(30.6001*(month+1))+day+1720981.5)+hour_jd#Hofmann

    return (jd)#julian date

def gps_time(jd):
    a = int(jd+.5)
    b = a+1537
    c = int((b-122.1)/365.25)
    e = int(365.25*c)
    f = int((b-e)/30.6001)
    d = b-e-int(30.6001*f)+((jd+0.5)-int(jd+0.5))
    day_of_week = (int(jd+0.5))%7
    # We add +1 as the GPS week starts at Saturday midnight
    sec_of_week = ((d%1)+day_of_week+1.0)*86400.0
    return (sec_of_week)#Second in the GPS week

def check_dt(dt_eph):#to account for beggining or end of week crossovers, Leick, pg 84
    if dt_eph>302400:
            dt_eph=dt_eph-604800
    if dt_eph<-302400:
            dt_eph=dt_eph+604800
    return (dt_eph)

def jd2GPSweek(jd):
    a=jd-44244-2400000.5
    week=int(a/7)
    weekday=int(a-week*7)
    return(week,weekday)

