# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
# Name:        GPS_file
# Purpose:     This software compute the satellite coordinate positiong and returns
#              the elevation and azimut form the coordinates of the station
# References: 
#               Leick,A. (2004): GPS satellite surveying
#               https://gssc.esa.int/navipedia/index.php/Main_Page
#-------------------------------------------------------------------------------
# 输入参数 n文件名称 n文件头文件行数 目标卫星历元时间(小时表示) 卫星编号 伪距观测值 接收机位置
def nav (fichero_nav,num_lin_cab_eph,hora_epoca,sv,observa_cion_str,matriz_approx_pos):
    import numpy
    import time_lib
    import coor_trans
    
    omegae_dot = 7.2921151467e-5#Earth rotation speed
    GM = 3.986005e14 # Earth's universal gravitational
    v_light=299792458 #Speed of light

    #INPUT DATA
    p2_obs=float(observa_cion_str)
    num_lin_cab=float(num_lin_cab_eph)
    sv=int(sv)
    hora_epoca=float(hora_epoca)
    x_aprox2=[]
    x_aprox2.append(float((matriz_approx_pos[0])[0]))
    x_aprox2.append(float((matriz_approx_pos[1])[0]))
    x_aprox2.append(float((matriz_approx_pos[2])[0]))
    naveg_obs=open(fichero_nav,'r')
    cont=0
    while cont<num_lin_cab:#We jump the header
        linea=naveg_obs.readline()
        linea=linea
        cont +=1

    while True:
        linea2=naveg_obs.readline()
        if not linea2:break

        obs_nav = linea2
        
        ###Search for the navigation information closest to the observation time 
        
        #To acount for observations at the begginig or end of the navigation file
        if hora_epoca==0 or hora_epoca==23:
            hora_top=2.1
        else:
            hora_top=1.1
            
        try:
            hora_ref=int(obs_nav[15:17])+float(obs_nav[17:20])/60+float(obs_nav[20:23])/3600-hora_epoca

            if (int(obs_nav[1:3])==sv) and abs(hora_ref)<hora_top:
                year_eph=int(obs_nav[4:8])
                month_eph=int(obs_nav[9:11])
                day_eph=int(obs_nav[12:14])
                linea3=naveg_obs.readline()
                linea4=naveg_obs.readline()
                linea5=naveg_obs.readline()
                linea6=naveg_obs.readline()
                linea7=naveg_obs.readline()
                
                #Emission time calculation
                jde=time_lib.jd(year_eph,month_eph,day_eph,hora_epoca)
                sec_of_week=time_lib.gps_time(jde)
    
                t_emision=sec_of_week-p2_obs/v_light
                
                toe=(float((linea5[0:23].replace("D","e"))))
                dt_eph=time_lib.check_dt(t_emision-toe)
                t_corr=(float((obs_nav[23:42].replace("D","e"))))+ \
                       (float((obs_nav[42:61].replace("D","e"))))*dt_eph+ \
                       (float((obs_nav[61:80].replace("D","e"))))*dt_eph*dt_eph                 
                dt_eph_corr=time_lib.check_dt((t_emision-t_corr)-toe)
                t_corr=(float((obs_nav[23:42].replace("D","e"))))+ \
                       (float((obs_nav[42:61].replace("D","e"))))*dt_eph_corr+ \
                       (float((obs_nav[61:80].replace("D","e"))))*dt_eph_corr*dt_eph_corr
                tx_GPS=time_lib.check_dt((t_emision-t_corr)-toe) 
                
                semi_eje=(float((linea4[61:80].replace("D","e"))))**2
                mean_motion=numpy.sqrt(GM/(semi_eje**3))+float((linea3[42:61].replace("D","e")))
                mean_anomaly=float((linea3[61:80].replace("D","e")))+mean_motion*tx_GPS
                ano_ecc=mean_anomaly
                ecc=float((linea4[23:42].replace("D","e")))
                for i_ano_ecc in range(10):
                    ano_ecc = mean_anomaly+ecc*numpy.sin(ano_ecc)
                true_ano = numpy.arctan2((numpy.sqrt(1-ecc**2)*numpy.sin(ano_ecc)), \
                          (numpy.cos(ano_ecc)-ecc))
                phi = true_ano+float((linea6[42:61].replace("D","e")))
                
                #parameters readed from navigation file
                cuc=float((linea4[0:23].replace("D","e")))
                cus=float((linea4[42:61].replace("D","e")))
                crc=float((linea6[23:42].replace("D","e")))
                crs=float((linea3[23:42].replace("D","e")))
                i0=float((linea6[0:23].replace("D","e")))
                idot=float((linea7[0:23].replace("D","e")))                
                cic=float((linea5[23:42].replace("D","e")))
                cis=float((linea5[61:80].replace("D","e")))
                omega0=float((linea5[42:61].replace("D","e")))
                omegadot=float((linea6[61:80].replace("D","e")))
                
                ### Satellite coordinates
                du=cuc*numpy.cos(2*phi)+cus*numpy.sin(2*phi)
                dr=crc*numpy.cos(2*phi)+crs*numpy.sin(2*phi)
                di=cic*numpy.cos(2*phi)+cis*numpy.sin(2*phi)
                u_eph=phi+du
                r_eph=semi_eje*(1-ecc*numpy.cos(ano_ecc))+dr
                i_eph=i0+idot*tx_GPS+di
                omega=omega0+(omegadot-omegae_dot)*tx_GPS-omegae_dot*toe
                x1=numpy.cos(u_eph)*r_eph
                y1=numpy.sin(u_eph)*r_eph
                sat_x=x1*numpy.cos(omega)-y1*numpy.cos(i_eph)*numpy.sin(omega)
                sat_y=x1*numpy.sin(omega)+y1*numpy.cos(i_eph)*numpy.cos(omega)
                sat_z=y1*numpy.sin(i_eph)
                # 可以用sp3文件直接导入卫星精密位置 omegae_dot是常数 后续函数不变 但是不够实时
                ###Rotated satellite coordinates
                omegatau=(p2_obs/v_light)*omegae_dot
                sat_x_rot=sat_x*numpy.cos(omegatau)+sat_y*numpy.sin(omegatau)
                sat_y_rot=-sat_x*numpy.sin(omegatau)+sat_y*numpy.cos(omegatau)
                sat_z_rot=sat_z
        
                (lat_ECEF,lon_ECEF,h_ECEF)=coor_trans.geod(x_aprox2)
                
                (el_sat,azi_sat)=coor_trans.horiz(lat_ECEF,lon_ECEF,h_ECEF,x_aprox2,sat_x_rot,sat_y_rot,sat_z_rot)
                el_sat=(el_sat*180/numpy.pi)
                azi_sat=(azi_sat*180/numpy.pi)
                naveg_obs.close()
                
                return(sat_x_rot,sat_y_rot,sat_z_rot,el_sat,azi_sat)
                                
                # break

        except:
            pass

# sp3文件的读取 9历元内插
def sp3 (fichero_sp3,num_lin_sp3_eph,sv,matriz_approx_pos,numepoh,satnum,obs_cab,sp3body,obsgpstime,gpssecond):
    import math
    import numpy as np
    import time_lib
    import coor_trans
    # 常数定义
    omegae_dot = 7.2921151467e-5#Earth rotation speed
    v_light=299792458 #Speed of light

    #对sp3文件进行差分处理
    # 添加存储矩阵 初始化
    sp3x=np.zeros(9)
    sp3y=np.zeros(9)
    sp3z=np.zeros(9)
    sp3t=np.zeros(9)
    sp3 = open(fichero_sp3,'r')
    iv=0
    firstGPS=gpssecond[0]
    gpstime=obsgpstime
    if  numepoh==96:
        iv = int(1 + (gpstime - firstGPS)/(15*60))
    if  numepoh==288 or numepoh==289:
        iv = int(1 + (gpstime - firstGPS)/(5*60))
    if  numepoh==2881 or numepoh==2880:
        iv = int(1 + (gpstime - firstGPS)/(30))
#  iv 在sp3中的目标历元数
# 将时间存储 后续进行读取
    sp3.close()
    trigpt = numepoh -5
    i1 = iv-4 
    i2 = iv+4 
    if (iv<5):
        i1 = 1
        i2 = 9
    if(iv>trigpt):
        i1=numepoh-8
        i2=numepoh
# 判定所需的9个历元顺序 sp3存储在sp3body

    numepoh=int(numepoh)
    satnum=int(satnum)
    sv=int(sv)
    x_aprox2=[]
    x_aprox2.append(float((matriz_approx_pos[0])[0]))
    x_aprox2.append(float((matriz_approx_pos[1])[0]))
    x_aprox2.append(float((matriz_approx_pos[2])[0]))

# 读取所需的9个历元 对特定一颗卫星 sv

    k=0
    for i in range(i1,i2+1,1):
        sp3t[k]=gpssecond[i-1]
        sp3x[k]=sp3body[i-1,sv-1,0]*1000
        sp3y[k]=sp3body[i-1,sv-1,1]*1000
        sp3z[k]=sp3body[i-1,sv-1,2]*1000
        k +=1
# 差分
    toffset=0.07
    for i in range(3):
        sat_x=dpolint(sp3t,sp3x,obsgpstime-toffset)
        sat_y=dpolint(sp3t,sp3y,obsgpstime-toffset)
        sat_z=dpolint(sp3t,sp3z,obsgpstime-toffset)
        omegatau=omegae_dot*toffset
        xnew=sat_x*math.cos(omegatau)+sat_y*math.sin(omegatau)
        ynew=-sat_x*math.sin(omegatau)+sat_y*math.cos(omegatau)
        znew=sat_z
        range1= (xnew-x_aprox2[0])**2 + (ynew-x_aprox2[1])**2 + (znew-x_aprox2[2])**2
        toffset = math.sqrt(range1)/v_light
        #         # 
    (lat_ECEF,lon_ECEF,h_ECEF)=coor_trans.geod(x_aprox2)
    (el_sat,azi_sat)=coor_trans.horiz(lat_ECEF,lon_ECEF,h_ECEF,x_aprox2,xnew,ynew,znew)
        #         # (el_sat,azi_sat)=coor_trans.horiz(lat_ECEF,lon_ECEF,h_ECEF,x_aprox2,sat_x,sat_y,sat_z)
    el_sat=(el_sat*180/math.pi)
    azi_sat=(azi_sat*180/math.pi)
    return(xnew,ynew,znew,el_sat,azi_sat)
    # return(0,0,0,0,0)                    

def readsp3body(filename,num_lin_sp3_eph,numepoh,satnum):
    import numpy as np
    import time_lib
    # 读取sp3的整体 按照历元 卫星 X Y Z 存储信息[noph][4]
    sp3 = open(filename,'r')
    cont=0
    gpssecond=[]
    while cont<num_lin_sp3_eph:
        linea=sp3.readline()
        linea=linea
        cont+=1
# 跳过表头
    sp3body = np.zeros((numepoh,satnum,3))
    # sp3body=[]
    epnum=0
    while True:
        linea1=sp3.readline()
        if not linea1:break
        if 'EOF' in linea1:
            break
        if '*' in linea1:
            sp3_eoph_time=linea1.split()
            jd=time_lib.jd(int(sp3_eoph_time[1]),int(sp3_eoph_time[2]),int(sp3_eoph_time[3]),int(sp3_eoph_time[4]))+(int(sp3_eoph_time[5])/(60*24))+(float(sp3_eoph_time[6])/(60*60*24))
            sec=time_lib.gps_time(jd)
            gpssecond.append(sec)
            epnum +=1
            linea1=sp3.readline()
        # sp3body[epnum-1,satnum,0]=int(linea1[2:4])
        if 'PG' in linea1:
            satnum1=int(linea1[2:4])-1
            sp3body[epnum-1,satnum1,0]=float(linea1[5:18])
            sp3body[epnum-1,satnum1,1]=float(linea1[19:32])
            sp3body[epnum-1,satnum1,2]=float(linea1[33:46])
    return(sp3body,gpssecond)

def dpolint(sp3t,sp3xyz,tc):
    # 差分过程 n表示次数 可能有点问题
    out=0
    n=9
    ns=0
    dif=abs(tc-sp3t[0])
    for i in range(9):
        dift=abs(tc-sp3t[i])
        if dift<dif :
            ns=i
            dif=dift
    # 数组赋值 拷贝
    c=sp3xyz.copy()
    d=sp3xyz.copy()
    out=sp3xyz[ns]
    if ns!=0:
        ns=ns-1
    for m in range(1,n):
        for i in range(n-m):
            ho=sp3t[i]-tc
            hp=sp3t[i+m]-tc
            w=c[i+1]-d[i]
            den=ho-hp
            if den==0:
                break
            den=w/den
            d[i]=hp*den
            c[i]=ho*den
        if 2*ns<(n-m):
            dy=c[ns+1]
        else:
            dy=d[ns]
            ns=ns-1
        out=out+dy
    return out
