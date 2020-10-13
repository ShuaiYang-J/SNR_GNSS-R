# -*- coding: utf-8 -*-

# observation and navigation RINEX files (version 3)
#epoch and satellite,id,SNR,elevation and azimut of the satellite, and year, month and day.
#-------------------------------------------------------------------------------

import os
import read_nav
import datetime
import time_lib
####USER INPUT PARAMETERS
path='test/'#输出 输入目录
# path='data/input/'
pseu='C1C'#观测伪距的频率
se='S1C'#SNR的信号频率
####END USER INPUT PARAMETERS

system='G'#G for GPS satellite constellation
if system=='G':
    nav_sys='n' 
    # fichero_sal=open(path+'output1/SNR-GPS.txt','w')
    fichero_sal=open(path+'SNR-GPS.txt','w')
    fichero_sal.write('epoch;epoch_time;sat;SNR;Elev;Azi;year;month;day \n')

###标记了文件生成的时间 path_in 输入文件目录
path_in=path
epoca=0
lista_dias=[]
# 读取文件夹 选择o文件并存储
filename = next(os.walk(path_in))[2]
for i in filename:
    if i[-1] =='o':
        fichero = (path+i)
        name=str(i[0:4])
        year_a=int(i[9:11])+2000
        doi=int(i[4:7])        
        d=datetime.date(year_a,1,1)+datetime.timedelta(doi-1)
        dl=str(year_a)+"-"+str(d.month)+'-'+str(d.day)
        lista_dias.append(str(dl))        
dates = [datetime.datetime.strptime(ts, "%Y-%m-%d") for ts in lista_dias]
dates.sort()
lista_diasb=[]
for i in dates:
     aa = i.strftime('%j')
     dl2=name+aa
     lista_diasb.append(dl2)
##########################

#####This block of the code read all observation files containing in the previous list
for i in lista_diasb:    
    for j in filename:
    
        if i in j and j[-1] =='o':
            fichero=path_in+j
            print ("working on file=", j)
            
            ####calculation of the number of lines contained in the header of the navigation file            
            fichero_nav=path_in+j[:-1]+nav_sys
            naveg_cab = open(fichero_nav,'r')
            num_lin_cab_eph=0
            for linea in naveg_cab.readlines():
                if 'RINEX VERSION / TYPE' in linea:
                    Version=int(linea[5:6])
                num_lin_cab_eph +=1            
                if 'END OF HEADER' in linea:
                    break                       
            naveg_cab.close()
            ####

            ####read observation file header information
            archivo_cab = open(fichero,'r')
            num_lin_cab=0
            while True:
                linea=archivo_cab.readline()
                num_lin_cab +=1
               
                if 'APPROX POSITION XYZ' in linea:
                    approx_pos = linea.split() 
                    matriz_approx_pos=[]
                    for i_approx_pos in range(3):
                        matriz_approx_pos.append([])
                        for j_approx_pos in range(1):
                            matriz_approx_pos[i_approx_pos].append(float(approx_pos[i_approx_pos]))

                if 'OBS TYPES' in linea and system in linea:
                    types_obs = linea.split() 
                    try:
                        num_types_obs = int(types_obs[1])

                        if num_types_obs>13:
                            pos_f_1=linea.split()
                            linea=archivo_cab.readline()
                            pos_f_2=linea.split()
                            pos=pos_f_1[0:15]+pos_f_2
                        else:
                            pos=linea.split()
                        pos_p2=pos.index(pseu)-2
                        pos_s2=pos.index(se)-2
                        #     pos_p2_1=linea.split()
                        #     pos_p2_2=pos_p2_2.split()
                        # pos_p2=linea.split().index(pseu)-2
                        # pos_s2=linea.split().index(se)-2
                    except:
                        pass
 
                if 'TIME OF FIRST OBS' in linea:
                    time_first_obs = linea.split()
                    year = int(time_first_obs[0])
                    mes = int(time_first_obs[1])
                    dia = int(time_first_obs[2])
                    # 计算GPS周 jd儒略日
                    jd=time_lib.jd(year,mes,dia,0)
                    GPSweek,GPSday=time_lib.jd2GPSweek(jd)
                if 'END OF HEADER' in linea:
                    break
            archivo_cab.close()  
            ####
            # 读取sp3文件的头部
            fichero_sp3=path_in+'igs'+str(GPSweek)+str(int(GPSday))+'.sp3'
            sp3_cab = open(fichero_sp3,'r')
            num_lin_sp3_eph=0
            for linea in sp3_cab.readlines():
                num_lin_sp3_eph +=1            
                # sp3的截至标志
                if num_lin_sp3_eph==1:
                    numepoh=int(linea[36:39])
                if num_lin_sp3_eph==3: 
                    satnum=int(linea[4:6])
                if '*'==linea[0:1]:
                    sp3linea=linea.split()
                    if sp3linea[1]==str(year) and sp3linea[2]==str(mes) and sp3linea[3]==str(dia):
                        num_lin_sp3_eph=num_lin_sp3_eph-1
                        break                        
            sp3_cab.close()
            sp3body,gpssceond=read_nav.readsp3body(fichero_sp3,num_lin_sp3_eph,numepoh,satnum)
            ####Read observation file data information and calculation of the satellite azimut and elevation
            archivo_obs = open(fichero,'r')
            cont=0            
            while cont<num_lin_cab:#jump the header
                linea=archivo_obs.readline()
                cont +=1
                
            num_sat_acum=0
            cont_ecu=0
            while True and Version==3:
                linea2=archivo_obs.readline()
                if not linea2:break
                epoca +=1
                obs_cab = linea2.split()
                if int(obs_cab[1])==year and int(obs_cab[2])==mes and int(obs_cab[3])==dia:
                       hora_epoca=float(obs_cab[4])+float(obs_cab[5])/60+float(obs_cab[6])/3600
                       obsjd=time_lib.jd(int(obs_cab[1]),int(obs_cab[2]),int(obs_cab[3]),int(obs_cab[4]))+(int(obs_cab[5])/(60*24))+(float(obs_cab[6])/(60*60*24))
                       obsgpstime=time_lib.gps_time(obsjd)
                       ok_flag=int(obs_cab[7])
                       num_sat=int(obs_cab[8])            
                       ##for every satellite of the epoch
                       for i_sat in range(num_sat):
                           linea3=archivo_obs.readline()
                           lect_sats=linea3.split()
                           if system in lect_sats[0] and ok_flag==0:
                               ini_p2=3+pos_p2*16
                               ini_s2=3+pos_s2*16
                               
                               id_sat=lect_sats[0][1:3]
                               try:
                                   C2=float(linea3[ini_p2:ini_p2+14])
                                   S2=float(linea3[ini_s2:ini_s2+14])
                                   if system=='G':
                                    # n文件实时性更好
                                    #    (x_sat,y_sat,z_sat,el_sat,azi_sat)=read_nav.nav(fichero_nav,num_lin_cab_eph,hora_epoca,id_sat,C2,matriz_approx_pos)
                                    # sp3作为验证数据
                                       (x_sat,y_sat,z_sat,el_sat,azi_sat)=read_nav.sp3(fichero_sp3,num_lin_sp3_eph,id_sat,matriz_approx_pos,numepoh,satnum,obs_cab,sp3body,obsgpstime,gpssceond)
                                   ### Only store satellite information from elavation less than 30 degrees andhigher than 5 degrees
                                   if el_sat<30 and el_sat>5 :
                                       nueva_fil=epoca,hora_epoca,int(id_sat),S2,el_sat,azi_sat,year,mes,dia
                                       nueva_fil=str(nueva_fil).replace('(','').replace(')','').replace(',',';')
                                       fichero_sal.write(nueva_fil+'\n')
                               except:
                                  pass
fichero_sal.close()                       
