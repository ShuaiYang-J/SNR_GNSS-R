# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        SNR_list_sat
# Purpose:     This software read the file SNR-GPS.txt
#               and creates a sructure of directories, one per GPS satellite.
#               Each directory cointains the plots for the process of every
#               single rise or fall satellite track and a .txt file containing
#               the information of epoch, valid track (yes/no), month, day, id_sat, initial azimut of 
#               the track (in degrees), final azimuth of the track (in degrees),
#               initial elevtion of the traack (in degrees)
#               final elevation of the track (in degrees), amplitude of the adjusted
#               wave, standard deviation of the adjusted amplitude, phase of
#               adjusted wave (in degrees), standard deviation of the adjusted wave
#               (in degrees), computed reflector height (in m)
#               and satellite direction (falling or rising).
#               A second .txt file is also generating wiht the same information
#               but only for those tracks that pass control conditions

#-------------------------------------------------------------------------------

import os
from astropy.stats import LombScargle # to astropy V2
#from astropy.timeseries import LombSacrgle # from astropy V3
import numpy
import math
from math import pi
import matplotlib.pyplot as plt 

####USER INPUT PARAMETERS
#ocation of the output directory (THIS FOULDER SHOULD EXIST)
figPath='test/output/'
#INPUT FILE
archivo='test/SNR-GPS.txt'#Input file

##PARAMETERS
#observation data interval (from the GPS observation) s 秒
int_tiempo=5
#number of empty epochs to consider the readings as the same satellite track,
#if the readings interval are greather that this parameter, the track is divided in two
#different tracks
epocas_offset=7 ##35 for 1 second interval, 7 for a 5 seconds interval and 1 for 30 seconds interval
#mininum number of degrees in elevation to cover by a satellite to be considered
inc_eleva=10#
#minimun elevation of the satellite be considered
el_min=5
#Maximum elevation of the satellite to be considered
el_max=30
#minimum azimuth of the satellite from the station to be considered
azi_min=0 
#maximum azimuth of the satellite from the station to be considered
azi_max=360 
#antenna height measured in field 天线高
h0=1.80
#GPS L1 or L2 wavelength of the signals. L1=0.1904, L2=0.2443 波长
land=0.1904 

#Creating the directories structure (one directory per satellite)
#These folders must not exist, delete them or save them elsewhere before running the software again, 
#otherwise you will have an error in the software
# 单GPS
for i in range(32):
    dir_dib=figPath+str(i+1)+'/'
    os.mkdir(dir_dib)

    

h_ant=h0

#Creation of empty lists
el_lis=[] # Elevatio
SNR_lis=[] # SNR
SNRi_lis=[] # indirect SNR
e_lis=[] # epoch
sat_lis=[] # Satellite
az_lis=[] # Azimuth
SNRv_lis=[] # SNR in volts
h_lis=[] # Hour
e_lis1=[] # lists with the epochs for every track
#creation of empyt lists of lists
el_lis2=[]
SNR_lis2=[]
e_lis2=[]
SNRv_lis2=[]
az_lis2=[]
h_lis2=[]

### generation of the list of satellites ###
cont_2=0
num_lin_cab2 = 1
# 打开记录的txt文件
matriz2=open(archivo,'r') 

# jump the header
while cont_2<num_lin_cab2:
    linea3=matriz2.readline()
    cont_2 +=1

while True:
    linea4=matriz2.readline()
    if not linea4:break
# split() 可选隔断符号
    datos2 = linea4.split(';') 

    sat2=datos2[2] # id of the satellite
    if sat2 not in sat_lis: # if the satellite is not in the list, is added on it
        sat_lis.append(sat2) 
        
matriz2.close()

## For every satellite on the list ##

for q in range(len(sat_lis)):
    sat=sat_lis[q]
    print ("working on satellite:",str(sat))

   # empty the list every satellite start the process
    el_lis=[] 
    SNR_lis=[] 
    SNRi_lis=[] 
    SNRi_lis2=[]
    e_lis=[] 
    el_lis2=[]
    SNR_lis2=[]
    e_lis2=[]
    e_lis1=[]
    e_lis2=[]
    SNRv_lis=[]
    SNRv_lis2=[]
    az_lis=[]
    az_lis2=[]
    h_lis=[]
    h_lis2=[]    
    mes_lis=[]
    mes_lis2=[]
    dia_lis=[]
    dia_lis2=[]
    #frecuencia=[]
    frecuencia2=[]
    power2=[]

    # reading the file again, jumping the header
    cont=0
    num_lin_cab = 1        
    matriz=open(archivo,'r') 
    while cont<num_lin_cab:
        linea=matriz.readline()
        cont +=1

    while True:
        linea2=matriz.readline()
        if not linea2:break

        datos = linea2.split(';') 
        
        #Conditions for azimuth and alevation for the satellites
        if sat in datos[2] and float(datos[4])>el_min and float(datos[4])<el_max and float(datos[5])>azi_min and float(datos[5])<azi_max: 
            
            #Generation of the output files
            # 满足设定条件的卫星
            salida = open(figPath+str(sat.strip())+'/'+'final_values.txt', 'w') 
            #Header line
            salida.write('epoch good_track epoch_time month day id_sat Azi_ini aZi_end ele_ini ele_end A desvA Phi desvPhi ref_height SatelliteDirection\n') 
            salida2 = open(figPath+str(sat.strip())+'/'+'final_values_good.txt', 'w') 
            #Header line
            salida2.write('epoch good_track epoch_time month day id_sat Azi_ini aZi_end ele_ini ele_end A desvA Phi desvPhi ref_height SatelliteDirection\n') 
                                    
            e=float(datos[0]) # read the epocj
            e_lis.append(e) # generate the list with the epochs
            
            try:
                # are two epochs consecutive?
                e2= (e_lis[-2])
                if (e-e2)<epocas_offset:
                
                    # read data
                    az=float(datos[5])
                    el=float(datos[4])
                    SNR=float(datos[3])
                    e1=float(datos[0])
                    h=str(datos[1])
                    anio=datos[-3]
                    mes=datos[-2]
                    dia=datos[-1].strip()
                
                    SNR_v=10**(SNR/20) # SNR in volts

                    # append readed values to the lists
                    az_lis.append(az)
                    el_lis.append(el)
                    SNR_lis.append(SNR)
                    SNRv_lis.append(SNR_v)
                    e_lis1.append(e1)
                    h_lis.append(h)
                    mes_lis.append(mes)
                    dia_lis.append(dia)
                
                else: # if a jump between epoch is readen in the input file
                
                    if len(el_lis)!=0 and len(e_lis1)!=0: # if the lists are not empty
                    # store the list in a "list of lists" and empty the list
                    # Obtention of lists with as mach lists are tracks without interruptions exists
                        el_lis2.append(el_lis)
                        el_lis=[]
                        SNR_lis2.append(SNR_lis)
                        SNR_lis=[]
                        SNRv_lis2.append(SNRv_lis)
                        SNRv_lis=[]
                        az_lis2.append(az_lis)
                        az_lis=[]
                        e_lis2.append(e_lis1)
                        e_lis1=[]
                        h_lis2.append(h_lis)
                        h_lis=[]                        
                        mes_lis2.append(mes_lis)
                        mes_lis=[]                        
                        dia_lis2.append(dia_lis)
                        dia_lis=[]

            except:
                pass
    
    if len(el_lis)!=0 and len(e_lis1)!=0: # if the lists are not empty
        # adding the last list to the list of lists
        el_lis2.append(el_lis)
        SNR_lis2.append(SNR_lis)
        SNRv_lis2.append(SNRv_lis)
        az_lis2.append(az_lis)
        e_lis2.append(e_lis1)
        h_lis2.append(h_lis)
        mes_lis2.append(mes_lis)
        dia_lis2.append(dia_lis)

    lon=len(el_lis2) # longitude of the list of lists (valid tracks)
    for i in range(0,lon): # for every satellite track
    
        # load the lists for every track
        el_lis=el_lis2[i]
        SNR_lis=SNR_lis2[i]
        SNRv_lis=SNRv_lis2[i]
        az_lis=az_lis2[i]
        e_lis1=e_lis2[i]
        h_lis=h_lis2[i]
        mes_lis=mes_lis2[i]
        dia_lis=dia_lis2[i]
        obs = len(el_lis) 
        
        SNRi_lis=[] # empty list with the indirect SNR
        # adjustment of a second order polinomial
        SNRi_cal=[]
        coef = numpy.polynomial.polynomial.polyfit(el_lis,SNRv_lis,2)
            
        for j in range(len(el_lis)): # for every element of the list    
           # apply the computed coefficients
            SNR_cal=coef[2]*el_lis[j]*el_lis[j]+coef[1]*el_lis[j]+coef[0]
            SNR_v=SNRv_lis[j]       
            #indirect SNR= (observed SNR - Computed SNR)
            SNRi = SNR_v-SNR_cal        
            SNRi_lis.append(SNRi) # list with the indirect SNR values in volts            
            SNRi_cal.append(SNR_cal)
            
        SNRi_lis2.append(SNRi_lis)

        # Empty time list generation
        ts_lis=[]
        contador=0
 
        # initial value of amplitde and phase for the adjustment
        a0=1.0
        p0=0.0

        # empty lists for the adjustment
        Alis=[]
        Klis=[]
        Plis=[]
        R=[]#residuals

       # initial weights (every observation has a one) 
        for r in range(len(el_lis)):
              paux=1
              Plis.append(paux)

       # 1o iterations in the adjustment
        for it in range(10):
            
            Alis=[]
            Klis=[]

            # for every element of the list (the track we are working on)
            for g in range(len(el_lis)):
                # elevation in radians and indirect SNR in volts
                elev=math.radians(el_lis[g])
                SNRo=SNRi_lis[g]

                # partial derivatives for the design matrix A
                a1=numpy.cos(((4*pi*h0)/land)*numpy.sin(elev)+p0)
                a2=-a0*numpy.sin((((4*pi*h0)/land))*numpy.sin(elev)+p0)

                # Computed SNR 计算公式
                SNR_cal=a0*numpy.cos(((4*pi*h0)/land)*numpy.sin(elev)+p0)
    
                # auxilir lists
                Aaux=[a1,a2]
                SNRaux=[SNRo-SNR_cal]
                
                # Adding to the lists
                Alis.append(Aaux)
                Klis.append(SNRaux)

            # from list to matrix
            A=numpy.matrix(Alis)
            K=numpy.matrix(Klis)

            # weight matrix, weights are fixed to 1 for all observations
            dim = len (Plis) 
            diagonal = numpy.matrix(numpy.eye(dim,dim)) 
            cont_p = 0
            for q in range(dim):
                diagonal[cont_p,cont_p]= Plis[q] 
                cont_p +=1    
            P = diagonal
    
            # least squeare solution
            At=numpy.transpose(A)
            try:
                X=numpy.linalg.inv(At*P*A)*At*P*K 
                Qxx=numpy.linalg.inv(At*P*A)
                
            # Update of the values
                a01=X[0,0]
                p1=X[1,0]
    
                a0=a0+a01
                p0=p0+p1
                R=A*X-K
                gdl=(obs-2)
            #variance-covariance matrix
                var_pos=float(((numpy.transpose(R))*P*R)/(gdl))
    
                Sxx=var_pos*Qxx # 
                desv_A=numpy.sqrt(Sxx[0,0]) # standard deviation for amplitude
                desv_p=numpy.sqrt(Sxx[1,1]) # standard deviation for phase
            except:
                pass

        # list with the values incremented according with the sample observatio data
        # with the same dimention as the track we are working on
        for k in range(len(el_lis)):
            ts=contador*int_tiempo
            ts_lis.append(ts)
            contador=contador+1
            
        # only tracks containing observations for more than a half an hour and covering more than inc_eleva degress
        el_lis_sen=[]    
        if len(el_lis)>(30*60/int_tiempo) and (max(el_lis)-min(el_lis))>inc_eleva:
        
            az1=az_lis[0] # first azimuth of the track
            lon_az=len(az_lis)
            az2=az_lis[lon_az-1] # last azimuth of the track
            
            ep1=e_lis1[0] # first elevation of the track
            lon_e=len(e_lis1)
            ep2=e_lis1[lon_e-1] # last elevation of the track
            
            # LombScargle computation based on the sin of elevation angles
            for i in el_lis:
                el_lis_sen.append(numpy.sin(i*numpy.pi/180))
                
            (frec_sen,power_sen)=LombScargle(el_lis_sen,SNRi_lis).autopower(nyquist_factor=2)                
            frec2_sen=[]
            line_frec2_sen=[]
            aux=0
            for i in frec_sen:                
                
                if i*land/2 <2.5: ##maximum H allowed is of 2.5 meters  最大高度限制
                    frec2_sen.append(i*land/2)
                    line_sec=(power_sen[aux])
                    line_frec2_sen.append(line_sec)
                    aux=aux+1            

            rh=frec2_sen[line_frec2_sen.index(max(line_frec2_sen))] #Reflector height 反射高
            
            vector_el=[]
            vector_azi=[]
            for s in range(len(e_lis1)):
                vector_el.append(round(el_lis[s],1))
                vector_azi.append(round(az_lis[s],1))
            
            #Check if the track is rising or falling            
            ini_ele=vector_el[0]
            fin_ele=vector_el[-1]
            ini_azi=vector_azi[0]
            fin_azi=vector_azi[-1]            
            dir_el=vector_el[-1]-vector_el[0]
            if dir_el>=0:
                direction="satellite-rising"
            else:
                direction="satellite-falling"
            
            #####CONDITIONS to select the "good" tracks            
            
            period_cond=max(line_frec2_sen)/numpy.mean(line_frec2_sen)
           
            x_range=numpy.linspace(5,30,100)
            y_onda=[]
            y2_onda=[]
            y3_onda=[]
            for i in range(len(x_range)):
                y_line=a0*numpy.cos(((4*pi*h0)/land)*numpy.sin(x_range[i]*pi/180)+p0)
                y_onda.append(y_line)
                for j in range(len(el_lis)-1):
                    if (x_range[i]>el_lis[j] and x_range[i]<el_lis[j+1]) or (x_range[i]<el_lis[j] and x_range[i]>el_lis[j+1]):                        
                        y2_onda.append(SNRi_lis[j])
                        y3_onda.append(y_onda[i])
                                                          
            res_onda=[]
            y_onda_red=y3_onda
            for k in range(len(y_onda_red)):                    
                line_res_onda= y_onda_red[k]-y2_onda[k]
                res_onda.append(line_res_onda)
                       
            wave_mean = numpy.mean(res_onda)
            wave_std = numpy.std(res_onda)
            
            if period_cond>6 and abs(wave_mean)<=1.3 and wave_std<25.0 and len(y_onda_red)>85 and len(y2_onda)<110:
                buena="yes"
            else:
                buena="no"
            
            ####END CONDITIONS
            
            ## Save resutls in the output files            
            for s in range(len(e_lis1)):
                if s==0:
                    phase=(p0*180.0/numpy.pi)%360.0                     
                    salida.write(str(e_lis1[s])+' '+buena+' '+str(round(period_cond,2))+' '+str(h_lis[s].strip())+' '+str(mes_lis[s].strip())+' '+str(dia_lis[s].strip())+' '+str(sat.strip())+' '+str(ini_azi)+' '+str(fin_azi)+' '+str(ini_ele)+' '+str(fin_ele)+' '+str(abs(a0))+' '+str(desv_A)+' '+str(phase)+' '+str(desv_p)+' '+str(round(rh,3))+' '+str(direction)+'\n')
                    if buena=="yes":
                        salida2.write(str(e_lis1[s])+' '+buena+' '+str(h_lis[s].strip())+' '+str(mes_lis[s].strip())+' '+str(dia_lis[s].strip())+' '+str(sat.strip())+' '+str(ini_azi)+' '+str(fin_azi)+' '+str(ini_ele)+' '+str(fin_ele)+' '+str(abs(a0))+' '+str(desv_A)+' '+str(phase)+' '+str(desv_p)+' '+str(round(rh,3))+' '+str(direction)+'\n')
                    epoch=str(e_lis1[s])
          
            # Plots generation
            fechi=str(anio)+'-'+str(mes_lis[s])+'-'+str(dia_lis[s].strip()) 
            primer_az=round(az1,1)
            final_az=round(az2,1)
            font_size=8
            nombre_dib=epoch+'-'+str(fechi.strip())+'-'+str(sat.strip())+'-'+str(primer_az)+'-'+str(final_az)

            a_x=[]
            cont=0
            cont_fig=1
            for i in range(len(SNR_lis)):
                a_x.append(cont*int_tiempo)
                cont=cont+1
               
            plt.subplot(221)
            plt.plot(a_x,SNR_lis,'k')
            plt.title('(a) Complete Signal',fontsize=10,fontweight='bold')
            plt.xlabel('Number of Seconds',fontsize=font_size,fontweight='bold')
            plt.ylabel('SNR (dB-Hz)',fontsize=font_size,fontweight='bold')
            plt.xticks(rotation=90)
            plt.xticks(fontsize=font_size)
            plt.yticks(fontsize=font_size)
            
            plt.subplot(222)
            plt.plot(el_lis,SNRi_lis,'k')
            plt.title('(b) Reflected Signal ',fontsize=10,fontweight='bold')
            plt.xlabel('Elevation (deg.)',fontsize=font_size,fontweight='bold')
            plt.ylabel('SNR (Volts/volts)',fontsize=font_size,fontweight='bold')
            plt.xlim([5,30])                
            plt.xticks(fontsize=font_size)
            plt.yticks(fontsize=font_size)
            
            plt.subplot(223)
            plt.plot(frec2_sen,line_frec2_sen,'k')
            plt.title('(c) Lomb-Scargle Periodogram',fontsize=10,fontweight='bold')
            plt.xlabel('Reflector Height (m)',fontsize=font_size,fontweight='bold')
            plt.ylabel('Spectral Amplitude',fontsize=font_size,fontweight='bold')
            plt.xlim([0.000,2.5]) 
            plt.ylim([0.0,0.7])
            plt.xticks(fontsize=font_size)
            plt.yticks(fontsize=font_size)
 
            plt.subplot(224)
            x_range=numpy.linspace(5,30,100)
            plt.plot(el_lis,SNRi_lis,'k',label='Reflected Signal')
            plt.title('(d) Adjusted Wave',fontsize=10,fontweight='bold')
            plt.plot(x_range,(a0*numpy.cos(((4*pi*h0)/land)*numpy.sin(x_range*pi/180)+p0)),'grey', label='Adjusted Wave')
            plt.xlabel('Elevation (deg.)',fontsize=font_size,fontweight='bold')
            plt.ylabel('SNR (Volts/Volts)',fontsize=font_size,fontweight='bold')
            plt.xlim([5,30])
            plt.xticks(fontsize=font_size)
            plt.yticks(fontsize=font_size)
                        
            plt.tight_layout(pad=2.5, w_pad=2.5, h_pad=2.0)
                        
            fig_name=figPath+str(sat.strip())+'/'+str(nombre_dib)
            plt.savefig(fig_name+'.png', dpi=600)
            
            plt.close()
            
            # empty the lists for the next track
            frecuencia2=[]
            power2=[]

            # initialize max and medium values for the next track
            maximo=0
            media=0
            
    salida.close()
    salida2.close()
    matriz.close() 
    
