#读取SNR文件 画出每颗卫星的原始SNR值趋势 筛选合适的卫星和高度角
import matplotlib.pyplot as plt 
import numpy as np

save_path='C:/Data/SNRfigure/'
archivo='C:/Data/SNR_SanXia/SNR-GPS.txt'#Input file

# archivo='C:/Users/ShuaiYang/Desktop/SNR/SNRpy/data/input/SNR-GPS.txt'


cont_2=0
num_lin_cab2 = 1
matriz2=open(archivo,'r') 
Sat_list=[]

# 先筛选存在的卫星编号
while cont_2<num_lin_cab2:
    linea=matriz2.readline()
    cont_2 +=1
while True:
    linea=matriz2.readline()
    if not linea:break
    datos=linea.split(';')
    sat=int(datos[2])
    if sat not in Sat_list:
        Sat_list.append(sat) 
matriz2.close()
#  
for i in Sat_list:
    cont_2=0
    SNR=[] 
    ele=[]
    az=[]
    matriz2=open(archivo,'r') 
    while cont_2<num_lin_cab2:
        linea=matriz2.readline()
        cont_2 +=1
    while True:
        linea2=matriz2.readline()
        if not linea2:break
        datos2 = linea2.split(';')
        sat1=int(datos2[2]) # id of the satellite
        if sat1==i:
            if float(datos2[4])>0:
                SNR.append(float(datos2[3])) #存储数据类型要主要 数字 而不是str
                ele.append(float(datos2[4]))
                az.append(float(datos2[5])) 
            # time.append(datos2[4])
        # rh.append(float(datos2[15]))  # 存储数据类型要主要 数字 而不是str
        # time.append(float(datos2[3]))
    matriz2.close()
    ele=np.array(ele)
    SNR=np.array(SNR)
    # num=np.argmax(ele)
    # index = np.argsort(time) #返回排序后的索引号
    # time = time[index]
    # rh = rh[index]
    if len(ele)!=0:
        plt.subplot(311)
        plt.plot(SNR)
        plt.subplot(312)
        plt.plot(ele)
        plt.subplot(313)
        plt.plot(az)
        fig_name=save_path+str(i)
        plt.savefig(fig_name+'.png')
        plt.close()
