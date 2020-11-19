#读取SNR文件 画出每颗卫星的原始SNR值趋势 筛选合适的卫星和高度角
import matplotlib.pyplot as plt 
import numpy as np

save_path='C:/Data/RHfigure/'
archivo='C:/Data/SNR_SanXia/out/final_values.txt'#Input file
# archivo='C:/Data/out/final_values.txt'
# archivo='C:/Users/ShuaiYang/Desktop/final_values.txt'
# archivo='C:/Users/ShuaiYang/Desktop/test.txt'
# archivo='C:/Data/ZAOU201908.txt'

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
    datos=linea.split(' ')
    sat=int(datos[6])
    if sat not in Sat_list:
        Sat_list.append(sat) 
matriz2.close()
#  
for i in Sat_list:
    cont_2=0
    RH=[] 
    time=[]
    matriz2=open(archivo,'r') 
    while cont_2<num_lin_cab2:
        linea=matriz2.readline()
        cont_2 +=1
    while True:
        linea2=matriz2.readline()
        if not linea2:break
        datos2 = linea2.split(' ')
        sat1=int(datos2[6]) # id of the satellite
        if sat1==i:
            RH.append(float(datos2[15])) #存储数据类型要主要 数字 而不是str
            time.append(float(datos2[3]))
            # time.append(datos2[4])
        # rh.append(float(datos2[15]))  # 存储数据类型要主要 数字 而不是str
        # time.append(float(datos2[3]))
    matriz2.close()
    time=np.array(time)*60
    RH=np.array(RH)
    if len(time) > 2:
        plt.plot(time, RH)
        fig_name = save_path+str(i)
        plt.savefig(fig_name+'.png')
        plt.close()
