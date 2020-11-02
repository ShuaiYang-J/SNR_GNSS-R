# -*- coding: utf-8 -*-
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time
import re
from tqdm import tqdm
# import pandas as pd #对多个台站数据进行下载时需要

# 用户自定义需要下载数据的时间 长时间数据加for循环 日期选择要更改
Daystar=214 # 开始时间
Dayend=214 # 选择要下载的天数 截至日
Stationname=['MAYG'] #下载台站名称 多个台站时选择列表或xlx文件读取 pd包
Datatype='d' #下载类型 d是观测文件 n为单GPS的导航文件(_GN) p为混合(_MN)
# 设置需要下载测站的文件路径 多个台站时 选择文件进行存储名称
Strorpath='C:\\Data\\Down' #存储位置 双线\\
# 用户名和密码
username = "Yangshuai"
password = "Yang1997"

input=[]
for i in range(Daystar,Dayend+1,1):
    Day=str(i)
    input.append('2019'+'-'+Day.rjust(3,'0'))

# 需要ChromeDriver 谷歌真棒
# 定义驱动相关参数
options = webdriver.ChromeOptions()
## 禁用下载弹窗， 设置下载路径
prefs = {'profile.default_content_settings.popups': 0, 'download.default_directory': Strorpath}
options.add_experimental_option('prefs', prefs)
# 无窗口运行
options.headless = True

# 设置浏览器驱动
driver = webdriver.Chrome(chrome_options=options)
# 登录
def get_user(driver, username, password):
    driver.get("https://cddis.nasa.gov/archive/gnss/data/highrate/")
    time.sleep(2)
    # 模拟登录
    print("[Info] logging in ...")
    driver.find_element_by_id("username").send_keys(username)
    driver.find_element_by_id("password").send_keys(password)
    driver.find_element_by_name("commit").click()
    print("[Info] Finish logging!")
    time.sleep(5)

    return driver

def get_source_code(driver, input_num, hour):
    
    driver.get("https://cddis.nasa.gov/archive/gnss/data/highrate/")
    time.sleep(2)
    driver.find_element_by_id(input_num[:4]).click()
    time.sleep(2)
    driver.find_element_by_id(input_num[5:8]).click()
    time.sleep(2)
    driver.find_element_by_id(input_num[2:4] + Datatype).click()
    time.sleep(2)
    driver.find_element_by_id(hour).click() #高频数据 分时间段 00 - 23
    time.sleep(2)

    source = driver.page_source

    return source

def load_sta(source):
    # 通过正则表达式匹配网页源码中的站点数据压缩包名称
    pattern = re.compile('<a.*?id="(.*?)"\stitle="DataFile"')
    items = re.findall(pattern, source)
    # ----------------------------------
    # 读取需要下载数据的站点名称；仅读取第一列,即站名
    # df = pd.read_excel(bds_stas_dir, sheet_name='BDS', usecols=[0])
    # df.values为 numpy数组
    # df_sta = df.values.tolist()
    # bds_sta_names = []
    # for sta in df_sta:
        # bds_sta_names.append(sta[0])
    #--------------------------------------------------

    # 匹配一下需要下载的站点和网站中提供下载的站点，提取能够下载的站点数据压缩包名称
    download_stas = []
    for item in items:
        if item[:4] in Stationname:
            download_stas.append(item)

    return download_stas

if __name__ == '__main__':
    driver=get_user(driver, username, password)
    for i in range(len(input)):
        input_num=input[i]
        for j in range(24):
            hour=str(j).rjust(2,'0')     
            source = get_source_code(driver, input_num, hour)
            download_stas = load_sta(source)
            print("[Info] Downloading...")
            for id in tqdm(download_stas):
                driver.find_element_by_id(id).click()
                time.sleep(1)
        print("[Info] Finish!")
    driver.close()
