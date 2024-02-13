import os
import numpy as np
from app01.Data import DA
import pdb
import re
import time

def read_edi(filename, convert):
    "read file .edi"
    data = DA()
    f = open(filename,'r')
    line = f.readline()
    print (['Reading ',filename])
    location=np.zeros((3,1))
    ctime=np.zeros((2,3))

    while line:
        if line.isspace():
            line = f.readline()
            line = line.rstrip()
        if line.find('DATAID')>=0:
            ns=line.find('=')
            #global sitename
            sitename = eval(line[ns+1:])

        if line.find('ACQDATE') >= 0:
            ctemp = line.replace('ACQDATE', '')
            ctemp = ctemp.replace('=', '')
            ctemp = [float(k) for k in (ctemp.replace('/', ' ')).split(' ')]
            ctime[0, 0] = ctemp[0]
            ctime[0, 1] = ctemp[1]
            ctime[0, 2] = ctemp[2]

            # 下一行
            line = f.readline()
            ctemp = line.replace('FILEDATE', '')
            ctemp = ctemp.replace('=', '')
            ctemp = [float(k) for k in (ctemp.replace('/', ' ')).split(' ')]
            ctime[1, 0] = ctemp[0]
            ctime[1, 1] = ctemp[1]
            ctime[1, 2] = ctemp[2]

        if line.find('REFLAT')>=0:
            ctemp = line.replace('REFLAT','')
            ctemp = ctemp.replace('=','')
            ctemp = [float(k) for k in (ctemp.replace(':',' ')).split(' ')]
            if line.find('-')>=0:
                location[0] = ctemp[2]/3600+ctemp[1]/60+ctemp[0]
            else:
                location[0] = ctemp[2]/3600+ctemp[1]/60+ctemp[0]
        
            # 下一行
            line = f.readline()
            ctemp = line.replace('REFLONG','')
            ctemp = ctemp.replace('=','')
            ctemp = [float(k) for k in (ctemp.replace(':',' ')).split(' ')]
            if line.find('-')>=0:
                location[1] = ctemp[2]/3600+ctemp[1]/60+ctemp[0]
            else:
                location[1] = -(abs(ctemp[2])/3600+abs(ctemp[1])/60+abs(ctemp[0]))

            # 下一行
            line = f.readline()
            ctemp = line.replace('REFELEV','')
            ctemp = ctemp.replace('=','')
            ctemp = [float(k) for k in (ctemp.replace(':',' ')).split(' ')]
            location[2] = ctemp
        # elif line.find('SPECTRA')>=0:
        #     print('this edi might be a spectra edi')
        #     print('please convert it to impedance format first')
        elif line.find('>FREQ') >= 0:
            ns = line.find('//')+2
            ne = len(line)
            data.nfreq = int (line[ns:ne]) # 字符长度-1: 最后一个字符
            line = f.readline().rstrip()
            while line.find('>!****')<0:
                data.freq.extend(re.split(' -|  ',line)) # 两个空格

                str_index = f.tell()-1 # 标记所在行的位置
                line = f.readline().rstrip()
            data.freq = np.array(data.freq)
            data.freq = data.freq.astype(float) #字符型数组 转换为 float
            data.z = np.zeros((data.nfreq,18))
           
            print (' frequencies list found...')
            print(line[ns:ne-1],'frequencies in total')

            # 返回上一行
            f.seek(str_index)
            line = f.readline()

        

        # 未完待续
        elif line.find('>ZXXR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),0] = np.array(data.temp)
            print(' ZXXR found...')

            f.seek(str_index)
            line = f.readline()
            
        elif line.find('>ZXXI')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip()  
            data.z[0:len(data.temp),1] = np.array(data.temp)
            print(' ZXXI found...')

            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>ZXX.VAR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip()             
            data.z[0:len(data.temp),2] = np.array(data.temp)
            print(' ZVAR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>ZXYR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),3] = np.array(data.temp)
            print(' ZXYR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>ZXYI')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),4] = np.array(data.temp)
            print(' ZXYI found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>ZXY.VAR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),5] = np.array(data.temp)
            print(' ZXY.VAR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>ZYXR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),6] = np.array(data.temp)
            print(' ZYXR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>ZYXI')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),7] = np.array(data.temp)
            print(' ZYXR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>ZYX.VAR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),8] = np.array(data.temp)
            print(' ZYX.VAR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>ZYYR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),9] = np.array(data.temp)
            print(' ZYYR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>ZYYI')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),10] = np.array(data.temp)
            print(' ZYYI found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>ZYY.VAR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),11] = np.array(data.temp)
            print(' ZYY.VAR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>RHOXY')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.rho = np.zeros((data.nfreq,4))
            data.rho[0:len(data.temp),0] = np.array(data.temp)
            print(' RHOXY found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>RHOXY.ERR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.rho[0:len(data.temp),1] = np.array(data.temp)
            print(' RHOXY.ERR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>RHOYX')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip()
            data.rho[0:len(data.temp),2] = np.array(data.temp)
            print(' RHOYX found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>RHOYX.ERR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.rho[0:len(data.temp),3] = np.array(data.temp)
            print(' RHOYX.ERR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>PHSXY')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.phs = np.zeros((data.nfreq,4))
            data.phs[0:len(data.temp),0] = np.array(data.temp)
            print(' PHSXY found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>PHSXY.ERR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.phs[0:len(data.temp),1] = np.array(data.temp)
            print(' PHSXY.ERR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>PHSYX')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.phs[0:len(data.temp),2] = np.array(data.temp)
            print(' PHSYX found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>PHSYX.ERR')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.phs[0:len(data.temp),3] = np.array(data.temp)
            print(' RHOXY.ERR found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()

        # adding tipper here
        elif line.find('>TXR.EXP')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),12] = np.array(data.temp)
            print(' TXR.EXP found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>TXI.EXP')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),13] = np.array(data.temp)
            print(' TXI.EXP found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>TXVAR.EXP')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),14] = np.array(data.temp)
            print(' TXVAR.EXP found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>TYR.EXP')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            data.z[0:len(data.temp),15] = np.array(data.temp)
            print(' TYR.EXP found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>TYI.EXP')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0:
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            # print(data.temp)
            data.z[0:len(data.temp),16] = np.array(data.temp)
            print(' TYI.EXP found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()
        elif line.find('>TYVAR.EXP')>=0:
            data.temp=[]
            line = f.readline().rstrip()
            while line.find('>')<0 and line.strip()!='':
                data.temp.extend(re.split(' -|  ',line)) # 两个空格
                str_index = f.tell()-1
                line = f.readline().rstrip() 
            # print(data.temp)
            data.z[0:len(data.temp),17] = np.array(data.temp)
            print(' TYVAR.EXP found...')
            # 返回上一行
            f.seek(str_index)
            line = f.readline()


        line = f.readline()

    f.close()
    return data, location, sitename, ctime