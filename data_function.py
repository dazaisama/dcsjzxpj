import os
from urllib import request

import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib

matplotlib.use("Agg")

def median(s):
    n=len(s)                   #计算列表内元素数量
    if n==1:                   #这个要非常注意，当元素只有一个的时候，直接取值
        return s[0]         
    elif n%2!=0:               #如果元素数量为奇数
        m=sorted(s)            #排序一下
        mid=m[(n-1)//2]        #中间值等于元素总数量减一以后除以2，记得要用//
        return mid
    else:
        m=sorted(s)               
        mid=float(m[n//2-1]+m[n//2])/2
        return mid

def shujuchulir(data,k=5):
    # a simple script to plot apparent resistivities and phases for impedance
    
    #  usage: plot_Z(data,sitename,opt,prt)
    #  data: data structure of station
    #  opt: xxyy, xyyx or txty
    #  prt: options to print the plot made
    #       could be on, off, silent or overlap
    freq = data.freq
#    zxxr = data.z[:,0]
#    zxxi = -data.z[:,1]
#    zxxvar = data.z[:,2]
    zxyr=data.z[:,3]
    zxyi=data.z[:,4]
#    zxyvar=data.z[:,5]
    zyxr=data.z[:,6]
    zyxi=data.z[:,7]
#    zyxvar=data.z[:,8]
#    zyyr=data.z[:,9]
#    zyyi=data.z[:,10]
#    zyyvar=data.z[:,11]
        
    rhoxy = (zxyr**2 + zxyi**2)/freq/5
#    phsxy = np.arctan2(zxyi, zxyr) *180/math.pi
#    rhoxye = 2*zxyvar/((zxyr**2+zxyi**2)**0.5)
#    phsxye = zxyvar/((zxyr**2+zxyi**2)**0.5)
#    rhoxy = np.log10(rhoxy)
#    phsxye = phsxye*180/math.pi
    rhoyx = (zyxr**2 + zyxi**2)/freq/5
#    phsyx = np.arctan2(zyxi, zyxr) *180/math.pi
#    rhoyxe = 2*zyxvar/((zyxr**2+zyxi**2)**0.5)
#    phsyxe = zyxvar/((zyxr**2+zyxi**2)**0.5)
#    rhoyx = np.log10(rhoyx)
#    phsyxe = phsyxe*180/math.pi
    L=len(rhoxy)
    A=np.zeros([2*L,15])
    A[0:L,0]=range(1,L+1,1)
    A[L:2*L,0]=range(1,L+1,1)
    A[0:L,1]=freq
    A[L:L*2,1]=freq
    A[0:L,2]=rhoxy
    A[L:L*2,2]=rhoyx
    for i in range(0,L*2):
        j=i//L;
        m=median(A[L*j+1:L*(j+1),2]);
        x1 = (np.sort(A[L*j:L*(j+1),2])[-k]+np.sort(A[L*j:L*(j+1),2])[k])/2;
        #x1=sum(A[0:72,2])/L;
        #std1=np.std(A[L*j:L*(j+1),2]);
        #std1=np.std(A[0:72,2]);        
        A[i,7]=A[i,2]/x1; 
        if i%L==0:
            A[i,4]=1
            A[i,5]=0
            A[i,6]=1
        elif i%L==1:
            A[i,4]=((A[i-1,2]+A[i+1,2])/2)/A[i,2]
            A[i,5]=abs(A[i-1,2]-A[i,2])
            A[i,6]=1
        elif i%L==L-2:
            A[i,4]=((A[i-1,2]+A[i+1,2])/2)/A[i,2]
            A[i,5]=abs(A[i-1,2]-A[i,2])
            A[i,6]=1
        elif i%L==L-1:
            A[i,4]=1
            A[i,5]=abs(A[i-1,2]-A[i,2])
            A[i,6]=1
        else:
            A[i,4]=((A[i-1,2]+A[i+1,2])/2)/A[i,2]
            A[i,5]=abs(A[i-1,2]-A[i,2])
            A[i,6]=((A[i-2,2]+A[i+2,2])/2)/A[i,2]        
        A[i,11]=A[i,2]-x1;
        A[i,13]=A[i,2]-m;
        A[i,14]=A[i,2]/m;

    for i in range(0,L*2):
        j=i//L;        
        #std2=np.std(A[L*j:L*(j+1),4]);
        #std2=np.std(A[0:72,4]);
        #std3=np.std(A[L*j:L*(j+1),5]);
        #std4=np.std(A[L*j:L*(j+1),6]);
        #std5=np.std(A[L*j:L*(j+1),11]);
        x2 = (np.sort(A[L*j:L*(j+1),4])[-k]+np.sort(A[L*j:L*(j+1),4])[k])/2;
        #x2=sum(A[L*j:L*(j+1),4])/L;        
        #x2=sum(A[0:72,4])/72;
        x3 = (np.sort(A[L*j:L*(j+1),5])[-k]+np.sort(A[L*j:L*(j+1),5])[k])/2;
        #x3=sum(A[L*j:L*(j+1),5])/L;
        x4 = (np.sort(A[L*j:L*(j+1),6])[-k]+np.sort(A[L*j:L*(j+1),6])[k])/2;
        #x4=sum(A[L*j:L*(j+1),6])/L;
        x5 = (np.sort(A[L*j:L*(j+1),11])[-k]+np.sort(A[L*j:L*(j+1),11])[k])/2;
        x5=sum(A[L*j:L*(j+1),11])/L;
        A[i,8]=A[i,4]/x2;
        A[i,9]=A[i,5]/x3;
        A[i,10]=A[i,6]/x4;        
        A[i,12]=A[i,11]/x5;


    data=pd.DataFrame(A,columns=['','','','y','x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11'])
    return A, data
    
    
def shujuchulip(data,k=5):
        # a simple script to plot apparent resistivities and phases for impedance
        
        #  usage: plot_Z(data,sitename,opt,prt)
        #  data: data structure of station 
        #  opt: xxyy, xyyx or txty
        #  prt: options to print the plot made 
        #       could be on, off, silent or overlap
    freq = data.freq
#    zxxr = data.z[:,0]
#    zxxi = -data.z[:,1]
#    zxxvar = data.z[:,2]
    zxyr=data.z[:,3]
    zxyi=data.z[:,4]
#    zxyvar=data.z[:,5]
    zyxr=data.z[:,6]
    zyxi=data.z[:,7]
#    zyxvar=data.z[:,8]
#    zyyr=data.z[:,9]
#    zyyi=data.z[:,10]
#    zyyvar=data.z[:,11]
    
    #    rhoxy = (zxyr**2 + zxyi**2)/freq/5
    phsxy = np.arctan2(zxyi, zxyr) *180/math.pi
    #    rhoxye = 2*zxyvar/((zxyr**2+zxyi**2)**0.5)
    #    phsxye = zxyvar/((zxyr**2+zxyi**2)**0.5)
    #    rhoxy = np.log10(rhoxy)
    #    phsxye = phsxye*180/math.pi
    #    rhoyx = (zyxr**2 + zyxi**2)/freq/5
    phsyx = np.arctan2(zyxi, zyxr) *180/math.pi
    #    rhoyxe = 2*zyxvar/((zyxr**2+zyxi**2)**0.5)
    #    phsyxe = zyxvar/((zyxr**2+zyxi**2)**0.5)
    #    rhoyx = np.log10(rhoyx)
    #    phsyxe = phsyxe*180/math.pi
    L=len(phsxy)
    A=np.zeros([2*L,15])
    A[0:L,0]=range(1,L+1,1)
    A[L:2*L,0]=range(1,L+1,1)
    A[0:L,1]=freq
    A[L:L*2,1]=freq
    A[0:L,2]=phsxy
    A[L:L*2,2]=phsyx
    for i in range(0,L*2):
        j=i//L;
        m=np.median(A[L*j:L*(j+1),2]);
        x1 = (np.sort(A[L*j:L*(j+1),2])[-k]+np.sort(A[L*j:L*(j+1),2])[k])/2;        
        #x1=sum(A[0:72,2])/L;
        #std1=np.std(A[L*j:L*(j+1),2]);
        #std1=np.std(A[0:72,2]);        
        A[i,7]=A[i,2]/x1; 
        if i%L==0:
            A[i,4]=1
            A[i,5]=0
            A[i,6]=1
        elif i%L==1:
            A[i,4]=((A[i-1,2]+A[i+1,2])/2)/A[i,2]
            A[i,5]=abs(A[i-1,2]-A[i,2])
            A[i,6]=1
        elif i%L==L-2:
            A[i,4]=((A[i-1,2]+A[i+1,2])/2)/A[i,2]
            A[i,5]=abs(A[i-1,2]-A[i,2])
            A[i,6]=1
        elif i%L==L-1:
            A[i,4]=1
            A[i,5]=abs(A[i-1,2]-A[i,2])
            A[i,6]=1
        else:
            A[i,4]=((A[i-1,2]+A[i+1,2])/2)/A[i,2]
            A[i,5]=abs(A[i-1,2]-A[i,2])
            A[i,6]=((A[i-2,2]+A[i+2,2])/2)/A[i,2]        
        A[i,11]=A[i,2]-x1;
        A[i,13]=A[i,2]-m;
        A[i,14]=A[i,2]/m;

    for i in range(0,L*2):
        j=i//L;        
        #std2=np.std(A[L*j:L*(j+1),4]);
        #std2=np.std(A[0:72,4]);
        #std3=np.std(A[L*j:L*(j+1),5]);
        #std4=np.std(A[L*j:L*(j+1),6]);
        #std5=np.std(A[L*j:L*(j+1),11]);
        x2 = (np.sort(A[L*j:L*(j+1),4])[-k]+np.sort(A[L*j:L*(j+1),4])[k])/2;
        #x2=sum(A[L*j:L*(j+1),4])/L;        
        #x2=sum(A[0:72,4])/72;
        x3 = (np.sort(A[L*j:L*(j+1),5])[-k]+np.sort(A[L*j:L*(j+1),5])[k])/2;
        #x3=sum(A[L*j:L*(j+1),5])/L;
        x4 = (np.sort(A[L*j:L*(j+1),6])[-k]+np.sort(A[L*j:L*(j+1),6])[k])/2;
        #x4=sum(A[L*j:L*(j+1),6])/L;
        x5 = (np.sort(A[L*j:L*(j+1),11])[-k]+np.sort(A[L*j:L*(j+1),11])[k])/2;
        x5=sum(A[L*j:L*(j+1),11])/L;
        A[i,8]=A[i,4]/x2;
        A[i,9]=A[i,5]/x3;
        A[i,10]=A[i,6]/x4;        
        A[i,12]=A[i,11]/x5; 

                
    data=pd.DataFrame(A,columns=['','','','y','x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11'])
    return A, data

def file_name(file_dir):
    for root,dirs,files in os.walk(file_dir):
        print(files)
    cfilename = []
    for fname in files:
        cfilename.append(root + '/' + fname)
        
    return cfilename


def plot_Z(data, sitename, opt, prt,id):
    # a simple script to plot apparent resistivities and phases for impedance
    
    #  usage: plot_Z(data,sitename,opt,prt)
    #  data: data structure of station 
    #  opt: xxyy, xyyx or txty
    #  prt: options to print the plot made 
    #       could be on, off, silent or overlap
    freq = data.freq
    zxxr = data.z[:,0]
    zxxi = -data.z[:,1]
    zxxvar = data.z[:,2]
    zxyr=data.z[:,3]
    zxyi=data.z[:,4]
    zxyvar=data.z[:,5]
    zyxr=data.z[:,6]
    zyxi=data.z[:,7]
    zyxvar=data.z[:,8]
    zyyr=data.z[:,9]
    zyyi=data.z[:,10]
    zyyvar=data.z[:,11]


    # 覆盖现有画布
#    if prt == 'overlap':
#        plt.gcf
#        plt.figure(figsize=(10,6))
#    else:
#        # 新建画布
#        plt.figure(num = sitename,figsize=(10,6))


    if opt == 'xxyy':
        # %=========calculate rho and phase from data=======%
        # %please note that the errorbars of resistivities and phases here are
        # %calculated using Gary Egbert's formula in Z files.   
        rhoxx = (zxxr**2 + zxxi**2)/freq/5
        phsxx = np.arctan2(zxxi, zxxr) *180/math.pi
        rhoxxe = 2*zxxvar/((zxxr**2+zxxi**2)**0.5)
        phsxxe = zxxvar/((zxxr**2+zxxi**2)**0.5)
        rhoxx = np.log10(rhoxx)
        phsxxe = phsxxe*180/math.pi
        rhoyy = (zyyr**2 + zyyi**2)/freq/5
        phsyy = np.arctan2(zyyi, zyyr) *180/math.pi
        rhoyye = 2*zyyvar/((zyyr**2+zyyi**2)**0.5)
        phsyye = zyyvar/((zyyr**2+zyyi**2)**0.5)
        rhoyy = np.log10(rhoyy)
        phsyye = phsyye*180/math.pi

        minrho = min(min(rhoxx),min(rhoyy))
        minrho = math.floor(minrho - 0.5)
        maxrho = max(max(rhoxx),max(rhoyy))
        maxrho = round(maxrho + 1)

        # 开始画图
        ax1 = plt.subplot(2,1,1)
        ax1.errorbar(freq,rhoxx,yerr=rhoxxe,fmt='xr',ecolor = 'r',capsize = 0)
        ax1.errorbar(freq,rhoyy,yerr=rhoyye,fmt='+c',ecolor = 'c',capsize = 0)
        ax1.set_ylim(minrho,maxrho)


        ax2 = plt.subplot(2,1,2)
        ax2.errorbar(freq,phsxx,yerr=phsxxe,fmt='xr',ecolor = 'r',capsize = 0)
        ax2.errorbar(freq,phsyy,yerr=phsyye,fmt='+c',ecolor = 'c',capsize = 0)

        #修改之前 ax2.set(ylim=(-180,180), yticks=[-180,-135,-90,-45,0,45,90,135,180], size=12)

        ax2.set(ylim=(-180, 180), yticks=[-180, -135, -90, -45, 0, 45, 90, 135, 180])
        ax2.tick_params(axis='both', labelsize=12)

        #修改之前 ax2.set_yticks([-180,-135,-90,-45,0,45,90,135,180], size=12)

        ax2.set_yticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
        ax2.tick_params(axis='y', labelsize=12)


    elif opt == 'xyyx':
        # %=========calculate rho and phase from data=======%
        rhoxy = (zxyr**2 + zxyi**2)/freq/5
        phsxy = np.arctan2(zxyi, zxyr) *180/math.pi
        rhoxye = 2*zxyvar/((zxyr**2+zxyi**2)**0.5)
        phsxye = zxyvar/((zxyr**2+zxyi**2)**0.5)
        rhoxy = np.log10(rhoxy)
        phsxye = phsxye*180/math.pi
        rhoyx = (zyxr**2 + zyxi**2)/freq/5
        phsyx = np.arctan2(zyxi, zyxr) *180/math.pi
        rhoyxe = 2*zyxvar/((zyxr**2+zyxi**2)**0.5)
        phsyxe = zyxvar/((zyxr**2+zyxi**2)**0.5)
        rhoyx = np.log10(rhoyx)
        phsyxe = phsyxe*180/math.pi
        
        minrho = min(min(rhoyx),min(rhoxy))
        minrho = math.floor(minrho - 0.5)
        maxrho = max(max(rhoyx),max(rhoxy))
        maxrho = round(maxrho + 1)

        # 开始画图
        fig,ax1 = plt.subplots(figsize=(5, 3))
        ax1.errorbar(freq, rhoxy, yerr=rhoxye, fmt='xr', ecolor='r', capsize=0)
        ax1.errorbar(freq, rhoyx, yerr=rhoyxe, fmt='.c', ecolor='c', capsize=0)
        ax1.set_ylim(0, 3)

        ax1.grid(linestyle='-.', which='both', axis='both', color='#e8e8e8')  # 浅灰色'#e8e8e8'
        ax1.set(xlim=[min(freq) / 3, max(freq) * 3], xscale='log', ylabel=r'$log_{10}$' + 'app. resistivity (Ohm*m)')


        ax1.invert_xaxis()
        if prt == 'on' or prt == 'silent':
            # plt.savefig(r"kayabezprint\\"+eval(sitename)+ opt +'.eps',dpi=600,format='eps')
            # plt.savefig(r"kayabezprint\\"+eval(sitename)+ opt +'.png',dpi=600,format='png')
#            plt.savefig(r"E:\大创\汇报\123\\"+eval(sitename)+ opt +'or'+'.eps',dpi=200,format='eps')
#            plt.savefig(r"E:\大创\汇报\123\\"+eval(sitename)+ opt +'or'+'.png',dpi=200,format='png')
#        if prt != 'silent':
            # print("))00")
            # plt.savefig(r"media/outputPictures/" + sitename + id + opt + 'or' + '.eps', dpi=200, format='eps')
            plt.savefig(r"media/outputPictures/" + id + sitename + opt + 'or' + '.png', dpi=200, format='png')
            plt.savefig(r"app01/static/png/" + id + sitename + opt + 'or' + '.png', dpi=200, format='png')
        # if prt != 'silent':


        fig, ax2 = plt.subplots(figsize=(5, 3))
        ax2.errorbar(freq, phsxy, yerr=phsxye, fmt='xr', ecolor='r', capsize=0)
        ax2.errorbar(freq, phsyx, yerr=phsyxe, fmt='.c', ecolor='c', capsize=0)
        ax2.set_ylim(0, 90)
        ax2.set_yticks([0, 15, 30, 45, 60, 75, 90])

        ax2.grid(linestyle='-.', which='both', axis='both', color='#e8e8e8')
        ax2.set(xlim=[min(freq) / 3, max(freq) * 3], xscale='log', xlabel='frequency (Hz)', ylabel='phase(degree)')
        ax2.invert_xaxis()
        if prt == 'on' or prt == 'silent':
            # plt.savefig(r"kayabezprint\\"+eval(sitename)+ opt +'.eps',dpi=600,format='eps')
            # plt.savefig(r"kayabezprint\\"+eval(sitename)+ opt +'.png',dpi=600,format='png')
#            plt.savefig(r"E:\大创\汇报\123\\"+eval(sitename)+ opt +'op'+'.eps',dpi=200,format='eps')
#            plt.savefig(r"E:\大创\汇报\123\\"+eval(sitename)+ opt +'op'+'.png',dpi=200,format='png')
#        if prt != 'silent':
#            plt.show()
#             plt.savefig(r"media/outputPictures/" + sitename + id + opt + 'op' + '.eps', dpi=200, format='eps')
            plt.savefig(r"media/outputPictures/" + id + sitename  + opt + 'op' + '.png', dpi=200, format='png')
            plt.savefig(r"app01/static/png/" + id + sitename + opt + 'op' + '.png', dpi=200, format='png')
        plt.close()
        # if prt != 'silent':
        #     plt.show()
        Filename1 =id + sitename + opt + 'or' + '.png'
        Filename2 =id + sitename + opt + 'op' + '.png'
        return Filename1,Filename2

def plot_R(data, predict, sitename,id):

        freq = data.freq
        zxxr = data.z[:, 0]
        zxxi = -data.z[:, 1]
        zxxvar = data.z[:, 2]
        zxyr = data.z[:, 3]
        zxyi = data.z[:, 4]
        zxyvar = data.z[:, 5]
        zyxr = data.z[:, 6]
        zyxi = data.z[:, 7]
        zyxvar = data.z[:, 8]
        zyyr = data.z[:, 9]
        zyyi = data.z[:, 10]
        zyyvar = data.z[:, 11]
        rhoxy = (zxyr ** 2 + zxyi ** 2) / freq / 5
        phsxy = np.arctan2(zxyi, zxyr) * 180 / math.pi
        rhoxye = 2 * zxyvar / ((zxyr ** 2 + zxyi ** 2) ** 0.5)
        phsxye = zxyvar / ((zxyr ** 2 + zxyi ** 2) ** 0.5)
        phsxye = phsxye * 180 / math.pi
        rhoyx = (zyxr ** 2 + zyxi ** 2) / freq / 5
        phsyx = np.arctan2(zyxi, zyxr) * 180 / math.pi
        rhoyxe = 2 * zyxvar / ((zyxr ** 2 + zyxi ** 2) ** 0.5)
        phsyxe = zyxvar / ((zyxr ** 2 + zyxi ** 2) ** 0.5)
        phsyxe = phsyxe * 180 / math.pi
        rhoyx = np.log10(rhoyx)
        rhoxy = np.log10(rhoxy)

        d = 0;
        e = 0;
        f = 0;
        L = int(len(predict) / 2)

        freq12 = np.zeros(L);
        freq10 = np.zeros(L)
        rhoxy12 = np.zeros(L);
        rhoxy10 = np.zeros(L);
        rhoxye12 = np.zeros(L);
        rhoxye10 = np.zeros(L);
        for c in range(0, L):
            if predict[c] == 2:
                freq12[d] = freq[c];
                rhoxy12[d] = rhoxy[c];
                rhoxye12[d] = rhoxye[c];
                d = d + 1;
            else:
                freq10[f] = freq[c];
                rhoxy10[f] = rhoxy[c];
                rhoxye10[f] = rhoxye[c];
                f = f + 1;
        freq12.resize(d)
        freq10.resize(f)
        rhoxy12.resize(d)
        rhoxy10.resize(f)
        rhoxye12.resize(d)
        rhoxye10.resize(f)

        d1 = 0;
        e1 = 0;
        f1 = 0;
        freq22 = np.zeros(L);
        freq20 = np.zeros(L);
        rhoyx22 = np.zeros(L);
        rhoyx20 = np.zeros(L);
        rhoyxe22 = np.zeros(L);
        rhoyxe20 = np.zeros(L);
        for c in range(L, 2 * L):
            if predict[c] == 2:
                freq22[d1] = freq[c - L];
                rhoyx22[d1] = rhoyx[c - L];
                rhoyxe22[d1] = rhoyxe[c - L];
                d1 = d1 + 1;
            else:
                freq20[f1] = freq[c - L];
                rhoyx20[f1] = rhoyx[c - L];
                rhoyxe20[f1] = rhoyxe[c - L];
                f1 = f1 + 1;

        freq22.resize(d1)
        freq20.resize(f1)
        rhoyx22.resize(d1)
        rhoyx20.resize(f1)
        rhoyxe22.resize(d1)
        rhoyxe20.resize(f1)
        # print(freq20)

        fig, ax1 = plt.subplots(figsize=(5, 3))

        ax1.errorbar(freq12, rhoxy12, yerr=rhoxye12, fmt='xr', ecolor='w', capsize=0)
        ax1.errorbar(freq10, rhoxy10, yerr=rhoxye10, fmt='xk', ecolor='w', capsize=0)
        ax1.errorbar(freq22, rhoyx22, yerr=rhoyxe22, fmt='.c', ecolor='w', capsize=0)
        ax1.errorbar(freq20, rhoyx20, yerr=rhoyxe20, fmt='.k', ecolor='w', capsize=0)
        ax1.set_ylim(0, 3)

        #        ax2 = plt.subplot(2,1,2)
        #        ax2.errorbar(freq,phsxy,yerr=phsxye,fmt='xr',ecolor = 'r',capsize = 0)
        #        ax2.errorbar(freq,phsyx,yerr=phsyxe,fmt='.c',ecolor = 'c',capsize = 0)
        #        ax2.set_ylim(0,90)
        #        ax2.set_yticks([0, 15, 30, 45, 60, 75, 90])

        ax1.grid(linestyle='-.', which='both', axis='both', color='#e8e8e8')  # 浅灰色'#e8e8e8'
        #    ax2.grid(linestyle='-.',which = 'both',axis='both',color='#e8e8e8')
        ax1.set(xlim=[min(freq) / 3, max(freq) * 3], xscale='log', ylabel=r'$log_{10}$' + 'app. resistivity (Ohm*m)')
        #    ax2.set(xlim=[min(freq)/3,max(freq)*3],xscale = 'log',xlabel='frequency (Hz)', ylabel = 'phase(degree)')
        ax1.invert_xaxis()
        #    ax2.invert_xaxis()

        # plt.savefig(r"kayabezprint\\"+eval(sitename)+ opt +'.eps',dpi=600,format='eps')
        # plt.savefig(r"kayabezprint\\"+eval(sitename)+ opt +'.png',dpi=600,format='png')
        # plt.savefig(r"media/outputPictures/" + sitename + id + 'r' + '.eps', dpi=200, format='eps')
        plt.savefig(r"media/outputPictures/" + id + sitename + 'r' + '.png', dpi=200, format='png')
        plt.savefig(r"app01/static/png/" + id + sitename + 'r' + '.png', dpi=200, format='png')
        plt.close()
#        plt.savefig(r"E:\大创\汇报\123\\" +'r'+'.eps',dpi=200,format='eps')
#        plt.savefig(r"E:\大创\汇报\123\\" +'r'+'.png',dpi=200,format='png')
        Filename =id + sitename + 'r' + '.png'
        return Filename

def plot_P(data, predict, sitename,id):
    freq = data.freq
    zxxr = data.z[:, 0]
    zxxi = -data.z[:, 1]
    zxxvar = data.z[:, 2]
    zxyr = data.z[:, 3]
    zxyi = data.z[:, 4]
    zxyvar = data.z[:, 5]
    zyxr = data.z[:, 6]
    zyxi = data.z[:, 7]
    zyxvar = data.z[:, 8]
    zyyr = data.z[:, 9]
    zyyi = data.z[:, 10]
    zyyvar = data.z[:, 11]
    rhoxy = (zxyr ** 2 + zxyi ** 2) / freq / 5
    phsxy = np.arctan2(zxyi, zxyr) * 180 / math.pi
    rhoxye = 2 * zxyvar / ((zxyr ** 2 + zxyi ** 2) ** 0.5)
    phsxye = zxyvar / ((zxyr ** 2 + zxyi ** 2) ** 0.5)
    phsxye = phsxye * 180 / math.pi
    rhoyx = (zyxr ** 2 + zyxi ** 2) / freq / 5
    phsyx = np.arctan2(zyxi, zyxr) * 180 / math.pi
    rhoyxe = 2 * zyxvar / ((zyxr ** 2 + zyxi ** 2) ** 0.5)
    phsyxe = zyxvar / ((zyxr ** 2 + zyxi ** 2) ** 0.5)
    phsyxe = phsyxe * 180 / math.pi
    rhoyx = np.log10(rhoyx)
    rhoxy = np.log10(rhoxy)

    d = 0;
    f = 0;
    L = int(len(predict) / 2)

    freq12 = np.zeros(L);
    freq10 = np.zeros(L)
    phsxy12 = np.zeros(L);
    phsxy10 = np.zeros(L);
    phsxye12 = np.zeros(L);
    phsxye10 = np.zeros(L);
    for c in range(0, L):
        if predict[c] == 2:
            freq12[d] = freq[c];
            phsxy12[d] = phsxy[c];
            phsxye12[d] = phsxye[c];
            d = d + 1;
        else:
            freq10[f] = freq[c];
            phsxy10[f] = phsxy[c];
            phsxye10[f] = phsxye[c];
            f = f + 1;
    freq12.resize(d)
    freq10.resize(f)
    phsxy12.resize(d)
    phsxy10.resize(f)
    phsxye12.resize(d)
    phsxye10.resize(f)

    d1 = 0;
    f1 = 0;
    freq22 = np.zeros(L);
    freq20 = np.zeros(L);
    phsyx22 = np.zeros(L);
    phsyx20 = np.zeros(L);
    phsyxe22 = np.zeros(L);
    phsyxe20 = np.zeros(L);
    for c in range(L, 2 * L):
        if predict[c] == 2:
            freq22[d1] = freq[c - L];
            phsyx22[d1] = phsyx[c - L];
            phsyxe22[d1] = phsyxe[c - L];
            d1 = d1 + 1;
        else:
            freq20[f1] = freq[c - L];
            phsyx20[f1] = phsyx[c - L];
            phsyxe20[f1] = phsyxe[c - L];
            f1 = f1 + 1;

    freq22.resize(d1)
    freq20.resize(f1)
    phsyx22.resize(d1)
    phsyx20.resize(f1)
    phsyxe22.resize(d1)
    phsyxe20.resize(f1)
    # print(freq20)

    # 开始画图
    fig, ax1 = plt.subplots(figsize=(5, 3))
    ax1.errorbar(freq12, phsxy12, yerr=phsxye12, fmt='xr', ecolor='w', capsize=0)
    ax1.errorbar(freq10, phsxy10, yerr=phsxye10, fmt='xk', ecolor='w', capsize=0)
    ax1.errorbar(freq22, phsyx22, yerr=phsyxe22, fmt='.c', ecolor='w', capsize=0)
    ax1.errorbar(freq20, phsyx20, yerr=phsyxe20, fmt='.k', ecolor='w', capsize=0)
    ax1.set_ylim(0, 90)
    ax1.set_yticks([0, 15, 30, 45, 60, 75, 90])

    #        ax2 = plt.subplot(2,1,2)
    #        ax2.errorbar(freq,phsxy,yerr=phsxye,fmt='xr',ecolor = 'r',capsize = 0)
    #        ax2.errorbar(freq,phsyx,yerr=phsyxe,fmt='.c',ecolor = 'c',capsize = 0)
    #        ax2.set_ylim(0,90)
    #        ax2.set_yticks([0, 15, 30, 45, 60, 75, 90])

    ax1.grid(linestyle='-.', which='both', axis='both', color='#e8e8e8')  # 浅灰色'#e8e8e8'
    #    ax2.grid(linestyle='-.',which = 'both',axis='both',color='#e8e8e8')
    ax1.set(xlim=[min(freq) / 3, max(freq) * 3], xscale='log', xlabel='frequency (Hz)', ylabel='phase(degree)')
    #    ax2.set(xlim=[min(freq)/3,max(freq)*3],xscale = 'log',xlabel='frequency (Hz)', ylabel = 'phase(degree)')
    ax1.invert_xaxis()
    #    ax2.invert_xaxis()

    # plt.savefig(r"media/outputPictures/" + sitename + id + 'p' + '.eps', dpi=200, format='eps')
    plt.savefig(r"media/outputPictures/" + id + sitename + 'p' + '.png', dpi=200, format='png')
    plt.savefig(r"app01/static/png/" + id + sitename + 'p' + '.png', dpi=200, format='png')
    plt.close()
#    plt.savefig(r"E:\大创\汇报\123\\" +'p'+'.eps',dpi=200,format='eps')
#    plt.savefig(r"E:\大创\汇报\123\\" +'p'+'.png',dpi=200,format='png') 
    Filename =id + sitename + 'p' + '.png'
    return Filename


def plot_ZZ(data1,data2, sitename, opt, prt):
    # a simple script to plot apparent resistivities and phases for impedance
    
    #  usage: plot_Z(data,sitename,opt,prt)
    #  data: data structure of station 
    #  opt: xxyy, xyyx or txty
    #  prt: options to print the plot made 
    #       could be on, off, silent or overlap
    freq1 = data1.freq
    zxxr1 = data1.z[:,0]
    zxxi1 = -data1.z[:,1]
    zxxvar1 = data1.z[:,2]
    zxyr1=data1.z[:,3]
    zxyi1=data1.z[:,4]
    zxyvar1=data1.z[:,5]
    zyxr1=data1.z[:,6]
    zyxi1=data1.z[:,7]
    zyxvar1=data1.z[:,8]
    zyyr1=data1.z[:,9]
    zyyi1=data1.z[:,10]
    zyyvar1=data1.z[:,11]
    
    freq2 = data2.freq
    zxxr2 = data2.z[:,0]
    zxxi2 = -data2.z[:,1]
    zxxvar2 = data2.z[:,2]
    zxyr2=data2.z[:,3]
    zxyi2=data2.z[:,4]
    zxyvar2=data2.z[:,5]
    zyxr2=data2.z[:,6]
    zyxi2=data2.z[:,7]
    zyxvar2=data2.z[:,8]
    zyyr2=data2.z[:,9]
    zyyi2=data2.z[:,10]
    zyyvar2=data2.z[:,11]


    # 覆盖现有画布
    if prt == 'overlap':
        plt.gcf
        plt.figure(figsize=(10,6))
    else:
        # 新建画布
        plt.figure(num = sitename,figsize=(10,6))


    # if opt == 'xxyy':
    #     # %=========calculate rho and phase from data=======%
    #     # %please note that the errorbars of resistivities and phases here are
    #     # %calculated using Gary Egbert's formula in Z files.   
    #     rhoxx = (zxxr**2 + zxxi**2)/freq/5
    #     phsxx = np.arctan2(zxxi, zxxr) *180/math.pi
    #     rhoxxe = 2*zxxvar/((zxxr**2+zxxi**2)**0.5)
    #     phsxxe = zxxvar/((zxxr**2+zxxi**2)**0.5)
    #     rhoxx = np.log10(rhoxx)
    #     phsxxe = phsxxe*180/math.pi
    #     rhoyy = (zyyr**2 + zyyi**2)/freq/5
    #     phsyy = np.arctan2(zyyi, zyyr) *180/math.pi
    #     rhoyye = 2*zyyvar/((zyyr**2+zyyi**2)**0.5)
    #     phsyye = zyyvar/((zyyr**2+zyyi**2)**0.5)
    #     rhoyy = np.log10(rhoyy)
    #     phsyye = phsyye*180/math.pi

    #     minrho = min(min(rhoxx),min(rhoyy))
    #     minrho = math.floor(minrho - 0.5)
    #     maxrho = max(max(rhoxx),max(rhoyy))
    #     maxrho = round(maxrho + 1)

    #     # 开始画图
    #     ax1 = plt.subplot(2,1,1)
    #     ax1.errorbar(freq,rhoxx,yerr=rhoxxe,fmt='xr',ecolor = 'r',capsize = 0)
    #     ax1.errorbar(freq,rhoyy,yerr=rhoyye,fmt='+c',ecolor = 'c',capsize = 0)
    #     ax1.set_ylim(minrho,maxrho)


    #     ax2 = plt.subplot(2,1,2)
    #     ax2.errorbar(freq,phsxx,yerr=phsxxe,fmt='xr',ecolor = 'r',capsize = 0)
    #     ax2.errorbar(freq,phsyy,yerr=phsyye,fmt='+c',ecolor = 'c',capsize = 0)
    #     ax2.set(ylim=(-180,180), yticks=[-180,-135,-90,-45,0,45,90,135,180], size=12)
    #     ax2.set_yticks([-180,-135,-90,-45,0,45,90,135,180], size=12)
        

    if opt == 'xyyx':
        # %=========calculate rho and phase from data=======%
        rhoxy1 = (zxyr1**2 + zxyi1**2)/freq1/5
        phsxy1 = np.arctan2(zxyi1, zxyr1) *180/math.pi
        rhoxye1 = 2*zxyvar1/((zxyr1**2+zxyi1**2)**0.5)
        phsxye1 = zxyvar1/((zxyr1**2+zxyi1**2)**0.5)
        rhoxy1 = np.log10(rhoxy1)
        phsxye1 = phsxye1*180/math.pi
        rhoyx1 = (zyxr1**2 + zyxi1**2)/freq1/5
        phsyx1 = np.arctan2(zyxi1, zyxr1) *180/math.pi
        rhoyxe1 = 2*zyxvar1/((zyxr1**2+zyxi1**2)**0.5)
        phsyxe1 = zyxvar1/((zyxr1**2+zyxi1**2)**0.5)
        rhoyx1 = np.log10(rhoyx1)
        phsyxe1 = phsyxe1*180/math.pi
        
        minrho1 = min(min(rhoyx1),min(rhoxy1))
        minrho1 = math.floor(minrho1 - 0.5)
        maxrho1 = max(max(rhoyx1),max(rhoxy1))
        maxrho1 = round(maxrho1 + 1)
        
        rhoxy2 = (zxyr2**2 + zxyi2**2)/freq2/5
        phsxy2 = np.arctan2(zxyi2, zxyr2) *180/math.pi
        rhoxye2 = 2*zxyvar2/((zxyr2**2+zxyi2**2)**0.5)
        phsxye2 = zxyvar2/((zxyr2**2+zxyi2**2)**0.5)
        rhoxy2 = np.log10(rhoxy2)
        phsxye2 = phsxye2*180/math.pi
        rhoyx2 = (zyxr2**2 + zyxi2**2)/freq2/5
        phsyx2 = np.arctan2(zyxi2, zyxr2) *180/math.pi
        rhoyxe2 = 2*zyxvar2/((zyxr2**2+zyxi2**2)**0.5)
        phsyxe2 = zyxvar2/((zyxr2**2+zyxi2**2)**0.5)
        rhoyx2 = np.log10(rhoyx2)
        phsyxe2 = phsyxe2*180/math.pi
        
        minrho2 = min(min(rhoyx2),min(rhoxy2))
        minrho2 = math.floor(minrho2 - 0.5)
        maxrho2 = max(max(rhoyx2),max(rhoxy2))
        maxrho2 = round(maxrho2 + 1)



        # 开始画图
        ax1 = plt.subplot(2,2,1) 

        ax1.errorbar(freq1,rhoxy1,yerr=rhoxye1,fmt='xr',ecolor = 'r',capsize = 0)
        ax1.errorbar(freq1,rhoyx1,yerr=rhoyxe1,fmt='.c',ecolor = 'c',capsize = 0)
        ax1.set_ylim(0,3)
        ax1.set_title('poor')
        
        ax2 = plt.subplot(2,2,2) 

        ax2.errorbar(freq2,rhoxy2,yerr=rhoxye2,fmt='xr',ecolor = 'r',capsize = 0)
        ax2.errorbar(freq2,rhoyx2,yerr=rhoyxe2,fmt='.c',ecolor = 'c',capsize = 0)
        ax2.set_ylim(0,3)
        ax2.set_title('good')

        ax3 = plt.subplot(2,2,3)
        ax3.errorbar(freq1,phsxy1,yerr=phsxye1,fmt='xr',ecolor = 'r',capsize = 0)
        ax3.errorbar(freq1,phsyx1,yerr=phsyxe1,fmt='.c',ecolor = 'c',capsize = 0)
        ax3.set_ylim(0,90)
        ax3.set_yticks([0, 15, 30, 45, 60, 75, 90])
       
        
        ax4 = plt.subplot(2,2,4)
        ax4.errorbar(freq2,phsxy2,yerr=phsxye2,fmt='xr',ecolor = 'r',capsize = 0)
        ax4.errorbar(freq2,phsyx2,yerr=phsyxe2,fmt='.c',ecolor = 'c',capsize = 0)
        ax4.set_ylim(0,90)
        ax4.set_yticks([0, 15, 30, 45, 60, 75, 90])
        
    
    ax1.grid(linestyle='-.',which = 'both',axis='both',color='#e8e8e8') #浅灰色'#e8e8e8'
    ax2.grid(linestyle='-.',which = 'both',axis='both',color='#e8e8e8')
    ax3.grid(linestyle='-.',which = 'both',axis='both',color='#e8e8e8') #浅灰色'#e8e8e8'
    ax4.grid(linestyle='-.',which = 'both',axis='both',color='#e8e8e8')
    ax1.set(xlim=[min(freq1)/3,max(freq1)*3],xscale = 'log', ylabel = r'$log_{10}$'+'app. resistivity (Ohm*m)')
    ax3.set(xlim=[min(freq2)/3,max(freq2)*3],xscale = 'log',xlabel='frequency (Hz)', ylabel = 'phase(degree)')
    ax2.set(xlim=[min(freq1)/3,max(freq1)*3],xscale = 'log', ylabel = r'$log_{10}$'+'app. resistivity (Ohm*m)')
    ax4.set(xlim=[min(freq2)/3,max(freq2)*3],xscale = 'log',xlabel='frequency (Hz)', ylabel = 'phase(degree)')

    ax1.invert_xaxis()
    ax2.invert_xaxis()
    ax3.invert_xaxis()
    ax4.invert_xaxis()


    # if prt == 'on' or prt == 'silent':
    #     # plt.savefig(r"kayabezprint\\"+eval(sitename)+ opt +'.eps',dpi=600,format='eps')
    #     # plt.savefig(r"kayabezprint\\"+eval(sitename)+ opt +'.png',dpi=600,format='png')
    #     plt.savefig(r"D:\svm\liaoxi_code\liaoxiprint\\"+eval(sitename)+ opt +'.eps',dpi=2000,format='eps')
    #     plt.savefig(r"D:\svm\liaoxi_code\liaoxiprint\\"+eval(sitename)+ opt +'.png',dpi=2000,format='png')
    # if prt != 'silent':
    #     plt.show()
    
    plt.savefig('poor and good',dpi=2000)
    