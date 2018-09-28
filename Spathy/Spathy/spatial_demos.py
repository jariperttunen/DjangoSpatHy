#!/home/jarip/spathy/bin/python
# -*- coding: utf-8 -*-
"""
Spatial demonstrations of Spathy

Created on Tue Sep 19 13:21:26 2017

@author: slauniai
"""
#This trick is not to use the default X-windows with matplotlib
#on some Linux servers. 
import matplotlib
matplotlib.use('Agg')
import os
import spotpy
import numpy as np
from scipy import stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from spathy_4 import spathy_driver
import spathypath

eps = np.finfo(float).eps

setupfile=spathypath.setupfile
#setupfile = 'c:\pyspace\spathy\ini\spathy4_default.ini'
# setupfile = 'c:\pyspace\spathy\ini\spathy4_ch_1default.ini'
print('setupfile spatial demos:',setupfile)
spa, ncf = spathy_driver(setupfile, catchment_id='3', ncf=True, cpy_outputs=True, bu_outputs=True, top_outputs=True)
print('spathy_driver end')
gis = spa.GisData
LAIc = gis['LAI_conif']
LAId = gis['LAI_decid']
soil = gis['soil']

peat_ix = np.where(soil == 4)
med_ix = np.where(soil == 2)
coarse_ix = np.where(soil == 1)

# cpy results
cres = spa.cpy.results
bres = spa.bu.results
tres = spa.top.results

TR = np.array(cres['Transpi'])  # tranpi
EF = np.array(cres['Efloor'])  # floor evap
IE = np.array(cres['Interc'])  # interception evap

ET = TR + EF + IE

# snow, seek maximum timing
SWE = np.array(cres['SWE'])  # SWE
a = np.nansum(SWE, axis=1)
a = np.nansum(a, axis=1)
swe_max_ix = int(np.where(a == np.nanmax(a))[0])

# rootzone water storage; seek maximum and minimum
W = np.array(bres['Wliq'])
a = np.nansum(W, axis=1)
a = np.nansum(a, axis=1)
ix_slow = int(np.where(a == np.nanmin(a))[0])
ix_shi = int(np.where(a == np.nanmax(a))[0])

DR = np.array(bres['Drain'])*1e3

P = np.sum(spa.FORC['Prec'])*spa.dt  # precip

# plot figures of annual ratios
#sns.color_palette('muted')
#plt.figure()
#plt.subplot(221)
#plt.imshow(LAIc + LAId); plt.colorbar(); plt.title('LAI (m$^2$m,^{-2}$)')
#
#plt.subplot(222)
#plt.imshow(np.sum(ET, axis=0)/P); plt.colorbar(); plt.title('ET/P (-)')
#
#plt.subplot(223)
#plt.imshow(np.sum(TR, axis=0)/P); plt.colorbar(); plt.title('T$_r$/P (-)')
#
#plt.subplot(224)
#plt.imshow(np.sum(IE, axis=0)/P); plt.colorbar(); plt.title('E/P (-)')

# sns.color_palette('muted')
print("Plot figure")
pdf = PdfPages('SpathyFigures.pdf')

fig=plt.figure()
plt.subplot(321)
sns.heatmap(LAIc + LAId, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('LAI (m$^2$m$^{-2}$)')

plt.subplot(322)
sns.heatmap(LAId /(LAIc + LAId), cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('LAI$_d$ fraction (-)')

plt.subplot(323)
sns.heatmap(np.sum(ET, axis=0)/P, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('ET/P (-)')

plt.subplot(324)
sns.heatmap(np.sum(IE, axis=0)/P, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('E/P (-)')

plt.subplot(325)
sns.heatmap(np.sum(TR, axis=0)/P, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('T$_r$/P (-)')

plt.subplot(326)
sns.heatmap(np.sum(EF, axis=0)/P, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('E$_f$/P (-)')
pdf.savefig()

plt.figure()
sns.heatmap(SWE[swe_max_ix,:,:], cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('SWE (mm)')
plt.savefig('Figure1.png')
pdf.savefig()

plt.figure()
Wliq = np.array(bres['Wliq'])
ss = np.array(bres['Wliq'] / spa.bu.poros)
S = np.array(tres['S'])
s_hi = spa.top.local_s(S[ix_shi])
s_hi[s_hi<0] = 0.0
s_low = spa.top.local_s(S[ix_slow])
s_low[s_low<0] = 0.0

plt.subplot(321)
sns.heatmap(Wliq[ix_slow,:,:], cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('Wliq (m3/m3)')
plt.subplot(322)
sns.heatmap(Wliq[ix_shi,:,:], cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('Wliq (m3/m3)')
plt.subplot(323)
sns.heatmap(ss[ix_slow,:,:], cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('root zone sat. ratio (-)')
plt.subplot(324)
sns.heatmap(ss[ix_shi,:,:], cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('root zone sat. ratio (-)')
plt.subplot(325)
sns.heatmap(s_low, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('topmodel sat deficit (m)')
plt.subplot(326)
sns.heatmap(s_hi, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); plt.title('topmodel sat deficit (m)')
plt.savefig('Figure2.png')
pdf.savefig()
# plot LAI -relations

plt.figure()
x = LAId + LAIc
print('Before with sns.color_palette')
i=3
with sns.color_palette('muted'):
    plt.subplot(221)
    plt.plot(x[peat_ix], np.sum(ET, axis=0)[peat_ix]/P, 'o', alpha=0.3, label='peat')
    plt.plot(x[med_ix], np.sum(ET, axis=0)[med_ix]/P, 'o', alpha=0.3, label='medium textured')
    plt.plot(x[coarse_ix], np.sum(ET, axis=0)[coarse_ix]/P, 'o', alpha=0.3, label='coarse textured')
    plt.ylabel('ET/P'); plt.xlabel('LAI (m$^2$m$^{-2}$)')
    plt.legend(loc=4)

    plt.subplot(222)
    plt.plot(x[peat_ix], np.sum(TR, axis=0)[peat_ix]/P, 'o', alpha=0.3, label='peat')
    plt.plot(x[med_ix], np.sum(TR, axis=0)[med_ix]/P, 'o', alpha=0.3, label='medium textured')
    plt.plot(x[coarse_ix], np.sum(TR, axis=0)[coarse_ix]/P, 'o', alpha=0.3, label='coarse textured')
    plt.ylabel('T$_r$/P'); plt.xlabel('LAI (m$^2$m$^{-2}$)')
    # plt.legend()

    plt.subplot(223)
    plt.plot(x[peat_ix], np.sum(IE, axis=0)[peat_ix]/P, 'o', alpha=0.3, label='peat')
    plt.plot(x[med_ix], np.sum(IE, axis=0)[med_ix]/P, 'o', alpha=0.3, label='medium textured')
    plt.plot(x[coarse_ix], np.sum(IE, axis=0)[coarse_ix]/P, 'o', alpha=0.3, label='coarse textured')
    plt.ylabel('E/P'); plt.xlabel('LAI (m$^2$m$^{-2}$)')
    # plt.legend()
    
    plt.subplot(224)
    plt.plot(x[peat_ix], np.sum(EF, axis=0)[peat_ix]/P, 'o', alpha=0.3, label='peat')
    plt.plot(x[med_ix], np.sum(EF, axis=0)[med_ix]/P, 'o', alpha=0.3, label='medium textured')
    plt.plot(x[coarse_ix], np.sum(EF, axis=0)[coarse_ix]/P, 'o', alpha=0.3, label='coarse textured')
    plt.ylabel('E$_f$/P'); plt.xlabel('LAI (m$^2$m$^{-2}$)')
    # plt.legend()
    plt.savefig('Figure'+str(i)+'.png')
    i=i+1
    pdf.savefig()
    plt.figure()
    y = LAId / (LAIc + LAId)
    
    plt.plot(x[peat_ix], y[peat_ix], 'o', alpha=0.3, label='peat')
    plt.plot(x[med_ix], y[med_ix], 'o', alpha=0.3, label='peat')
    plt.plot(x[coarse_ix], y[coarse_ix], 'o', alpha=0.3, label='peat')
    plt.savefig('Figure'+str(i)+'.png')
    i=i+1
    pdf.savefig()
    #plt.close()
#r = np.size(LAI)
#x = np.reshape(LAIc + LAId,r)
#tr = np.reshape(np.sum(TR, axis=0) / P, r)
#et = np.reshape(np.sum(TR + EF + IE, axis=0) / P, r)
#ef = np.reshape(np.sum(EF, axis=0) / P, r)
#ie = np.reshape(np.sum(IE, axis=0) / P, r)
#c = np.reshape(LAId / (LAIc + LAId), r)
#s = np.reshape(soil, r)
#
#
#dat = {'LAI': x, 'dfrac': c, 'TR': tr, 'ET': et, 'IE': ef,
#       'IE': ie, 'soil': s}
#       
#data = pd.DataFrame(dat, index=np.arange(r))
#data = data.dropna()
#
#plt.figure()
#sns.stripplot(x='LAI', y='ET', data=data, hue='soil')
#plt.show()
#plt.savefig('Spathy.pdf')
# plot Qt: lisää tämä, piirtää mallinnetun valunnan aikasarjan
Qt = 1e3 * np.array(tres['Qt'])  # mm/d
plt.figure()
plt.plot(Qt)
plt.ylim([0, 20])
plt.savefig('Figure'+str(i)+'.png')
#plt.close()
pdf.savefig()
del Qt
pdf.close()
print("Demos done, end of program")
