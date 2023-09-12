#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from pylab import *
import cmath

isMD=True

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.subplot.bottom'] = 0.15

d1_1=pd.read_csv('../CorrelationFunction/Chi3mom_a3tuned_meff_R1.0.dat',delimiter='\t',comment="#")
d1_3=pd.read_csv('../CorrelationFunction/Chi3mom_a3tuned_meff_R3.0.dat',delimiter='\t',comment="#")
d1_5=pd.read_csv('../CorrelationFunction/Chi3mom_a3tuned_meff_R5.0.dat',delimiter='\t',comment="#")

d2_1=pd.read_csv('../CorrelationFunction/LY4_a3tuned_meff_R1.0.dat',delimiter='\t',comment="#")
d2_3=pd.read_csv('../CorrelationFunction/LY4_a3tuned_meff_R3.0.dat',delimiter='\t',comment="#")
d2_5=pd.read_csv('../CorrelationFunction/LY4_a3tuned_meff_R5.0.dat',delimiter='\t',comment="#")

fig=plt.figure()
#subplots_adjust(hspace=0.0,wspace=0.0,left=0.2,right=0.85)
subplots_adjust(hspace=0.0,wspace=0.0,top=0.97,left=0.24,right=0.96)
ax = subplot(1,1,1)

ax.plot(d1_1["q(MeV/c)"],d1_1["CF"],label=r"Chi3 R=1 fm",linewidth=2,color='r',linestyle='-')
ax.plot(d1_3["q(MeV/c)"],d1_3["CF"],label=r"Chi3 R=3 fm",linewidth=2,color='r',linestyle='--')
ax.plot(d1_5["q(MeV/c)"],d1_5["CF"],label=r"Chi3 R=5 fm",linewidth=2,color='r',linestyle=':')

ax.plot(d2_1["q(MeV/c)"],d2_1["CF"],label=r"LY-IV R=1 fm",linewidth=2,color='k',linestyle='-')
ax.plot(d2_3["q(MeV/c)"],d2_3["CF"],label=r"LY-IV R=3 fm",linewidth=2,color='k',linestyle='--')
ax.plot(d2_5["q(MeV/c)"],d2_5["CF"],label=r"LY-IV R=5 fm",linewidth=2,color='k',linestyle=':')


ax.legend(loc='best',frameon=0,numpoints=1,fontsize=16)
#ax.set_xlim(0,10)
#ax.set_ylim(-1.1,1.1)
#ax.set_xticks(np.arange(0,2.1,0.5),fontsize=14)
#plt.yticks(arange(-40,30.1,10), fontsize=14)
ax.set_ylabel(r'C(q)',fontsize=16)
ax.set_xlabel(r'q (MeV/c)',fontsize=16)

ax.tick_params(axis='x',top='true',bottom='true',direction='in')
ax.tick_params(axis='y',which='both',left='true',right='true',direction='in')
ax.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig("CF.pdf",dpi=300)
plt.savefig("CF.png",dpi=300)
#plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。