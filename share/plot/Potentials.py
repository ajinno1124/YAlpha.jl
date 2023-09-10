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

d1_dense=pd.read_csv('./Potentials/Potential_Chi2mom_a3tuned_meff.dat',delimiter='\t',comment="#")
d2_dense=pd.read_csv('./Potentials/Potential_Chi3mom_a3tuned_meff.dat',delimiter='\t',comment="#")
d3_dense=pd.read_csv('./Potentials/Potential_LY4_a3tuned_meff.dat',delimiter='\t',comment="#")

fig=plt.figure()
#subplots_adjust(hspace=0.0,wspace=0.0,left=0.2,right=0.85)
subplots_adjust(hspace=0.0,wspace=0.0,top=0.97,left=0.24,right=0.96)
ax = subplot(1,1,1)

ax.plot(d1_dense["r(fm)"],d1_dense["U_local(MeV)"],label=r"Chi2 U",linewidth=2,color='b',linestyle='-')
ax.plot(d2_dense["r(fm)"],d2_dense["U_local(MeV)"],label=r"Chi3 U",linewidth=2,color='r',linestyle='-')
ax.plot(d3_dense["r(fm)"],d3_dense["U_local(MeV)"],label=r"LY-IV U",linewidth=2,color='k',linestyle='-')

ax.plot(d1_dense["r(fm)"],d1_dense["U_m(MeV)"],label=r"Chi2 Um",linewidth=2,color='b',linestyle='--')
ax.plot(d2_dense["r(fm)"],d2_dense["U_m(MeV)"],label=r"Chi3 Um",linewidth=2,color='r',linestyle='--')
ax.plot(d3_dense["r(fm)"],d3_dense["U_m(MeV)"],label=r"LY-IV Um",linewidth=2,color='k',linestyle='--')


ax.legend(loc='best',frameon=0,numpoints=1,fontsize=16)
ax.set_xlim(0,5)
ax.set_ylim(-35,10)
#ax.set_xticks(np.arange(0,2.1,0.5),fontsize=14)
#plt.yticks(arange(-40,30.1,10), fontsize=14)
ax.set_ylabel(r'V (MeV)',fontsize=16)
ax.set_xlabel(r'r (fm)',fontsize=16)

ax.tick_params(axis='x',top='true',bottom='true',direction='in')
ax.tick_params(axis='y',which='both',left='true',right='true',direction='in')
ax.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig("Potentials.pdf",dpi=300)
plt.savefig("Potentials.png",dpi=300)
#plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。