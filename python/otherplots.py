#!/usr/bin/env python2

import ROOT
import numpy as np

data = np.loadtxt('data.dat')
tm, tr, pm, pr, dt, pvx0, ip, k, eta, L, ms, fvxt, fvyt, fvt, fvxp, fvyp, fvp = data.T

print ms
print ip

ms_vs_ip = ROOT.TGraph()
i = 0
for m,n in zip(ms,ip):
    ms_vs_ip.SetPoint(i,m,n)
    i = i + 1
ms_vs_ip.SetMarkerStyle(7)
ms_vs_ip.Draw('AP')

raw_input('')
