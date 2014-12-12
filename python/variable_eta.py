#!/usr/bin/env python2

import ROOT
import numpy as np
import sys

data = np.loadtxt(sys.argv[1])
tm, tr, pm, pr, dt, pvx0, tvx0, ip, k, eta, L, ms, fvxt, fvyt, fvt, fvxp, fvyp, fvp, tts = data.T

dx_vs_ip = ROOT.TGraph()
ms_vs_ip = ROOT.TGraph()
i = 0
for m,n in zip(ms,ip):
    dx_vs_ip.SetPoint(i,n,tr[0]+pr[0]-m)
    ms_vs_ip.SetPoint(i,n,m)
    i = i + 1

dx_vs_ip.SetMarkerStyle(7)
ms_vs_ip.SetMarkerStyle(7)
ms_vs_ip.SetMarkerColor(ROOT.kRed)

mg = ROOT.TMultiGraph()
mg.Add(dx_vs_ip)

mg.Draw('AP')

raw_input('')
