#!/usr/bin/env python2

import ROOT
import numpy as np
import sys

data = np.loadtxt(sys.argv[1])
tm, tr, pm, pr, dt, pvx0, tvx0, ip, k, eta, L, ms, fvxt, fvyt, fvt, fvxp, fvyp, fvp, tts = data.T

dx_vs_L = ROOT.TGraph()
ms_vs_L = ROOT.TGraph()
i = 0
for m,n in zip(ms,L):
    dx_vs_L.SetPoint(i,tr[0]+pr[0]-m,n)
    ms_vs_L.SetPoint(i,m,n)
    i = i + 1

dx_vs_L.SetMarkerStyle(7)
ms_vs_L.SetMarkerStyle(7)
ms_vs_L.SetMarkerColor(ROOT.kRed)

mg = ROOT.TMultiGraph()
#mg.Add(ms_vs_L)
mg.Add(dx_vs_L)

mg.Draw('AP')

raw_input('')
