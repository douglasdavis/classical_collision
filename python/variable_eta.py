#!/usr/bin/env python2

import ROOT
import numpy as np
import sys

data = np.loadtxt(sys.argv[1])
tm, tr, pm, pr, dt, pvx0, tvx0, ip, k, eta, L, ms, fvxt, fvyt, fvt, fvxp, fvyp, fvp, tts = data.T

dx_vs_eta = ROOT.TGraph()
ms_vs_eta = ROOT.TGraph()
i = 0
for m,n in zip(ms,eta):
    dx_vs_eta.SetPoint(i,n,tr[0]+pr[0]-m)
    ms_vs_eta.SetPoint(i,n,m)
    i = i + 1

dx_vs_eta.SetMarkerStyle(7)
ms_vs_eta.SetMarkerStyle(7)
ms_vs_eta.SetMarkerColor(ROOT.kRed)

mg = ROOT.TMultiGraph()
mg.Add(dx_vs_eta)

mg.Draw('AP')

raw_input('')
