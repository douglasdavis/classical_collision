#!/usr/bin/env python2

import ROOT
import numpy as np
import sys

data = np.loadtxt(sys.argv[1])
tm, tr, pm, pr, dt, pvx0, tvx0, ip, k, eta, L, ms, fvxt, fvyt, fvt, fvxp, fvyp, fvp = data.T

print ms
print eta

ms_vs_eta = ROOT.TGraph()
i = 0
for m,n in zip(ms,eta):
    ms_vs_eta.SetPoint(i,m,n)
    i = i + 1
ms_vs_eta.SetMarkerStyle(7)
ms_vs_eta.Draw('AP')

raw_input('')
