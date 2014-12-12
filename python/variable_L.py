#!/usr/bin/env python2

import ROOT
import numpy as np
import sys
from styleit import gstyle

gstyle()

data = np.loadtxt(sys.argv[1])
tm, tr, pm, pr, dt, pvx0, tvx0, ip, k, eta, L, ms, fvxt, fvyt, fvt, fvxp, fvyp, fvp, tts = data.T

dx_vs_L  = ROOT.TGraph()
vfp_vs_L = ROOT.TGraph()
vft_vs_L = ROOT.TGraph()
el_vs_L  = ROOT.TGraph()

dx_vs_L.SetTitle(';#lambda;#Delta x_{max}')
vfp_vs_L.SetTitle(';#lambda;Projectile v_{f}')
vft_vs_L.SetTitle(';#lambda;Target v_{f}')
el_vs_L.SetTitle(';#lambda;Energy Loss')

i = 0
for iL, ipvx0, ims, ifvt, ifvp in zip(L,pvx0,ms,fvt,fvp):
    dx_vs_L.SetPoint(i,iL,tr[0]+pr[0]-ims)
    vfp_vs_L.SetPoint(i,iL,ifvp)
    vft_vs_L.SetPoint(i,iL,ifvt)
    el_vs_L.SetPoint(i,iL,0.5*pm[0]*ipvx0*ipvx0-(0.5*pm[0]*ifvp*ifvp+0.5*tm[0]*ifvt*ifvt))
    i = i + 1
                     
dx_vs_L.SetMarkerStyle(7)
el_vs_L.SetMarkerStyle(7)
vft_vs_L.SetMarkerStyle(7)
vft_vs_L.SetMarkerColor(ROOT.kBlue)
vfp_vs_L.SetMarkerStyle(7)
vfp_vs_L.SetMarkerColor(ROOT.kRed)

c1 = ROOT.TCanvas()
dx_vs_L.Draw('AP')
c2 = ROOT.TCanvas()
el_vs_L.Draw('AP')
c3 = ROOT.TCanvas()
vfp_vs_L.Draw('AP')
c4 = ROOT.TCanvas()
vft_vs_L.Draw('AP')
raw_input('')
