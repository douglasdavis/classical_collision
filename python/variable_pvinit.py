#!/usr/bin/env python2

import ROOT
import numpy as np
import sys
from styleit import gstyle

gstyle()

data = np.loadtxt(sys.argv[1])
tm, tr, pm, pr, dt, pvx0, tvx0, ip, k, eta, L, ms, fvxt, fvyt, fvt, fvxp, fvyp, fvp, tts = data.T

dx_vs_pvinit  = ROOT.TGraph()
vfp_vs_pvinit = ROOT.TGraph()
vft_vs_pvinit = ROOT.TGraph()
el_vs_pvinit  = ROOT.TGraph()

dx_vs_pvinit.SetTitle(';Projectile v_{0};#Delta x_{max}')
vfp_vs_pvinit.SetTitle(';Projectile v_{0};Projectile v_{f}')
vft_vs_pvinit.SetTitle(';Projectile v_{0};Target v_{f}')
el_vs_pvinit.SetTitle(';Projectile v_{0};Energy Loss')

i = 0
for ipvinit, ipvx0, ims, ifvt, ifvp in zip(pvx0,pvx0,ms,fvt,fvp):
    dx_vs_pvinit.SetPoint(i,ipvinit,tr[0]+pr[0]-ims)
    vfp_vs_pvinit.SetPoint(i,ipvinit,ifvp)
    vft_vs_pvinit.SetPoint(i,ipvinit,ifvt)
    el_vs_pvinit.SetPoint(i,ipvinit,0.5*pm[0]*ipvx0*ipvx0-(0.5*pm[0]*ifvp*ifvp+0.5*tm[0]*ifvt*ifvt))
    i = i + 1
                     
dx_vs_pvinit.SetMarkerStyle(7)
el_vs_pvinit.SetMarkerStyle(7)
vft_vs_pvinit.SetMarkerStyle(7)
vft_vs_pvinit.SetMarkerColor(ROOT.kBlue)
vfp_vs_pvinit.SetMarkerStyle(7)
vfp_vs_pvinit.SetMarkerColor(ROOT.kRed)

c1 = ROOT.TCanvas()
dx_vs_pvinit.Draw('AP')
c2 = ROOT.TCanvas()
el_vs_pvinit.Draw('AP')
c3 = ROOT.TCanvas()
vfp_vs_pvinit.Draw('AP')
c4 = ROOT.TCanvas()
vft_vs_pvinit.Draw('AP')
raw_input('')
