#!/usr/bin/env python2

import ROOT
import numpy as np
import sys
from styleit import gstyle

gstyle()

data = np.loadtxt(sys.argv[1])
tm, tr, pm, pr, dt, pvx0, tvx0, ip, k, eta, L, ms, fvxt, fvyt, fvt, fvxp, fvyp, fvp, tts = data.T

dx_vs_eta  = ROOT.TGraph()
vfp_vs_eta = ROOT.TGraph()
vft_vs_eta = ROOT.TGraph()
el_vs_eta  = ROOT.TGraph()

dx_vs_eta.SetTitle(';#eta;#Delta x_{max}/r_{t}')
vfp_vs_eta.SetTitle(';#eta;Projectile v_{f}')
vft_vs_eta.SetTitle(';#eta;Target v_{f}')
el_vs_eta.SetTitle(';#eta;(E_{0}-E_{f})/E_{0}')

i = 0
for ieta, ipvx0, ims, ifvt, ifvp in zip(eta,pvx0,ms,fvt,fvp):
    dx_vs_eta.SetPoint(i,ieta,(tr[0]+pr[0]-ims)/tr[0])
    vfp_vs_eta.SetPoint(i,ieta,ifvp)
    vft_vs_eta.SetPoint(i,ieta,ifvt)
    el_vs_eta.SetPoint(i,ieta,(0.5*pm[0]*ipvx0*ipvx0-(0.5*pm[0]*ifvp*ifvp+0.5*tm[0]*ifvt*ifvt))/(0.5*pm[0]*ipvx0*ipvx0))
    i = i + 1
                     
dx_vs_eta.SetMarkerStyle(7)
el_vs_eta.SetMarkerStyle(7)
vft_vs_eta.SetMarkerStyle(7)
vft_vs_eta.SetMarkerColor(ROOT.kBlue)
vfp_vs_eta.SetMarkerStyle(7)
vfp_vs_eta.SetMarkerColor(ROOT.kRed)

c1 = ROOT.TCanvas()
dx_vs_eta.Draw('AP')
c2 = ROOT.TCanvas()
el_vs_eta.Draw('AP')
c3 = ROOT.TCanvas()
vfp_vs_eta.Draw('AP')
c4 = ROOT.TCanvas()
vft_vs_eta.Draw('AP')
raw_input('')
