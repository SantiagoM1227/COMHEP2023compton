import sys, ROOT, argparse, pickle, os, subprocess, pathlib, math
from importlib import import_module
import numpy as np
from time import asctime

abspath = pathlib.Path('data/E7GEV/LEFT1.MC').absolute()   #Pz=-1
print(str(abspath))
with open(abspath, 'rb') as fp:
  MC = pickle.load(fp)
HXe = MC['Xe']
HEp = MC['Ep']
HXe.SetStats(0); HXe.SetTitle('Electrons: MC')
HEp.SetStats(0); HEp.SetTitle('Photons: MC')

abspath = pathlib.Path('data/E7GEV/RIGHT1.MC').absolute()   #Pz=1
print(str(abspath))
with open(abspath, 'rb') as fp:
  MC2 = pickle.load(fp)
HXe2 = MC2['Xe']
HEp2 = MC2['Ep']
HXe2.SetStats(0);
HEp2.SetStats(0);





HXeAsym = ROOT.TH1D('Asym',  'electrons asymetry', HXe.GetNbinsX(), HXe.GetXaxis().GetXmin(), HXe.GetXaxis().GetXmax())
HEpAsym = ROOT.TH1D('Asym',  'photons asymetry', HEp.GetNbinsX(), HEp.GetXaxis().GetXmin(), HEp.GetXaxis().GetXmax())

v1 = np.array(HXe)[1:-1].astype(float)
v2 = np.array(HXe2)[1:-1].astype(float)
errAsym = np.divide(2*np.real(np.emath.sqrt(v2**2*v1+v2*v1**2)),(v1+v2)**2,where=(v1+v2)>0,out=np.zeros_like(v1))
for k in range(0,v1.size):
  if ((v1[k]+v2[k])>0):
    HXeAsym.SetBinContent(k+1,(v1[k]-v2[k])/(v1[k]+v2[k]))
    HXeAsym.SetBinError(k+1,errAsym[k])
    
v1 = (np.array(HEp)[1:-1]).astype(float)
v2 = (np.array(HEp2)[1:-1]).astype(float)
errAsym = np.divide(2*np.real(np.emath.sqrt(v2**2*v1+v2*v1**2)),(v1+v2)**2,where=(v1+v2)>0,out=np.zeros_like(v1))
for k in range(0,v1.size):
  if ((v1[k]+v2[k])>0):
    HEpAsym.SetBinContent(k+1,(v1[k]-v2[k])/(v1[k]+v2[k]))
    HEpAsym.SetBinError(k+1,errAsym[k])
    
HEp_fit_func = ROOT.TF1("fit_func", "pol{}".format(5),HEp.GetXaxis().GetXmin(), HEp.GetXaxis().GetXmax())

HEp_fit = HEpAsym.Fit(HEp_fit_func, 'R')


cv = ROOT.TCanvas('cv','cv',0,0,1600,1200);
HXeAsym.Draw('p')
ROOT.gStyle.SetOptStat(0)
cv.Print('plotAsym_elec.png')

HEpAsym.Draw('P')
ROOT.gStyle.SetOptStat(0)
cv.Print('plotAsym_phot.png')

 


