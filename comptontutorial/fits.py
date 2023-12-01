#!/usr/bin/env python3

import sys, ROOT, argparse, pickle, os, subprocess, pathlib, math
from importlib import import_module
from parameters import Constants
from crosssection import pCALO , eCherenkov
from time import asctime

def cpuinfo():
  if (os.path.exists('/proc/cpuinfo')):
    with open('/proc/cpuinfo', 'r') as f: info = f.readlines()
    model = [el for el in info if 'model name' in el]
    return model[0].strip().split(': ')[1]
  else:
    model = subprocess.check_output('sysctl -n machdep.cpu.brand_string',shell=True).decode('utf-8')
    return model.strip()
    
def Roun(x, dx): #rounding
    if dx==None:
      n = 1
      if x == 0.0: a = x
      else:
        a = round(abs(x), n)
        while  a == 0.0:
          n += 1
          a = round(abs(x), n)
          if a > 0.0: break
        if abs(x) > 1: x = round(x, 1)
        else: x = round(x, n+1)
    else:
      if dx == 0.0: dx = x
      else:
        n = 1
        a = round(abs(dx), n)
        while  a == 0.0:
          n += 1
          a = round(abs(dx), n)
          if a > 0.0: break
        if abs(dx) > 1:  x = round(x, 1)
        else:  x = round(x, n+1)
    return x
  
def Roune(x, dx):
    x = Roun(x, dx)
    if x == 0.0: x = ' fixed'
    else: x = ' ± ' + str(x)
    return x

class POLARIMETER:
  def __init__(self, filename=None):
    if filename:
      abspath = pathlib.Path(filename).absolute()
      print(str(abspath))
      with open(abspath, 'rb') as fp:
        self.MC = pickle.load(fp)
    else:
      print('Please specify filename!'); exit()
    # we now save the hardware setting in a uuid filename
    sfname = filename.replace('.MC','_hardsave.py')
    with open(sfname, 'w') as fp:  fp.writelines(self.MC['HW'])
    if (sfname[:2]=='./'):
      sfname = sfname[2:]
    elif (sfname[1:]=='/'):
      sfname = sfname[1:]
    hardsave = import_module(sfname.replace('/','.')[:-3])
    self.Laser, self.Beam, self.Det = hardsave.Laser(), hardsave.Beam(), hardsave.Detector()
    self.HEp = self.MC['Ep']
    self.DEp = self.HEp.Clone()
    self.HXe = self.MC['Xe']
    self.DXe = self.HXe.Clone()
    self.HXe.SetStats(0); self.HXe.SetTitle('Electrons: MC')
    self.HEp.SetStats(0); self.HEp.SetTitle('Photons: MC')
    self.DXe.SetStats(0); self.DXe.SetTitle('Electrons: residuald')
    self.DEp.SetStats(0); self.DEp.SetTitle('Photons: residuals')
    #

  def ParametersMC(self):
    self.ParametersTable = ROOT.TPaveText(.05, .05, .95, .96)
    self.ParametersTable.AddText('Monte-Carlo Parameters:')
    self.ParametersTable.AddText('Laser #lambda_{0} = %5.3f um' %       (1e+4*self.Laser.λo))
    self.ParametersTable.AddText('Electron E_{0} = %6.3f GeV' %         (1e-9*self.Beam.Eo))
    self.ParametersTable.AddText('Electron #gamma = %5.3f#times10^{3}'% (1e-3*self.Beam.γ))
    self.ParametersTable.AddText('Compton x_{#c} = %5.3f'      %        (     self.Beam.xc))
    self.ParametersTable.AddText('Bend: #gamma#theta_{0} = %5.3f'  %    (     self.Beam.γ*self.Beam.θo))
    self.ParametersTable.AddText('PC = %5.3f' % (self.Laser.PC))
    self.ParametersTable.AddText('Pz = %5.3f' % (self.Beam.Pz))

  def PhotonsSetFunction(self):
    self.pfit = pCALO(Det = self.Det, setup = self.Beam)
    self.pEY = ROOT.TF1('pEY', self.pfit , self.Det.E_min, self.Det.E_max, 5)
    self.pEY.SetParName(0, 'PCPz');   self.pEY.SetParameter(0, 0.); #self.pEY.SetParLimits(0,-1.1, 1.1)
    self.pEY.SetParName(1, 'resolA');   self.pEY.SetParameter(1, 0.03); self.pEY.SetParLimits(1,0., 1.)
    self.pEY.SetParName(2, 'resolB');   self.pEY.SetParameter(2, 0.00); self.pEY.SetParLimits(2,0., 1.)
    self.pEY.SetParName(3, 'resolC');   self.pEY.SetParameter(3, 0.01); self.pEY.SetParLimits(3,0., 1.)
    self.pEY.SetParName(4,' norm');    self.pEY.SetParameter(4, 1)  # amplitude
    self.pEY.SetNpx(self.Det.E_npix)
    self.pEY.SetTitle('Photons: Fit'); #self.pEY.GetXaxis().SetTitle('E [eV]')

  def FitPhotons(self):
    fixed_parameters = [1,2,3] #[2,3] / []
    for p in fixed_parameters: self.pEY.FixParameter(p, self.pEY.GetParameter(p))
    self.FitPhotonsResult =  self.HEp.Fit(self.pEY, 'SVNMEQ')
    return not self.FitPhotonsResult.Status()

  def PhotonsResiduals(self):
    self.pZeros = 0
    for bin in range(1, self.Det.E_npix + 1):
      H = self.DEp.GetBinContent(bin)
      if H:
        F = self.pEY.Eval(self.Det.E_min+(bin-0.5)*self.Det.E_pix)
        err = self.DEp.GetBinError(bin)
        if err:
          self.DEp.SetBinContent(bin, (F - H)/err**0.5)
      if H<1: self.pZeros += 1
    self.DEp.SetTitle('Photons: (Fit - MC)/(MC)^{1/2}')

  def ElectronsSetFunction(self):
    self.efit = eCherenkov(Det = self.Det, setup = self.Beam)
    self.eX = ROOT.TF1('eX', self.efit ,self.Det.X_beam, self.Det.X_beam+self.Det.X_size, 4)
    L1, L2 = self.Beam.leip_L, self.Beam.spec_L
    θo, ymax  = self.Beam.θo    , self.Beam.ymax
    Eo = self.Beam.Eo
    u = ymax/(1-ymax)
    σX = 0.1
    self.eX.SetParName(0,  'X1');      self.eX.SetParameter(0,   0.0 ); #self.eX.SetParLimits(0,  -10, 10); self.eX.SetParError(0,0.1)# beam position x, mm
    self.eX.SetParName(1,  'σX');      self.eX.SetParameter(1,  σX ) ; self.eX.SetParLimits(1, 0., 100.); # resolution
    self.eX.SetParName(2,  'PCPz');   self.eX.SetParameter(2, 0.0 ); #self.eX.SetParLimits(2,  -1.1, 1.1)
    self.eX.SetParName(3, 'norm');    self.eX.SetParameter(3, 1)
    self.eX.SetNpx(self.Det.X_npix)
    self.eX.SetTitle('Electrons: Fit'); self.eX.GetXaxis().SetTitle('X, mm')

  def FitElectrons(self):
    fixed_parameters = [] #
    for p in fixed_parameters: self.eX.FixParameter(p, self.eX.GetParameter(p))
    self.FitElectronsResult =  self.HXe.Fit(self.eX, 'SVNMEQ') # Use observed errors instead of expected
    return not self.FitElectronsResult.Status()

  def ElectronsResiduals(self):
    for   binx in range(1, self.Det.X_npix+1):
      H = self.DXe.GetBinContent(binx)
      if H:
        F = self.eX.Eval(self.Det.X_beam + (binx-0.5)*self.Det.X_pix)
        err = self.DXe.GetBinError(binx)
        if err:
          self.DXe.SetBinContent(binx, (F-H)/err**0.5)
    self.DXe.SetTitle('Electrons: (Fit - MC)/(MC)^{1/2}')

  def ElectronsResults(self):
    chi2 = self.FitElectronsResult.Chi2()
    NDF  = self.FitElectronsResult.Ndf() ; prob = self.FitElectronsResult.Prob()
    print('χ² = {:9.1f}, NDF = {:7d}, p-value:{:7.5f}'.format(chi2, NDF, prob))
    X1   = self.eX.GetParameter(0);       dX1 = self.eX.GetParError(0)
    σX   = self.eX.GetParameter(1);       dσX = self.eX.GetParError(1)
    self.ElectronsTable = ROOT.TPaveText(.05, .05, .95, .96)
    self.ElectronsTable.AddText(cpuinfo())
    self.ElectronsTable.AddText('Electrons fit: t = {:.3f} s (CPU {:.3f} s)'.format(ROOT.gBenchmark.GetRealTime('EFit'), ROOT.gBenchmark.GetCpuTime('EFit')))
    self.ElectronsTable.AddText('#chi^{2}/NDF = %.1f/%d | Prob = %.4f' % (chi2, NDF, prob))
    self.ElectronsTable.AddText('X_{1} = %7.4f #pm %5.3f mm' % (X1, dX1))
    self.ElectronsTable.AddText('#sigma_{X} = %7.3f #pm %5.3f mm' % (σX, dσX))
    self.ElectronsTable.AddText('PCPz = %05.3f #pm %5.3f'  % (self.eX.GetParameter(2), self.eX.GetParError(2)))
  
  def PhotonsResults(self):
    chi2 = self.FitPhotonsResult.Chi2(); print('Chi2: %f'        % (chi2))
    NDF  = self.FitPhotonsResult.Ndf();  print('NDF: %d'         % (NDF) )
    prob = self.FitPhotonsResult.Prob(); print('Probability: %f' % (prob))
    self.PhotonsTable = ROOT.TPaveText(.05, .05, .95, .96)
    self.PhotonsTable.AddText(cpuinfo())
    self.PhotonsTable.AddText('Photons fit: t = {:.0f} s (CPU {:.0f} s)'.format(ROOT.gBenchmark.GetRealTime('PFit'), ROOT.gBenchmark.GetCpuTime('PFit')))
    self.PhotonsTable.AddText('#chi^{2}/NDF = %.1f/%d | Prob = %.4f' % (chi2,NDF,prob))
    self.PhotonsTable.AddText('PCPz = %05.3f #pm %5.3f'  % (self.pEY.GetParameter(0), self.pEY.GetParError(0)))

  def ParametersT(self, PA):
    name1 = 'Monte-Carlo Parameters:'
    Laser_lambda0 = 1e+4*self.Laser.λo
    Electron_E0 = 1e-9*self.Beam.Eo
    Electron_gamma = 1e-3*self.Beam.γ
    Compton_xc = self.Beam.xc
    Bend_gamma_theta0 = self.Beam.γ*self.Beam.θo
    PC = self.Laser.PC
    Pz = self.Beam.Pz
    PA[0] = name1 
    PA[1] = 'Laser λo = ' + str(Roun(Laser_lambda0, dx = None)) + ' um'
    PA[2] = 'Electron Eo = ' + str(Roun(Electron_E0, dx = None)) + ' GeV'
    PA[3] = 'Electron γ = ' + str(Roun(Electron_gamma, dx = None)) + '*10^3'
    PA[4] = 'Compton xc = ' + str(Roun(Compton_xc, dx = None))
    PA[5] = 'Bend: γθo = ' + str(Roun(Bend_gamma_theta0, dx = None))
    PA[6] = 'PC = ' + str(Roun(PC, dx = None))
    PA[7] = 'Pz = ' + str(Roun(Pz, dx = None))
    return PA
   
  def ElectronsResultsT(self, PA):
    chi2 = self.FitElectronsResult.Chi2()
    NDF  = self.FitElectronsResult.Ndf() ; prob = self.FitElectronsResult.Prob()
    X1   = self.eX.GetParameter(0);       dX1 = self.eX.GetParError(0)
    σX = self.eX.GetParameter(1);       dσX = self.eX.GetParError(1)
    PCPz = self.eX.GetParameter(0);     dPCPz  = self.eX.GetParError(0)
    PA[0] = 'Electrons fit: '
    PA[1] = 'chi2/NDF = ' + str(Roun(chi2, dx = None)) + '/' + str(Roun(NDF, dx = None)) + ' | Prob = ' + str(Roun(prob, dx = None))
    PA[2] = 'X1 = ' + str(Roun(X1, dX1)) + str(Roune(dX1, dX1)) + ' mm'
    PA[3] = '#sigma_{X} = ' + str(Roun(σX, dσX)) + str(Roune(dσX, dσX)) + ' mm'
    PA[4] = 'PCPz = ' + str(Roun(PCPz, dPCPz)) + str(Roune(dPCPz, dPCPz))
    return PA
  
  def PhotonsResultsT(self, PA):
    chi2 = self.FitPhotonsResult.Chi2()
    NDF  = self.FitPhotonsResult.Ndf()
    prob = self.FitPhotonsResult.Prob()
    PCPz   = self.pEY.GetParameter(0);     dPCPz  = self.pEY.GetParError(0)
    PA[0] = 'Photons fit:'
    PA[1] = 'chi2/NDF = ' + str(Roun(chi2, dx = None)) + '/' + str(Roun(NDF, dx = None)) + ' | Prob = ' + str(Roun(prob, dx = None))
    PA[2] = 'PCPz = ' + str(Roun(PCPz, dPCPz)) + str(Roune(dPCPz, dPCPz))
    return PA
    
class DISPLAY:
  def __init__(self):
    ROOT.gStyle.SetOptFit(1111);    ROOT.gStyle.SetOptStat('ne');    ROOT.gStyle.SetFitFormat("8.3g")
    ROOT.gROOT.SetStyle("Plain");   ROOT.gROOT.ForceStyle()      #   ROOT.gStyle.SetPalette(56)
    self.cv = ROOT.TCanvas('cv','cv',0,0,1600,1200);               self.cv.Divide(3,3)

  def ShowOnPad(self, nPad, entity, fname='', grid=False, goption=''):
    self.cv.cd(nPad) 
    if grid: self.cv.GetPad(nPad).SetGrid()
    entity.Draw(goption)
    self.cv.Modified();    self.cv.Update()
    if nPad==9:
      self.cv.Print(fname + '.png')
   

  def ShowOnCanvas(self, entity, grid=False, goption=''):
    self.cv1 = ROOT.TCanvas('cv1','cv1',10,10,1000,800)
    if grid: self.cv1.SetGrid()
    entity.Draw(goption)
    self.cv1.Modified();    self.cv1.Update()


def main(argv):
  parser = argparse.ArgumentParser(description='Process cli arguments.')
  parser.add_argument('-f','--folder', nargs='?', default='.', type=str, help='Name of the datafile.', const='.')
  parser.add_argument('-v','--verbose', action='store_true', help = 'ask questions to progress')
  parser.add_argument('filename', type=str, help='Name of the datafile.')
  args = parser.parse_args()
  
  fname = args.filename
  fname = fname.replace('.MC', '_Fit')
  foldername = args.folder
  filedirname = foldername+'/'+args.filename
  fitfiledirname =foldername+'/'+fname

  DATA = POLARIMETER(filename=filedirname)
  LOOK = DISPLAY()
  DATA.ParametersMC()
  LOOK.ShowOnPad(nPad=1, entity = DATA.HEp, fname = fitfiledirname, grid = True,goption='E')
  LOOK.ShowOnPad(nPad=2, entity = DATA.HXe, fname = fitfiledirname, grid = True,goption='E')
  LOOK.ShowOnPad(nPad=3, entity = DATA.ParametersTable, fname = fitfiledirname)
  LOOK.ShowOnCanvas(entity = DATA.HXe, grid = True, goption='E')
  '''ROOT.gBenchmark.Start('PFit')
  DATA.PhotonsSetFunction()
  DATA.FitPhotons(); print('photon fit results')
  DATA.PhotonsResiduals()
  ROOT.gBenchmark.Stop('PFit')
  LOOK.ShowOnPad(nPad=4, entity = DATA.HEp, fname = fitfiledirname, grid = True,goption='E')
  LOOK.ShowOnPad(nPad=4, entity = DATA.pEY, fname = fitfiledirname, grid = True,goption='SAME')
  LOOK.ShowOnPad(nPad=5, entity = DATA.DEp, fname = fitfiledirname, grid = True)
  DATA.PhotonsResults()
  LOOK.ShowOnPad(nPad=6, entity = DATA.PhotonsTable, fname = fitfiledirname)'''
  if args.verbose: arg = input('Fit electrons ?')
  else:            arg = False
  ROOT.gBenchmark.Start('EFit')
  DATA.ElectronsSetFunction()
  DATA.FitElectrons()
  DATA.ElectronsResiduals()
  ROOT.gBenchmark.Stop('EFit')
  LOOK.ShowOnPad(nPad=7, entity = DATA.HXe, fname = fitfiledirname, grid = True, goption='E')
  LOOK.ShowOnPad(nPad=7, entity = DATA.eX, fname = fitfiledirname, grid = True, goption='SAME')
  LOOK.ShowOnPad(nPad=8, entity = DATA.DXe, fname = fitfiledirname, grid = True, goption='COLZ')
  DATA.ElectronsResults()
  LOOK.ShowOnPad(nPad=9, entity = DATA.ElectronsTable, fname = fitfiledirname)
  
  C = {}; D = {}; E = {}
  I = DATA.ParametersT(PA = C); J = DATA.ElectronsResultsT(PA = D);
  #K = DATA.PhotonsResultsT(PA = E)
  f= open(foldername+'/'+fname + '.txt',"w+")
  for i in range(len(I)):
    f.write(I[i] + '\n')
  f.write( '\n'+ '\n')
  for i in range(len(J)):
    f.write(J[i] + '\n')
  f.write( '\n'+ '\n')
  for i in range(len(E)):
    f.write(E[i] + '\n')
  f.close()
  
  if args.verbose: argo = input()
  exit()

if __name__ == "__main__": main(sys.argv)



