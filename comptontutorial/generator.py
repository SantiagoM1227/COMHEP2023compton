#!/usr/bin/env python3
import ROOT, sys, pickle, copy, argparse, os
from ctypes import c_double
from time import asctime
from parameters import Laser, Beam, Constants, Detector as Det
from crosssection import xsec

class MonteCarlo(xsec):
  def __init__(self,seed):
    cDouble = c_double*1
    ROOT.gRandom.SetSeed(seed)
    xsec.__init__(self)
    self.rannor    = ROOT.gRandom.Rannor
    self.y, self.φ = cDouble(), cDouble()
    self.xe, self.ye = cDouble(), cDouble()
    self.se = cDouble()
    self.fmax = ROOT.TMath.Max(self.fXS.Eval(0.,0.),self.fXS.Eval(Beam.ymax,0.))
    self.xsectot = 2*Constants.π*Constants.re**2/Beam.xc*( (1-4/Beam.xc-8/Beam.xc**2)*ROOT.TMath.Log(1+Beam.xc)+0.5+8/Beam.xc-0.5/(1+Beam.xc)**2 + Laser.PC*Beam.Pz*((1+2/Beam.xc)*ROOT.TMath.Log(1+Beam.xc)-5/2.+1/(1+Beam.xc)-0.5/(1+Beam.xc)**2))
    self.lumi = Laser.E/(Laser.ωo*Constants.e)*Beam.Q/Constants.e/(2*Constants.π*ROOT.TMath.Sqrt(Laser.σy**2+Beam.σy**2)*ROOT.TMath.Sqrt(Laser.σx**2+Beam.σx**2)*ROOT.TMath.Sqrt(1+ROOT.TMath.Tan(Laser.θ0/2.)**2*(Laser.σz**2+Beam.σz**2)/(Laser.σx**2+Beam.σx**2)))
    self.lumi0 =Laser.E/(Laser.ωo*Constants.e)*Beam.Q/Constants.e/(2*Constants.π*ROOT.TMath.Sqrt(Laser.σy**2+Beam.σy**2)*ROOT.TMath.Sqrt(Laser.σx**2+Beam.σx**2))
    self.prob =self.xsectot*self.lumi*1e-2

  def generate(self):
    #self.fXS.GetRandom2(self.y, self.φ)       # get u,φ
    self.y=Beam.ymax*ROOT.gRandom.Rndm()
    self.φ=2*Constants.π*ROOT.gRandom.Rndm()
    self.func = self.fmax*ROOT.gRandom.Rndm()
    while(self.func>self.fXS.Eval(self.y,self.φ)):
      self.y=Beam.ymax*ROOT.gRandom.Rndm()
      self.φ=2*Constants.π*ROOT.gRandom.Rndm()
      self.func = self.fmax*ROOT.gRandom.Rndm()
    self.rannor(self.xe, self.ye)               # get x,y from normal distribution
    x   = self.xe[0] * Beam.σx                # electron betatron x-coordinate [mm]
    y    = self.ye[0] * Beam.σy                # electron betatron y-coordinate [mm]
    self.rannor(      self.xe, self.ye)         # Get x,y from normal distribution
    θex  = self.xe[0] * Beam.ηx                # electron betatron x-angle [rad]
    θey  = self.ye[0] * Beam.ηy                # electron betatron y-angle [rad]
    sin  = ROOT.TMath.Sin(self.φ)
    cos  = ROOT.TMath.Cos(self.φ)
    θp   = ROOT.TMath.ACos(1/Beam.β*(1.-Beam.xc*(1-self.y)/(2*Beam.γ**2*self.y)))
    θpx  = θex + θp*cos                 # photon   scattering x-angle [rad]
    θpy  = θey + θp*sin                 # photon   scattering y-angle [rad]
    θe   = ROOT.TMath.ATan(ROOT.TMath.Tan(θp)*self.y/(1-self.y))                     # electron scattering   angle [rad]
    θex  = θex - θe*cos                       # electron scattering x-angle [rad]
    θey  = θey - θe*sin                       # electron scattering y-angle [rad]
    uΔ   = self.y/(1-self.y)                # electron bending angle / θo
    Δ = ROOT.gRandom.Gaus(0,Beam.σE)
    uΔ += -Δ/(1+Δ)
    Egam = self.y*Beam.Eo*(1+Δ)
    resol = Egam*ROOT.TMath.Sqrt(Det.resolA**2/Egam/1e9+(Det.resolB/Egam/1e9)**2+Det.resolC**2)
    Egam = ROOT.gRandom.Gaus(Egam,resol)
    return x, θpx, θpy, θex, uΔ, Egam

def main(argv):
  parser = argparse.ArgumentParser(description='Process cli arguments.')
  parser.add_argument('-f','--folder', nargs='?', default='.', type=str, help='destination where to store the results.', const='.')
  parser.add_argument('-n','--nevents', nargs='?', default=1000, type=int, help='number of events to generate.', const=100)
  parser.add_argument('-v','--verbose', action='store_true', help = 'ask questions to progress')
  parser.add_argument('-s','--seed', nargs='?', default=1000, type=int, help='number of events to generate.', const=1000)
  def_name =  'PL_'+str(int(Laser.PL*10))+'_PC_'+str(int(Laser.PC*10))+'_PT_'+str(int(Beam.PT*10))+'_Pz_'+str(int(Beam.Pz*10)) 
  parser.add_argument('-N','--Name', nargs='?', default=def_name, type=str, help='Name of file.', const=def_name)
  
  args = parser.parse_args()
  foldername = args.folder
  name= args.Name
  
  
  print('Generating MC for %3d events' % args.nevents)
  
  print('Laser ωo:                %.3f eV'   % (Laser.ωo))
  print('Comptony xc-parameter:    %.4f'      % (Beam.xc))
  print('beam γ-factor:           %.3f'      % (Beam.γ))
  print('beam bending angle:      %.3f mrad' % (1000*Beam.θo))
  print('beam bending angle:      %.3f 1/γ'  % (Beam.γ*Beam.θo))
  print('γ-to-beam distance :     %.3f mm'   % (1000    *Beam.θo*Beam.spec_L))
  print('horizontal size:         %.3f mm'   % (Beam.σx))
  print('vertical size:           %.3f mm'   % (Beam.σy))
  print('horizontal spread add:   %.3f mm'   % (1000*Beam.ηx*Beam.leip_L))
  print('vertical spread add:     %.3f mm'   % (1000*Beam.ηy*Beam.leip_L))
  if args.verbose: arg = input('Continue, OK?')
  else:            arg = True

  Xemin, Xemax = Det.X_beam, Det.X_beam + Det.X_size
  Xe  = ROOT.TH1F('Xe',  'electrons  X', Det.X_npix, Xemin, Xemax)
  Xe.GetXaxis().SetTitle('X, mm'); Xe.GetYaxis().SetTitle('e^{#pm}/pixel');

  Ep  = ROOT.TH1F('Ep',  'photons  Energy', Det.E_npix, Det.E_min, Det.E_max)
  Ep.GetXaxis().SetTitle('E, [eV]');   Ep.GetYaxis().SetTitle('photon/bin')

  G = MonteCarlo(int(args.seed))

  cv = ROOT.TCanvas('cv','cv',0,0,1600,1200)
  cv.Divide(1,2)
  for p in range(2): cv.GetPad(p+1).SetGrid()
  L1     = 1000*Beam.leip_L                         # scattering base [mm]
  L2     = 1000*Beam.spec_L*ROOT.TMath.Tan(Beam.θo)    # tan(θo) * bending base [mm]
  icos   = 1./ROOT.TMath.Cos(Beam.θo)
  def cvscreen():
    cv.cd(1);   Ep.Draw();       cv.cd(2);    Xe.Draw()
    cv.Modified(); cv.Update()
    sys.stdout.flush()
  def process():
    x, θpx, θpy, θex, uΔ, Egam = G.generate()
    xe = icos*(x + L1*ROOT.TMath.Tan(θex)) + L2*uΔ
    #xe = L2*uΔ
    Ep.Fill(Egam);
    if xe > Det.X_beam:
      Xe.Fill(xe);
  
  for i in range(int(args.nevents*G.prob)):
    process()
    if (i%100000 == 0):
      print('\rProgress:  %3d%%' % int(i*100/int(args.nevents)), end=' ')
  cvscreen()
    
  print()
  
  c = {}
  c[0] = './'+foldername+'/'+ name + '.MC'
  c[1] = './'+foldername+'/'+ name +  '.png'
  c[2] = './'+foldername+'/'+ name +  'hw.txt'
  
  if args.verbose:
    arg = not ('n' or 'N') in input('Save data? (Yes/no)');

  if arg:
    with open('parameters.py', 'r') as fp:  content = fp.readlines()
    SAVE = {}
    SAVE['HW'] = content
    SAVE['Xe'] = Xe; SAVE['Ep'] = Ep
    with open(c[0], 'wb') as fp:  pickle.dump(SAVE, fp, -1)
    print('Data saved to %s' % c[0])
    cv.Print(c[1])
    with open('parameters.py') as fp:  acontent = fp.readlines()
    #print(acontent)
    with open(c[2], 'w') as f:
      for i in range(len(acontent)):
        f.write(acontent[i])
    print('Data saved to %s' % c[2])

  exit(0)

if __name__ == "__main__": main(sys.argv)



