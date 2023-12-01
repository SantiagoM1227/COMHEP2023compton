import ROOT
from parameters import Constants, Laser, Beam, Detector as Det
import numpy as np;
from scipy.ndimage     import convolve1d     as CONV
from scipy.ndimage     import gaussian_filter1d     as GAUSS
from scipy.special     import erf         as ERF

#double y(double omega, double beta, double theta, double Eei) {
#    return omega/Eei*(1+beta)/(1-beta*TMath::Cos(theta)+omega*(1+TMath::Cos(theta))/Eei);
#}

class xsec: # DIFFERENTIAL CROSS SECTIONS FOR MONTE-CARLO (dσ / du dφ) +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
  def __init__(self): # x is y_c, y is φ
    r = 'x/([0]*(1-x))'
    A = '1-x+1/(1-x)'
    F0 = '({}-4*{}*(1-{}))/[0]'.format(A,r,r)
    FL = '+[1]*4*{}*(1-{})*cos(2*(y-[2]))/[0]'.format(r,r)
    FCz = '+[3]*[6]*{}*(2-x)*(1-2*{})'.format(r,r)
    FCT = '+[3]*[5]*2*{}*sqrt((1+[0])*x*([0]/(1+[0])-x))*cos(y-[5])/[0]'.format(r)
    XS = '{}'.format(F0+FL+FCz+FCT); #    print(XS)
    self.fXS = ROOT.TF2('fXS', XS)
    self.fXS.SetRange( 0, 0, Beam.ymax, 2*Constants.π)
    self.fXS.SetParameter(0, Beam.xc)
    self.fXS.SetParameter(1, Laser.PL)
    self.fXS.SetParameter(2, Laser.φL)
    self.fXS.SetParameter(3, Laser.PC)
    self.fXS.SetParameter(4, Beam.PT)
    self.fXS.SetParameter(5, Beam.φe)
    self.fXS.SetParameter(6, Beam.Pz)
    self.fXS.SetNpx(1000); self.fXS.SetNpy(3600)

class Photons: # THIS IS THE CROSS SECTION FOR SCATTERED PHOTONS (dσ / dηx dηy) +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
  def __init__(self, setup):
    self.xc    = setup.xc
    self.ymax  = setup.ymax
    self.Eo    = setup.Eo
    self.xs    = [None for i in range(2)]
    self.Eg    =  None

  def Components(self):
    y  = self.Eg/self.Eo                          # y parameter
    r = y/(self.xc*(1-y))
    A = 1-y+1/(1-y)
    self.xs[0] = (y<self.ymax)*(A-4*r*(1-r))/self.xc # unpolarized xSection
    self.xs[1] = (y<self.ymax)*r*(2-y)*(1-2*r)         # longitudinal electron polarization

  def Total(self, parameters):
    Pz = parameters
    self.xst   = 1e5*(self.xs[0] + Pz*self.xs[1])
    
class pCALO: # THIS IS THE CALORIMETER FOR SCATTERED PHOTONS  +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
  def __init__(self, Det, setup):
    self.iteration, self.position, self.compton, self.resolution = 0, [], [], []
    self.E_npix = Det.E_npix
    self.E_min = Det.E_min
    self.E_pix = Det.E_pix
    self.Em = 1;
    self.E_npix *= self.Em
    self.E_pix = self.E_pix/self.Em
    self.xmid_n = np.zeros(self.E_npix, dtype=np.double)
    self.convol = np.zeros(self.E_npix, dtype=np.double)
    self.XS     = Photons(setup)
    for epix in range(self.E_npix): self.xmid_n[epix] = self.E_min + (epix+0.5)*self.E_pix
    self.XS.Eg = self.xmid_n
    self.XS.Components()
    
  def __call__(self, x, p):
    compton, resolution, change = p[0], [p[i] for i in (1,2,3)], False
    if (self.compton != compton) or change: # if Compton parameters changed
      self.compton = compton;            change = True
      self.XS.Total(compton)
    if (self.resolution != resolution) or change: # if resolution parameters changed
      self.resolution = resolution;            change = True
      sigma = self.xmid_n*np.sqrt(p[1]**2/(self.xmid_n/1e9)+(p[2]/(self.xmid_n/1e9))**2+p[3]**2)
      for n in range(self.E_npix):
        self.convol[n] = 0;
        for k in range(self.E_npix):
          # fails in the first bins where the reoslution is small wrt to the bin size.
          #gausterm = np.exp(-((self.xmid_n[n]-self.xmid_n[k])/sigma[k])**2/2)/np.sqrt(2*Constants.π*sigma[k]**2)
          # assumes that the sigma is constant over the bin size.
          gausterm = 0.5*(ERF((self.xmid_n[n]-self.xmid_n[k]-self.E_pix/2.)/sigma[k]/np.sqrt(2.))-ERF((self.xmid_n[n]-self.xmid_n[k]+self.E_pix/2.)/sigma[k]/np.sqrt(2.)))
          self.convol[n] += self.XS.xst[k]*gausterm
    if (change):
      self.iteration += 1;
      if np.random.rand()<0.001:  print('iteration: %d ' % (self.iteration))
    res = 0
    for ik in range(self.Em):
      #print(int((x[0]-self.E_min)/self.E_pix)+ik-int(self.Em/2.))
      res += self.convol[int((x[0]-self.E_min)/self.E_pix)+ik-int(self.Em/2.)]
    #res = self.convol[int((x[0]-self.E_min)/self.E_pix)]
    return p[4]*res


class Electrons: # THIS IS THE CROSS SECTION FOR SCATTERED ELECTRONS  (dσ / dx dy) +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

  def __init__(self, setup):
    self.xc    = setup.xc
    self.ymax  = setup.ymax
    self.θo    = setup.θo
    self.spec_L = setup.spec_L
    self.Eo    = setup.Eo
    self.xs    = [None for i in range(2)]
    self.x    =  None

  def Components(self):
    u = self.x/(1000*self.spec_L*ROOT.TMath.Tan(self.θo))
    Ee = self.Eo/(1+u)
    y  = 1-Ee/self.Eo                          # y parameter
    r = y/(self.xc*(1-y))
    A = 1-y+1/(1-y)
    self.xs[0] = (y<self.ymax)*(y>0)*(u>0)*(A-4*r*(1-r))/self.xc*(Ee/self.Eo)**2 # unpolarized xSection
    self.xs[1] = (y<self.ymax)*(y>0)*(u>0)*r*(2-y)*(1-2*r)*(Ee/self.Eo)**2     # longitudinal electron polarization
    
  def Total(self, parameters):
    Pz = parameters
    self.xst   = 1e5*(self.xs[0] + Pz*self.xs[1])/(1000*self.spec_L*ROOT.TMath.Tan(self.θo))

class eCherenkov: # THIS IS THE CHERENKOV FOR SCATTERED ELECTRONS  +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
  def __init__(self, Det, setup):
    self.iteration, self.ellipse, self.compton, self.resolution = 0, [], [], []
    self.X_npix, self.X_pix, self.X_beam = Det.X_npix, Det.X_pix, Det.X_beam
    self.Xm = 10;
    self.X_npix *= self.Xm
    self.X_pix = self.X_pix/self.Xm
    self.pixcx       = np.zeros( self.X_npix,                 dtype = np.double) # x pixel centers
    self.convol       = np.zeros( self.X_npix,                 dtype = np.double) # x pixel centers
    self.XS          = Electrons(setup)
          
  def __call__(self,x,p):
    # p are parameters of the function
    # x is variable (position in x and y)
    # used to fill a ROOT.TF1
    X0, σX, compton, change = p[0], p[1], p[2], False
    if self.ellipse != X0: # if geometry parameters changed
      self.ellipse = X0; change = True
      for xpix in range(self.X_npix): self.pixcx[xpix] = self.X_beam + X0 + (xpix+0.5)*self.X_pix
      self.XS.x = self.pixcx
      self.XS.Components()
    if self.compton != compton or change: # if Compton parameters changed
      self.compton = Pz = compton; change = True
      self.XS.Total(compton)
    if self.resolution != σX or change: # if Compton parameters changed
      self.resolution = σX; change = True
      self.convol = GAUSS(self.XS.xst,sigma=(σX/self.X_pix),truncate=6.0, mode='nearest')
    if change:
      self.iteration += 1;
      if np.random.rand()<0.005:
        print('iteration: %d ' % (self.iteration))
    res = 0
    for ik in range(self.Xm):
      res += self.convol[int((x[0]-self.X_beam)/self.X_pix)+ik-int(self.Xm/2.)]
    #if (index>-1 & index<self.X_npix):
    #  res = self.XS.xst[int((x[0]-self.X_beam)/self.X_pix)]
    return p[3]*res
