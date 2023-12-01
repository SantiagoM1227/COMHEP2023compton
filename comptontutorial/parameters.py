import numpy as np;
# Α Β Γ Δ Ε Ζ Η Θ Ι Κ Λ Μ Ν Ξ Ο Π Ρ Σ Τ Υ Φ Χ Ψ Ω
# α β γ δ ε ζ η θ ι κ λ μ ν ξ ο π ρ σ τ υ φ χ ψ ω

class Constants:
  c      = 2.99792458e+10              # speed of light in vacuum, cm/s
  me     = 0.510998928e+6              # electron rest energy, eV/c^2
  h      = 6.62607015e-34              # Plank constant, J*s
  e      = 1.602176634e-19             # electron charge magnitude, C
  re     = 2.8179403267e-13            # classical electron radius, cm
  barn   = 1.0e-24                     # one barn, m^2
  hc     = h*c/e                       # hc, eV*cm
  π      = np.pi                       # pi
  deg    = π/180.                      # 1 deg in rads

class Laser:
  λo     = 515e-7                    # λo, cm
  ωo     = Constants.hc/λo             # ωo, eV
  PL     = 0.0                         # laser linear polarization horizontal (ξ1=1) or vertical ξ1=-1
  φL     = 0.0                         # laser linear polarization (φ=π/4 ξ2=1) or vertical (φ=-π/4 ξ2=-1)
  PC     = (1 - PL**2)**0.5            # laser circular polarization degree ξ3 [-1:1]
  θ0     = 3*Constants.deg             # beams crossing angle
  E      = 1e-6                        #
  σx     = 0.5                         # size in x, mm
  σy     = 0.5                         # size in y, mm
  σz     = 1                           # size in z, mm

class Beam:
  Eo     = 45.0e+9                      # electron (beam) energy, eV
  σE     = 1.e-3                       # beam energy spread, a.u.
  γ      = Eo/Constants.me             # electron Lorentz factor
  β      = np.sqrt(1-1/(γ**2))
  xc     = 2*γ*Laser.ωo/Constants.me*(1+β*np.cos(Laser.θ0));   # xc parameter
  ymax   = xc/(1+xc)
  θo     = .0021341                    # bending angle, rad
  bend_l = 24.12                       # dipole length, m
  spec_L = 0.5*bend_l + 88.0           # spectrometer arc length, m
  leip_L = 0.5*bend_l + spec_L + 5.0   # laser electron i.p. distance, m
  εx     = 4.5e-9                      # emittance x, m*rad
  εy     = 4.5e-11                     # emittance y, m*rad
  βx     = 96.46                       # βx, m
  βy     = 127.09                      # βy, m
  σx     = 1000*(εx*βx)**0.5           # betatron size in x, mm
  σy     = 1000*(εy*βy)**0.5           # betatron size in y, mm
  ηx     = (εx/βx)**0.5                # betatron angular spread in x, rad
  ηy     = (εy/βy)**0.5                # betatron angular spread in y, rad
  PT     = 0.                          # transverse electron spin polariztion along y
  φe     = Constants.π/2.              # laser linear polarization (φ=π/4 ξ2=1) or vertical (φ=-π/4 ξ2=-1)
  Pz     = -1.0                         # longitudinal electron spin polariztion along z
  Q      = 10e-9
  σz     = 1                           # size in z, mm

class Detector:                        # Electron Pixel Detector
  X_beam = 10.0                        # beam-detector horizontal space, mm
  X_size = 170.0                        # detector X-size (horizontal)  , mm
  X_npix = 1200                        # number of pixels in X
  X_pix  = X_size/float(X_npix)        # one pixel in X, mm
  E_npix = 200                         # number of bins for photon detector
  E_min = 0.                           # min energy for photon detector
  E_max = Beam.ymax*Beam.Eo*1.1        # max energy for photon detector
  E_pix = (E_max-E_min)/float(E_npix)  # bin size in energy units for photon detector
  resolA = 0.01                        # resolution of the photon detector
  resolB = 0.00                        # resolution of the photon detector
  resolC = 0.00                        # resolution of the photon detector
