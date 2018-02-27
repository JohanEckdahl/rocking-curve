from math import pi, asin, cos, sin, sqrt
import matplotlib.pyplot as plt

def sign(n): 	
	return 1 if n >= 0 else -1


#---------------------#
#-- Input Variables --#
#---------------------#
e_charge = 1.602e-19
c_light = 2.99792458e8
hbar = 1.0545887e-34
keV = e_charge*1000.

E = 4.0*keV
phi_a = 0.0
d_h = 3.1355e-10
pol = 's'

chi_r0 =-1.0843e-4
chi_i0 = 1.5225e-5

abs_chi_rh = 5.7236e-05
abs_chi_ih = 1.0601e-05  
arg_crh_over_cih = -pi


#-----------------------#
#-- Calculated Values --#
#-----------------------#
lam = 2.*pi*hbar*c_light/E
theta_B = asin(0.5*lam/d_h)
psi_0 = 0.5*pi - theta_B + phi_a
psi_h = 0.5*pi + theta_B + phi_a

C_pol = 1.0 if pol == 's' else cos(2.*theta_B)
gamma_0 = cos(psi_0)
gamma_h = cos(psi_h)
gamma   = gamma_h/gamma_0

chi_rhr = abs_chi_rh * sin(0.5*arg_crh_over_cih); 
chi_rhi =-abs_chi_rh * cos(0.5*arg_crh_over_cih); 
chi_ihr =-abs_chi_ih * sin(0.5*arg_crh_over_cih);
chi_ihi =-abs_chi_ih * cos(0.5*arg_crh_over_cih);

chi_0    = chi_r0  + 1j*chi_i0 
chi_h    = chi_rhr + 1j*chi_rhi + 1j*(chi_ihr + 1j*chi_ihi)
chi_hbar = chi_rhr - 1j*chi_rhi + 1j*(chi_ihr - 1j*chi_ihi)

chih_times_chihbar = abs_chi_rh**2 - abs_chi_ih**2 + 2.0*1j*(chi_rhr*chi_ihr - chi_ihi*chi_rhi) 
abs_chih_over_chihbar = (abs_chi_rh**2 + abs_chi_ih**2 - 2.*(chi_rhr*chi_ihi - chi_rhi*chi_ihr))/(abs_chi_rh**2 + abs_chi_ih**2 + 2.*(chi_rhr*chi_ihi - chi_rhi*chi_ihr))

D_theta_0s = -0.5*chi_0*(1.0 - gamma)/sin(2.0*theta_B)


#---------------#
#-- Functions --#
#---------------#

# there is something wrong with the handling of complex numbers here
def eta(D_theta):
	return (D_theta - abs(D_theta_0s)) * sin(2.0*theta_B)/(abs(C_pol)*sqrt(abs(gamma))*sqrt(abs(chih_times_chihbar)))
	return (D_theta - D_theta_0s) * sin(2.0*theta_B)/(abs(C_pol)*sqrt(abs(gamma))*sqrt(chih_times_chihbar))

def R(D_theta):
	return abs_chih_over_chihbar * abs(eta(D_theta) - sign(eta(D_theta))*sqrt(eta(D_theta)**2 - 1.0))**2

#----------------#
#-- Print Data --#
#----------------#
print ('E:                   ' , E/e_charge, ' eV')
print ('lambda:              ' , lam*1e10, ' A')
print ('theta_B:             ' , theta_B/pi*180.)
print ('psi_0:               ' , psi_0/pi*180.)
print ('psi_h:               ' , psi_h/pi*180.)
print ('gamma_0:             ' , gamma_0)
print ('gamma_h:             ' , gamma_h)
print ('gamma:               ' , gamma)
print ('chi_h*chi_hbar:      ' , chih_times_chihbar, chi_h*chi_hbar)
print ('abs(chi_h/chi_hbar): ' , abs_chih_over_chihbar, abs(chi_h/chi_hbar))
print ('D_theta_0s:          ' , D_theta_0s)


#----------#
#-- Plot --#
#----------#

plt.plot([R(x/1000) for x in range(-1000,1000)])
plt.show()


'''
set log y


plot[-1e-3:1e-3][] R(x)

pause -1
'''
