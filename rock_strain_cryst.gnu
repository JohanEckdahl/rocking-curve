set term x11

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

lambda = 2.*pi*hbar*c_light/E
theta_B = asin(0.5*lambda/d_h)
psi_0 = 0.5*pi - theta_B + phi_a
psi_h = 0.5*pi + theta_B + phi_a

C_pol = (pol eq 's') ? 1.0 : cos(2.*theta_B)
gamma_0 = cos(psi_0)
gamma_h = cos(psi_h)
gamma   = gamma_h/gamma_0

chi_rhr = abs_chi_rh * sin(0.5*arg_crh_over_cih); 
chi_rhi =-abs_chi_rh * cos(0.5*arg_crh_over_cih); 
chi_ihr =-abs_chi_ih * sin(0.5*arg_crh_over_cih);
chi_ihi =-abs_chi_ih * cos(0.5*arg_crh_over_cih);

i = {0,1}
chi_0    = chi_r0  + i*chi_i0 
chi_h    = chi_rhr + i*chi_rhi + i*(chi_ihr + i*chi_ihi)
chi_hbar = chi_rhr - i*chi_rhi + i*(chi_ihr - i*chi_ihi)

chih_times_chihbar = abs_chi_rh**2 - abs_chi_ih**2 + 2.0*i*(chi_rhr*chi_ihr - chi_ihi*chi_rhi) 
abs_chih_over_chihbar = (abs_chi_rh**2 + abs_chi_ih**2 - 2.*(chi_rhr*chi_ihi - chi_rhi*chi_ihr))/(abs_chi_rh**2 + abs_chi_ih**2 + 2.*(chi_rhr*chi_ihi - chi_rhi*chi_ihr))

D_theta_0s = -0.5*chi_0*(1.0 - gamma)/sin(2.0*theta_B)

eta(D_theta) = (D_theta - D_theta_0s) * sin(2.0*theta_B)/(abs(C_pol)*sqrt(abs(gamma))*sqrt(chih_times_chihbar))

R(D_theta) = abs_chih_over_chihbar * abs(eta(D_theta) - sgn(eta(D_theta))*sqrt(eta(D_theta)**2 - 1.0))**2

print 'E:                   ' , E/e_charge, ' eV'
print 'lambda:              ' , lambda*1e10, ' A'
print 'theta_B:             ' , theta_B/pi*180.
print 'psi_0:               ' , psi_0/pi*180.
print 'psi_h:               ' , psi_h/pi*180.
print 'gamma_0:             ' , gamma_0
print 'gamma_h:             ' , gamma_h
print 'gamma:               ' , gamma
print 'chi_h*chi_hbar:      ' , chih_times_chihbar, chi_h*chi_hbar 
print 'abs(chi_h/chi_hbar): ' , abs_chih_over_chihbar, abs(chi_h/chi_hbar)
print 'D_theta_0s:          ' , D_theta_0s

set log y


plot[-1e-3:1e-3][] R(x)

pause -1