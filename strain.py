#!/usr/bin/env python3

#This script plots the rocking curve of a crystal with constant strain gradient

from scipy.special import pbdv
from math import sqrt

blah=1.
def sign(exp):
	return 1
plusminus = 1

#Independent Variables
c=blah
eta = blah
k=blah
chi_0 = blah
chi_h = blah
chi_h_bar = blah
gamma_0 = blah
gamma_h = blah
gamma = blah

#Calculated Values
beta_h = blah
B = blah
D_1v = blah
D_v = blah
D_ratio=abs(D_1v/D_v)**2

#Function argument q_0
epsilon_n = (chi_0 +sqrt(abs(gamma)chi_h*chi_h_bar)*abs(c))*(-eta+plusminus*sqrt(eta**2+sign(gamma))))/2
alpha = 2*epsilon_n
q_0 =(k/2)*((chi_0/gamma_0)+(chi_0-alpha)/abs(gamma_h))



#Rocking Curve Function
def R_n(q_0):
	return 1
