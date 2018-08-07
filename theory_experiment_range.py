# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 11:55:17 2018

@author: Scott
"""
import math
import numpy as np
import matplotlib.pyplot as plt

# Initialize Parameters
ni = 1.515  # RI of substrate
n1 = 1.366  # Average RI in sample
nm = 1.429  # RI of moving material (nucleosome)
NAi = 0.55  # Illimunation NA
NAc = 1.4  # Collection NA
gamma = 0.052232  # Fresnel coefficent
k = n1 * 10000000 * 2*math.pi/550  # Wavenumber in medium

phi = np.arange(0, 0.5, 0.001)  # Volume Fraction of moving material
lc = np.arange(0, 300, 0.1) * 0.0000001  # Correlation Lengths
phi = phi[None, :] #Add 2nd dimention for later operations
lc = lc[None, :] #Add 2nd dimention for later operations


#sigma_low = 0.00010889 * math.pi
#sigma_high = 0.00021542 * math.pi
sigma_low = 0.00010889
sigma_high = 0.00021542

prefactor =  ni*16* math.pi**2 * gamma**2 * (NAc/NAi)**2 * k**3 * ((nm-n1)/n1)**2
prefactor *= (2*math.pi)**(-5/2)

sigma_sq = prefactor * np.transpose(lc**3) * (phi * (1-phi))
sigma = sigma_sq**0.5

sigma_bounds = ((sigma < sigma_high) & (sigma > sigma_low)).choose(0, 1)

extra_factor = 1 / ((1 + 4 * k**2 * lc**2) * (1 + k**2 * lc**2 * (4 + NAi**2)))
# extra_factor = 1;


sigma = (sigma_sq * np.transpose(extra_factor))**0.5
sigma_bounds2 = ((sigma < sigma_high) & (sigma > sigma_low)).choose(0, 1)

plt.figure()
plt.subplot(221)
plt.imshow(sigma_bounds, aspect='auto', extent=[0, 0.5, 300, 0])
plt.xlabel('Phi')
plt.ylabel('Lc (nm)')
plt.title('Assumption: k*lc << 1')

plt.subplot(222)
plt.imshow(sigma_bounds2, aspect='auto', extent=[0, 0.5, 300, 0])
plt.xlabel('Phi')
plt.ylabel('Lc (nm)')
plt.title('No Assumptions')

plt.subplot(223)
plt.imshow(sigma_bounds, aspect='auto', extent=[0, 0.5, 300, 0])
plt.xlabel('Phi')
plt.ylabel('Lc (nm)')
plt.ylim(20, 0)
plt.title('Zoomed In, k*lc << 1')

plt.subplot(224)
plt.imshow(sigma_bounds2, aspect='auto', extent=[0, 0.5, 300, 0])
plt.xlabel('Phi')
plt.ylabel('Lc (nm)')
plt.ylim(20, 0)
plt.title('Zoomed In, No Assumption')
