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
n1 = 1.46  # Average RI in sample
nm = 1.5959  # RI of moving material (nucleosome)
NAi = 0.55  # Illimunation NA
NAc = 1.4  # Collection NA
gamma = 0.094957  # Fresnel coefficent
k = n1 * 10000000 * 2*math.pi/550  # Wavenumber in medium

#phi = np.arange(0, 0.01, 0.0001)  # Volume Fraction of moving material
phi = 0.003
lc = np.arange(0, 400, 0.1) * 0.0000001  # Correlation Lengths

#phi = phi[None, :]  #Add 2nd dimention for later operations
lc = lc[None, :]  #Add 2nd dimention for later operations

sigma_low = 0.0000325318
sigma_high = 0.0000906242


prefactor = ni * 16 * math.pi**2 * gamma**2 * (NAc/NAi)**2 * k**3 * ((nm-n1)/n1)**2
prefactor *= (2*math.pi)**(-5/2)

sigma_sq = prefactor * np.transpose(lc**3) * (phi * (1-phi))
sigma = sigma_sq**0.5

sigma_bounds = ((sigma < sigma_high) & (sigma > sigma_low)).choose(0, 1)

extra_factor = 1 / ((1 + 4 * k**2 * lc**2) * (1 + k**2 * lc**2 * (4 + NAi**2)))
# extra_factor = 1;

#sig_n_sq = ((nm-n1)/n1)**2 * (phi * (1-phi))

sigma2 = (sigma_sq * np.transpose(extra_factor))**0.5
#sigma_L = gamma**2 * sig_n_sq * 0.25 * (1 - 1/(1+(k*np.transpose(lc)*NAc)**2)**0.5)
sigma_bounds2 = ((sigma < sigma_high) & (sigma > sigma_low)).choose(0, 1)

#plt.imshow(sigma,  aspect='auto', extent=[0.0001, 0.01, 300, 0])

# Comparing sigma vs lc at set phi
plt.figure()
plt.subplot(121)
plt.plot(np.transpose(lc)/0.0000001, sigma)
plt.xlabel('Lc (nm)')
plt.ylabel('Sigma_t')
plt.title('Assumption: k*lc << 1')
plt.subplot(122)
plt.plot(np.transpose(lc)/0.0000001, sigma2)
plt.xlabel('Lc (nm)')
plt.ylabel('Sigma_t')
plt.title('No Assumptions')
#plt.subplot(133)
#plt.plot(np.transpose(lc)/0.0000001, (sigma2**2 + sigma_L)**0.5)
#plt.xlabel('Lc (nm)')
#plt.ylabel('Sigma_t')
#plt.title('No Assumptions')



plt.figure()
plt.subplot(221)
plt.imshow(sigma_bounds, aspect='auto', extent=[0.0001, 0.01, 100, 0])
plt.xlabel('Phi')
plt.ylabel('Lc (nm)')
plt.title('Assumption: k*lc << 1')

plt.subplot(222)
plt.imshow(sigma_bounds2, aspect='auto', extent=[0.0001, 0.01, 100, 0])
plt.xlabel('Phi')
plt.ylabel('Lc (nm)')
plt.title('No Assumptions')

plt.subplot(223)
plt.imshow(sigma_bounds, aspect='auto', extent=[0.0001, 0.01, 100, 0])
plt.xlabel('Phi')
plt.ylabel('Lc (nm)')
plt.ylim(20, 0)
plt.title('Zoomed In, k*lc << 1')

plt.subplot(224)
plt.imshow(sigma_bounds2, aspect='auto', extent=[0.0001, 0.01, 100, 0])
plt.xlabel('Phi')
plt.ylabel('Lc (nm)')
plt.ylim(20, 0)
plt.title('Zoomed In, No Assumption')
