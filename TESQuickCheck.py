
# coding: utf-8

# In[9]:

import numpy as np
import matplotlib.pyplot as plt


# In[19]:

#Import parameters from bosh. Labeled as state_y.py

R_b = 201   #20180701 Cool Down
R = 0.7*10.5*10**(-3)    #bosh state_y.py don't forget factor of 0.5~1
R_sh = 385*10**(-6)   #bosh state_y.py
L = 0.4*10**(-6)   #bosh state_y.py
G450 = 30*10**(-12)   #bosh state_y.py
C450 = 1.0*10**(-12)   #bosh state_y.py
betaG = 2.2   #bosh state_y.py
betaC = 3   #bosh state_y.py
T = 0.45

Vb = 1800/(2**(16))*2.5
deltaVb = 10/(2**(16))*2.5

I_b = Vb/R_b
I = I_b/(1+R/R_sh)
I_sh = I_b/(1+R_sh/R)
V = I_sh*R_sh

deltaI_b = deltaVb/R_b
deltaI = deltaI_b/(1+R/R_sh)
deltaI_sh = deltaI_b/(1+R_sh/R)
deltaV = deltaI_sh*R_sh

P_J = I*I*R
G = G450*(T/0.45)**(betaG)
C = C450*(T/0.45)**(betaC)


# In[20]:

#parameters to fit, give them intial values at present
beta_I = 0.01
alpha_I = 200  #bosh state_y.py


# In[21]:

#once we have fitted parameters, we can work out the following time constants
L_I = alpha_I*P_J/(G*T)
tau = (C450/G450)*(T/0.45)**(betaC-betaG)

M11 = (R_sh+R*(1+beta_I))/L
M12 = L_I*G/(I*L)
M21 = -I*R*(2+beta_I)/C
M22 = (1-L_I)/tau


Invpar = M11*M22-M21*M12
InvM11 = M22/Invpar
InvM12 = -M12/Invpar
InvM21 = -M21/Invpar
InvM22 = M11/Invpar

#for bias step deltaV is not zero but deltaP_opt is zero
#Ap = -deltaI*((G*tau)/(I*R*(2+beta_I))*np.sqrt((M11-M22)**2+4*M21*M12))**(-1) #homogenous solution
#Am = -Ap #homogenous solution
lamp = (M11+M22+np.sqrt((M11-M22)**2+4*M21*M12))/2
lamm = (M11+M22-np.sqrt((M11-M22)**2+4*M21*M12))/2
vp1 = ((1-L_I-lamp*tau)/(2+beta_I))*(G/(I*R))
vm1 = ((1-L_I-lamm*tau)/(2+beta_I))*(G/(I*R))
vp2 = 1
vm2 = 1
Ap = (deltaV/L)*(vm1*InvM21-InvM11)/(vp1-vm1) #inhomogenous solution
Am = -InvM21*deltaV/L-Ap #inhomogenous solution


# In[22]:

time_range = 5*10**(-3)
t = np.arange(0, time_range, 1*10**(-6))

dI1 = Ap*np.exp(-lamp*t)*vp1
dI2 = Am*np.exp(-lamm*t)*vm1
dI = dI1+dI2+InvM11*deltaV/L


plt.plot(10**3*t, 10**6*dI1, label="time constant t+ = %.4f ms" %(1/lamp*1000))
plt.legend()
plt.xlabel("time (ms)")
plt.ylabel("I (uA)")
plt.show()

plt.plot(10**3*t, 10**6*dI2, label="time constant t- = %.4f ms" %(1/lamm*1000))
plt.legend()
plt.xlabel("time (ms)")
plt.ylabel("I (uA)")
plt.show()

plt.plot(10**3*t, 10**6*dI1, label="time constant t+ = %.4f ms" %(1/lamp*1000))
plt.plot(10**3*t, 10**6*dI2, label="time constant t- = %.4f ms" %(1/lamm*1000))
plt.plot(10**3*t, 10**6*dI, label="total")
plt.legend()
plt.xlabel("time (ms)")
plt.ylabel("I (uA)")
plt.savefig('current change.png', format='png', dpi=1000)
plt.show()


# In[23]:

time_range = 5*10**(-3)
t = np.arange(0, time_range, 1*10**(-6))

dT1 = Ap*np.exp(-lamp*t)*vp2
dT2 = Am*np.exp(-lamm*t)*vm2
dT = dT1+dT2++InvM21*deltaV/L

plt.plot(10**3*t, 10**6*dT1, label="time constant t+ = %.4f ms" %(1/lamp*1000))
plt.legend()
plt.xlabel("time (ms)")
plt.ylabel("T (uK)")
plt.show()

plt.plot(10**3*t, 10**6*dT2, label="time constant t- = %.4f ms" %(1/lamm*1000))
plt.legend()
plt.xlabel("time (ms)")
plt.ylabel("T (uK)")
plt.show()

plt.plot(10**3*t, 10**6*dT1, label="time constant t+ = %.4f ms" %(1/lamp*1000))
plt.plot(10**3*t, 10**6*dT2, label="time constant t- = %.4f ms" %(1/lamm*1000))
plt.plot(10**3*t, 10**6*dT, label="total")
plt.legend()
plt.xlabel("time (ms)")
plt.ylabel("T (uK)")
plt.savefig('temperature change.png', format='png', dpi=1000)
plt.show()


# In[ ]:



