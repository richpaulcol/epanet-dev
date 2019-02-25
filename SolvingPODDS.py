import numpy as np
import pylab as pp

rho = 998
g = 9.81
mu = 1e-3

L = 500.			#m
roughness = 0.075		#mm
d = 0.25			#m
A = np.pi*d**2 / 4

H1 = 20				#m
dt = 60
times = np.arange(0,3*60*60,dt)
Q = np.ones(times.shape)*5	#lps
Q[10:] = 12   			#lps

Q = np.load('Flow.npy')
times = np.arange(0,Q.size*dt,dt)
Q = Q/1000.
V = Q/A

Re = rho * V * d / mu
lam = 0.25 / (np.log10((roughness/1000.) / (3.7*d) + 5.74/(Re**0.9))**2) 
hl = lam*L*V**2 / (2*g*d)

H2 = H1- hl

sf = hl/L
tau_a = rho * g * d * sf / 4

tau_init  = 0.157


dx = dt * max(V) 			#Spatial Discretisation stepsize (m)
x = np.arange(0,L,dx)	#Spatial Discretisation 
TMass_old = np.ones(x.shape)*1e-6	#Initial Turbidity (kg)
TMass_new = np.ones(x.shape)*1e-6
#####   Old School PODDS version


C_max = 100.
P = 0.005
k = -3.

C_init = tau_init *k + C_max

C = C_init

dN_PODDS_old = []
dC_PODDS_old = []
C_PODDS_old = []
tau_c_PODDS_old = []
T_end_old = []
for i in range(times.size):
	tau_c = (C - C_max) / k
	R = P/Q[i] *(tau_a[i] - tau_c)
	
	if R < 0:
		R = 0
	
	As = 4 * Q[i] * dt / d
	dN = R* As
	dC = R* dt*Q[i]
	
	C -= dC
	
	dN_PODDS_old.append(dN)
	dC_PODDS_old.append(dC)
	C_PODDS_old.append(C)
	tau_c_PODDS_old.append(tau_c)
	
	PODDSMass = dN
	lam = ( V[i]*dt / dx)
	#if Q>0:
	TMass_old[1:] = lam*TMass_old[:-1] + (1-lam)*TMass_old[1:] + PODDSMass
	TMass_old[0] = PODDSMass
	
	T_end_old.append(TMass_old[-1])
	
##### PODDS proper
alpha = -k
beta = P / (-k)

dN_PODDS_new = []
tau_c_PODDS_new = []
tau_c = tau_init
T_end_new = []
for i in range(times.size):
	
	excess = tau_a[i] - tau_c
	if excess <=0:
		excess = 0
	dtau_c = beta * (excess) *dt
	dN = (4 / d) * alpha * beta * (excess) * dt
	tau_c += dtau_c
	
	dN_PODDS_new.append(dN)
	tau_c_PODDS_new.append(tau_c)
	
	PODDSMass = dN
	lam = ( V[i]*dt / dx)
	#if Q>0:
	TMass_new[1:] = lam*TMass_new[:-1] + (1-lam)*TMass_new[1:] + PODDSMass
	TMass_new[0] = PODDSMass
	
	T_end_new.append(TMass_new[-1])
	
	
#pp.figure()
#pp.plot(times,dN_PODDS_old)
#pp.plot(times,dN_PODDS_new)
#pp.show()

pp.figure()
pp.plot(times/3600.,T_end_new)
pp.plot(times/3600.,T_end_old)
pp.show()

