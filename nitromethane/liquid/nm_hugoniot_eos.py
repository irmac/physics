#!/usr/bin/env python
from matplotlib import pylab
from pylab import *
from scipy import *


#This program calculates the P-T Hugoniot curve. The
#input is the volume and temp dependent specific heat
#capacity, the hugoniot PV data, and PV data from an
#isothermal compression.


#Iain MacLeod

#Units used: [M]: grams
#            [L]: cm
#            [T]: 10^-5 sec
#
#   [E] = [M][L]^2  =10^3 J
#            [T]^2
#
# so 1J = 10^-3 of our Energy units
# This is the reason for the 1E-3
# which appears in the specific
# heats, since they are computed
# in terms of J/gK


#Temp and specific vol at atmospheric pressure

T_0 = 298.0 #K     #Temp of isotherm
V_0 = 0.8913#g/cc  #specific Volume of nitromethane
                   #at 298K and 1 atm

#The final Compressed vol
V_end = 0.3913
dV = 0.01
N=int((V_0-V_end)/dV)

rho_0= float(1.0/V_0)
#constant appearing in Cv
m2 = 1E-3*1.714

#Isothermal PV Curve data from Murnaghan Fit
B_T0 = 1.32 #GPa
n    = 7.144

a1   = 1.248
a2   = 1.637
a3   = 5.117

#Sonic Velocity in NM (a la Lysne and Hardesty JCP 59:12:6512)
c0 = 1.318 #cm/10^-5 sec

#Shock Velocity
def U_s(UP): #U_s=shock wave velocity, U_p = particle velocity
    return float(1.248+(1-1.248)*exp(-1*5.117*UP) + 1.637*UP)

def V1(UP,target_vol):
	return V_0*( 1 - UP/(1.248 +(1-1.248)*exp(-1*5.117*UP) +1.637*UP)) - target_vol

#Volume and Temp Dependent Specific Heat Capacity
def Cv_T(v,T):
    mu  = float(V_0/v - 1)
    theta = 2326.0 + 102.2*mu
    x = theta/T
    return 1E-3*(1.146 +1.714*pow(x,2)*exp(x)/pow((exp(x)-1) ,2 ) )


#Volume Dedendent Specific Heat Capacity evaluated along T=298 isotherm
def Cv_T0(v):
    mu = float(V_0/v - 1)
    theta = 2326.0 + 102.2*mu
    x=theta/T_0
    return 1E-3*(1.146 +1.714*pow(x,2)*exp(x)/pow( (exp(x) -1) ,2 ) )

#Volume dependent Eintein temperature
def theta(v):
    mu = float(V_0/v - 1)
    return float(2326.0 + 102.2*mu)

#Volume derivative of Einstein temperature
def dtheta(v):
    return -1*102.2*V_0/(v*v)

#Build some arrays
P_H=zeros(N,typecode='f')
P_H2=zeros(N,typecode='f')
P_T0=zeros(N,typecode='f')
V=zeros(N,typecode='f')
V2=zeros(N,typecode='f')
UP=zeros(N,typecode='f')
PVh=open('PV_hugoniot.dat','w')
PVT0=open('PV_T0.dat','w')
dP=open('dP.dat','w')

for i in range(0,N):
   V[i] = V_0 - i*dV
   if i==0:
       UP[i]= optimize.fsolve(V1,0,args=(V[i]))
   else:
       UP[i]= optimize.fsolve(V1,UP[i-1],args=(V[i]))

   P_H[i]  = float((UP[i] * U_s(UP[i])/V_0)*c0*c0)
   P_T0[i] = float(B_T0/n * ( pow((V_0/V[i]),n) - 1 ))
   PVh.write("%f\t%f\n" %(V[i],P_H[i]))
   PVT0.write("%f\t%f\n" %(V[i],P_T0[i]))
   dP.write("%f\t%f\n" %(V[i],P_H[i]-P_T0[i]))


#    V[i]    = float(V_0*(1 - Up_c[i]/U_s(Up_c[i])))
#
#
#    PVh.write("%f\t%f\n" %(V[i],P_H[i]))

dP.close()

P_T0[0]=0.0
P_H[0]=0.0

#Compute dP/dV for the hugoniont curve.
#This is done by fitting a spline to
#the curve, and computing the derivative
#of the spline.

#We need to get the X data in accending
#order for the spline routine to work
Pfit=zeros(N,typecode='f')
Pfit1=zeros(N,typecode='f')

Vfit=zeros(N,typecode='f')
Vfit1=zeros(N,typecode='f')

dPdV1=zeros(N,typecode='f')
dPdV=zeros(N,typecode='f')
for i in range(0,N):
    Pfit[i]=P_H[N-1-i]
    Vfit[i]=V[N-1-i]
for i in range(0,N):
    Vfit1[i]=Vfit[N-1-i]
    Pfit1[i]=Pfit[N-1-i]

    print "%f\t%f"%(Pfit1[i],Vfit1[i])




#Ok, lets fit
#Cubic Spline Fit used.
tck = interpolate.splrep(Vfit,Pfit,s=0)
tck1 = interpolate.splrep(Pfit1,Vfit1,s=0)
print " \n\n"
print tck1

pnew = arange(0,21,1)
print pnew
Vnew = interpolate.splev(pnew,tck1,der=0)

dPdV1 = interpolate.splev(Vfit,tck,der=1)

#figure(1)
#plot(Vfit,dpdT)
#show()

#reverse the storage order of data points
#because we are looking at a compression
#process
f=open("dpdt.dat","w")
g=open("PV2.dat","w")
for i in range(0,N):
    dPdV[i]=dPdV1[N-1-i]
    f.write("%f\t%f\n"%(V[i],dPdV[i]))
    if(i<21):
	    g.write("%f\t%f\n"%(pnew[i],Vnew[i]))
f.close()
g.close()


#################################################
#     ccc    c    c    ccc     ccccccc
#    c      c c   c   c           c
#    c     ccccc  c   c           c
#     ccc c     c ccc  ccc        c
#
#For every volume point, compute b(v,T)=(dP/dT)_v
#This requires knowledge of the shock temperature.
#b(v,T) and T itself is then computed iteratively
#to self consistency.
#
#The steps are
#For each volume point:

############################################################################
#
#         / C(v,T) \                   ____________
# T_in -->          -->Calc T_out --> [T_out == T_in]-YES-->SAVE T and b(v,T)
# ^       \ b(v,T) /                         |
# |                                          NO
# |                                          |
# |                                      T_in=T_out
# |                                          |
# --------------------------------------------
############################################################################

#Compute Shock Temperature
Shock_T=zeros(N,typecode='f')
BVT=zeros(N,typecode='f')
CVT=zeros(N,typecode='f')
CVT0=zeros(N,typecode='f')
F=zeros(N,typecode='f')
G=zeros(N,typecode='f')
IF=zeros(N,typecode='f')
T2=zeros(N,typecode='f')
T1=zeros(N,typecode='f')
T3=zeros(N,typecode='f')
GAMMA=zeros(N,typecode='f')
#Initial Values
X0          = theta(V_0)/T_0
CVT[0]      = Cv_T(V_0,T_0)
CVT0[0]     = Cv_T0(V_0)
IF[0]       = 1.0
T3[0]= 0
T2[0]= T_0
T1[0]= T_0
bvt0 = 0.001637  + m2*(dtheta(V_0)/theta(V_0))*X0*X0*exp(X0)/pow((exp(X0)-1),2)#GPa/K
BVT[0] = bvt0   -m2*(dtheta(V_0)/theta(V_0))*X0*X0*exp(X0)/pow((exp(X0)-1),2)

print BVT[0]
print m2*(dtheta(V_0)/theta(V_0))*X0*X0*exp(X0)/pow((exp(X0)-1),2)

GAMMA[0]=BVT[0]/CVT[0]
#print"||"
#print X0
#print exp(X0)-1
Shock_T[0]  = T_0

print "CVT\tCVT0\tIF\tBVT\n"
#print("%f\t%f\t%f\t%f\n"%(CVT[0],CVT0[0],IF[0],BVT[0]))
f=open('cv.dat','w')

# F(v)=1/2*(p + (v0-v)*dPdV)
for i in range(0,N):
    F[i]= 0.5*(P_H[i] + (V_0-V[i])*(dPdV[i]))
    G[i]=exp(V[i]*0.001637)
    f.write("%f\t%f\n"%(V[i],Cv_T(V[i],298)))
f.close
#loop over volumes
for i in range(1,N):

    #if i == 0:
    T_in = Shock_T[i-1] + 1
    T_out=0.0

    #length of the neccessary arrays
    LEN=i+1

    cvT=zeros(LEN,typecode='f')
    cvT0=zeros(LEN,typecode='f')
    PV_H=zeros(LEN,typecode='f')
    PV_T0=zeros(LEN,typecode='f')
    bvt=zeros(LEN,typecode='f')
    v=zeros(LEN,typecode='f')
#    print len(cvT)
    intf=zeros(LEN,typecode='f')
    P=zeros(LEN,typecode='f')
    Q=zeros(LEN,typecode='f')
    Z=zeros(LEN,typecode='f')
    #fill the arrays with values already known,
    #ie the i-1 values already computed.
    for j in range(0,LEN):
        v[j]=V[j]

        cvT[j]  = CVT[j]
        cvT0[j] = Cv_T0(V[j])
        PV_H[j] = P_H[j]
        PV_T0[j]= P_T0[j]
        bvt[j]  = BVT[j]
        intf[j] = IF[j]


    #Self Consistency loop for T and BVT
    while abs(T_in - T_out) > 0.01:
        if(T_out != 0.0):
            T_in = T_out
        #build new bvt for the ith volume point
        k=LEN-1
        A  = (PV_H[k] - PV_T0[k])/(T_in - T_0)
        #print A
        B  = m2*dtheta(V[k])/theta(V[k])
       # print B
        C  = theta(V[k])/(T_in - T_0)
        x  = theta(V[k])/T_in
        x0 = theta(V[k])/T_0
       # print x0
        D  =( 1.0/(exp(x)-1) - 1.0/(exp(x0)-1) )
        E  = x0*x0*exp(x0)/pow((exp(x0) - 1 ),2)
       # print E
        cvT[k]= Cv_T(V[k],T_in)
        bvt[k]=  A + B*(C*D - E)  + (dtheta(V[k])/theta(V[k]))*(cvT0[k] - cvT[k])
        #MPa/K

        #Now Build the integrating factor value for this volume

        for j in range(0,LEN):
            P[j]   = bvt[j]/((cvT[j]))
#            print "%d\t%f\t%f\t%f"%(j,bvt[j],cvT[j],P[j])

        I1=integrate.trapz(P,v)
        intf[k] = exp(I1)


        print intf[i]
        #Now construct the F(v)*intfact(v)/Cv integrand
        for kk in range(0,LEN):
            Q[kk] = F[kk]*intf[kk]/cvT[kk]
            Z[kk]=G[kk]*F[kk]


        integral2=integrate.trapz(Q,v)
        z= (1.0/intf[i]) * integral2

        T_i = T_0/intf[i] + (1.0/intf[i])*integral2

#        for j in range(0,LEN):
#            print "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f"%(j,bvt[j],cvT[j],P[j],Q[j],F[k],intf[j],z,Shock_T[j])
#        print "===============================================================================================================================\n"
#        print T_i
        T_out = T_i


    #T2[i]=T_0*exp(1.637*(V_0-V[i])) + exp(-1.0*1.637*V[i]/1.2)*integrate.trapz(Z,v)
    T2[i]= (dtheta(V[k])/theta(V[k]))*(cvT0[k] - cvT[k])
    #Store the SelfCons calcd quantities
    BVT[i]=bvt[i]
    Shock_T[i]= T_i
    IF[i]=intf[i]
    CVT[i]=cvT[i]
    GAMMA[i]=BVT[i]/CVT[i]
    T1[i]=T_0/intf[i]
    T3[i]=A
    #print T_i



PTfit=interpolate.splrep(P_H,Shock_T,s=0)
pnew = arange(0,21,1)
tnew = interpolate.splev(pnew,PTfit,der=0)





hh=open("PT.dat","w")
for i in range(0,21):
	#write T(P) data to file
	hh.write("%f\t%f\n"%(pnew[i],tnew[20-i]))
	#obain T(V) data




hh.close()
f=open("TS.dat","w")
g=open("BVT.dat","w")
h=open("PP.dat","w")
o=open("GAMMA.dat","w")
p=open("cvt.dat","w")
q=open("QQ.dat","w")
for i in range(0,N):
    f.write("%f\t%f\t%f\t%f\t%f\t%f\n"%(V[i],Shock_T[i],T2[i],T1[i],T3[i],IF[i]))
    g.write("%f\t%f\n"%(V[i],BVT[i]))
    h.write("%f\t%f\n"%(P_H[i],Shock_T[i]))
    o.write("%f\t%f\n"%(V[i],GAMMA[i]))
    p.write("%f\t%f\t%f\n"%(V[i],CVT[i],Cv_T0(V[i])))
    q.write("%f\t%f\n"%(V[i],Q[i]))

f.close()
g.close()
h.close()
o.close()
p.close()
q.close()

