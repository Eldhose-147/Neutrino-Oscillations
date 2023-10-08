import numpy as np
from cmath import exp,sin,cos
import matplotlib.pyplot as plt

m21=complex(0,7.5*pow(10,-5))
m31=complex(0,2.5*pow(10,-3))
O12=0.58381263
O13=0.1504474
O23=0.7347836
pe=[]
pm=[]
pu=[]
cp_phi=0
E=1
def U(dO,t):
    phi=complex(0,dO)
    H=np.matrix([[1, 0, 0],
                [0, exp(-(2*m21*1.27*t)/(E)), 0],
                [0, 0, exp(-(2*m31*1.27*t)/(E))]])


    U1=np.matrix([[cos(O12)*cos(O13), sin(O12)*cos(O13), sin(O13)*exp(-phi)],
                [-sin(O12)*cos(O23)-cos(O12)*sin(O13)*sin(O23)*exp(phi), cos(O12)*cos(O23)-sin(O12)*sin(O13)*sin(O23)*exp(phi), cos(O13)*sin(O23)],
                [sin(O12)*sin(O23)-cos(O12)*sin(O13)*cos(O23)*exp(phi), -cos(O12)*sin(O23)-sin(O12)*sin(O13)*cos(O23)*exp(phi), cos(O13)*cos(O23)]])

    U1d=U1.getH()
    A=np.matmul(U1,H)
    return(np.matmul(A,U1d))
    
for t in range (0,1200,10):
    state=np.matrix([[0, 0, 0],
                    [0, 1, 0],
                    [0, 0, 0]])

    B=np.matmul(U(cp_phi,t),state)
    fstate=np.matmul(B,U(cp_phi,t).getH())
    pe.append(fstate[0,0])
    pm.append(fstate[1,1])
    pu.append(fstate[2,2])


dis=np.arange(0,1200,10)
plt.plot(dis,pe,color='g',label='electron neutrino')
plt.plot(dis,pm,color='r',label='muon neutrino')
plt.plot(dis,pu,color='b',label='tau neutrino')
plt.xlabel('L/E (km/GeV)')
plt.ylabel('Oscillation Probability')
plt.show()