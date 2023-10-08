import numpy as np
from cmath import exp,sin,cos
import matplotlib.pyplot as plt

m12v0=complex(0,7.5*pow(10,-5))
m13v0=complex(0,2.5*pow(10,-3))
O12v0=0.58381263
O13v0=0.1504474

m12v1=complex(0,9.36339*pow(10,-4))
m13v1=complex(0,2.4858*pow(10,-3))
O12v1=1.5341
O13v1=0.246599
O23=0.7347836

m12v2=complex(0,0.0000717308)
m13v2=complex(0,0.0024937)
O12v2=0.646648
O13v2=0.151048

E=1

pe1=[]
pe2=[]
pe3=[]
def U(t,m21,m31,O12,O13):
    phi=0
    H=np.matrix([[1, 0, 0],
                [0, exp(-(2*m21*1.27*t)/(E)), 0],
                [0, 0, exp(-(2*m31*1.27*t)/(E))]])


    U1=np.matrix([[cos(O12)*cos(O13), sin(O12)*cos(O13), sin(O13)*exp(-phi)],
                [-sin(O12)*cos(O23)-cos(O12)*sin(O13)*sin(O23)*exp(phi), cos(O12)*cos(O23)-sin(O12)*sin(O13)*sin(O23)*exp(phi), cos(O13)*sin(O23)],
                [sin(O12)*sin(O23)-cos(O12)*sin(O13)*cos(O23)*exp(phi), -cos(O12)*sin(O23)-sin(O12)*sin(O13)*cos(O23)*exp(phi), cos(O13)*cos(O23)]])

    U1d=U1.getH()
    A=np.matmul(U1,H)
    return(np.matmul(A,U1d))
    
for t in range (12000,18000,10):
    state=np.matrix([[0, 0, 0],
                    [0, 1, 0],
                    [0, 0, 0]])

    X=np.matmul(U(t,m12v0,m13v0,O12v0,O13v0),state)
    fstate1=np.matmul(X,U(t,m12v0,m13v0,O12v0,O13v0).getH())
    Y=np.matmul(U(t,m12v1,m13v1,O12v1,O13v1),state)
    fstate2=np.matmul(Y,U(t,m12v1,m13v1,O12v1,O13v1).getH())
    Z=np.matmul(U(t,m12v2,m13v2,O12v2,O13v2),state)
    fstate3=np.matmul(Z,U(t,m12v2,m13v2,O12v2,O13v2).getH())
    pe1.append(fstate1[1,1])
    pe2.append(fstate2[1,1])
    pe3.append(fstate3[1,1])


dis=np.arange(12000,18000,10)
plt.plot(dis,pe1,color='b',label='V=0')
plt.plot(dis,pe2,color='r',label='V=$10^{-3}$')
plt.plot(dis,pe3,color='g',label='V=10$^{-5}$')
plt.xlabel('L/E (km/GeV)')
plt.ylabel('Oscillation Probability')
plt.xlim(12000,18000)
plt.legend()
plt.savefig('matter_muon_muon.jpg')
plt.show()




