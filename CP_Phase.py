import numpy as np
from cmath import exp,sin,cos,pi
import matplotlib.pyplot as plt

m21=complex(0,7.42*pow(10,-5))
m31=complex(0,2.51*pow(10,-3))
O12=0.58381263
O13=0.1504474
O23=0.7347836
t=295

pe1=[]
pe2=[]
pe3=[]

def Mat(A,B,C):
    D=np.matmul(A,B)
    return(np.matmul(D,C))

def U(dO,E):
    phi=complex(0,dO)
    H=np.matrix([[1, 0, 0],
                [0, exp(-(2*m21*1.27*t)/E), 0],
                [0, 0, exp(-(2*m31*1.27*t)/E)]])


    U1=np.matrix([[cos(O12)*cos(O13), sin(O12)*cos(O13), sin(O13)*exp(-phi)],
                [-sin(O12)*cos(O23)-cos(O12)*sin(O13)*sin(O23)*exp(phi), cos(O12)*cos(O23)-sin(O12)*sin(O13)*sin(O23)*exp(phi), cos(O13)*sin(O23)],
                [sin(O12)*sin(O23)-cos(O12)*sin(O13)*cos(O23)*exp(phi), -cos(O12)*sin(O23)-sin(O12)*sin(O13)*cos(O23)*exp(phi), cos(O13)*cos(O23)]])

    U1d=U1.getH()
    return(Mat(U1,H,U1d))
    
n=1000 #data points
for E in range (1,n+1,1):
    state=np.matrix([[0, 0, 0],
                    [0, 1, 0],
                    [0, 0, 0]])

    X=np.matmul(U(0,E/n),state)
    fstate1=np.matmul(X,U(0,E/n).getH())
    Y=np.matmul(U(pi/2,E/n),state)
    fstate2=np.matmul(Y,U(pi/2,E/n).getH())
    Z=np.matmul(U(pi,E/n),state)
    fstate3=np.matmul(Z,U(pi,E/n).getH())



    pe1.append(fstate1[0,0])
    pe2.append(fstate2[0,0])
    pe3.append(fstate3[0,0])

    
dis=np.linspace(0.001,1,n)
plt.semilogx([dis],[pe1])
plt.semilogx([dis],[pe2])
plt.semilogx([dis],[pe3])
plt.xlim(0.075,1)
plt.ylim(0,0.2)
plt.plot(dis,pe1,color='g',label='CP phase=0')
plt.plot(dis,pe2,color='r',label='CP phase=pi/2')
plt.plot(dis,pe3,color='b',label='CP phasepi')
plt.xlabel("Energy (E)")
plt.ylabel("Probability")
plt.title("Muon to Electron")
plt.legend()
plt.savefig('Cp_violation.jpg')
plt.show()



