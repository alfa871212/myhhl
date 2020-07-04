from qiskit.circuit.library import QFT,XGate,RZZGate
from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit
from qiskit.quantum_info.states import Statevector, partial_trace
from qiskit.aqua.components.reciprocals import LookupRotation
from qiskit.quantum_info.states import state_fidelity
from qiskit.ignis.verification.tomography import state_tomography_circuits,StateTomographyFitter
from qiskit import execute,Aer
from copy import deepcopy
from qiskit.tools import job_monitor
from qiskit.extensions.simulator import snapshot
from qiskit.visualization import plot_state_city,circuit_drawer
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from qiskit.extensions.simulator import snapshot
from qiskit.providers.aer.extensions import SnapshotStatevector



def myRzz(theta):
    qr=QuantumRegister(1)
    qc=QuantumCircuit(qr)
    qc.u1(theta, qr)
    qc.x(qr)
    qc.u1(theta, qr)
    qc.x(qr)
    circuit_drawer(qc,output='mpl',filename='rzz.png')
    gate=qc.to_gate()
    gate.name=f'myRzz_{theta}'
    return gate
def Agate(numLis,inverse=False,plot_gate=False):
    if len(numLis)!=6:
        raise Exception("Parameter Lis error!")
    qr = QuantumRegister(2)
    qc = QuantumCircuit(qr)
    qc.cz(qr[0],qr[1])
    qc.rx(numLis[0],qr[1])
    sqrtXdg = XGate().power(0.5).inverse()
    qc.append(sqrtXdg,qargs=[qr[1]])
    
    qc.append(myRzz(numLis[1]).control(),qargs=qr[:])
    if plot_gate:
        qc.barrier()
    qc.rx(numLis[2],qr[0])
    qc.append(myRzz(numLis[3]).control(),qargs=[qr[1]]+[qr[0]])
    qc.cx(qr[0],qr[1])
    qc.rx(numLis[4],qr[0])
    qc.cx(qr[0],qr[1])
    qc.cz(qr[0],qr[1])
    #print(qc)
    if plot_gate:
        circuit_drawer(qc,output='mpl',filename='Agate.png')
        return
    ret = qc.to_gate()
    if inverse:
        ret.inverse()
    ret.name = f'gate_1/{numLis[5]}'
    return ret
def createQC(r):
    #parameter setting lis
    paraLis16=[ 1.57e+00,  5.24e-01, -9.92e-07,  5.24e-01, -7.89e-07,16]
    paraLis8=[ 3.93,  1.31, -0.79,  1.31,  0.79,8] 
    paraLis4=[1.18, 1.44, 1.18, 1.44, 1.96,4] 
    paraLis2=[-3.34,  0.46,  2.16,  0.46, -3.73,2] 
    #paraLis16 =[0.196,0.379,0.981,0.589,1.178,16]
    #paraLis8 =[1.963,1.115,1.963,2.615,1.178,8]
    #paraLis4 =[-0.785,1.017,3.927,2.517,2.356,4]
    #paraLis2 =[-9.014*10**(-9),-0.75,1.571,0.75,-1.571,2]
    #paraLis = [paraLis16,paraLis8,paraLis4,paraLis2]
    paraLis = [[1.5707963049051756, 0.5235986332337932, -4.82404012048593e-08, 0.5235989479710231, -1.3335823001363274e-09, 16], [3.1415773715120245, 1.0471999995628776, 1.5707796550096538, 1.0471959121890322, -1.57079788886257, 8], [-2.356204293590994, 1.3089929020124322, 5.497795996751179, 1.3090022941259296, 0.7853896298185202, 4], [1.178098234809345, -0.6544998552823452, 1.1780979112483747, -0.6544939363207527, -4.319687165797595, 2]]
    #Register and circuit 
    s = QuantumRegister(1,name='s')
    j = QuantumRegister(4,name='j')
    q = QuantumRegister(2,name='q')
    cr = ClassicalRegister(1,name='cr')
    crtmp= ClassicalRegister(2,name='crtmp')
    qc = QuantumCircuit(s,j,q,cr,crtmp)
    #Gate preparation
    qft = QFT(4,inverse=False)
    iqft = QFT(4,inverse=True)
    gateLis =[]
    gateinvLis=[]
    for i in range(4):
        toPut = Agate(paraLis[i]).control()
        gateLis.append(toPut)
    qc.h(j)
    qc.h(q)
    qc.barrier()
    for i in range(4):
        gate = gateLis[i]
        qc.append(gate,qargs=[j[i]]+q[:])
        qc.barrier()
    qc.append(iqft,qargs=j[:])
    qc.barrier()
    #qc.swap(j[1],j[3])
    for i in range(4):
        angle = 2**(3-i)*np.pi
        angle*=2**(-r+1)
        qc.cry(angle,j[i],s[0])
        qc.barrier()
    qc.append(qft,qargs=j[:])
    qc.barrier()
    for i in range(4):
        gate = gateLis[3-i].inverse()
        qc.append(gate,qargs=[j[3-i]]+q[:])
        qc.barrier()
    
    qc.h(j)

    qc.measure(s,cr)
    qc.barrier()
    qc.measure(q,crtmp)
    #qc.snapshot('res',qubits=q)
    #print(qc)
    style = {'fontsize':30,'figwidth':200

        }
    circuit_drawer(qc,scale=1.6,output='mpl',filename='qc.png',style=style)
    return qc
#myRzz(np.pi)
""" numLis=[]
for i in range(6):
    numLis.append(i)
Agate(numLis,False,True) """
createQC(5)