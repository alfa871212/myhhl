import typing
from qiskit.circuit.library import QFT,XGate,RZZGate
import qiskit
from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit
import numpy as np
import scipy.linalg as la


QubitType = typing.Tuple[qiskit.QuantumRegister, int]
def myRzz(theta):
    qr=QuantumRegister(1)
    qc=QuantumCircuit(qr)
    qc.u1(theta, qr)
    qc.x(qr)
    qc.u1(theta, qr)
    qc.x(qr)
    gate=qc.to_gate()
    gate.name=f'myRzz_{theta}'
    return gate

def Agate(numLis,inverse=False,plot_gate=False):
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
        print(qc)
        return
    ret = qc.to_gate()
    if inverse:
        ret.inverse()
    return ret
def swap(U):
    from copy import deepcopy
    cpy = deepcopy(U)
    cpy[[1, 2], :] = cpy[[2, 1], :]
    cpy[:, [1, 2]] = cpy[:, [2, 1]]
    return cpy
def error(paraLis,power):
    
    ancilla = QuantumRegister(1)
    b = QuantumRegister(2)
    cr = ClassicalRegister(1)
    qc = QuantumCircuit(ancilla,b,cr)
    gate =Agate(paraLis)
    invgate=Agate(paraLis).inverse()
    qc.append(invgate,qargs=b[:])
    unitary_sim = qiskit.Aer.get_backend('unitary_simulator')
    res = qiskit.execute(qc, unitary_sim).result()
    unitary = res.get_unitary()
    A = .25 * np.array(
            [[15, 9, 5, -3], [9, 15, 3, -5], [5, 3, 15, -9], [-3, -5, -9, 15]])
    t0 = 2 * np.pi
            
    expA = swap(la.expm(-1.j * A * t0 * (2 ** power/ 16)))
    unit = unitary[1::2, 1::2]
    np.set_printoptions(precision=2)
    err = la.norm(unit - expA)
    print("{: g}".format(err), end='\r', flush=True)
    return err



def main():
    # Optimise!
    import scipy.optimize as opt
    #paraLis16 =[0.1,0.3,0.9,0.5,1.2]
    #paraLis8 =[2.0,1.0,2.0,3.0,1.0]
    #paraLis4 =[-0.5,1.0,3.5,2.5,2.0]
    #paraLis2 =[0,-0.5,1.5,0.75,-2.0]
    paraLis16 =[0.196,0.379,0.981,0.589,1.178]
    paraLis8 =[1.963,1.115,1.963,2.615,1.178]
    paraLis4 =[-0.785,1.017,3.927,2.517,2.356]
    paraLis2 =[-9.014*10**(-9),-0.75,1.571,0.75,-1.571]
    lis = [paraLis16,paraLis8,paraLis4,paraLis2]
    reslis=[]
    for i in range(len(lis)):
        opt_res = opt.minimize(error,lis[i],args=(4-i))
        if opt_res.success:
            res = opt_res.x.tolist()
            
            res.append(2**(4-i))
           
            reslis.append(res)
        print(opt_res)
    print(reslis)
    


if __name__ == '__main__':
    main()