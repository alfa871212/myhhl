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
from qiskit.visualization import plot_state_city
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
        print(qc)
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
    for i in range(4):
        gate = gateLis[i]
        qc.append(gate,qargs=[j[i]]+q[:])
    qc.append(iqft,qargs=j[:])
    #qc.swap(j[1],j[3])
    for i in range(4):
        angle = 2**(3-i)*np.pi
        angle*=2**(-r+1)
        qc.cry(angle,j[i],s[0])
    qc.append(qft,qargs=j[:])

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
    
    return qc
def createQC2(r):
    #parameter setting lis
    paraLis16 =[0.196,0.379,0.981,0.589,1.178,16]
    paraLis8 =[1.963,1.115,1.963,2.615,1.178,8]
    paraLis4 =[-0.785,1.017,3.927,2.517,2.356,4]
    paraLis2 =[-9.014*10**(-9),-0.75,1.571,0.75,-1.571,2]
    paraLis = [paraLis16,paraLis8,paraLis4,paraLis2]
    #Register and circuit 
    s = QuantumRegister(1,name='s')
    j = QuantumRegister(4,name='j')
    q = QuantumRegister(2,name='q')
    #cr = ClassicalRegister(1,name='cr')
    #crtmp= ClassicalRegister(2,name='crtmp')
    qc = QuantumCircuit(s,j,q)
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
    for i in range(4):
        gate = gateLis[i]
        qc.append(gate,qargs=[j[i]]+q[:])
    qc.append(iqft,qargs=j[:])
    #qc.swap(j[1],j[3])
    for i in range(4):
        angle = 2**(3-i)*np.pi
        angle*=2**(-r+1)
        qc.cry(angle,j[i],s[0])
    qc.append(qft,qargs=j[:])

    for i in range(4):
        gate = gateLis[3-i].inverse()
        qc.append(gate,qargs=[j[3-i]]+q[:])
    qc.barrier()
    qc.h(j)

    #qc.measure(s,cr)
    qc.barrier()
    print(qc)
    tomo_circuits = state_tomography_circuits(qc,q)
    tomo_circuits_noanc = deepcopy(tomo_circuits)
    to_put = ClassicalRegister(1)
    for i in tomo_circuits:
        i.add_register(to_put)
        i.measure(s,to_put[0])
        print(i)
    
    backend = Aer.get_backend('qasm_simulator')
    results = execute(tomo_circuits,backend,shots=10**5).result()
    print(results.get_counts())
    probs = []
    for circ in tomo_circuits:
        counts = results.get_counts(circ)
        s, f = 0, 0
        for k, v in counts.items():
            if k[0] == "1":
                s += v
            else:
                f += v
        probs.append(s / (f + s))
    #results_noanc = tomo_postselect(results)
    data = StateTomographyFitter(results,tomo_circuits)
    print(data)
    #omo_data = StateTomographyFitter(results_noanc,tomo_circuits_noanc)
    rho_fit = data.fit()
    print(rho_fit)
    '''
    vec = np.sqrt(np.diag(rho_fit))
    hhl_results(vec)
    '''
    
def tomo_postselect(results):
    new_results = deepcopy(results)

    for resultidx, _ in enumerate(results.results):
        old_counts = results.get_counts(resultidx)
        new_counts = {}

        # change the size of the classical register
        new_results.results[resultidx].header.creg_sizes = [new_results.results[resultidx].header.creg_sizes[0]]
        new_results.results[resultidx].header.clbit_labels = new_results.results[resultidx].header.clbit_labels[0:-1]
        new_results.results[resultidx].header.memory_slots = new_results.results[resultidx].header.memory_slots - 1

        for reg_key in old_counts:
            reg_bits = reg_key.split(' ')
            if reg_bits[0] == '1':
                new_counts[reg_bits[1]] = old_counts[reg_key]

        data_counts = new_results.results[resultidx].data.counts
        new_results.results[resultidx].data.counts = new_counts if isinstance(data_counts, dict) else data_counts.from_dict(new_counts)

    return new_results
def error(arr1,arr2):
    sum=0
    for i in range(4):
        sum+=(arr1[i]-arr2[i])**2
def fidelity(vec1,vec2):
    sv1 = Statevector(vec1)
    sv2 = Statevector(vec2)
    return state_fidelity(sv1,sv2)
def normalize(vec):
    
    norm_vec = norm(vec)
    retlis=[]
    #print(norm_vec)
    for i in range(len(vec)):
        retlis.append(vec[i]/norm_vec)
    #print(retlis)
    ret = np.asarray(retlis)
    return ret


    
def hhl_results(vec):
    res_vec =vec
    vector = [1/2,1/2,1/2,1/2]
    vector = np.asarray(vector)
    matrix = [[15, 9, 5, -3],
            [9, 15, 3, -5],
            [5, 3, 15, -9],
            [-3, -5, -9, 15]]
    matrix = np.asarray((matrix))
    # Rescaling the output vector to the real solution vector
    tmp_vec = matrix.dot(res_vec)
    f1 = np.linalg.norm(vector) / np.linalg.norm(tmp_vec)
    # "-1+1" to fix angle error for -0.-0.j
    f2 = sum(np.angle(vector * tmp_vec.conj() - 1 + 1)) / (np.log2(matrix.shape[0]))
    print(f1 * res_vec * np.exp(-1j * f2))

#createQC2(5)

fid_lis=[]
idx_lis=[]
prob_lis=[]
i = 1.00
matrix = [[15, 9, 5, -3],
            [9, 15, 3, -5],
            [5, 3, 15, -9],
            [-3, -5, -9, 15]]
matrix = np.asarray((matrix))
   #vector = [1, 0]
vector = [1/2,1/2,1/2,1/2]
vector = np.asarray(vector)
while i<=10.0:
    print('-'*50)
    i = np.round(i,decimals=3)
    print(f'Simulating with r={i}...')
    #backend = Aer.get_backend('qasm_simulator')
    backend = Aer.get_backend('qasm_simulator')
    qc = createQC(i)
    idx_lis.append(i)
    i+=1
    
    res_state = execute(qc,backend,shots=10**5).result()
    resdic=res_state.get_counts()

    reslis=list(resdic.items())
    print(reslis)

    prob=0
    occur=np.zeros(4)
    for j in range(len(reslis)):
        if reslis[j][0][3]=='1':
            idx_occur = int(reslis[j][0][:2],2)
            #print(idx_occur)
            occur[idx_occur]=(reslis[j][1])**(1/2)
            prob+=(reslis[j][1]/10**5)
    prob_lis.append(prob)       
    
    #snapshots = res_state.data()['snapshots']['statevector']['res']
    #print(snapshots)
    #print(reslis)
    #print(occur)
    print('Our simulation...')
    #occur=np.round(occur/occur[1],decimals=4)
    occur=normalize(occur)
    occur[0]=-1*occur[0]
    print(occur)
    
    print("Ans should be...")
    sol = np.round(np.array([-1,7,11,13]),decimals=4)
    ret=normalize(sol)
    print(ret)
    fid =fidelity(occur,ret)
    fid_lis.append(fid)
    print(f"Fidelity={fid}")
    print(f"Probability={prob}")
    print('-'*50)
    
print('*'*50)
print(f'Max fidelity is {max(fid_lis)}, when r={1+fid_lis.index(max(fid_lis))}')
print(f'Min fidelity is {min(fid_lis)}, when r={1+fid_lis.index(min(fid_lis))}')
fidmax=max(fid_lis)
idx=fid_lis.index(fidmax)
prob_ave=sum(prob_lis)/len(prob_lis)
if prob_lis[idx]>prob_ave:
    print(f"Choosing proper r={1+0.01*idx}, with fidelity={fid_lis[idx]} and probability={prob_lis[idx]}.")
'''

fid_ave=sum(fid_lis)/len(fid_lis)
prob_ave=sum(prob_lis)/len(prob_lis)
for i in range(len(fid_lis)):
    if fid_lis[i]>fid_ave and prob_lis[i]>prob_ave:
        print(f"Choosing proper r={1+0.1*i}, with fidelity={fid_lis[i]} and probability={prob_lis[i]}.")
'''

x_ax = idx_lis
y_ax = np.asarray(fid_lis)
y_ax_p = np.asarray(prob_lis)
'''
weighted_ax = y_ax+y_ax_p
weighted_max=max(weighted_ax)
print(weighted_max)
'''
plt.plot(x_ax,y_ax,label='fidelity')
plt.plot(x_ax,y_ax_p,label='probability')
#plt.plot(x_ax,weighted_ax,label='Weighted')
plt.legend()
plt.savefig('./r_tmp2.png')


