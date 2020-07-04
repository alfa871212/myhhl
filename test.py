from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit
from qiskit.extensions.simulator import snapshot
import simulation as sim

qc = QuantumCircuit(3) 
qc.h(0) 
qc.h(1) 
qc.snapshot("1", qubits=[1])
qc.snapshot("2", qubits=[2])
print(qc)
qc.measure_all()


res_state = sim.mySim(qc)
snapshots = res_state.data()['snapshots']['statevector']
print(snapshots)