import numpy as np
from qiskit.quantum_info import Operator
from qiskit.extensions import HamiltonianGate
from qiskit.aqua.operators.legacy import MatrixOperator
from qiskit.aqua import QuantumInstance
from qiskit.aqua.algorithms import QPE
from qiskit import Aer
from qiskit.visualization import circuit_drawer
from qiskit.aqua.components.initial_states import Custom
A = .25 * np.array([[15, 9, 5, -3], [9, 15, 3, -5], [5, 3, 15, -9], [-3, -5, -9, 15]])
t0=2*np.pi
#expA = la.expm(-1.j * A * t0)
op=MatrixOperator(A)
vector = [1/2,1/2,1/2,1/2]
init_state = Custom(2, state_vector=vector)
algo=QPE(operator=op,state_in=init_state)


result = algo.run(QuantumInstance(Aer.get_backend('qasm_simulator')))
print(result)
