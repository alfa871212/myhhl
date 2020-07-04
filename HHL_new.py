from qiskit import Aer,IBMQ
IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q-hub-ntu')
from qiskit.circuit.library import QFT
from qiskit.aqua import QuantumInstance, aqua_globals
from qiskit.quantum_info import state_fidelity
from qiskit.aqua.algorithms import HHL, NumPyLSsolver
from qiskit.aqua.components.eigs import EigsQPE
from qiskit.aqua.components.reciprocals import LookupRotation
from qiskit.aqua.operators import MatrixOperator
from qiskit.aqua.components.initial_states import Custom
from qiskit.visualization import circuit_drawer
import numpy as np


def create_eigs(matrix, num_ancillae, num_time_slices, negative_evals):
    ne_qfts = [None, None]
    if negative_evals:
        num_ancillae += 1
        ne_qfts = [QFT(num_ancillae - 1), QFT(num_ancillae - 1).inverse()]
    ret = EigsQPE(MatrixOperator(matrix=matrix),
                   QFT(num_ancillae).inverse(),
                   num_time_slices=num_time_slices,
                   num_ancillae=num_ancillae,
                   expansion_mode='suzuki',
                   expansion_order=2,
                   evo_time=None,  # This is t, can set to: np.pi*3/4
                   negative_evals=negative_evals,
                   ne_qfts=ne_qfts)
    #print(ret.construct_circuit(mode='circuit'))
    return ret

def fidelity(hhl, ref):
    solution_hhl_normed = hhl / np.linalg.norm(hhl)
    solution_ref_normed = ref / np.linalg.norm(ref)
    fidelity = state_fidelity(solution_hhl_normed, solution_ref_normed)
    print("Fidelity:\t\t %f" % fidelity)

def main():
   #matrix = [[1, -1/3], [-1/3, 1]]
   matrix = [[15, 9, 5, -3],
            [9, 15, 3, -5],
            [5, 3, 15, -9],
            [-3, -5, -9, 15]]
   #vector = [1, 0]
   vector = [1/2,1/2,1/2,1/2]
   orig_size = len(vector)
   matrix, vector, truncate_powerdim, truncate_hermitian = HHL.matrix_resize(matrix, vector)
# Initialize eigenvalue finding module
   eigs = create_eigs(matrix, 3, 50, False)
   num_q, num_a = eigs.get_register_sizes()

# Initialize initial state module
   init_state = Custom(num_q, state_vector=vector)

# Initialize reciprocal rotation module
   reci = LookupRotation(negative_evals=eigs._negative_evals, evo_time=eigs._evo_time)
   
   algo = HHL(matrix, vector, truncate_powerdim, truncate_hermitian, eigs,
           init_state, reci, num_q, num_a, orig_size)
   circuit_drawer(algo.construct_circuit(),output='mpl',filename='./hhlnew.png')
   result = algo.run(QuantumInstance(Aer.get_backend('statevector_simulator')))

   print("Solution:\t\t", np.round(result['solution'], 5))

   result_ref = NumPyLSsolver(matrix, vector).run()
   print("Classical Solution:\t", np.round(result_ref['solution'], 5))

   print("Probability:\t\t %f" % result['probability_result'])
   fidelity(result['solution'], result_ref['solution'])

if __name__ == '__main__':
  main()
｀｀


