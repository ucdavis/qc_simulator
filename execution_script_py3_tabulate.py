from quantum_simulator_py3_tabulate import *

# write code to execute here
quantum_state = create_state(3,[0])
quantum_state = apply_unitary(X,[2,1,0],quantum_state)
print_me(quantum_state, 'probabilities')
