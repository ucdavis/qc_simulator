from quantum_simulator_py2_no_tabulate import *

#Write code here
quantum_state = create_state(3,[0])
quantum_state = apply_unitary(X,[0,1,2],quantum_state)
print_me(quantum_state, 'full')
