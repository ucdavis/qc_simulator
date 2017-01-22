from quantum_simulator import *
'''
mynewstate = create_state(2, [0,complex(1,1),complex(1,1),complex(1,1)])

print_me(mynewstate)

for i in range(10):
    mynewstate = grover_iteration(mynewstate, [0])
    print_me(mynewstate, 'probabilities')
    print is_normalised(mynewstate)

'''
#################### Execution

quantum_state = create_state(3,[0])
print Rx(-math.pi/4)
quantum_state = apply_unitary(Y, [1,2], quantum_state)
#project_on_blochsphere(quantum_state)

'''
quantum_state = create_state(4,[15])
print_me(quantum_state, 'full')

quantum_state = apply_total_unitary(X, [0,1,2], quantum_state)

print_me(quantum_state, 'full')
#project_on_blochsphere(state)
measure(quantum_state, 4)
'''


'''
    ##### CODE SNIPPET SWAP TEST #####

# initialize initial quantum state |000>
# 0th qubit is used as control
# 1st and 2nd qubit are to be swapped
quantum_state = create_state(3,[0])

# Put the first qubit into state |+> = 0.707*|0> + 0.7071*|1>
# Inner product is <q1|q2> = 0.7071
quantum_state = apply_total_unitary(H,[1],quantum_state)

# print the quantum state before SWAP test
print_me(quantum_state, 'full')

### START: MAIN SWAP TEST ALGORITHM

# Put control qubit into uniform superposition
quantum_state = apply_total_unitary(H,[0],quantum_state)

# FREDKIN Gate (controlled SWAP)
# decomposed into three Toffoli gates
quantum_state = apply_total_unitary(X,[0,1,2],quantum_state)
quantum_state = apply_total_unitary(X,[0,2,1],quantum_state)
quantum_state = apply_total_unitary(X,[0,1,2],quantum_state)

# interfere the two states with another Hadamard on the control qubit
quantum_state = apply_total_unitary(H, [0], quantum_state)

### END: MAIN SWAP TEST ALGORITHM

# print the final quantum state
print_me(quantum_state, 'full')

# Measure and gather statistics
measure(quantum_state,10)
'''
