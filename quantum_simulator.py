#!/usr/bin/env python
import numpy as np # for numerics
import random # for random number generation
from math import sqrt,pi,e # some common function
import math
import qutip as qp # qutip library for Bloch sphere visualisations
import cmath # library for complex numbers

#################### State functions
def create_state(num_qubits, lis):
    state = np.zeros(num_qubits)
    # if the list has less entries than amplitudes (2**n),
    # interpret the entries as positions
    if len(lis) < 2**num_qubits:
        # check if all entries are valid positions
        if any(isinstance(item, complex) or item > 2**num_qubits \
                   or not isinstance(item, int) for item in lis) :
            raise StandardError('Cannot interpret input of State() creator.'\
                                    ' Please enter a list of valid positions.')
        # initialise state
        state = np.array([1./sqrt(len(lis)) if i in lis else 0 \
                                 for i in range(2**num_qubits)])
    # else if the list has as many entries as amplitudes (2**n),
    # interpret the entries as amplitudes
    elif len(lis) == 2**num_qubits:
        state = np.array(lis)
        if not is_normalised(state):
            state = renormalise(state)
            print 'Note thate the state you generated was normalised ' \
                    'automatically '

    else:
        raise StandardError('Cannot interpret input of State() creator.'\
                                ' Please enter a list of valid amplitudes or positions.')
    return state



def renormalise(state): # Renormalise the amplitude vector to unit length
        normalis_factor = float(np.sqrt(np.vdot(state,state)))
        state = state/normalis_factor
        return state

def is_normalised(state): #Check if a state is normalised
    if np.isclose(float(np.vdot(state,state)), float(1.0)):
        return True
    else:
        return False

def measure(state, runs = 1, output = 'outcomes'): #perform measurements on the state
# simulates a repeated generation and projective measurement of the state
# Options: 'outcomes' prints every result while 'stats' prints an overview
    results = np.random.choice(len(state), runs, p=[abs(el)**2 for el in state])
    if output == 'outcomes':
        print "Measurement Results"
        print "Index  Basis state "
        print "-----   ----------- "
        for el_res in results:
            print "{0:04}".format(el_res),'   |', "".join(( "{0:0", \
                                str(int(np.log2(len(state)))),"b}")).format(el_res),'>'
    #if output == 'stats':
    #still to do
    return None

def print_me(state, style = None): # print out current state.
# Options:
# None/empty - simple output of state information
# 'slim' - As in 'None' but without zero amplitude entries
# 'amplitudes' - only array of amplitudes
    np.set_printoptions(precision=3, suppress = True)
    print
    if style == None: # print all nonzero amplitudes
        print "Index Probability Amplitude Basis state "
        print "----- ----------- --------- ----------- "
        for i in range(len(state)):
            if not np.isclose(state[i],0.0):
                basis_string = "".join(( "{0:0", \
                                str(int(np.log2(len(state)))),"b}"))
                print '', "{0:04}".format(i), '    ', \
                        "{0:.3f}".format(abs(state[i])**2), '   ', \
                        "{0:.3f}".format(state[i]), \
                        "".join(('  |',basis_string.format(i) , '>' ))
    if style == 'full': # print all amplitudes
        print "Index Probability Amplitude Basis state "
        print "----- ----------- --------- ----------- "
        for i in range(len(state)):
            basis_string = "".join(( "{0:0", str(int(np.log2(len(state)))),"b}"))
            print '', "{0:04}".format(i), '    ', \
                        "{0:.3f}".format(abs(state[i])**2), '   ', \
                        "{0:.3f}".format(state[i]), \
                        "".join(('  |',basis_string.format(i) , '>' ))

    if style == 'amplitudes':
        print "Amplitudes: ", state

    print
    return None

def project_on_blochsphere(state):
    if len(state) == 2:
        alpha = state[0]
        beta = state[1]
        bloch = qp.Bloch() # initialise Bloch sphere
        # Define the x,y and z axes for the Bloch sphere
        #x_basis = (qp.basis(2,0)+(1+0j)*qp.basis(2,1)).unit()
        #y_basis = (qp.basis(2,0)+(0+1j)*qp.basis(2,1)).unit()
        #z_basis = (qp.basis(2,0)+(0+0j)*qp.basis(2,1)).unit()
        #bloch.add_states([x_basis,y_basis,z_basis]) # add the axes to the Bloch sphere
        bloch.vector_color = ['g'] # Bloch vector colour
        bloch.vector_width = 3 # define Bloch vector width
        #onestate = [[0.,0.,-1.]] # define |1> state vector
        #bloch.add_vectors(onestate) # add |1> state vector

        # Find and eliminate global phase
        angle_alpha = cmath.phase(alpha)
        angle_beta = cmath.phase(beta)

        if angle_beta < 0:
            if angle_alpha < angle_beta:
                alpha_new = alpha/cmath.exp(1j*angle_beta)
                beta_new = beta/cmath.exp(1j*angle_beta)
            else:
                alpha_new = alpha/cmath.exp(1j*angle_alpha)
                beta_new = beta/cmath.exp(1j*angle_alpha)
        else:
            if angle_alpha > angle_beta:
                alpha_new = alpha/cmath.exp(1j*angle_beta)
                beta_new = beta/cmath.exp(1j*angle_beta)
            else:
                alpha_new = alpha/cmath.exp(1j*angle_alpha)
                beta_new = beta/cmath.exp(1j*angle_alpha)

        if abs(alpha) == 0 or abs(beta) == 0:
            if alpha == 0:
                bloch.clear()
                down = [0,0,-1]
                bloch.add_vectors(down)
            else:
                bloch.clear()
                up = [0,0,1]
                bloch.add_vectors(up)
        else:
            # Compute theta and phi from alpha and beta
            theta = 2*cmath.acos(alpha_new)
            phi = -1j*cmath.log(beta_new/cmath.sin(theta/2))

            # Compute the cartesian coordinates
            x = cmath.sin(theta)*cmath.cos(phi)
            y = cmath.sin(theta)*cmath.sin(phi)
            z = cmath.cos(theta)

            # Create the new state vector and plot it onto the Bloch sphere
            new_vec = [x.real,y.real,z.real]
            bloch.add_vectors(new_vec)

        bloch.show()
    else:
        raise StandardError('Bloch projection is only supported'\
                                'for single qubit states.')


#################### Gate functions

# define some elementary gates
i_ = np.complex(0,1)
H = 1./sqrt(2)*np.array([[1, 1],[1, -1]])
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -i_],[i_, 0]])
Z = np.array([[1,0],[0,-1]])
eye = np.eye(2,2)
S=np.array([[1,0],[0,i_]])
Sdagger=np.array([[1,0],[0,-i_]])
T=np.array([[1,0],[0, e**(i_*pi/4.)]])
Tdagger=np.array([[1,0],[0, e**(-i_*pi/4.)]])

def apply_total_unitary(gate_matrix, qubit_pos, quantum_state):

    if (len(qubit_pos) == 1):

        num_qubits = int(math.log(len(quantum_state),2))

        print (qubit_pos[0] > num_qubits)
        if (qubit_pos[0] < 0) or (qubit_pos[0] > num_qubits) :
            raise StandardError('Your selected qubit position is out of range.'\
                            ' Please choose a valid qubit position.')
        else:
            if cmp(gate_matrix.shape, (2,2)) != 0:
                raise StandardError('Cannot create new Gate(). '\
                                'Input matrix must be 2x2.')
            unitary_list = [gate_matrix if qubit_pos[0] == i else np.eye(2,2) for i in range(num_qubits)]

            # calculate the Kronecker tensor product with (n-1) identity matrices
            # to obtain a 2**n x 2**n matrix that can be acted on n qubits
            u_new = unitary_list[0]
            for k in range(num_qubits-1):
                u_new = np.kron(u_new, unitary_list[k+1])
            gate = u_new

            # apply the 2**n x 2**n matrix to the current quantum_state
            quantum_state = np.dot(gate,quantum_state)

            return quantum_state

    elif (len(qubit_pos) == 2):
        test

    else:
        raise StandardError('Too many qubits specified.'\
                                ' Please enter a maximum of 2 valid positions.')

    return None


#################### Execution

state = create_state(1,[0.5, 0.5])
print_me(state, 'full')
state = apply_total_unitary(X, [0], state)
print_me(state, 'full')
project_on_blochsphere(state)
measure(state, 4)
