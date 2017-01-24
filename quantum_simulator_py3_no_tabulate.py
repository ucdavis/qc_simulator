#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np # for numerics
import random # for random number generation
from math import sqrt,pi,e # some common function
import math
#import qutip as qp # qutip library for Bloch sphere visualisations
import cmath # library for complex numbers
from collections import Counter
#from tabulate import tabulate # for nice printing

#################### Defining some elementary gates
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

def Rx(angle):
    angle = float(angle)
    return np.array([[cmath.cos(angle/2),-i_*cmath.sin(angle/2)],[-i_*cmath.sin(angle/2), cmath.cos(angle/2)]])

def Ry(angle):
    angle = float(angle)
    return np.array([[cmath.cos(angle/2),-cmath.sin(angle/2)],[cmath.sin(angle/2),cmath.cos(angle/2)]])

def Rz(angle):
    angle = float(angle)
    return np.array([[cmath.exp(-i_*angle/2),0],[0,cmath.exp(i_*angle/2)]])

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
            print( 'Note thate the state you generated was normalised ' \
                    'automatically ')

    else:
        raise StandardError('Cannot interpret input of State() creator.'\
                                ' Please enter a list of valid amplitudes or positions.')
    return state



def renormalise(state): # Renormalise the amplitude vector to unit length
        normalis_factor = float(np.sqrt(np.vdot(state,state)))
        state = state/normalis_factor
        return state

def is_normalised(state): #Check if a state is normalised
    if np.isclose(float(np.vdot(state,state)), float(1.0),rtol = 1e-03):
        return True
    else:
        return False

def measure(state, runs = 1, output = 'outcomes'): #perform measurements on the state
# simulates a repeated generation and projective measurement of the state
# Options: 'outcomes' prints every result while 'stats' prints an overview
    results = np.random.choice(len(state), runs, p=[abs(el)**2 for el in state])
    if output == 'outcomes':
        print ("Measurement Results")
        print ("Index  Basis state ")
        print ("-----   ----------- ")
        for el_res in results:
            print(el_res,'       |', "".join(( "{0:0", \
                            str(int(np.log2(len(state)))),"b}")).format(el_res),'>')

    if output == 'stats':
        hist_dict = Counter(results)
        indices = list(hist_dict.keys())
        occurences = [value/float(runs) for value in list(hist_dict.values())]
        print ("\n Measurement Statistics:")
        print ("rel. occ.   Index   Basis state")
        print ("---------   ------  -----------")
        for i in range(len(indices)):
            print(occurences[i], "          ",indices[i],'         |', "".join(( "{0:0", \
                            str(int(np.log2(len(state)))),"b}")).format(results[i]),'>')
        #print tabulate(printdata, headers = ['rel. occ.', 'Index', 'Basis state'])


    return None

def print_me(state, style = None): # print out current state.
# Options:
# None/empty - simple output of state information
# 'slim' - As in 'None' but without zero amplitude entries
# 'amplitudes' - only array of amplitudes
    np.set_printoptions(precision=3, suppress = True)
    print
    if style == None: # print all nonzero amplitudes
        print("Index Probability Amplitude Basis state ")
        print("----- ----------- --------- ----------- ")
        for i in range(len(state)):
            if not np.isclose(state[i],0.0):
                basis_string = "".join(( "{0:0", \
                                str(int(np.log2(len(state)))),"b}"))
                print('', "{0:04}".format(i), '    ', \
                        "{0:.3f}".format(abs(state[i])**2), '   ', \
                        "{0:.3f}".format(state[i]), \
                        "".join(('  |',basis_string.format(i) , '>' )))
    if style == 'full': # print all amplitudes
        print("Index Probability Amplitude Basis state ")
        print("----- ----------- --------- ----------- ")
        for i in range(len(state)):
            basis_string = "".join(( "{0:0", str(int(np.log2(len(state)))),"b}"))
            print('', "{0:04}".format(i), '    ', \
                        "{0:.3f}".format(abs(state[i])**2), '   ', \
                        "{0:.3f}".format(state[i]), \
                        "".join(('  |',basis_string.format(i) , '>' )))

    if style == 'amplitudes':
        print("Amplitudes: ", state)

    if style == 'probabilities':
        print("Probabilities:\n ", ["{0:.3f}".format(np.abs(item)**2) \
                                    for item in state])

    print
    return None

def grover_iteration(state, marked_pos):
# performs a Grover iteration on a quantum state
    # check if list is of desired format
    if any(item > len(state) for item in marked_pos)\
       or any( not isinstance(item, int) for item in marked_pos):
        raise StandardError('Cannot interpret the list of marked positions'\
                                    ' in grover_iteration()')

    marked_state = [- el if i in marked_pos else el \
                    for i,el in enumerate(state)]
    rotated_state = [-el + 2*np.mean(marked_state) for el in marked_state]
    return rotated_state

'''
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
        #onestate = [[0.,0.,-1.]] # define |1> state vector
        #bloch.add_vectors(onestate) # add |1> state vector

        bloch.vector_color = ['g'] # Bloch vector colour
        bloch.vector_width = 3 # define Bloch vector width

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
                                ' for single qubit states.')

'''
#################### Gate functions

def apply_unitary(gate_matrix, qubit_pos, quantum_state):
    num_qubits = int(math.log(len(quantum_state),2))

    # check if input matrix is a 2x2 matrix
    if (gate_matrix.shape > (2,2)) - (gate_matrix.shape < (2,2)) != 0:
        raise StandardError('Cannot create total unitary. '\
                        'Input matrix must be 2x2.')

    # check if input matrix is unitary
    if np.allclose(np.linalg.inv(gate_matrix),gate_matrix.conjugate().transpose()) == False:
        raise StandardError('Cannot create total unitary.'\
                            ' Input matrix must be unitary.')

    if any(item > num_qubits-1 or item < 0 for item in qubit_pos):
        raise StandardError('Cannot apply quantum gate.'\
                            ' Qubit position is not valid.')

    if (len(qubit_pos) == 1):

        # check if qubit positions are valid
        if (qubit_pos[0] < 0) or (qubit_pos[0] > num_qubits) :
            raise StandardError('Your selected qubit position is out of range.'\
                            ' Please choose a valid qubit position.')
        else:
            # create a list of gates representing the required tensor product
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

        # isolate control and target qubits
        control = qubit_pos[0]
        target = qubit_pos[1]

        if control==target:
            raise StandardError('Target and control are the same. '\
                                'Please choose different target and '\
                                'control qubits.')

        # CASE 1: ADJACENT QUBITS ARE CHOSEN AS CONTROL AND TARGET
        if (control == target-1) or (target == control-1):
            checker = False
            if (target == control-1):
                checker = True
                save_control = control
                control = target
                target = save_control

            # initialize empty 4 x 4 matrix for controlled gate
            cgate = np.zeros((4,4),dtype=np.complex128)

            # if control position is reached:
            # perform the outer product |1><1| and, thereafter, the tensor product with the unitary that shall be controlled
            cgate += np.kron(np.matrix(create_state(1,[0,1])).transpose()*np.matrix(create_state(1,[0,1])), gate_matrix)

            # perform the outer product |0><0| and, thereafter, the tensor product with the identity matrix
            cgate += np.kron(np.matrix(create_state(1,[1,0])).transpose()*np.matrix(create_state(1,[1,0])), eye)
            # convert to array
            cgate = np.array(cgate)

            if num_qubits > 2:
                # perform the tensor products with identity matrices
                for k in range(num_qubits):
                    # pre-multiply identities
                    if all([(k<control),(k<target)]):
                        cgate = np.kron(eye,cgate)
                    # post-multiply identities
                    elif all([(k>control),(k>target)]):
                        cgate = np.kron(cgate, eye)

                if checker:
                    save_control = control
                    control = target
                    target = save_control

                # use the Hadamard trick to reverse the direction of the CNOT gate
                if control > target:
                    quantum_state = apply_unitary(H,[control],quantum_state)
                    quantum_state = apply_unitary(H,[target],quantum_state)
                    quantum_state = np.dot(cgate,quantum_state)
                    quantum_state = apply_unitary(H,[control],quantum_state)
                    quantum_state = apply_unitary(H,[target],quantum_state)

                    return quantum_state

            if checker:
                save_control = control
                control = target
                target = save_control


            # use the Hadamard trick to reverse the direction of the CNOT gate
            if control > target:
                quantum_state = apply_unitary(H,[control],quantum_state)
                quantum_state = apply_unitary(H,[target],quantum_state)
                quantum_state = np.dot(cgate,quantum_state)
                quantum_state = apply_unitary(H,[control],quantum_state)
                quantum_state = apply_unitary(H,[target],quantum_state)
            else:
                # apply the 2**n x 2**n matrix to the quantum state
                quantum_state = np.dot(cgate,quantum_state)

            return quantum_state

        else:
            # obtain the respective gate matrix with the
            # create_controlledGate function
            cgate = create_controlledGate(gate_matrix, qubit_pos, len(quantum_state), num_qubits)

            # apply the 2**n x 2**n matrix to the quantum state
            quantum_state = np.dot(cgate,quantum_state)

            return quantum_state

    # Controlled controlled case: currently only allows for Toffoli
    elif (len(qubit_pos) == 3):

        control1 = qubit_pos[0]
        control2 = qubit_pos[1]
        target = qubit_pos[2]

        # check if input gate is X > only Toffoli allowed for now
        if (gate_matrix==X).all() == False:
            raise StandardError('Cannot create the controlled controlled U gate. '\
                                'Only Toffoli supported so far. '\
                                'Input matrix must be the X gate.')

        quantum_state = apply_unitary(H,[target],quantum_state)
        quantum_state = apply_unitary(X,[control1,target],quantum_state)
        quantum_state = apply_unitary(Tdagger,[target],quantum_state)
        quantum_state = apply_unitary(X,[control2,target],quantum_state)
        quantum_state = apply_unitary(T,[target],quantum_state)
        quantum_state = apply_unitary(X,[control1,target],quantum_state)
        quantum_state = apply_unitary(Tdagger,[target],quantum_state)
        quantum_state = apply_unitary(X,[control2,target],quantum_state)
        quantum_state = apply_unitary(T,[target],quantum_state)
        quantum_state = apply_unitary(T,[control1],quantum_state)
        quantum_state = apply_unitary(X,[control2,control1],quantum_state)
        quantum_state = apply_unitary(H,[target],quantum_state)
        quantum_state = apply_unitary(T,[control2],quantum_state)
        quantum_state = apply_unitary(Tdagger,[control1],quantum_state)
        quantum_state = apply_unitary(X,[control2,control1],quantum_state)

        return quantum_state

    else:
        raise StandardError('Too many qubits specified.'\
                                ' Please enter a maximum of 2 valid positions.')

    return None

def create_controlledGate(gate_matrix, qubit_pos, num_amplitudes, num_qubits):
    control = qubit_pos[0]
    target = qubit_pos[1]

    if ((control-target) == -(num_qubits-1)):
        cgate = np.eye(num_amplitudes,num_amplitudes)

        iteration_list = np.array(int(num_amplitudes/2))
        print(iteration_list)
        value_save = int(num_amplitudes/2)
        for k in range(int(num_amplitudes/4-1)):
            iteration_list = np.append(iteration_list,value_save+2)
            value_save = value_save+2

        for m in iteration_list:
            cgate[np.array([m,m+1])]=cgate[np.array([m+1,m])]
        return cgate

    elif ((control-target) == num_qubits-1):
        cgate = np.eye(num_amplitudes,num_amplitudes)

        iteration_list = np.array(1)
        value_save = 1
        for k in range(int(num_amplitudes/4-1)):
            iteration_list = np.append(iteration_list,value_save+2)
            value_save = value_save+2

        for m in iteration_list:
            cgate[np.array([m,m+int(num_amplitudes/2)])]=cgate[np.array([m+int(num_amplitudes/2),m])]

        return cgate

    elif (control <= num_qubits-2) and (target <= num_qubits-2):
        pre_cgate = create_controlledGate(gate_matrix, qubit_pos, int(num_amplitudes/2), num_qubits-1)
        cgate = np.kron(pre_cgate,eye)

        return cgate

    elif (control == num_qubits-1) or (target == num_qubits-2):
        pre_cgate = create_controlledGate(gate_matrix, [qubit_pos[0]-1, qubit_pos[1]-1], int(num_amplitudes/2), num_qubits-1)
        cgate = np.kron(eye,pre_cgate)

        return cgate

    else:
        qubit_pos = [ x-1 for x in qubit_pos]
        pre_cgate = create_controlledGate(gate_matrix, qubit_pos, int(num_amplitudes/2), num_qubits-1)
        cgate = np.kron(eye,pre_cgate)

        return cgate
