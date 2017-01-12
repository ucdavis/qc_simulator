#!/usr/bin/env python

import numpy as np # for numerics
import random # for random number generation
from tabulate import tabulate # for nice printing
from math import sqrt,pi,e # some common function
import qutip as qp # qutip library for Bloch sphere visualisations
import cmath # library for complex numbers

####
## States
####

class State:
# call with State(number_of_qubits, list_of_basisstatepositions_OR_amplitudes)

    # ground state as default instance of a State
    def __init__(self, lis = [0]):

        # if the list has less entries than amplitudes 2**n,
        # interpret the entries as positions
        if len(lis) < 2**num_qubits:
            # check if all entries are valid positions
            if any(isinstance(item, complex) \
                   or item > 2**num_qubits \
                   or not isinstance(item, int) \
                   for item in lis) :
                raise StandardError('Cannot interpret input of State() creator.'\
                                ' Please enter a list of valid positions.')
            # initialise state
            self.state = np.array([1./sqrt(len(lis)) if i in lis else 0 \
                                   for i in range(2**num_qubits)])
        # else if the list has as many entries as amplitudes 2**n,
        # interpret the entries as amplitudes
        elif len(lis) == 2**num_qubits:
            #TODO: check if normalised
            self.state = np.array(lis)

        else:
            raise StandardError('Cannot interpret input of State() creator.'\
                                ' Please enter a list of valid amplitudes or positions.')

    def print_amplitudes(self):
        np.set_printoptions(precision=3, suppress = True)
        print self.state

    def print_basisstates(self):
        np.set_printoptions(precision=3, suppress = True)
        print self.state

    def print_full(self):
        np.set_printoptions(precision=3, suppress = True)
        for i in range(len(self.state)):
            basis_string = "".join(( "{0:0", str(int(np.log2(len(self.state)))),"b}"))
            print round(self.state[i],2), "".join(('  |',basis_string.format(i) , '>' ))

    def print_full_pretty(self):
        np.set_printoptions(precision=3, suppress = True)
        data = []
        for i in range(len(self.state)):
            basis_string = "".join(( "{0:0", str(int(np.log2(len(self.state)))),"b}"))
            row = [ self.state[i],  "".join(('  |',basis_string.format(i) , '>' )) ]
            data.append(row)
        print
        print tabulate(data, headers=['Amplitude',' Basis state'] )
        print

    def project_on_blochsphere(self):
        if num_qubits == 1:
            alpha = self.state[0]
            beta = self.state[1]
            bloch = qp.Bloch() # initialise Bloch sphere
            # Define the x,y and z axes for the Bloch sphere
            x_basis = (qp.basis(2,0)+(1+0j)*qp.basis(2,1)).unit()
            y_basis = (qp.basis(2,0)+(0+1j)*qp.basis(2,1)).unit()
            z_basis = (qp.basis(2,0)+(0+0j)*qp.basis(2,1)).unit()
            bloch.add_states([x_basis,y_basis,z_basis]) # add the axes to the Bloch sphere
            bloch.vector_color = ['k','k','k','k','g'] # Bloch vector colours
            bloch.vector_width = 3 # define Bloch vector width
            onestate = [[0.,0.,-1.]] # define |1> state vector
            bloch.add_vectors(onestate) # add |1> state vector

            # Find and eliminate global phase
            angle_alpha = cmath.phase(alpha)
            #print "angle_alpha:", angle_alpha
            angle_beta = cmath.phase(beta)
            #print "angle_beta:", angle_beta

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
                    #alpha_new = alpha/cmath.exp(1j*angle_beta)
                    #beta_new = beta/cmath.exp(1j*angle_beta)
                else:
                    alpha_new = alpha/cmath.exp(1j*angle_alpha)
                    beta_new = beta/cmath.exp(1j*angle_alpha)
            #print "alpha_new:", alpha_new
            #print "beta_new:", beta_new

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
        	    #print "theta:", theta
        	    #print "phi:", phi

        	    # Compute the cartesian coordinates
        	    x = cmath.sin(theta)*cmath.cos(phi)
        	    y = cmath.sin(theta)*cmath.sin(phi)
        	    z = cmath.cos(theta)
        	    #print "x:", x
        	    #print "y:", y
        	    #print "z:", z

        	    # Create the new state vector and plot it onto the Bloch sphere
        	    new_vec = [x.real,y.real,z.real]
                bloch.add_vectors(new_vec)

            bloch.show()

    #zero_state = np.array('1; 0')
    #one_state = np.array('0; 1')
'''
####
## States2
####

class State2:
# call with State(number_of_qubits, list_of_basisstatepositions_OR_amplitudes)

    # ground state as default instance of a State
    def __init__(self, qubits, lis = [0]):

        # if the list has less entries than amplitudes 2**n,
        # interpret the entries as positions
        if len(lis) < 2**qubits:
            # check if all entries are valid positions
            if any(isinstance(item, complex) \
                   or item > 2**qubits \
                   or not isinstance(item, int) \
                   for item in lis) :
                raise StandardError('Cannot interpret input of State() creator.'\
                                ' Please enter a list of valid positions.')
            # initialise state
            temp = np.matrix([1./sqrt(len(lis)) if i in lis else 0 \
                                   for i in range(2**qubits)])
            self.state = temp.transpose()
        # else if the list has as many entries as amplitudes 2**n,
        # interpret the entries as amplitudes
        elif len(lis) == 2**qubits:
            #TODO: check if normalised
            self.state = np.matrix(lis).transpose()
        else:
            raise StandardError('Cannot interpret input of State() creator.'\
                                ' Please enter a list of valid amplitudes or positions.')

    def print_amplitudes(self):
        np.set_printoptions(precision=3, suppress = True)
        print self.state

    def print_basisstates(self):
        np.set_printoptions(precision=3, suppress = True)
        print self.state

    def print_full(self):
        np.set_printoptions(precision=3, suppress = True)
        for i in range(len(self.state)):
            basis_string = "".join(( "{0:0", str(int(np.log2(len(self.state)))),"b}"))
            print round(self.state[i],2), "".join(('  |',basis_string.format(i) , '>' ))

    def print_full_pretty(self):
        np.set_printoptions(precision=3, suppress = True)
        data = []
        for i in range(len(self.state)):
            basis_string = "".join(( "{0:0", str(int(np.log2(len(self.state)))),"b}"))
            row = [ self.state[i],  "".join(('  |',basis_string.format(i) , '>' )) ]
            data.append(row)
        print
        print tabulate(data, headers=['Amplitude',' Basis state'] )
        print

    def project_on_blochsphere(self):
        alpha = self.state[0]
        beta = self.state[1]
        bloch = qp.Bloch() # initialise Bloch sphere
        # Define the x,y and z axes for the Bloch sphere
        x_basis = (qp.basis(2,0)+(1+0j)*qp.basis(2,1)).unit()
        y_basis = (qp.basis(2,0)+(0+1j)*qp.basis(2,1)).unit()
        z_basis = (qp.basis(2,0)+(0+0j)*qp.basis(2,1)).unit()
        bloch.add_states([x_basis,y_basis,z_basis]) # add the axes to the Bloch sphere
        bloch.vector_color = ['k','k','k','k','g'] # Bloch vector colours
        bloch.vector_width = 3 # define Bloch vector width
        onestate = [[0.,0.,-1.]] # define |1> state vector
        bloch.add_vectors(onestate) # add |1> state vector

        # Find and eliminate global phase
        angle_alpha = cmath.phase(alpha)
        #print "angle_alpha:", angle_alpha
        angle_beta = cmath.phase(beta)
        #print "angle_beta:", angle_beta

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
        #print "alpha_new:", alpha_new
        #print "beta_new:", beta_new

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
        	#print "theta:", theta
        	#print "phi:", phi

        	# Compute the cartesian coordinates
        	x = cmath.sin(theta)*cmath.cos(phi)
        	y = cmath.sin(theta)*cmath.sin(phi)
        	z = cmath.cos(theta)
        	#print "x:", x
        	#print "y:", y
        	#print "z:", z

        	# Create the new state vector and plot it onto the Bloch sphere
        	new_vec = [x.real,y.real,z.real]
        	bloch.add_vectors(new_vec)

        bloch.show()

    zero_state = np.matrix('1; 0')
    one_state = np.matrix('0; 1')
'''
####
## Gates
####

class Gate:

    # identity as default instance of a Gate
    def __init__(self, unitary, opt angle, qubit_pos = -1):

        # check if input matrix is unitary
        if np.allclose(np.linalg.inv(unitary),unitary.conjugate().transpose()) == False:
            raise StandardError('Cannot create new Gate().'\
                                'Input matrix must be unitary.')

        # if gate is applied to all n qubits check if its dimension are 2**n x 2**n
        if qubit_pos == -1:
            if cmp(unitary.shape, (2**num_qubits,2**num_qubits)) != 0:
                raise StandardError('Cannot create new Gate().'\
                   'Input matrix must be 2^n x 2^n.')
            self.gate = np.array(unitary)

        # if gate is applied to a single qubit check if matrix has dimension 2x2
        if 0 <= qubit_pos < num_qubits:
            if cmp(unitary.shape, (2,2)) != 0:
                raise StandardError('Cannot create new Gate().'\
                                'Input matrix must be 2x2.')
            unitary_list = [unitary if qubit_pos == i else Gate.eye for i in range(num_qubits)]

            # calculate the Kronecker tensor product with (n-1) identity matrices
            # to obtain a 2**n x 2**n matrix that can be acted on n qubits
            #u_new = np.kron(unitary_list[0], unitary_list[1])
            u_new = unitary_list[0]
            for k in range(num_qubits-1):
                 u_new = np.kron(u_new, unitary_list[k+1])
            self.gate = u_new

    # define some elementary gates
    i_ = np.complex(0,1)
    H = 1./sqrt(2)*np.array([[1, 1],[1, -1]])
    X = np.matrix('0 1; 1 0')
    Y = np.matrix([[0, -i_],[i_, 0]])
    Z = np.matrix([[1,0],[0,-1]])
    eye = np.eye(2,2)
    S=np.matrix([[1,0],[0,i_]])
    Sdagger=np.matrix([[1,0],[0,-i_]])
    T=np.matrix([[1,0],[0, e**(i_*pi/4.)]])
    Tdagger=np.matrix([[1,0],[0, e**(-i_*pi/4.)]])
    Rx(x) = np.array([angle])

    #TODO: CNOT
class ControlledGate:

    def __init__(self, unitary, control, target):

        # CASE 1: ADJACENT QUBITS ARE CHOSEN AS CONTROL AND TARGET
        if (control == target-1) or (target == control-1):

            #control_list = [0 if control == i else 1 for i in range(2)]
            #print control_list
            #identity_list = [0 if control_list[j] == 1 else 1 for j in range(2)]
            #print identity_list
            # initialize empty 4 x 4 matrix for CNOT
            self.gate = np.zeros((4,4))
            for k in range(2):
                if k == control:
                    # if control position is reached:
                    # perform the outer product |1><1| and, thereafter, the tensor product with the unitary that shall be controlled
                    self.gate += np.kron(State2(1,[0,1]).state*State2(1,[0,1]).state.transpose(), unitary)
                else:
                    # perform the outer product |0><0| and, thereafter, the tensor product with the identity matrix
                    self.gate += np.kron(State2(1,[1,0]).state*State2(1,[1,0]).state.transpose(), Gate.eye)

            if control > target:
                self.gate = np.kron(Gate.H,Gate.H)*self.gate*np.kron(Gate.H,Gate.H)
            if num_qubits > 2:
                for k in range(num_qubits-2):
                    self.gate = np.kron(self.gate, Gate.eye)
                '''
                order_list = [0]*(num_qubits-1)
                marker = False
                for i in range(num_qubits-1):
                    if ((control == i) or (target == i)) and (marker == False):
                        order_list[i] = self.gate
                        marker = True
                    else:
                        order_list[i] = Gate.eye
                        #print "test", order_list[i]
                        #order_list = [self.gate; i += 1; if (control == i) or (target == i) else Gate.eye for i in range(num_qubits)]
                u_new = order_list[0]
                for k in range(num_qubits-2):
                    u_new = np.kron(u_new, order_list[k+1])
                    self.gate = u_new
                    print self.gate
                '''
        # OTHER CASES:
        else:

            # initialize zero matrix
            CNOT_initial = np.zeros((2**n,2**n))
            for i in range(2**n-2):
                CNOT_initial[i:i+2,i:i+2] = Gate.eye



#SCRIPT

# DO NOT CHANGE THE GLOBAL DEFINITION OF num_qubits!
global num_qubits;
# define the number of qubits for this simulation
num_qubits = 1


myState1 = State([0])
myState1.print_full_pretty()
newGate = Gate(Gate.X, 0)
Gate(Gate.Rx(pi/2),0)
newGate2 = Gate(Gate.H, 0)
#myState1.print_full_pretty()
#myState2 = State2.one_state
print myState1.state
#print myState1.state
#myState1 = myState1.apply_gate(newGate)
print newGate2.gate
myState1.state = np.dot(newGate2.gate,myState1.state)
#print myState1.state
#myState1.state = newGate.gate*myState1.state
#print myState1.state
#newGate = Gate(Gate.X, 1)
#print newGate.gate
#myState1.state = newGate.gate*myState1.state

#print type(myState1.state)
myState1.print_full_pretty()
#print np.matrix(myState1.state).transpose()
#print myState1.state
CX = ControlledGate(Gate.X,1,0)
print type(CX.gate)
#print CX.gate
#myState1.state = np.dot(myState1.state,CX.gate)
#print CX.gate
#myState1.print_full_pretty()
myState1.project_on_blochsphere()
