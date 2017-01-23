# Quantum Computing Simulator (qc_simulator)

[![PyPI](https://img.shields.io/pypi/pyversions/Django.svg)]()

A simple educational **quantum computing simulator** written in Python. The library was specifically designed for the **28th Chris Engelbrecht Summer School of the South African National Institute for Theoretical Physics**.

This repository contains:

1. [The quantum computing simulator library for Python 2.7](quantum_simulator_py2.py).
2. [The quantum computing simulator library for Python 3.4, 3.5 and 3.6](quantum_simulator_py3.py).
3. [An execution script for Python 2.7](execution_script_py2.py) to use the quantum computing simulator library.
4. For people using Python 3.4, 3.5 or 3.6, please use [the execution script for these versions](execution_script_py3.py).
5. [The README](README.md) - containing the full documentation for the library.


## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
  - [Quantum States](#quantum-states)
  - [Quantum Gates](#quantum-gates)
  - [Measurement & Stats](#measurement-&-stats)
- [Contribute](#contribute)
- [License](#license)

## Background

This library was designed and written by [@mariaschuld](https://github.com/mariaschuld) and [@markf94](https://github.com/markf94) from the [Centre for Quantum Technology](http://quantum.ukzn.ac.za/) at the [University of KwaZulu-Natal](www.ukzn.ac.za). It is meant as an educational tool to simulate quantum computations. The main aim was to come up with a syntax that is easy to read & understand for beginners. Additionally, the developed library is scalable to any number of qubits unless your RAM allows it. The library includes initialisation of arbitrary quantum states and amplitude distributions, Bloch sphere projections, a variety of single-qubit gates, all-to-all connected CNOT and Toffoli (CCNOT) gates. With easy syntax and without visible for-loops the simulation can be repeated many times and the resulting statistics are displayed.

## Install

This project uses [Python](http://python.org/) and is compatible with versions 2.7, 3.4, 3.5 and 3.6. The following packages are required: [numpy](http://www.numpy.org), [QuTiP](http://qutip.org), [random](https://docs.python.org/2/library/random.html), [collections](https://docs.python.org/2/library/collections) and [tabulate](https://pypi.python.org/pypi/tabulate). **Go check them out if you don't have them locally installed.**

Just download this repository as a ZIP file or clone it via your favourite shell. Unzip the folder and you are ready to go. No additional installation required.

## Usage



### Quantum states

To*******************
CREATE A NEW STATE:

create_state(num_qubits, lis)

INPUTS:
-> num_qubits:
				number of qubits of the quantum system
-> lis:
				if lis is a list of integers with length < num_qubits, it is
				interpreted as a list of positions of nonzero amplitudes in a
				uniform superposition

				if lis is a list of numbers with length = num_qubits, it is
				interpreted as an amplitude vector

OUTPUTS:
-> numpy array of amplitudes representing a quantum state


-----Examples:

mynewstate = create_state(2, [0,3])
print mynewstate

>>> [ 0.707  0.     0.     0.707]

mynewstate = create_state(2, [1/sqrt(2),1/sqrt(2), 0,0])
print mynewstate

>>> [ 0.707  0.707  0.     0.   ]

mynewstate = create_state(2, [1, 2, 3, 4])
print mynewstate

>>> Note thate the state you generated was normalised automatically
>>> [ 0.183  0.365  0.548  0.73 ]

mynewstate = create_state(2, [1,2,3,4,5])
print mynewstate

>>> StandardError: Cannot interpret input of State() creator. Please enter a list of valid amplitudes or positions.

******************
RENORMALISE VECTOR:

renormalise(state)

INPUTS:
			-> state: a numpy array representing the amplitude vector of
								a non-normalised quantum state
OUTPUTS:
			-> the same amplitude vector but normalised to unit length


-----Examples:

my_normalised_state = normalise(state)
print my_normalised_state

>>> [ 0.707  0.     0.     0.707]

******************
FIND OUT IF VECTOR IS NORMALISED:

is_normalised(state)

INPUTS:
			-> state: a numpy array representing the amplitude vector of
								a non-normalised quantum state
OUTPUTS:
			-> True or False (boolean), depending on whether the state is normalised
					(error tolerance is set to 1e-03 = 0.001)

-----Examples:

print is_normalised([ 0.707  0.     0.     0.707])

>>> True

print is_normalised([ 1, 2, 3, 4])

>>> False

***********************
MEASURE A STATE

measure(state, runs, output)

INPUTS:
			-> state: a numpy array representing the amplitude vector of
								a quantum state
			-> runs:  number of times the measurement has to be repeated
								(simulating a repeated preparation of the state before
								each measurement)
			-> output: Type of printout
								'outcomes' - prints a table of the (stochastic) measurement
														outcomes
								'stats' - prints a histogram of the (stochastic) measurement
													outcomes

OUTPUTS:
			-> prints the measurement results


-----Examples:

measure(create_state(4, [0,2,4,6]),2, 'outcomes')

>>> Measurement Results:
>>>   Run    Index  Basis state
>>> -----  -------  -------------
>>>     1        4  |0100>
>>>     2        0  |0000>



measure(create_state(4, [0,2,4,6]),4, 'stats')

>>> Measurement Statistics:
>>>   rel. occ.    Index  Basis state
>>> -----------  -------  -------------
>>>        0.25        2  |0010>
>>>        0.5         4  |0100>
>>>        0.25        6  |0110>





****************************
PRINT FUNCTION

print_me(state, style)

INPUT:
			-> state: a numpy array representing the amplitude vector of
								a quantum state
			-> style: How to print
								-- 'None' or no entry: prints a table of only
										the nonzero basis states
								-- 'full': prints all states
								-- 'amplitudes': prints the amplitude vector only
OUTPUT:
			-> printout of measurements

----Examples:

print print_me(create_state(2, [0]), None)

>>> Quantum State:
>>>   Index    Probability    Amplitude  Basis state
>>> -------  -------------  -----------  -------------
>>>       0            0.5        0.707  |00>
>>>       2            0.5        0.707  |10>



print_me(create_state(2, [0,2]), 'full')

>>>  Quantum State:
>>>   Index    Probability    Amplitude  Basis state
>>> -------  -------------  -----------  -------------
>>>       0            0.5        0.707  |00>
>>>       1            0          0      |01>
>>>       2            0.5        0.707  |10>
>>>       3            0          0      |11>




*******************************
GROVER ITERATION:

grover_iteration(state, marked_pos)

INPUT:
			-> state
			-> marked_pos: list of integers of the marked positions


### Quantum gates

If your README is compliant with Standard-Readme and you're on GitHub, it would be great if you could add the badge. This allows people to link back to this Spec, and helps adoption of the README. The badge is **not required**.

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

To add in Markdown format, use this code:

```
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)
```

### Measurement & stats

To see how the specification has been applied, see the [example-readmes](example-readmes/).

## Contribute

Feel free to dive in! [Open an issue](https://github.com/mariaschuld/qc_simulator/issues/new)

Or send [Maria](mailto:mariaschuld@gmail.com) or [Mark](markfingerhuth@protonmail.com) an email.

## License

GNU (c) Maria Schuld & Mark Fingerhuth
