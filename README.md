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
	- [Initialising quantum states](#initialising-quantum-states)
	- [Amplitude normalisation](#amplitude-normalisation)
	- [Printing quantum states](#printing-quantum-states)
	- [Quantum gates](#quantum-gates)
	- [Special functions](#special-functions)
	- [Measurement & statistics](#measurement-and-statistics)
- [Contribute](#contribute)
- [License](#license)

## Background

This library was designed and written by [@mariaschuld](https://github.com/mariaschuld) and [@markf94](https://github.com/markf94) from the [Centre for Quantum Technology](http://quantum.ukzn.ac.za/) at the [University of KwaZulu-Natal](www.ukzn.ac.za). It is meant as an educational tool to simulate quantum computations. The main aim was to come up with a syntax that is easy to read & understand for beginners. Additionally, the developed library is scalable to any number of qubits only limited by the amount of available RAM in your computer. The library includes initialisation of arbitrary quantum states and amplitude distributions, Bloch sphere projections, a variety of single-qubit gates, all-to-all connected controlled-U and Toffoli (CCNOT) gates. With easy syntax and without visible for-loops the simulation can be repeated many times with clear output of the resulting statistics.

## Install

This project uses [Python](http://python.org/) and is compatible with versions 2.7, 3.4, 3.5 and 3.6. The following packages are required: [numpy](http://www.numpy.org), [QuTiP](http://qutip.org), [random](https://docs.python.org/2/library/random.html), [collections](https://docs.python.org/2/library/collections) and [tabulate](https://pypi.python.org/pypi/tabulate). **Go check them out if you don't have them locally installed.**

Just download this repository as a ZIP file or clone it via your favourite shell. Unzip the folder and you are ready to go. No additional installation required.

## Usage



### Initialising quantum states

General syntax:

```
create_state(num_qubits, lis)
```

INPUTS:

-> `num_qubits`: number of qubits of the quantum system.

-> `lis`: if `lis` is a list of integers with length < num_qubits, it is
				interpreted as a list of positions of nonzero amplitudes in a
				uniform superposition. If `lis` is a list of numbers with
				length = num_qubits, it is interpreted as an amplitude vector.

OUTPUT:

-> numpy array of amplitudes representing a quantum state.

--------------------------------------------------------------------------------

**Example 1:**

Suppose you want to initialise a uniform superposition over the 0th and 3rd
two-qubit states |00> and |11>. Therefore, we interpret the second input parameter
`lis` as a list of positions of nonzero amplitudes in a
uniform superposition. The syntax is as follows:

```
mynewstate = create_state(2, [0,3])
print mynewstate
```

This outputs the desired amplitude vector `[ 0.707  0.     0.     0.707]`.

--------------------------------------------------------------------------------

**Example 2:**

Suppose you want a superposition over the first two two-qubit states |00> and |01>.
Instead of using the input variable `lis` as a list of qubit positions we can fill
it with four amplitudes corresponding to the four two-qubit states.
To initialise a uniform superpositon over |00> and |01>, the following syntax can
be used:
```
mynewstate = create_state(2, [1/sqrt(2),1/sqrt(2), 0,0])
print mynewstate
```
This outputs the desired amplitude vector `[ 0.707  0.707  0.     0.   ]`.

--------------------------------------------------------------------------------

**Example 3:**

```
mynewstate = create_state(2, [1, 2, 3, 4])
print mynewstate
```

Console output:
```
Note thate the state you generated was normalised automatically
[ 0.183  0.365  0.548  0.73 ]
```

--------------------------------------------------------------------------------

**Example 4:**

```
mynewstate = create_state(2, [1,2,3,4,5])
print mynewstate
```

```
StandardError: Cannot interpret input of State() creator.
Please enter a list of valid amplitudes or positions.
```

--------------------------------------------------------------------------------

--------------------------------------------------------------------------------

### Amplitude normalisation

**Renormalising an amplitude vector**

```
renormalise(state)
```

INPUTS:

-> `state`: a numpy array representing the amplitude vector of a non-normalised quantum state

OUTPUT:

-> the same amplitude vector but normalised to unit length

--------------------------------------------------------------------------------

**Example 1:**

```
my_normalised_state = normalise(state)
print my_normalised_state
```
This outputs the normalised amplitude vector `[ 0.707  0.     0.     0.707]`.

--------------------------------------------------------------------------------

**Finding out if an amplitude vector is normalised**

```
is_normalised(state)
```

INPUTS:

-> `state`: a numpy array representing the amplitude vector of
								a non-normalised quantum state

OUTPUT:

-> `True` or `False` (boolean), depending on whether the state is normalised
					(error tolerance is set to 1e-03 = 0.001)

**Example 1:**

```
print is_normalised([ 0.707  0.     0.     0.707])
```
This outputs `True`.


```
print is_normalised([ 1, 2, 3, 4])
```
This outputs `False`.

--------------------------------------------------------------------------------

--------------------------------------------------------------------------------

### Printing quantum states

```
print_me(state, style)
```

INPUT:

-> `state`: a numpy array representing the amplitude vector of a quantum state

-> `style`: How to print

	> -- `None` or no entry: prints a table of only the nonzero basis states

	> -- `full`: prints all states

	> -- `amplitudes`: prints the amplitude vector only

OUTPUT:

-> printout of measurements

--------------------------------------------------------------------------------

**Example 1:**

```
print print_me(create_state(2, [0]), None)
```

```
Quantum State:
  Index    Probability    Amplitude  Basis state
-------  -------------  -----------  -------------
      0            0.5        0.707  |00>
      2            0.5        0.707  |10>
```

```
print_me(create_state(2, [0,2]), 'full')
```

```
Quantum State:
Index    Probability    Amplitude  Basis state
------  -------------  -----------  -------------
    0            0.5        0.707  |00>
    1            0          0      |01>
    2            0.5        0.707  |10>
    3            0          0      |11>
```

### Quantum gates



### Special functions

#### Grover iteration

```
grover_iteration(state, marked_pos)
```

INPUT:

-> `state`: a numpy array representing the amplitude vector of a quantum state

-> `marked_pos`: list of integers of the marked positions

### Measurement and statistics

```
measure(state, runs, output)
```

INPUTS:

-> `state`: a numpy array representing the amplitude vector of a quantum state

-> `runs`:  number of times the measurement has to be repeated
(simulating a repeated preparation of the state before each measurement)

-> `output`: Type of printout:

-> `outcomes` - prints a table of the (stochastic) measurement outcomes <-

-> `stats` - prints a histogram of the (stochastic) measurement outcomes <-

OUTPUT:

			-> prints the measurement results

--------------------------------------------------------------------------------

**Example 1:**

```
measure(create_state(4, [0,2,4,6]),2, 'outcomes')
```

```
Measurement Results:
  Run    Index  Basis state
-----  -------  -------------
    1        4  |0100>
    2        0  |0000>
```
--------------------------------------------------------------------------------

**Example 2:**

```
measure(create_state(4, [0,2,4,6]),4, 'stats')
```

```
Measurement Statistics:
  rel. occ.    Index  Basis state
-----------  -------  -------------
       0.25        2  |0010>
       0.5         4  |0100>
       0.25        6  |0110>
```

## Contribute

Feel free to dive in! [Open an issue](https://github.com/mariaschuld/qc_simulator/issues/new)

Or send [Maria](mailto:mariaschuld@gmail.com) or [Mark](markfingerhuth@protonmail.com) an email.

## License

GNU (c) Maria Schuld & Mark Fingerhuth
