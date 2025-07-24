⚠️ This wiki is still under construction. Original documentation for the Hibridon package can be found [here](https://www2.chem.umd.edu/groups/alexander/hibridon/hib43/)


![image](./resources/hibridon-logo.png)


1. [What is Hibridon?](#what-is-Hibridon)
2. [How do I get the Hibridon source code?](#how-do-I-get-the-Hibridon-source-code)
3. [How does it work?](#how-does-it-work)
4. [How do I use it?](#how-do-i-use-it)
5. [I want to contribute](#i-want-to-contribute)

## What is Hibridon?

Hibridon is a program package to solve the close-coupled equations which occur in the quantum treatment of inelastic atomic and molecular collisions.

Gas-phase scattering, photodissociation, collisions of atoms and/or molecules with flat surfaces, and bound states of weakly-bound complexes can be treated.s

Although the documentation present on this Wiki contains some general formal material, we assume here that you have some background with the quantum theory of inelastic scattering.

#### History of Hibridon
The Hibridon package was written by M. H. Alexander, D. E. Manolopoulos, H.-J. Werner, B. Follmeg and P. J. Dagdigian, with contributions by P. F. Vohralik, D. Lemoine, G. Corey, R. Gordon, B. Johnson, T. Orlikowski, A. Berning, A. Degli-Esposti, C. Rist, B. Pouilly, G. van der Sanden, M. Yang, F. de Weerd, S. Gregurick, J. Klos and F. Lique.

Support for its development was provided by grants (to MHA) from the U. S. National Science Foundation, the U. S. Army Research Office, the Office of Scientific Research of the U.S. Air Force, and (to HJW) from the German Fonds der Chemische Industrie.


## How do I get the Hibridon source code?

The Hibridon package source code and the instructions to compile it can be found on the homepage of the [Hibridon GitHub repository](https://github.com/hibridon/hibridon).

## How does it work?

Execution of the Hibridon package is controlled by the Hibridon driver.

- Description of the commands to the driver can be found [here](Commands)

- The scattering calculation is controlled by:
  - system independent [parameters](Parameters)
  - system dependent [parameters](Parameters)
  - logical [flags](Flags)
  - input and output [files](Files), and job [names](Files)

- At present Hibridon can treat a number of different types of collision [systems](Systems). This is made possible by the presence of multiple system specific subroutines:
  - [BASIS](basis) subroutines which define the Hamiltonian matrix
  - [SYSDAT](SYSDAT) subroutines which control the input of system specific data
  - System specific [POT](POT) subroutines which determine the distance dependence of the potential matrix [V(R)](Close-coupled-equations).

## How do I use it?
A number of [examples](Examples) are available to help you learn how to use the Hibridon package.

## I want to contribute
Hibridon v5 is an open-source software under the [GNU General Public Licence v3](https://www.gnu.org/licenses/gpl-3.0.html).
Everyone is free to use, distribute, modify, and distribute modified versions of this software.

The GitHub repository is public. We encourage contributors to submit their modifications by opening Pull Requests on this repository. Several useful guides on how to collaborate using GitHub can be found [here](https://lab.github.com).
