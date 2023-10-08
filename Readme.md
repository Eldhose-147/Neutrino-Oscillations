<a name="br1"></a> 

1 Introduction to Quantum Computing and NMR

1\.1 Quntum bit

In classical computer a bit is the basic unit of classical information.A bit can have two

states either 0 or 1. A qubit on the other hand is the basic unit of quantum information.

The two possible states of a qubit are |0⟩ and |1⟩ but unlike a bit a qubit can exist as a

linear combination of states called a superposition.

|ϕ⟩ = α |0⟩ + β |1⟩

where α and β are complex numbers.

When we measure a qubit we get either the result 0, with probability α<sup>2</sup>, or the result 1,

with probability β<sup>2</sup>. Naturally, α<sup>2</sup> + β<sup>2</sup> = 1, since the probabilities must sum to one. We

can rewrite |ϕ⟩ as :

|ϕ⟩ = e<sup>iγ</sup>(cos(θ/2) |0⟩ + e<sup>iψ</sup>sin(θ/2) |1⟩

Now since e<sup>iγ</sup> is a global phase and can be ignored. Therefore any qubit can be written

in the format:

|ϕ⟩ = (cos(θ/2) |0⟩ + e<sup>iψ</sup>sin(θ/2) |1⟩

In case of a multi qubit system the state can be written as :

|ϕ⟩ = Σ α |ϕ<sup>1</sup>⟩ ⊗ |ϕ<sup>2</sup>⟩ .... ⊗ |ϕ<sup>n</sup>⟩

j

j

j

j

j

where |ϕ<sup>i</sup> ⟩ ϵ(|0⟩ , |1⟩)

j

1\.2 Density Matrix Representation

The density matrix formulation is very useful in describing the state of an ensemble quan-

tum system such as an ensemble of spins in NMR. For an ensemble with the probability

p to be in |ψ ⟩ state, the density operator is given as:

i

i

ρ = Σ p |ψ ⟩ ⟨ψ |

i

i

i

i

where Σ p = 1. If the system has only state |ψ⟩ then the desity operator is given by:

i

i

ρ = |ψ⟩ ⟨ψ|

A density matrix has th following properties:

\1) ρ is hermitian i\.e\. ρ = ρ†

\2) ρ is positive

3)T r[ρ] = 1

\4) No global phase ambiguity

For a pure state state ensemble the T r[ρ<sup>2</sup>] = 1 but for a mixed ensemble the Tr[ρ<sup>2</sup>] < 1.

This can be used to distinguish between a pure and mixed state.

1\.3 Quantum Gates

Quantum Gates are the building blocks of quantum circuits. Basically any unitary matrix

forms a quantum gate. Unlike the gates of a classical computer a quantum gate is

reversible since it is a unitary matrix. It has been shown that a set of gates that consists

of all one-qubit quantum gates and the two-qubit exclusive-OR gate is universal in the

5



<a name="br2"></a> 

sense that all unitary operations can be expressed as compositions of such gates. One

set of such gates is the one dimensional rotation gates R , R , R and the CNOT gate.

X

Y

Z

Some of the basic quantum gates in matrix form are:

\1) Hadamard gate

ꢀ

ꢁ

1

1

1

√

<sub>2</sub> 1 −1

2)Pauli X gate

ꢀ

ꢁ

0 1

1 0

3)Pauli Y gate

ꢀ

ꢁ

0 −i

i

0

4)Pauli Z gate

ꢀ

ꢁ

1

0

0 −1

5)CNOT gate





1 0 0 0

0 1 0 0

0 0 0 1

0 0 1 0













1\.4 Nuclear Magnetic Resonance

Nuclear magnetic resonance (NMR) describes a phenomenon wherein, an ensemble of

nuclear spins precessing in a static magnetic ﬁeld, absorb and emit radiation in the ra-

diofrequency range in resonance with their Larmor frequencies. The total spin of the

nuclei is described by the operator I and the component of the spin along the z-axis is

described by the operator I<sub>Z</sub>.

I |I, m⟩ = I(I + 1) |I, m⟩

I |I, m⟩ = m |I, m⟩

Z

Atomic nuclei with non-zero spin also possess a magnetic dipole moment µ which is

given as :

µ = γ h¯I where γ is the gyromagnetic ratio

n

n

When an atomic nuclei with magnetic dipole moment µ is placed in an external static

magnetic ﬁeld B , the nuclear state |I, m⟩ assume energy state from m=-I,-I+1,...+I.

0

This spliting of the energy is called Zeeman eﬀect and is described bt the hamiltonian:

H = −µ.B = −µ B = −γ h¯B I = −h¯ω I

Z

0

z

0

n

0

Z

L

Z

here ω is the magnitude of lamor frequency which is the frequency with which the nuclei

L

precess about the axis of the magnetic ﬁeld.

ω = γ B

0

L

n

From the above equation we see that the lamor frequency depends on the gyro mag-

netic ratio γ<sub>n</sub>. This quantity is diﬀerent for diﬀerent nuclei and hence the lamor frequency.

6



<a name="br3"></a> 

Under the hamiltonian < I<sub>X</sub> > and < I > oscillates with frequency ω , whereas I

Y

L

Z

remains stationary. The eigenvalue of the Hamiltonian is given by:

E<sub>m</sub> = −mh¯ω<sub>L</sub>

1\.4.1 Resonance Phenomenon

The Larmor frequencies of the nuclear spins in a static magnetic ﬁeld of a few Tesla are

of the order of MHz. The transition between the spin states can be realized by means of

selective RF(radio frequency) pulses applied perpendicular to the quantization axis.

B (t) = 2B cos(Ωt + Φ)i

1

1

where Ω is the frequency of the magnetic ﬁeld. The Hamiltonian is given by :

H<sub>RF</sub> = −µ.B (t) = −γ h¯I [2B cos(Ωt + Φ)

1

n

x

1

The RF Hamiltonian can be treated as a perturbation to the main Zeeman Hamiltonian,

considering that the magnitude of B (t), is much smaller than that of B .Therefore, the

1

0

dominating role is still played by H<sub>Z</sub> and the eﬀect of H<sub>RF</sub> can be determined using stan-

dard time-dependent perturbation theory. Brieﬂy, the result is that, when the frequency

of the RF ﬁeld is close to the Larmor frequency (Ω ≈ ω ), i.e., on resonance, there is

L

transitions between the eigenstates of H given by probability:

Z

P<sub>m−→n</sub> ∝ γ<sup>2</sup>h¯ B | ⟨m| I |n⟩ |

2

2

2

n

1

X

1\.5 NMR Quantum Computing

Nuclear magnetic resonance quantum computing is one of the several proposed approaches

for constructing a quantum computer, that uses the spin states of nuclei within molecules

7



<a name="br4"></a> 

as qubits. In 1997, D. G. Cory and I. L. Chuang independently proposed a NMR quantum

com- puter that can be programmed much like a quantum computer. Their computational

model uses an ensemble quantum computer wherein the results of a measurement are the

expectation values of the observable. This computational model can be realized by NMR

spectroscopy on macroscopic ensembles of nuclear spins. NMR Quantum computing

using liquid sample was introduced known as the liquid state NMR. This uses molecules

having non zero spin like CHCl ( with C<sup>13</sup> and H<sup>1</sup>), Si(CH ) (with one Si<sup>29</sup>, C<sup>13</sup> and

3

3 4

one H<sup>1</sup>)..etc.

1\.5.1 NMR Qubit

In case of an NMR quantum computer a qubit is realised by a spin <sup>1</sup> nuclei. The

2

hamiltonian is given by:

H = −h¯ω I

Z

L

z

Now in case of N qubit sytem we have two nuclei of non zero spin <sup>1</sup> and also a J coupling

2

interaction between them. The J -coupling (also called indirect or scalar coupling) is also

an interaction between the nuclear magnetic dipole moments of neighbor nuclei, but in

this case the interaction is not direct, being mediated by the electron cloud involved in

the chemical bonds between the corresponding atoms.

The Hamiltonian is given by:

H = −Σ<sup>n</sup> h¯ω I<sup>i</sup> + 2πΣ<sub>i<j</sub>h¯JI<sub>Z</sub><sup>i</sup> I<sub>Z</sub><sup>j</sup>

i=1

1

Z

1\.5.2 Initialisation

In NMR quantum computing the initialisation begins with setting the system to a pure

state. For an ensemble of identical nuclei in thermal equilibrium, the population of each

energy level is given by the Boltzmann distribution. For a two-level system

I =

ratio

<sup>1</sup>, with the population

2

n<sub>−</sub>

and

n<sub>+</sub>

of the m=<sup>1</sup> and m=~~−~~<sup>1</sup> levels, respectively the

2

2

n

−

n

\+

<sub>= e</sub>−¯hω<sub>L</sub>/k<sub>B</sub>T

where k<sub>B</sub> is the Boltzmann constant and T is the absolute temperature of the ensemble.

All the energy levels are almost equally populated and the initial state is mixed. Under

the high temperature approximation the initial state of the system is given by:

1

ρ ≈ (I + ϵ∆ρ<sub>eq</sub>)

2<sup>N</sup>

where I is identity ϵ is the purity factor and ρ is the deviation density matrix. The

eq

problem of pure states in NMR can be overcome by preparing a pseudopure state which

is isomorphic to a pure state. An ensemble of a pure state is given by ρ<sub>pure</sub> = |psi⟩ ⟨ψ|

and the corresponding pseudopure state is given by:

1~~−~~ϵ

N

ρ =

eq

I + ϵ |ψ⟩ ⟨ψ|

2

A pseudopure state in NMR can be prepared by several methods such as spatial averag-

ing, temporal averaging and logical labelling.

8



<a name="br5"></a> 

2 Neutrino Oscillations

2\.1 Introduction

Neutrinos are fermionic particles ﬁrst theorized by Wolfgang Pauli and later discovered in

the Cowan–Reines neutrino experiment.They belong to the neutral lepton family, which

makes it incredibly diﬃcult to detect(because they only interact via the weak force) even

though they are the most abundant subatomic particle in the universe.Three ﬂovors of

neutrinos were predicted by the standard model, electron muon and tau neutrinos. The

electron neutrino was discovered by Reines and Cowan by observing electron antineutrinos

released from a nuclear reactor in South Carolina, 26 years after Pauli’s hypothesis.The

second kind of neutrino, muon neutrino, was discovered by Lederman, Schwartz and

Steinberger in 1962 with the ﬁrst accelerator neutrino beam at Brookhaven National

Laboratory. The tau neutrino was only discovered in 2000 in the Fermi laboratory.

The idea of neutrino oscillation ﬁrst put forward by Pontecorvo in 1957, proposed that

neutrino-antineutrino transition may occur in analogy with the neutral kaon mixing.

Later Pontecorvo, Maki, Nakagawa and Sakata speculated that neutrino could change

ﬂavor while in ﬂight, called “neutrino oscillation”, if neutrinos have mass and mixing

exists between ﬂavor and mass eigenstates. In 1960 Raymond Davis observed a deﬁcit in

the electron neutrino produced in the nuclear fusion in sun from the value predicted by

the standard solar model. This could only have been possible if neutrinos changed ﬂavors

as the travel to the earth. Later other experiments showed that the deﬁcit was 70% from

the predicted value.

This shortfall of solar neutrinos was explained by experiments led by this year’s laure-

ates, Arthur McDonald of Queen’s University in Kingston, Canada, and Takaaki Kajita

of the University of Tokyo. The work revealed that the three ﬂavors of neutrino can

interconvert as the particles stream through space. Kajita headed a team working at the

Super-Kamiokande detector in Japan. The detector consists of a tank of 50,000 tons of

water buried underground to shield it from cosmic rays. It detects electron and muon

neutrinos coming not from the Sun but from the Earth’s atmosphere, where neutrinos

are produced by collisions of cosmic rays with atmospheric atoms. Very rarely, these

neutrinos will collide with atomic nuclei in the detector’s water molecules and generate

ﬂashes of light.

Later in 2001 and 2002, McDonald conﬁrmed the same result with the detector called

the Sudbury Neutrino Observatory (SNO) in Sunbury, Canada. These developments

shared the 2015 Nobel Prize in Physics for the discovery of neutrino oscillations, which

shows that neutrinos have mass.

2\.2 PMNS Theory

Pontecorvo-Maki- Nakagawa-Sakata speculated that if neutrinos have mass mixing ex-

ists between the ﬂovor and mass eigenstates.PMNS(Pontecorvo-Maki- Nakagawa-Sakata)

theory states that each of the neutrino ﬂavor basis( electron, muon and tau neutrinos)

can be written as a superposition of three mass eigenstates(m1,m2 and m3). The trans-

formation from neutrino mass basis to ﬂavor basis is achieved by the PMNS matrix.

9



<a name="br6"></a> 

|ψ ⟩ = U

|ψ ⟩

α

PMNS

mi

where α=e,µ,τ labelling the neutrino ﬂavors and i=1,2,3 for the neutrino mass states

The PMNS matrix is characterised by three "mixing angles" θ12, θ13, θ23 and one

CP-violating phase δcp.The PMNS matrix is given by:



 

 



1

0

0

c<sub>13</sub>

0

−s<sub>13</sub>e<sup>δ</sup><sub>cp</sub>

0 0s<sub>13</sub>e<sup>−iδcp</sup>

c<sub>12</sub> s<sub>12</sub>

0

U<sub>PMNS</sub> = 0 c<sub>23</sub> s  

1

0

0

c<sub>13</sub>

 −s

c<sub>12</sub> 0

23

23

12

0 −s<sub>23</sub>

c

0

0

1

The dynamics of neutrino oscillation in vacuum is governed by the Hamiltonian written

in the mass eigenbasis as:





0

0

0

0

H = 0 δm<sup>2</sup> /2E



0

12

0

0

δm<sup>2</sup><sub>13</sub>/2E

where δm<sup>2</sup> = m<sup>2</sup> − m<sup>2</sup>(≈ 7.5 × 10−5eV <sup>2</sup>) , δm<sup>2</sup> = m<sup>2</sup> − m<sup>2</sup>(≈ 2.5 × 10−3eV <sup>2</sup>) and E is

12

2

1

13

3

1

the neutrino energy (E ≈ p(momentum)).

From an initial state |v(0)⟩ the neutrino evolves in time t in matrix form as

|v(t)⟩ = e<sup>−iHt</sup> |v(0)⟩

where H is the Hamiltonian in the ﬂavor basis.

Now we know the hamiltonian in the mass basis and can be converted to the ﬂavor basis

using the the PMNS matrix.

†

H = U<sub>PMNS</sub>H<sub>0</sub>U

PMNS

HU<sup>†</sup>

Thus |v(t)⟩ = e−iU

|v(0)⟩

1

P MNS

P MNS





0

0

0

also written as |v(t)⟩ = U<sub>PMNS</sub> 0 e−iδm<sup>2</sup> t/2E

 U†

|v(0)⟩

12

PMNS

0

0

e<sup>−iδm213t/2E</sup>

Now in case of matter interaction the Hamiltonian is given by: H = H + H

m

1

where H is given by : H = U<sub>PMNS</sub>H<sub>0</sub>U†

PMNS





V (x) 0 0

and H =  0

0 0

1

0

0 0

√

Here V ( x ) is the neutrino weak interaction potential energy V = 2G N (where G

F

e

F

is the Fermi constant, N<sub>e</sub> is the electron number density in the medium,). Here we are

only taking the case of constant matter interaction i.e.

V(x)=const

10



<a name="br7"></a> 

The hamiltonian H is given by;



m











0

0

0

0

V (x) 0 0

H<sub>m</sub> = U<sub>PMNS</sub> 0 δm<sup>2</sup> /2E

 U

†

\+

0

0

0 0

0 0

12

0

PMNS

0

δm<sup>2</sup><sub>13</sub>/2E

2\.3 Why Quantum Simulation ?

The idea of quantum simulations is the use of quantum systems to simulate other quantum

systems in programmable fashion. It is impossible to truly simulate a quantum mechanical

event using a classical universal device. According to the lecture by Feyman in 1981, the

only way to simulate events in nature is by using a system that is quantum mechanical.

Oscillating neutrino beams exhibit quantum coherence over distances of thousands

of kilometers. Precise measurements of parameters in the PMNS framework might lead

to new physics beyond the Standard Model. However it is very diﬃcult to determine

this in neutrino oscillation experiments. Quantum simulation are the best alternatives

to study this phenomenology. Thus if we could encode the state of neutrinos in some n

qubit hilbert space and ﬁnd unitary gates to evolve the quantum state, we could simulate

neutrino oscillation on a quantum computer. However with the present level of technology

quantum simulations are noisy. Thus noise reducing algorithms are also implemented to

get accurate results.

At present many companies like IBM, Microsoft, Google etc .. have made their quan-

tum computers public. Also other quantum computers like NMR quantum computer,

photonic quantum computer facilities are also available.

3 Numerical Simulation

The ﬁrst task of the project was to do the numerical simulation of of neutrino oscillation.

The PMNS matrix has 4 parameters to be speciﬁed.The parameters that are reasonably

well measured are the solar mixing angle θ<sub>12</sub> ≈ 34◦ the reactor mixing angle θ<sub>13</sub> ≈ 8.5◦

and the solar mass splitting δm<sup>2</sup> ≈ 7.5×10−<sup>5</sup>eV <sup>2</sup>. Parameters with well-determined par-

12

tial information are the atmospheric mixing angle θ ≈ 45◦ and the atmospheric mass

23

splitting δm<sup>2</sup> ≈ 2.5×10−<sup>3</sup>eV <sup>2</sup>. However, the phase δ in CP-violation still needs to be dis-

cp

13

covered with signiﬁcant certainties. Thus numerical simulatons for δ = 0, π/2, π, −π/2

was done. The numerical simulation are used to determine

\1) The probability of ﬁnding the various neutrino ﬂavors\.

\2) The distance over which the probability of ﬁnding the a neutrino ﬂavor is maximum

\3) Compare it with the experimental result

\4) The period of neutrino oscillation for various energy E

\5) The eﬀect of matter interaction on neutrino oscillation

3\.1 Vacuum

Survival probability when the initial states were electron,muon and tau neutrinos with

energy E=1GeV as a function of L/E (Where L is the distance given in kilometer and E

is in GeV) were determined.

The initial state of the sytem was taken to be |v(0)⟩ where α = e, µ, τ labelling the

α

neutrino ﬂavors electron, muon and tau respectively. The system was evolved according

to the hamiltonian

11



<a name="br8"></a> 





0

0

0

0

H = U<sub>PMNS</sub> 0 δm<sup>2</sup> /2E

 U

†

12

0

PMNS

0

δm<sup>2</sup> /2E

13

Thus the state of the system at any particular time ’t’ is given by:

|v(t)⟩ = e<sup>−iHt</sup> |v(0)⟩

Now the probability of state being in various ﬂavors is determined by:

| ⟨v | |v ⟩ |<sup>2</sup>

α

t

The graph was plotted for probability vs L/E with α = e, µ, τ.The following are the

results for the three neutrino ﬂavors as initial states.

Figure 1: Oscillation probability P −→ P , α = e, µ, τ for δ = 0

e

α

cp

Figure 2: Oscillation probability P −→ P , α = e, µ, τ for δ = 0

µ

α

cp

12



<a name="br9"></a> 

Figure 3: Oscillation probability P −→ P , α = e, µ, τ for δ = 0

τ

α

cp

Figure 4: Probability of v −→ v for diﬀerent E at L= 295 km

µ

e

3\.2 Matter Interaction

In case of matter interaction the hamiltonian













0

0

0

0

V (x) 0 0

H<sub>m</sub> = U<sub>PMNS</sub> 0 δm<sup>2</sup> /2E

 U

†

\+

0

0

0 0

0 0

12

0

PMNS

0

δm<sup>2</sup><sub>13</sub>/2E

and the state of the system at any particular time ’t’ is given by: |v(t)⟩ = e−iH <sup>t</sup> |v(0)⟩ .

m

Here V(x) is taken to be a constant. The numerical simulation for V (x) = 0, 10−<sup>3</sup>, 10−5

is done with the initial state being v<sub>µ</sub>(muon nutrino)

13



<a name="br10"></a> 

Figure 5: Oscillation probability P −→ P for δ = 0

µ

e

cp

Figure 6: Oscillation Probability P −→ P for δ = 0

µ

µ

cp

14



<a name="br11"></a> 

4 Quantum Simulations

In order to do quantum simulations of neutrino oscillations, the system needs to be

represented in a two-qubit hilbert space. Each of the orthogonal states |0, 0⟩ , |0, 1⟩ , |1, 0⟩

represents the three neutrino ﬂavors - electron, muon, tau respectively. The forth state

|1, 1⟩ is taken to be a sterile neutrino (Sterile neutrinos are hypothetical particles that are

believed to interact only via gravity and not via any of the other fundamental interactions

of the Standard Model) that is considered decoupled to the other three states.

The PMNS matrix in the two-qubit hilbert space can be represented as :

U<sub>PMNS</sub> = R (θ ).R (θ , δ ).R (θ12)

23 23

13 13 cp

12





1

0

0

0

0

0

0

1





c<sub>23</sub> s<sub>23</sub>

R (θ ) = 



23 23





0 −s

c<sub>23</sub>

0

23

0

0







<sup>−iδ</sup>cp

c<sub>13</sub>

0

0 s e

13

0

0

0

1



1

0

0

0

c<sub>13</sub>

0

R (θ , δ ) = 

12



13 13 cp



δ<sub>cp</sub>



−s e

13

0







c

s<sub>12</sub> 0 0

c<sub>12</sub> 0 0



−s

R (θ12 = 

12



12





0

0

0

0

1 0

0 1

where c and s stands for cos and sin respectively.

The hamiltonian in two-qubit hilbert space can be written as :









0

0

0

0

<sup>0 δm2</sup>12<sup>/2E</sup>

0

0

0

0

δm /2E 0

H = 



0



2

13



0

0

0

1

4\.1 Building the Quantum Circuit

To make the quantum circuit, the PMNS matrix needs to be written in terms of the basic

quantum gates such as the controlled-u rotation gates, C-NOT gate and the Pauli X gate.

The controlled-u rotation gate when the target qubit is the LSB(Least signiﬁcant bit):





1

0

0

0

1

0

0





0 cos(θ/2)e<sup>−i(ϕ+λ)</sup>

0

0

−

sin(θ/2)e<sup>−i(ϕ−λ)</sup>

CU<sup>0</sup>(θ, λ, ϕ) = 







0

0

sinθ/2)e<sup>i(ϕ−λ)</sup>

cos(θ/2)e<sup>i(ϕ+λ)</sup>

The controlled-u rotation gate when the target qubit is the MSB (Most signiﬁcant bit):





1 0

0 1

0

0

0

0





CU<sup>1</sup>(θ, λ, ϕ) = 





<sub>0 0 cos(θ/2)e</sub>−i(ϕ+λ) −sin(θ/2)e

−i(ϕ−λ) 

0 0 sinθ/2)e<sup>i(ϕ−λ)</sup>

cos(θ/2)e

i(ϕ+λ)

15



<a name="br12"></a> 

The C-NOT gate when the target is the LSB:

1 0 0 0









0 1 0 0

0 0 0 1

C<sup>0</sup> = 



X

0 0 1 0

The C-NOT gate when the target is the MSB:









1 0 0 0

0 0 0 1

C<sup>1</sup> = 



0 0 1 0

X

0 1 0 0

The Pauli X gate is:

0 1

1 0

ꢀ

ꢁ

X =

Now each of the gates R (θ ), R (θ , δ )andR (θ12) is expressed in terms of

23 23

13 13 cp

12

CU<sup>0</sup>(θ, λ, ϕ), CU<sup>1</sup>(θ, λ, ϕ), C<sup>0</sup> , C<sup>1</sup> and X gates as the follows (note that q is the MSB

X

X

0

and q<sub>1</sub> is the LSB):

\1) R (θ12) = (I ⊗ X)\.CU<sup>1</sup>(−2θ , 0, 0)\.(I ⊗ X)

12

12

\2) R (θ13, δ ) = (X ⊗ I)\.CU<sup>0</sup>(−2θ , δ , −δ<sub>CP</sub> )\.(X ⊗ I)

13

cp

13 CP

\3) R (θ23) = CU<sup>0</sup> \.CU<sup>1</sup>(−2θ , 0, 0)\.CU<sup>0</sup> )

23

X

23

X

16



<a name="br13"></a> 









1

0

2

12

0

0

0

0

x

0 e<sup>−iδm</sup>

t/2E

0

e<sup>−iδm213t/2E</sup>

0

The matrix e−<sup>iH t</sup> = 



(x can be any value) can be

0





0

0

0

0

ꢀ

ꢁ

1

0

written in terms of the P(λ)(where P(λ) =

) gate as follows:

0 e<sup>iλ</sup>

e<sup>iH</sup><sub>0</sub><sup>t</sup> = P(−δm<sup>2</sup> t/2E) ⊗ P(−δm<sup>2</sup> t/2E)

12

13

Therefore the circuit for the unitary matrix U = U<sub>PMNS</sub>.e−<sup>iH</sup><sub>0</sub><sup>t</sup>.U†

is as follows:

PMNS

<sup>−iH</sup>0<sup>t</sup>

.R ( θ ).R ( θ , δ ).R ( θ )

−

−

−

U = R (θ ).R (θ , δ ).R (θ 2).e

23 23

13 13 cp

12

1

12 12 13 13 cp 23 23

Figure 7: Unitary Matrix U

4\.2 IBM Qiskit

IBM has developed some of the most advanced quantum computers in the world.Moreover

they have given free cloud access to their quantum computers through Qiskit framework.

Qiskit is an open-source software development kit for working with quantum computers

at the level of circuits, pulses, and algorithms and running them on prototype quantum

devices on IBM Quantum Experience or on simulators on a local computer. The primary

version of Qiskit uses Python programming language and is well supported on the Jupyter

Notebook environment.

4\.2.1 QasmSimulator

The QasmSimulator backend is designed to mimic an actual IBM quantum computer.

The quantum measurements in a real IBM quantum computer has a lot of noise in the

result. But the QasmSimulator can simulate quantum circuits both ideally and subject to

noise modeling. The circuit in ﬁgure 5 was run in the QasmSimulator to verify the circuit

and compare it with the numerical simulations. The following results were obtained.

17



<a name="br14"></a> 

Figure 8: Oscillation probability P −→ P , α = e, µ, τ for δ = 0

e

α

cp

Figure 9: Oscillation probability P −→ P , α = e, µ, τ for δ = 0

µ

α

cp

18



<a name="br15"></a> 

Figure 10: Oscillation probability P −→ P , α = e, µ, τ for δ = 0

τ

α

cp

Figure 11: Oscillation probability P −→ P , α = e, µ, τ for δ = π

µ

α

cp

19



<a name="br16"></a> 

The QasmSimulator measurements matches exactly with the theoritical results. This

proves that the quantum circuit built is accurate.

4\.2.2 Real Quantum Simulations

The quantum circuit was run on the "ibmq-manila" backend. The results obtained had

a lot of noise which made the measurements deviate from thetheoritical results. Thus

the "CompleteMeasFitter" and "complete-meas-cal" functions were imported from the

qiskit.ignis.mitigation package. This reduced the error by signiﬁcant amount. The fol-

lowing results were obtained from the real quantum computer simulations.

Figure 12: Oscillation probability P −→ P , α = e, µ, τ for δ = 0

e

α

cp

20



<a name="br17"></a> 

Figure 13: Oscillation probability P −→ P , α = e, µ, τ for δ = 0

µ

α

cp

Figure 14: Oscillation probability P −→ P , α = e, µ, τ for δ = 0

τ

α

cp

21



<a name="br18"></a> 

Figure 15: Oscillation probability P −→ P , α = e, µ, τ for δ = π

τ

α

cp

5 Conclusion

During my 2-month internship, I gained valuable experience in the ﬁeld of quantum

simulation. Neutrino oscillation is a quantum phenomenon that has been studied for

many decades. But the cost and time required to do neutrino oscillation experiments is

substantial. Quantum computers oﬀer a promising solution to this problem. Quantum

computers are able to represent and manipulate quantum states with exponentially more

precision than classical computers. This means that quantum computers could be used

to simulate neutrino oscillations with much greater accuracy than classical computers.

I learned about the theoretical foundations of neutrino oscillations and the challenges

involved in simulating them on a quantum computer. I was able to make quantum circuit

to simulate the three ﬂavor neutrino oscillation. The quantum circuit was based on the

PMNS theory that explains neutrino oscillations. The simulations were ﬁrst done on

Qiskit, the framework to access the IBM quantum computer. This helped me in learn

about building quantum circuits for any quantum system’s. I was also able to learn Qiskit

in detail and learnt the various packages in it, especially for error correction.

The most valuable part of my internship was to work in the NMR quantum computing

facility available in the lab. I learnt the working principle NMR and about how a qubit

can be realised in an NMR quantum computer. I was also able to help in making quantum

circuits for the neutrino oscillation simulation in the NMR quantum computer. Working

in the lab also helped me realise limitations in NMR quantum computer and the methods

used to overcome the limitations in real-time.

22



<a name="br19"></a> 

6 References

[1] Ioannisian, A. and Pokorski, S. (2018) ‘Three neutrino oscillations in matter’, Physics

Letters B, 782, pp. 641–645. doi:10.1016/j.physletb.2018.06.001.

[2] Freund, M. (2001) ‘Analytic approximations for three neutrino oscillation parameters

and probabilities in matter’, Physical Review D, 64(5). doi:10.1103/physrevd.64.053003.

[3] Giganti, C., Lavignac, S. and Zito, M. (2018) ‘Neutrino oscillations: The rise of the

pmns paradigm’, Progress in Particle and Nuclear Physics, 98, pp. 1–54. doi:10.1016/j.ppnp.2017.10.001

[4] Argüelles, C.A. and Jones, B.J. (2019) ‘Neutrino oscillations in a quantum processor’,

Physical Review Research, 1(3). doi:10.1103/physrevresearch.1.033176.

[5] Molewski, M.J. and Jones, B.J.P. (2022) ‘Scalable qubit representations of neutrino

mixing matrices’, Physical Review D, 105(5). doi:10.1103/physrevd.105.056024.

[6] Nguyen, H.C. et al. (2023) ‘Simulating neutrino oscillations on a superconducting

qutrit’, Physical Review D, 108(2). doi:10.1103/physrevd.108.023013.

[7] Nielsen, M.A. and Chuang, I.L. (2022) Quantum Computation and Quantum Infor-

mation. Cambridge: Cambridge University Press.

[8] Oliveira, I.S. (2007) NMR Quantum Information Processing. Amsterdam: Elsevier.

23


