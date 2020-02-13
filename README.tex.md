# Master-Thesis

**[Master Thesis](https://amslaurea.unibo.it/19642/1/Nicola%20Franco.pdf)**

## Abstract 

A Distributed observer design is described for estimating the state of a continuous-time, input free, linear system. This thesis explains how to construct the local estimators, which comprise the observer inputs and outputs, and it is shown which are the requirements to deal with this structure. Every agent senses an output signal from the system and distributes it
across a fixed-time network to its neighbors. The information flow increases the capability of each agent to estimate the state of the system and uses collaboration to improve the quality of data.
The proposed solution has several positive features compared to recent results in the literature, which include milder assumptions on the network connectivity and the maximum dimension of the state of each observer does not exceed the order of the plant. The conditions are reduced to certain detectability requirements for each cluster of agents in the network, where a cluster is identified as a subset of agents that satisfy specific properties. Instead, the dimension of each observer is reduced to the number of possible observable states of the system, collected by the agent and by the neighbors.

## System Example

Simulations are presented in order to show the effectiveness of the proposed solution. Firstly a generic system with sinusoidal behavior, in which agents cooperate to estimate the plant dynamics. Then in the second section, the properties for the convergence of the observers are showed in a use case relates to the search and rescue problem. Matlab is the tool used
to test the algorithm and to present the results.

Specifically, it is posted a case in which one of the agents is without output measures. Furthermore, the graph is not strongly connected, in fact, one agent is only able to provide his information to the rest of the network and not to receive. Let’s consider the following continuous-time system, 

\begin{equation}
\dot{x} = 
\left[\begin{array}{cccc}0 & 1 & 0 & 0 \\-1 & 0 & 0 & 0 \\0 & 0 & 0 & -2 \\0 & 0 & 2 & 0\end{array}\right]
x,
\end{equation}
