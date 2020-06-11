# rsiena  <img src="inst/rsienalogo.png" align="right" width="150"/>

![CRAN/METACRAN](https://img.shields.io/cran/l/RSiena)
![CRAN/METACRAN](https://img.shields.io/cran/v/RSiena)
![GitHub issues](https://img.shields.io/github/issues-raw/snlab-nl/rsiena)
![GitHub All Releases](https://img.shields.io/github/downloads/snlab-nl/rsiena/total)

## About

SIENA is a program for the statistical analysis of network data, with the focus on social networks.
Networks here are understood as entire (complete) networks, not as personal (egocentered) networks: it is assumed that a set of nodes (social actors) is given, and all ties (links) between these nodes are known - except perhaps for a moderate amount of missing data.
SIENA is designed for analyzing various types of data as dependent variables:

### Longitudinal network data:
This refers to repeated measures of networks on a given node set (although it is allowed that there are some changes in the node set). Models can be specified with actor-oriented as well as tie-oriented dynamics; but mainly the former.

### Longitudinal data of networks and behavior:
This is like longitudinal network data, but in addition there are one or more changing nodal variables that are also treated as dependent variables, and referred to as behavior. The network will influence the dynamics of the behavior, and the behavior will influence the dynamics of the network. In other words, this is about the co-evolution of networks and behavior.

### Cross-sectional network data.
'Cross-sectional' means that only one observation is available. This method uses exponential random graph models ('ERGMs'), also called p* models.
The ERG model is implemented in SIENA version 3, but not any more in version 4 (RSiena). It is still available, but no longer maintained.
The name SIENA stands for Simulation Investigation for Empirical Network Analysis.

## Migration in progress...

We are migrating RSiena development and releases to this repository.

Please stay tuned...
