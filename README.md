# rsiena  <img src="https://raw.githubusercontent.com/stocnet/rsiena/main/inst/rsienalogo.png" align="right" width="150"/>

![CRAN/METACRAN](https://img.shields.io/cran/l/RSiena)
![CRAN/METACRAN](https://img.shields.io/cran/v/RSiena)
![GitHub R package version](https://img.shields.io/github/r-package/v/stocnet/rsiena)
![GitHub issues](https://img.shields.io/github/issues-raw/stocnet/rsiena)
![GitHub All Releases](https://img.shields.io/github/downloads/stocnet/rsiena/total)
![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)


SIENA is a program for the statistical analysis of network data, with the focus on social networks.
Networks here are understood as entire (complete) networks, not as personal (egocentered) networks: 
it is assumed that a set of nodes (social actors) is given, and all ties (links) between these nodes are known - 
except perhaps for a moderate amount of missing data.
The name SIENA stands for Simulation Investigation for Empirical Network Analysis.
The R package is called RSiena.

## Installation

For most people, the best way to install RSiena is to install the latest version from CRAN:

```r
install.packages("RSiena")
```

The latest binary release on GitHub will have newer features:

```r
# On Windows:
install.packages("https://github.com/stocnet/rsiena/releases/latest/download/RSiena.zip", repos = NULL)

# On Linux
install.packages("https://github.com/stocnet/rsiena/releases/latest/download/RSiena.tar.gz", repos = NULL)

# On Mac
install.packages("https://github.com/stocnet/rsiena/releases/latest/download/RSiena.tgz", repos = NULL)
```

To install the source version from GitHub install the `{remotes}` package and then run the following. NB: this requires compilation of `C++` source files so it may take some time.

```r
# latest version
remotes::install_github("stocnet/rsiena@main")

# development version
remotes::install_github("stocnet/rsiena@develop")
```

## Data types

SIENA is designed for analyzing various types of data as dependent variables:

### Longitudinal network data:
This refers to repeated measures of networks on a given node set (although it is allowed that there are some changes in the node set). Models can be specified with actor-oriented as well as tie-oriented dynamics; but mainly the former.

Practical restrictions are that the number of actors should not be too large; a few hundred already is pretty large.

### Longitudinal data of networks and behavior:
This is like longitudinal network data, but in addition there are one or more changing nodal variables that are also treated as dependent variables, and referred to as behavior. The network will influence the dynamics of the behavior, and the behavior will influence the dynamics of the network. In other words, this is about the co-evolution of networks and behavior.

### Multivariate and two-mode networks:
Network data sets can be multivariate, i.e., be composed of multiple networks on the same node set.
Some or all of these networks can be two-mode networks. The restriction is that the first mode must be the same for all networks; the first mode is defined as the set of actors. The second mode node sets are allowed to differ across the various networks in a given data set. For such multivariate data sets, the model again is about the co-evolution of several networks; and this may be combined with behavior. 

## Manual:  
There is an extensive [manual](https://www.stats.ox.ac.uk/~snijders/siena/RSiena_Manual.pdf) which is complementary to the help pages in the package.

## Migration in progress...

We are migrating RSiena development and releases to this repository.

The main website is still [here](http://www.stats.ox.ac.uk/~snijders/siena/) for the time being,
however we are currently migrating many resources to [this website](http://stocnet.github.io/rsiena/),
and you can find [a wiki here](https://github.com/stocnet/rsiena/wiki) that holds much of the information on the original website,
including background on SAOMs and RSiena, and links to teaching materials, literature, and contributing people and projects.

## Installation

### From binary

Perhaps the easiest way to install RSiena is by installing a compiled binary.
Binaries for all major OSes -- Windows, Mac, and Linux -- 
can be found by clicking on the latest release for your OS [here](https://github.com/stocnet/rsiena/releases/latest).
For Windows you should use the `RSiena.zip`, for macOS it should be `RSiena.tgz`, and for Linux `RSiena.tar.gz`.

Once the file has been downloaded, install the binary appropriate for your Operating System like so:

`install.packages("~/Downloads/RSiena.zip", repos = NULL)`

amending the file suffix as necessary.

### From source

To install from source the latest main version of RSiena from Github, 
please install the `{remotes}` package from CRAN and then enter into the console:

`remotes::install_github("stocnet/rsiena", ref = "main")`

The development version of RSiena can be similarly installed as:

`remotes::install_github("stocnet/rsiena@develop")`

## Citation

To cite the RSiena package in publications use:

> Ruth M. Ripley, Tom A. B. Snijders, Zsofia Boda, Andras Voros, and Paulina Preciado (2023). Manual
> for Siena version 4.0. R package version 1.3.22.
> https://www.cran.r-project.org/web/packages/RSiena/.

A BibTeX entry for LaTeX users is

```bib
@TechReport{,
  title = {Manual for {Siena} version 4.0},
  author = {Ruth M. Ripley and Tom A. B. Snijders and Zsofia B'{o}da and Andr'{a}s V"{o}r"{o}s and Paulina Preciado},
  year = {2023},
  institution = {Oxford: University of Oxford, Department of Statistics; Nuffield College},
  note = {R package version 1.3.22. https://www.cran.r-project.org/web/packages/RSiena/},
}
```

For more references, see https://www.stats.ox.ac.uk/~snijders/siena/. 
