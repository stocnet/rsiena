# Submitting the package to CRAN.   
  
This text contains the main things to be noted and remembered.   


## Informative websites

- [CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html) 
is a compact description of what should be done for submitting a package.

- [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html) 
is an extensive description of requirements for R packages.

-  [RSiena website at CRAN](https://cran.r-project.org/web/packages/RSiena/)    
There will be questions about how any issues arising here were addressed.   
In particular,
[RSiena check results](https://CRAN.R-project.org/web/checks/check_results_RSiena.html)

## What you need on your machine (in addition to the package itself)
- R-release and R-devel should be installed, with all dependencies. 
- For the last step, also the reverse dependencies (not many),
  with all their dependencies, should be installed for R-devel.

## Checks
1. Of course checks should be done on your own machine, with R-devel.
2. For other operating systems, use package rhub.   
   When checking the tarball created at my own machine (Windows) from
   my local Github directory, I got warning messages   
   "Found the following Makefile(s) with CR or CRLF line endings: src/Makevars   
    Some Unix 'make' programs require LF line endings."    
   These were not obtained any more when I had copied files from 
   https://github.com/snlab-nl/rsiena/archive/refs/tags/v1.2.34.zip
   and made a tarball from those (after making some additional changes). 
   Apparently, there is something with line endings that is not correct
   on my own Windows setup for GitHub.
3. When all checks passed well (Windows, Mac, Unix), the reverse dependencies
   should be checked.