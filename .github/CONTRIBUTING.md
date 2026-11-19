# Contributing

Contributions to `rsiena`, whether in the form of issue identification, bug fixes, new code or documentation are encouraged and welcome:

* [Submit an issue](#issues)
* [Fix a bug or implement new features](#adding-new-code)
* [Document existing code](#documentation)
* [Make a release or submit to CRAN](#releases-code-maintainer-only) (code maintainer only)

Please note that the `rsiena` project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

## Issues

Please use the issue tracker on GitHub to identify problems or suggest new functionality, before submitting changes to the code.
We use issues to identify bugs and tasks, discuss feature requests, and to track implementation of changes.

When submitting an issue, please provide at least a 'Type' label that best describes what the issue is about.
The most useful issues are ones that precisely identify a bug, 
or propose a test that should pass but instead fails.

## Adding new code
Independent or assigned code contributions are most welcome.
When you write new R code, please follow the [tidyverse style guide](https://style.tidyverse.org/index.html).
The packages `lintr` and `goodpractice` help to find problems, and `styler` fixes some of the formatting:

```r
# basic lintr checking
lintr::lint_package(path = ".")

# goodpractice checks. Exclude line length 80
goodpractice::gp(path = ".", checks = goodpractice::all_checks()[-c(8)])

# styler fixes some of the styling issues
styler::style_file("filePath")
```

If you develop new code in `C++`, please follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
The C++ code is in `src/`. The package uses `configure` (made from `configure.ac` with autoconf) to compile it.

### Branches
`rsiena` code development is on GitHub, at <https://github.com/stocnet/rsiena>.
You can learn more about Git and GitHub in [R Packages](https://r-pkgs.org/software-development-practices.html).

We use two **main branches** in this project:

1. The `main` branch holds the released versions of the package.
Each merge into `main` makes a new release (see [Releases](#releases)).
Do not push to `main` directly. Use a pull request from `develop` or from a hotfix branch.

2. The `develop` branch holds the latest development stage.
Contributors can submit minor changes to this branch that improve existing functionality.
Put new features that can break existing functionality in a feature branch.

We use two types of **supporting branches**:

3. *Feature branches* are for new functionality. Start them from `develop`.
When the feature is complete, merge the branch into `develop` with a pull request.
Delete the branch if the feature is abandoned.

4. *Hotfix branches* are for fixes to severe bugs in a release. Start them from `main`.
Branch names start with `hotfix-`.
A hotfix lets the code maintainer fix a release without the (possibly unstable) changes in `develop`.

This branching model is based on <https://nvie.com/posts/a-successful-git-branching-model/>.

### Develop branch (minor changes to existing functionality)
1. Switch to the **develop** branch: `git switch develop`.
2. Get the latest version: `git pull origin develop`.
3. Make your changes.
4. Add a short description of the changes to `NEWS.md` (see [NEWS](#news)).
5. Stage your changes: `git add PATH-TO-CHANGED-FILE`.
6. Commit your changes [with an appropriate message](#commit-messages): `git commit -m "DESCRIPTION"`.
7. Push your commit: `git push origin develop`.

### Feature branches (new functionality)
1. Make the branch: `git switch -c myfeature develop`.
2. Push it: `git push -u origin myfeature`.
3. Open a pull request from `myfeature` to `develop` on GitHub.
4. After the merge, delete the branch: `git branch -d myfeature`.

### Hotfix branches (critical bugs in a release)
1. Make the branch: `git switch -c hotfix-1.6.16 main`.
2. Fix the bug, add a test, and raise the patch digit of `Version` in `DESCRIPTION` (see [Versioning](#versioning)).
3. Open a pull request from the hotfix branch to `main`. The merge makes the release.
4. Merge `main` back into `develop`: `git switch develop` and then `git merge main`.
5. Delete the hotfix branch: `git branch -d hotfix-1.6.16`.

### Releases (code maintainer only)
The GitHub Actions workflows make a release automatically:

- `prchecks.yml` runs on each pull request to `main`.
  It builds and checks the package on macOS, Windows and Ubuntu (R-release), and on Ubuntu (R-devel).
  It also fails if `README.Rmd` and `README.md` do not change together.
- `pushrelease.yml` runs on each push to `main`.
  It checks the package again, and then tags the commit as `v` plus the `Version` in `DESCRIPTION` (for example `v1.6.15`).
  It makes a GitHub release with the binaries and `RSiena_Manual.pdf`, and then deploys the website to the `gh-pages` branch.

To make a release:

1. On `develop`, raise `Version` and `Date` in `DESCRIPTION`, and in the table in `man/RSiena-package.Rd`.
2. Make sure that the top entry in `NEWS.md` has the same version.
3. If you changed `README.Rmd`, knit it again (see [Documentation](#documentation)).
4. Open a pull request from `develop` to `main`, and wait for the checks to pass.
5. Merge the pull request.

Do not make the tag by hand.
If the tag already exists, the release job stops with an error.
Then raise `Version` and push again.

### Commit messages
Commits that relate to existing issues should reference the updated status of those issues, and mention the issue number (preceded by a hash symbol: #) in the commit description:

``` Resolved #31 by adding a new function that does things, also updated documentation ```

Where the issue hash (i.e. #31) is preceded by `resolve`, `resolves`, `resolved`, `close`, `closes`, `closed`, `fix`, `fixes`, or `fixed` (capitalised or not), the status of the issue(s) mentioned is updated automatically.
Our current syntactical standard is to mention the issue first and then provide a short description of what the committed changes do in relation to that issue.
Any ancillary changes can be mentioned after a comma.

It should all be written in a single line, like so: `{verb} {issue} {describe main action/changes}, {additional actions/changes}`.

### Testing
The tests are plain R scripts in [tests/](../tests).
`R CMD check` runs each script.
A script fails if it gives an error, so use `stopifnot()` to test a result.
If a script has a saved output file (for example `newparallel.Rout.save`),
`R CMD check` also compares the output with that file.

Make sure that all tests pass before you commit changes to existing code.
To run all tests, build and check the package:

```sh
R CMD build .
R CMD check RSiena_1.6.15.tar.gz
```

Use the version number from `DESCRIPTION` in the file name.
Use `R CMD check --as-cran` before a release.

When you fix a bug or write a new function, add a test.
Put a test for one topic in one file, named after that topic (for example `tests/set_effect.R`).
Make the test fail on the old code and pass on the new code.
Use the data sets in the package (for example `s501`, `s502`, `s503`, `s50a`) when this is possible.

## Documentation
The help pages are hand-written Rd files in [man/](../man).
The package does not use roxygen2, so edit the Rd files directly.
If you add or change an effect, also update the effects section of the manual in `docs/manual/RSiena_Manual.tex`.

The website at <https://stocnet.github.io/rsiena/> is made with pkgdown from `pkgdown/_pkgdown.yml`.
Add a new help page to the `reference:` section of that file.
The articles are R Markdown files in `vignettes/articles/`.
They are not part of the package on CRAN.
To build the website locally:

```r
pkgdown::build_site(new_process = FALSE)
```

The site goes to `pkgdown/site/`, which Git ignores.

`README.md` is made from `README.Rmd`. Do not edit `README.md` directly.
After you change `README.Rmd`, knit it again and commit both files:

```r
rmarkdown::render("README.Rmd", output_options = list(html_preview = FALSE))
```

Please contact us directly with proposals for updating the documentation, or submit an issue if existing documentation is unclear.

## NEWS
Add each change to the top entry of `NEWS.md`.
Each version has a heading such as `# RSiena 1.6.15 (2026-11-19)`.
Put the changes under these second-level headings, in this sequence, and use only the headings that you need:
`## Effects`, `## Functionality`, `## Bug corrections`, `## Messages`, `## Documentation`, `## Coding`, `## Tests`, `## Package`.

## Versioning
The `rsiena` package uses [semantic versioning](https://semver.org/).
The version has the format Major.Minor.Patch:

- Raise the major digit for changes that are not compatible with earlier code.
- Raise the minor digit for new functionality that is compatible with earlier code.
- Raise the patch digit for bug corrections that are compatible with earlier code.

## Submitting to CRAN (code maintainer only)

### Informative websites

- [CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html)
  is a compact description of what should be done for submitting a package.
- [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
  is an extensive description of requirements for R packages.
- The [RSiena page at CRAN](https://cran.r-project.org/web/packages/RSiena/),
  in particular the
  [RSiena check results](https://CRAN.R-project.org/web/checks/check_results_RSiena.html).
  CRAN will ask how any issues reported there were addressed.

### What you need on your machine

- R-release and R-devel, each with all dependencies of the package.
- For the last step, the reverse dependencies (there are not many),
  with all their dependencies, installed for R-devel.

### Checks

1. Check the package on your own machine with R-devel.
2. Check the package on other operating systems with the `rhub` package.
3. When the checks pass on Windows, macOS and Linux,
   check the reverse dependencies.

On Windows, a tarball built from a local GitHub clone can give the warning
"Found the following Makefile(s) with CR or CRLF line endings: src/Makevars".
This comes from the line endings of the local Git setup.
A tarball built from the files of a GitHub release archive
(for example `https://github.com/stocnet/rsiena/archive/refs/tags/v1.2.34.zip`)
does not give this warning.

## For developers using macOS
To compile the package from source on macOS, you need:

- The Xcode Command Line Tools (`xcode-select --install`).
- The gfortran compiler that CRAN recommends for your version of R (see <https://mac.r-project.org/tools/>).
- autoconf and automake, for example from Homebrew (`brew install autoconf automake`).

Then install the package from the repository folder with `R CMD INSTALL .`.

If the installation fails with an error such as "symbol not found" after you change C++ code or switch branches,
old object files are the cause. Remove them and install again with `R CMD INSTALL --preclean .`.

If the compiler cannot find a library (for example `-lgfortran` or `-lintl`),
check the paths in `~/.R/Makevars` and in the `Makeconf` file of your R installation.
For `-lintl`, install gettext (`brew install gettext`).
