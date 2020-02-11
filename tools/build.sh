#!/bin/bash -x

[ -z "${R_HOME}" ] && R_HOME=`R RHOME`

## ensure that NAMESPACE contains load directive
echo -e "# Generated by roxygen2: do not edit by hand\nuseDynLib(OncoBayes2, .registration = TRUE)\n" > NAMESPACE

[ -z "$R_VERSION" ] && R_VERSION=`${R_HOME}/bin/R --version | head -n 1 | cut -d" " -f 3 | sed 's#.[0-9]$##'`

## compile dll
if [ "$R_VERSION" == "3.4" ];
then
   "${R_HOME}/bin/R" --slave -e 'library(devtools); devtools::compile_dll()' ## use devtools in old R environments
else
    ## use pkgbuild in new R environments
   "${R_HOME}/bin/R" --slave -e 'library(pkgbuild); pkgbuild::compile_dll()' > /dev/null 2>&1
fi

## create internal data-sets
"${R_HOME}/bin/R" --slave --file=tools/make-ds.R

## create documentation Rd files
#"${R_HOME}/bin/R" --slave -e 'library(devtools); devtools::load_all(); library(roxygen2); roxygen2::roxygenize()'
"${R_HOME}/bin/R" --slave -e 'library(roxygen2); roxygen2::roxygenize()'

##"${R_HOME}/bin/R" --slave -e 'library(devtools); devtools::load_all(); devtools::document()'

