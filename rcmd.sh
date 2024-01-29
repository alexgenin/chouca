#!/bin/bash
# 
# This builds the package from git export
# 

PKG="$(basename $(pwd))"
TMPDIR="/tmp/${PKG}build${RANDOM}"
ACTION="$(echo $1 | tr '[:upper:]' '[:lower:]')" 
shift

if [ "$ACTION" == "" ]; then 
  echo "No action provided" 
  exit 1 
fi

git clone "." "$TMPDIR/$PKG"

if [ "$ACTION" == "build" ]; then 
  R CMD build "$TMPDIR/$PKG"
fi

if [ "$ACTION" == "check" ]; then 
  cd "$TMPDIR"
  R CMD build "$PKG" 
  R CMD check --as-cran ${PKG}*.gz
fi

if [ "$ACTION" == "install" ]; then 
  cd "$TMPDIR"
  R CMD build "$PKG" 
  R -e "install.packages(dir(pattern = '^$PKG.*gz$'), repos = NULL)"
fi

rm -rf "$TMPDIR"

