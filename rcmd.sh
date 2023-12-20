#!/bin/bash
# 
# This builds the package from git export
# 

TMPDIR="/tmp/choucabuild$RANDOM"
ACTION="$(echo $1 | tr '[:upper:]' '[:lower:]')" 
shift

if [ "$ACTION" == "" ]; then 
  echo "No action provided" 
  exit 1 
fi
  
git clone "." "$TMPDIR/chouca"

if [ "$ACTION" == "build" ]; then 
  R CMD build "$TMPDIR/chouca"
fi

if [ "$ACTION" == "check" ]; then 
  cd "$TMPDIR"
  R CMD build "chouca" 
  R CMD check --as-cran chouca*.tar.gz
fi

rm -rf "$TMPDIR"

