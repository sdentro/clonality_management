#!/bin/bash
outdir=$1

cp *.py ${outdir}/
hash=`git log --pretty=format:'%h' -n 1`
touch ${outdir}/generators_version_${hash}
