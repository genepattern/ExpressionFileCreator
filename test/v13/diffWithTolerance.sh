#!/bin/sh
diffdir=`mktemp -d diffdir.XXXXXX`
mkdir $diffdir/in1
mkdir $diffdir/in2
base1=`basename $1`
base2=`basename $2`

ln $1 $diffdir/in1/$base1
ln $2 $diffdir/in2/$base2

cd $diffdir
Rscript --no-save --quiet --slave --no-restore ~/Documents/workspace/GenePattern/modules/DiffDatasets/src/run_diff_datasets_gpunit.R --first.input.file="in1/$base1" --second.input.file="in2/$base2" --round.method=round --round.digits=6
status=$?
cd ..
rm -rf $diffdir
exit $status