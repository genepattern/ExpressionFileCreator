#!/bin/sh
diffdir=`mktemp -d diffdir.XXXXXX`
mkdir $diffdir/in1
mkdir $diffdir/in2
base1=`basename $1`
base2=`basename $2`

ln -s $1 $diffdir/in1/$base1
ln -s $2 $diffdir/in2/$base2

cd $diffdir
/Library/Frameworks/R.framework/Versions/2.15/Resources/bin/Rscript --no-save --quiet --slave --no-restore ~/Documents/workspace/GenePattern/modules/DiffDatasets/src/run_diff_datasets.R ~/Documents/workspace/GenePattern/modules/DiffDatasets/src/ /Applications/GenePatternServer /Applications/GenePatternServer/patches/ --first.input.file="in1/$base1" --second.input.file="in2/$base2" --round.method=round --round.digits=6
status=$?
cd ..
rm -rf $diffdir
exit $status