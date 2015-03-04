#!/bin/sh
diffdir=`mktemp -d diffdir.XXXXXX`
base1=`basename $1`
base2=`basename $2`

ln -s $1 $diffdir/$base1
ln -s $2 $diffdir/$base2

cd $diffdir
/Library/Frameworks/R.framework/Versions/2.15/Resources/bin/Rscript --no-save --quiet --slave --no-restore ~/Documents/workspace/GenePattern/modules/DiffDatasets/src/run_diff_datasets.R ~/Documents/workspace/GenePattern/modules/DiffDatasets/src/ /Applications/GenePatternServer /Applications/GenePatternServer/patches/ --first.input.file="$base1" --second.input.file="$base2" --round.method=round --round.digits=6
cd ..
rm -rf $diffdir