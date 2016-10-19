#!/bin/bash

inputdir=$1;shift
annotated_output=$1;shift
summary_output=$1;shift
infofile=$1;shift
pair=$1;

perl $SCRIPT_PATH/hotspot.pl -d $inputdir \
-i $infofile \
-o $annotated_output \
-s $summary_output \
-p $pair

exit 0