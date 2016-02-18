#!/bin/bash

output=$1;shift
refg=$2
input=${9}

command -v table_annovar.pl >/dev/null 2>&1 || {
      echo "ERROR : table_annovar.pl not found. Add annovar scripts to your galaxy path !" ;
      return 1 ;
}

mkdir out
name=${input##*/}
name=${name%%.*}

perl $SCRIPT_PATH/mutspecAnnot.pl \
	--outfile out \
	--pathAVDBList $SCRIPT_PATH \
	--temp "./temp" \
	$* 2>&1

ls out/Mutational_Analysis/Annovar/
cp out/Mutational_Analysis/Annovar/${name}.${refg}_multianno.txt $output

exit 0

