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
	$*

ls out/Mutational_Analysis/Annovar/

if [ -e "out/Mutational_Analysis/Annovar/${name}.${refg}_multianno.txt" ]; then
	cp out/Mutational_Analysis/Annovar/${name}.${refg}_multianno.txt $output
fi

exit 0

