#!/bin/bash

newsign=$1
html=$2
ref=$3

output_dir=${html%%.*}_files

matrix=${newsign%.*}_files/NMF/Files/MatrixW-Normto100.txt

mkdir $output_dir

Rscript --no-save $SCRIPT_PATH/R/compareSignature_Galaxy.r $ref $matrix $output_dir 2>&1

# Convert the image into png format
cd $output_dir

echo "<html><body>" >> $html
echo "<center> <h2> Cosine similarity comparison </h2> </center>" >> $html

echo "<table>" >> $html
echo "<tr> <td> <center> <br/> <a href="Similarity_Matrix.txt">Similarity_Matrix.txt</a> </center> </td> </tr>" >> $html
echo "<tr>" >> $html
echo "<td><a href="Similarity_Matrix.png">" >> $html
echo "<img width="1000" src="Similarity_Matrix.png" /></a></td>" >> $html
echo "</tr>" >> $html
echo "</table>" >> $html

exit 0
