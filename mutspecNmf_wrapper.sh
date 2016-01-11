#!/bin/bash

#########################################
###     SPECIFY THE NUMBER OF CPU     ###
#########################################
cpu=8




html=$1;shift
parameters=$1;shift
source=$1;shift
input=$1

if [[ $source == "html" ]] 
then input=${input%%.*}_files/Mutational_Analysis/Figures/Input_NMF/Input_NMF_Count.txt
fi

output_dir=${html%%.*}_files
mkdir $output_dir

Rscript $SCRIPT_PATH/R/somaticSignature_Galaxy.r $parameters --cpu $cpu --input $input --output $output_dir 2>&1


## Test the existence of the files and graphs produced by NMF
if [[ ! -e "$output_dir/NMF/Files/MatrixW-Normto100.txt" ]]; then
	>&2 echo "error"
	exit
fi


echo "<html><body>" >> $html
echo "<center> <h2> NMF Mutational signatures analysis </h2> </center>" >> $html


echo "<table>" >> $html
echo "<tr> <br/> <th><h3>Heatmap of the mixture coefficient matrix</h3></th> </tr>" >> $html
echo "<tr> <td> <center> <br/> <a href="NMF/Files/Cluster_MixtureCoeff.txt">Cluster_MixtureCoeff.txt</a> </center> </td> </tr>" >> $html
echo "<tr>" >> $html

if [[ ! -e "$output_dir/NMF/Figures/Heatmap_MixtureCoeff.png" ]]; then
	echo "WARNING: NMF package can't plot the heatmap when the samples size is above 300. <br/>" >> $html
else
	echo "<td> <center> <a href="NMF/Figures/Heatmap_MixtureCoeff.png">" >> $html
	echo "<img src="NMF/Figures/Heatmap_MixtureCoeff.png" /></a> <center> </td>" >> $html
fi	
echo "</tr>" >> $html
echo "</table>" >> $html

echo "<br/><br/>" >> $html

echo "<table>" >> $html
echo "<tr>" >> $html
echo "<th><h3>Signature composition</h3></th>" >> $html
echo "</tr>" >> $html
echo "<tr><td> <center> <a href="NMF/Files/MatrixW-Normto100.txt">Composition somatic mutation (input matrix for the tool MutSpec-Compare)</a><center></td></tr>" >> $html
echo "<tr>" >> $html
echo "<td><a href="NMF/Figures/CompositionSomaticMutation.png">" >> $html
echo "<img width="1000" src="NMF/Figures/CompositionSomaticMutation.png" /></a></td>" >> $html
echo "</tr>	" >> $html
echo "</table>" >> $html
echo "<br/><br/>" >> $html

echo "<table>" >> $html
echo "<tr>" >> $html
echo "<th><h3>Sample contribution to signatures</h3></th>" >> $html
echo "</tr>" >> $html
echo "<tr><td> <center> <a href="NMF/Files/MatrixH-Inputggplot2.txt">Contribution mutation signature matrix</a></center></td></tr>" >> $html
echo "<tr>" >> $html
echo "<td><a href="NMF/Figures/ContributionMutationSignature.png">" >> $html
echo "<img width="700" src="NMF/Figures/ContributionMutationSignature.png" /></a></td>" >> $html
echo "</tr>	" >> $html
echo "</table>" >> $html
echo "<br/><br/>" >> $html


echo "<table>" >> $html
echo "<tr>" >> $html
echo "<th><h3>Average contributions of each signatures in each cluster</h3></th>" >> $html
echo "</tr>" >> $html
echo "<tr><td> <center> <a href="NMF/Files/Average_ContriByCluster.txt">Average contributions</a></center></td></tr>" >> $html
echo "<tr>" >> $html
echo "<td><a href="NMF/Figures/Average_ContriByCluster.png">" >> $html
echo "<img width="700" src="NMF/Figures/Average_ContriByCluster.png" /></a></td>" >> $html
echo "</tr>	" >> $html
echo "</table>" >> $html
echo "<br/><br/>" >> $html

echo "<br/><br/><br/><br/>" >> $html



exit 0
