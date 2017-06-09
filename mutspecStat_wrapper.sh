#!/bin/bash

#########################################
###     SPECIFY THE NUMBER OF CPU     ###
#########################################
cpu=8





#########################################
###        Recover the arguments      ###
#########################################
html=$1;shift
len_file_path=$1;shift
estimSign=$1;shift
parameters=$1;shift
working_dir=`pwd`



mkdir in
cd in

names=$(sed 's/\s/_/g' <<< $*)
names=$(sed 's/_\// \//g' <<< $names)
names=$(sed 's/_annotated//g' <<< $names)
names=$(sed 's/_filtered//g' <<< $names)
names=$(sed 's/\.txt_/_/' <<< $names)

for name in ${names}
do
  file=$(sed 's/=/ /' <<< $name);
  echo $file
  ln -s $file
done
cd ..

output_dir=${html%%.*}_files


#########################################
### Calculates the statistics         ###
#########################################

perl $SCRIPT_PATH/mutspecStat.pl --outfile $output_dir \
	--temp "$working_dir/temp" \
	--pathRscript $SCRIPT_PATH \
	$parameters \
	$working_dir/in


#########################################
### Estimate the number of signatures ###
#########################################
if [[ $estimSign > 0 ]]; then
  Rscript $SCRIPT_PATH/R/estimateSign_Galaxy.r --input $output_dir/Mutational_Analysis/Figures/Input_NMF/Input_NMF_Count.txt --stop $estimSign --cpu $cpu --output $output_dir/Mutational_Analysis/Figures/Estimate_Number_Signatures.png 2>&1
fi


#########################################
### Create css                          #
#########################################
css=$output_dir/Mutational_Analysis/style.css
echo ".legend{position:relative}.legend .legend-hidden{display:none;position:absolute;background-color:#fff;border:3px solid #03F;padding:3px;color:#000;font-size:1em;border-radius:10px;margin-top:-40px}.legend:hover .legend-hidden{display:block}" > $css


#########################################
### Create an archive with all results  #
#########################################
cd $output_dir
zip -r "$output_dir/Mutational_Analysis.zip" "Mutational_Analysis"



# HMTL page for the result of the tool
echo "<html>" >> $html
echo "<body>" >> $html

if [ -d $output_dir/Mutational_Analysis/Figures ]; then

echo "<center> <h2>Mutational spectra report summary</h2> </center>" >> $html

echo "<br/> Download the results" >> $html
echo "<br/><a href="Mutational_Analysis.zip">Mutational_Analysis.zip</a><br/>" >> $html

echo "<br/> Download the full report in Excel" >> $html

## One report with all the samples. Specify the full path
if [[ -e "$output_dir/Mutational_Analysis/Report_Mutation_Spectra.xls" ]]
then
	# Interpreted by Galaxy so don't need the full path
	echo "<br/><a href="Mutational_Analysis/Report_Mutation_Spectra.xls">Report_Mutation_Spectra.xls</a>" >> $html
fi
## One report for each samples
for file in $names
do
  name=$(echo ${file}| cut -d"=" -f2)
  name=${name%.*}

  # One report for each samples
  if [[ -e "$output_dir/Mutational_Analysis/Report_Mutation_Spectra-$name.xls" ]]
  then
    echo "<br/><a href="Mutational_Analysis/Report_Mutation_Spectra-$name.xls">Report_Mutation_Spectra-$name.xls</a>" >> $html
  fi
done
## One report for each samples: Pool_Data
if [[ $parameters =~ "--pooldata" ]]; then
  if [[ -e "$output_dir/Mutational_Analysis/Report_Mutation_Spectra-Pool_Data.xls" ]]; then
    echo "<br/><a href="Mutational_Analysis/Report_Mutation_Spectra-Pool_Data.xls">Report_Mutation_Spectra-Pool_Data.xls</a>" >> $html
  fi
fi


## Input file for NMF
if [[ -e "$output_dir/Mutational_Analysis/Figures/Input_NMF/Input_NMF_Count.txt" ]]
then
  # Interpreted by Galaxy so don't need the full path
  echo "<br/><br/> Download the input file for the tool mutSpec-NMF" >> $html
  echo "<br/><a href="Mutational_Analysis/Figures/Input_NMF/Input_NMF_Count.txt">Input_NMF_Count.txt</a><br/>" >> $html
fi

## Computed statistics for estimating the number of signatures
if [[ $estimSign > 0 ]]; then
  echo "<br/> Link to the computed statistics for estimating the number of signatures <br/>" >> $html
  if [[ -e "$output_dir/Mutational_Analysis/Figures/Estimate_Number_Signatures.png" ]]; then
    outEstimateSign="$output_dir/Mutational_Analysis/EstimatingSignatures.html"
    touch $outEstimateSign
    echo "<a href='Mutational_Analysis/EstimatingSignatures.html'>Estimating the number of signatures</a><br/>" >> $html
    echo "<br/> <center> <h2>Computed statistics for estimating the number of signatures</h2> </center> <br/>" >> $outEstimateSign
    echo "Several approaches have been proposed to choose the optimal number of signatures to extract with NMF. <br/>
          Brunet et al. 2004, proposed to take the first number of signature for which the cophenetic coefficient starts decreasing, <br/>
          Hutchins et al. 2008, suggested to choose the first value where the RSS curve presents an inflection point. <br/>
          Frigyesi et al. 2008, considered the smallest value at which the decrease in the RSS is lower than the decrease of the RSS obtained from random data. <br/><br/>
          The estimation are based on Brunet’s algorithm computed from 50 runs for each value of signature to estimate. <br/> <br/>
          The original data are shuffled for comparing the quality measures obtained with our data (Data x) and from randomized data (Data y). The curves for the actual data are in solid line, those for the randomized data are in dashed line. <br/> <br/>" >> $outEstimateSign
    echo "<img src="Figures/Estimate_Number_Signatures.png width="1000""/><br/></td>" >> $outEstimateSign
  else
    echo "<br/>There is not enough mutations for estimating the number of signatures <br/>" >> $html
    echo "Read the tool standard output for more detail<br/>"  >> $html
  fi
fi


## HMTL Link to the samples
echo "<br/> Link to individual samples <br/>" >> $html


## Consider only samples with at least one mutation
for name in `ls $output_dir/Mutational_Analysis/Figures/Impact_protein_sequence`
do

  ## Pool Data is handle separately
  if [ $name = "Pool_Data" ]; then
  	continue
  fi

  outfile="$output_dir/Mutational_Analysis/$name.html"
  touch $outfile # Create an empty file named $outfile
  echo "<a href='Mutational_Analysis/$name.html'>$name</a><br/>" >> $html

#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#                                                     			INDIVIDUAL SAMPLES                                                                        #
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
echo "<br/> <center> <h2>Mutational Spectra report for $name</h2> </center> <br/>" >> $outfile

echo "<html>" >> $outfile

echo "<head>" >> $outfile
echo "<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">" >> $outfile
# Link to the css style file for having a legend when we pass the mouse on the figures
echo "<link rel="stylesheet" href="style.css" />" >> $outfile
echo "</head>" >> $outfile

echo "<body>" >> $outfile

echo "<table>" >> $outfile
echo "<tr>" >> $outfile
echo "<th><h3>Overall mutation distribution</h3></th>" >> $outfile
echo "<th><h3>Impact on protein sequence</h3></th>" >> $outfile
echo "</tr><tr>" >> $outfile
echo "<td> <center> <a href="Figures/Overall_mutation_distribution/$name/$name-OverallMutationDistribution.txt">$name-OverallMutationDistribution.txt</a> </center> </td>" >> $outfile
echo "<td> <center> <a href="Figures/Impact_protein_sequence/$name/$name-DistributionExoFunc.txt">$name-DistributionExoFunc.txt</a> </center> </td>" >> $outfile
echo "</tr><tr>" >> $outfile

echo "<td>" >> $outfile
echo "<span class="legend"><img src="Figures/Overall_mutation_distribution/$name/$name-OverallMutationDistribution.png width="280""/>" >> $outfile
echo "<span class="legend-hidden">" >> $outfile
echo "<center><B>Overall Mutation Distribution</center></B><br/>Proportion of all mutation types (total count are indicated in parenthesis). For indels the counts are based on annotation retrieved from the database ExonicFunc.refGene<br/>" >> $outfile
echo "</td>" >> $outfile
echo "<td>" >> $outfile
echo "<span class="legend"><img src="Figures/Impact_protein_sequence/$name/$name-DistributionExoFunc.png width="400""/>" >> $outfile
echo "<span class="legend-hidden">" >> $outfile
echo "<center><B>Graph 1. Impact on protein sequence</center></B><br/>Impact of all mutations (SBS and Indel) on the protein sequence based on the ExonicFunc.refGene annotation<br/>" >> $outfile
echo "</td>" >> $outfile

echo "</tr>" >> $outfile
echo "</table>" >> $outfile

echo "<br/><br/>"  >> $outfile


echo "<table>" >> $outfile
echo "<tr>" >> $outfile
echo "<th><h3>SBS distribution</h3></th>" >> $outfile
echo "<th><h3>Stranded distribution of SBS</h3></th>" >> $outfile
echo "</tr><tr>" >> $outfile
echo "<td> <center> <a href="Figures/SBS_distribution/$name/$name-SBS_distribution.txt">$name-SBS_distribution.txt</a> </center> </td>" >> $outfile
echo "<td> <center> <a href="Figures/Stranded_Analysis/$name/$name-StrandBias.txt">$name-StrandBias.txt</a> </center> </td>" >> $outfile
echo "</tr><tr>" >> $outfile

echo "<td>" >> $outfile
echo "<span class="legend"><img src="Figures/SBS_distribution/$name/$name-SBS_distribution.png width="550""/>" >> $outfile
echo "<span class="legend-hidden">" >> $outfile
echo "<center><B>Graph 2. SBS distribution</center></B><br/>Proportion of each type of single base substitution (SBS)<br/>" >> $outfile
echo "</td>" >> $outfile
echo "<td>" >> $outfile
echo "<span class="legend"><img src="Figures/Stranded_Analysis/$name/$name-StrandBias.png width="400""/>" >> $outfile
echo "<span class="legend-hidden">" >> $outfile
echo "<center><B>Graph 3. Stranded distribution of SBS</center></B><br/>Count of the six substitution types on the transcribed and non-transcribed strand<br/>" >> $outfile
echo "</td>" >> $outfile

echo "</tr>" >> $outfile
echo "</table>" >> $outfile


echo "<br/><br/>"  >> $outfile


######################################################
#	Trinucleotide sequence context of SBS on genomic #
######################################################
echo "<table>" >> $outfile
echo "<h3>Trinucleotide sequence context of SBS on the genomic sequence</h3>" >> $outfile
echo "<tr>" >> $outfile
echo "<td> <center> <a href="Figures/Trinucleotide_Sequence_Context/$name/$name-MutationSpectraPercent-Genomic.txt">$name-MutationSpectraPercent.txt</a> </center> </td>" >> $outfile
echo "<td> <center> <a href="Figures/Trinucleotide_Sequence_Context/$name/$name-HeatmapPercent-Genomic.txt">$name-HeatmapPercent-Genomic.txt</a> </center> </td>" >> $outfile
echo "</tr><tr>" >> $outfile

echo "<td>" >> $outfile
echo "<span class="legend"><img src="Figures/Trinucleotide_Sequence_Context/$name/$name-MutationSpectraPercent-Genomic.png width="1000""/>" >> $outfile
echo "<span class="legend-hidden">" >> $outfile
echo "<center><B>Panel 1. Trinucleotide sequence context of SBS on the genomic sequence</center></B><br/>Proportion of the six substitution types with their trinucleotide sequence context (total number of mutation is shown in parenthesis)<br/>" >> $outfile
echo "</td>" >> $outfile
echo "<td>" >> $outfile
echo "<span class="legend"><img src="Figures/Trinucleotide_Sequence_Context/$name/$name-HeatmapPercent-Genomic.png width="250""/>" >> $outfile
echo "<span class="legend-hidden">" >> $outfile
echo "<center><B>Panel 1. Trinucleotide sequence context of SBS on the genomic sequence</center></B><br/>Proportion of the six substitution types with their trinucleotide sequence context<br/>" >> $outfile
echo "</td>" >> $outfile

echo "</tr>" >> $outfile
echo "</table>" >> $outfile


echo "<br/><br/>"  >> $outfile


##############################################################
#	Trinucleotide sequence context of SBS on coding sequence #
##############################################################
echo "<table>" >> $outfile
echo "<h3>Stranded analysis of trinucleotide sequence context of SBS</h3>" >> $outfile
echo "<tr>" >> $outfile
echo "<td> <center> <a href="Figures/Stranded_Analysis/$name/$name-StrandedSignaturePercent.txt">$name-StrandedSignaturePercent.txt</a> </center> </td>" >> $outfile
echo "</tr><tr>" >> $outfile

echo "<td>" >> $outfile
echo "<span class="legend"><img src="Figures/Stranded_Analysis/$name/$name-StrandedSignaturePercent.png width="1300""/>" >> $outfile
echo "<span class="legend-hidden">" >> $outfile
echo "<center><B>Panel 2. Stranded analysis of trinucleotide sequence context of SBS</center></B><br/>Proportion of SBS with their trinucleotide context considering the non-transcribed and transcribed strand<br/>" >> $outfile
echo "</td>" >> $outfile
echo "</tr>" >> $outfile
echo "</table>" >> $outfile

echo "<br/><br/>"  >> $outfile


#############################################
#	Sequence logo generated with Weblogo3 	#
#############################################
echo "<table>" >> $outfile
echo "<h3>Wider sequence context with Weblogo3</h3>" >> $outfile
# Legende de la figure : Panel 3. Wider sequence context on genomic strand generated with Weblogo3

# C>A
echo "<tr>" >> $outfile
if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/$name/$name-CA-Probability.png" ]]; then
  echo "<td>WARNING: No sequence for C>A </br> </td>" >> $outfile
else
  echo "<td><a href="Figures/WebLogo/$name/$name-CA.fa">$name-CA.fa</a><br/>" >> $outfile
  echo "<img src="Figures/WebLogo/$name/$name-CA-Probability.png width="600" "/><br/></td>" >> $outfile
fi
# C>G
if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/$name/$name-CG-Probability.png" ]]; then
  echo "<td> WARNING: No sequence for C>G </br> </td>" >> $outfile
else
  echo "<td><a href="Figures/WebLogo/$name/$name-CG.fa">$name-CG.fa</a><br/>" >> $outfile
  echo "<img src="Figures/WebLogo/$name/$name-CG-Probability.png width="600" "/><br/></td>" >> $outfile
fi
# C>T
if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/$name/$name-CT-Probability.png" ]]; then
  echo "<td> WARNING: No sequence for C>T </br> </td>" >> $outfile
else
 echo "<td><a href="Figures/WebLogo/$name/$name-CT.fa">$name-CT.fa</a><br/>" >> $outfile
 echo "<img src="Figures/WebLogo/$name/$name-CT-Probability.png width="600" "/><br/></td>" >> $outfile
fi
echo "</tr>" >> $outfile

# T>A
echo "<tr>" >> $outfile
if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/$name/$name-TA-Probability.png" ]]; then
  echo "<td>WARNING: No sequence for T>A </br> </td>" >> $outfile
else
  echo "<td><a href="Figures/WebLogo/$name/$name-TA.fa">$name-TA.fa</a><br/>" >> $outfile
  echo "<img src="Figures/WebLogo/$name/$name-TA-Probability.png width="600" "/><br/></td>" >> $outfile
fi
# T>C
if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/$name/$name-TC-Probability.png" ]]; then
  echo "<td>WARNING: No sequence for T>C </br> </td>" >> $outfile
else
  echo "<td><a href="Figures/WebLogo/$name/$name-TC.fa">$name-TC.fa</a><br/>" >> $outfile
  echo "<img src="Figures/WebLogo/$name/$name-TC-Probability.png width="600" "/><br/></td>" >> $outfile
fi
# T>G
if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/$name/$name-TG-Probability.png" ]]; then
  echo "<td>WARNING: No sequence for T>G </br> </td>" >> $outfile
else
  echo "<td><a href="Figures/WebLogo/$name/$name-TG.fa">$name-TG.fa</a><br/>" >> $outfile
  echo "<img src="Figures/WebLogo/$name/$name-TG-Probability.png width="600" "/><br/></td>" >> $outfile
fi
echo "</tr>" >> $outfile

echo "</table>" >> $outfile

echo "</body></html>" >> $outfile

done

#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#                                                                       POOL DATA                                                                     #
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
## HMTL Link to Pool_Data
if [[ $parameters =~ "--pooldata" ]]; then
  outfilePoolData="$output_dir/Mutational_Analysis/Pool_Data.html"
  touch $outfilePoolData # Create an empty file named $outfile
  echo "<a href='Mutational_Analysis/Pool_Data.html'>Pool_Data</a><br/>" >> $html


  echo "<br/> <center> <h2>Mutational Spectra report for Pool_Data</h2> </center> <br/>" >> $outfilePoolData
  echo "<html>" >> $outfilePoolData

  echo "<head>" >> $outfilePoolData
  echo "<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">" >> $outfilePoolData
  # Link to the css style file for having a legend when we pass the mouse on the figures
  echo "<link rel="stylesheet" href="style.css" />" >> $outfilePoolData
  echo "</head>" >> $outfilePoolData

  echo "<body>" >> $outfilePoolData

  echo "<table>" >> $outfilePoolData
  echo "<tr>" >> $outfilePoolData
  echo "<th><h3>Overall mutation distribution</h3></th>" >> $outfilePoolData
  echo "<th><h3>Impact on protein sequence</h3></th>" >> $outfilePoolData
  echo "</tr><tr>" >> $outfilePoolData
  echo "<td> <center> <a href="Figures/Overall_mutation_distribution/Pool_Data/Pool_Data-OverallMutationDistribution.txt">Pool_Data-OverallMutationDistribution.txt</a> </center> </td>" >> $outfilePoolData
  echo "<td> <center> <a href="Figures/Impact_protein_sequence/Pool_Data/Pool_Data-DistributionExoFunc.txt">Pool_Data-DistributionExoFunc.txt</a> </center> </td>" >> $outfilePoolData
  echo "</tr><tr>" >> $outfilePoolData

  echo "<td>" >> $outfilePoolData
  echo "<span class="legend"><img src="Figures/Overall_mutation_distribution/Pool_Data/Pool_Data-OverallMutationDistribution.png width="280""/>" >> $outfilePoolData
  echo "<span class="legend-hidden">" >> $outfilePoolData
  echo "<center><B>Overall Mutation Distribution</center></B><br/>Proportion of all mutation types (total count are indicated in parenthesis). For indels the counts are based on annotation retrieved from the database ExonicFunc.refGene<br/>" >> $outfilePoolData
  echo "</td>" >> $outfilePoolData
  echo "<td>" >> $outfilePoolData
  echo "<span class="legend"><img src="Figures/Impact_protein_sequence/Pool_Data/Pool_Data-DistributionExoFunc.png width="400""/>" >> $outfilePoolData
  echo "<span class="legend-hidden">" >> $outfilePoolData
  echo "<center><B>Graph 1. Impact on protein sequence</center></B><br/>Impact of all mutations (SBS and Indel) on the protein sequence based on the ExonicFunc.refGene annotation<br/>" >> $outfilePoolData
  echo "</td>" >> $outfilePoolData

  echo "</tr>" >> $outfilePoolData
  echo "</table>" >> $outfilePoolData

  echo "<br/><br/>"  >> $outfilePoolData


  echo "<table>" >> $outfilePoolData
  echo "<tr>" >> $outfilePoolData
  echo "<th><h3>SBS distribution</h3></th>" >> $outfilePoolData
  echo "<th><h3>Stranded distribution of SBS</h3></th>" >> $outfilePoolData
  echo "</tr><tr>" >> $outfilePoolData
  echo "<td> <center> <a href="Figures/SBS_distribution/Pool_Data/Pool_Data-SBS_distribution.txt">Pool_Data-SBS_distribution.txt</a> </center> </td>" >> $outfilePoolData
  echo "<td> <center> <a href="Figures/Stranded_Analysis/Pool_Data/Pool_Data-StrandBias.txt">Pool_Data-StrandBias.txt</a> </center> </td>" >> $outfilePoolData
  echo "</tr><tr>" >> $outfilePoolData

  echo "<td>" >> $outfilePoolData
  echo "<span class="legend"><img src="Figures/SBS_distribution/Pool_Data/Pool_Data-SBS_distribution.png width="550""/>" >> $outfilePoolData
  echo "<span class="legend-hidden">" >> $outfilePoolData
  echo "<center><B>Graph 2. SBS distribution</center></B><br/>Proportion of each type of single base substitution (SBS)<br/>" >> $outfilePoolData
  echo "</td>" >> $outfilePoolData
  echo "<td>" >> $outfilePoolData
  echo "<span class="legend"><img src="Figures/Stranded_Analysis/Pool_Data/Pool_Data-StrandBias.png width="400""/>" >> $outfilePoolData
  echo "<span class="legend-hidden">" >> $outfilePoolData
  echo "<center><B>Graph 3. Stranded distribution of SBS</center></B><br/>Count of the six substitution types on the transcribed and non-transcribed strand<br/>" >> $outfilePoolData
  echo "</td>" >> $outfilePoolData

  echo "</tr>" >> $outfilePoolData
  echo "</table>" >> $outfilePoolData


  echo "<br/><br/>"  >> $outfilePoolData


  ##########################################################
  # Trinucleotide sequence context of SBS on genomic: Pool #
  ##########################################################
  echo "<table>" >> $outfilePoolData
  echo "<h3>Trinucleotide sequence context of SBS on the genomic sequence</h3>" >> $outfilePoolData
  echo "<tr>" >> $outfilePoolData
  echo "<td> <center> <a href="Figures/Trinucleotide_Sequence_Context/Pool_Data/Pool_Data-MutationSpectraPercent-Genomic.txt">Pool_Data-MutationSpectraPercent.txt</a> </center> </td>" >> $outfilePoolData
  echo "<td> <center> <a href="Figures/Trinucleotide_Sequence_Context/Pool_Data/Pool_Data-HeatmapPercent-Genomic.txt">Pool_Data-HeatmapPercent-Genomic.txt</a> </center> </td>" >> $outfilePoolData
  echo "</tr><tr>" >> $outfilePoolData

  echo "<td>" >> $outfilePoolData
  echo "<span class="legend"><img src="Figures/Trinucleotide_Sequence_Context/Pool_Data/Pool_Data-MutationSpectraPercent-Genomic.png width="1000""/>" >> $outfilePoolData
  echo "<span class="legend-hidden">" >> $outfilePoolData
  echo "<center><B>Panel 1. Trinucleotide sequence context of SBS on the genomic sequence</center></B><br/>Proportion of the six substitution types with their trinucleotide sequence context (total number of mutation is shown in parenthesis)<br/>" >> $outfilePoolData
  echo "</td>" >> $outfilePoolData
  echo "<td>" >> $outfilePoolData
  echo "<span class="legend"><img src="Figures/Trinucleotide_Sequence_Context/Pool_Data/Pool_Data-HeatmapPercent-Genomic.png width="250""/>" >> $outfilePoolData
  echo "<span class="legend-hidden">" >> $outfilePoolData
  echo "<center><B>Panel 1. Trinucleotide sequence context of SBS on the genomic sequence</center></B><br/>Proportion of the six substitution types with their trinucleotide sequence context<br/>" >> $outfilePoolData
  echo "</td>" >> $outfilePoolData

  echo "</tr>" >> $outfilePoolData
  echo "</table>" >> $outfilePoolData


  echo "<br/><br/>"  >> $outfilePoolData


  ##################################################################
  # Trinucleotide sequence context of SBS on coding sequence: Pool #
  ##################################################################
  echo "<table>" >> $outfilePoolData
  echo "<h3>Stranded analysis of trinucleotide sequence context of SBS</h3>" >> $outfilePoolData
  echo "<tr>" >> $outfilePoolData
  echo "<td> <center> <a href="Figures/Stranded_Analysis/Pool_Data/Pool_Data-StrandedSignaturePercent.txt">Pool_Data-StrandedSignaturePercent.txt</a> </center> </td>" >> $outfilePoolData
  echo "</tr><tr>" >> $outfilePoolData

  echo "<td>" >> $outfilePoolData
  echo "<span class="legend"><img src="Figures/Stranded_Analysis/Pool_Data/Pool_Data-StrandedSignaturePercent.png width="1300""/>" >> $outfilePoolData
  echo "<span class="legend-hidden">" >> $outfilePoolData
  echo "<center><B>Panel 2. Stranded analysis of trinucleotide sequence context of SBS</center></B><br/>Proportion of SBS with their trinucleotide context considering the non-transcribed and transcribed strand<br/>" >> $outfilePoolData
  echo "</td>" >> $outfilePoolData
  echo "</tr>" >> $outfilePoolData
  echo "</table>" >> $outfilePoolData

  echo "<br/><br/>"  >> $outfilePoolData

  #####################################################
  # Sequence logo generated with Weblogo3: Pool   #
  #####################################################
  echo "<table>" >> $outfilePoolData
  echo "<h3>Sequence logo generated with Weblogo3</h3>" >> $outfilePoolData
  # C>A
  echo "<tr>" >> $outfilePoolData
  if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/Pool_Data/Pool_Data-CA-Probability.png" ]]; then
    echo "<td>WARNING: No sequence for C>A </br> </td>" >> $outfilePoolData
  else
    echo "<td><a href="Figures/WebLogo/Pool_Data/Pool_Data-CA.fa">Pool_Data-CA.fa</a><br/>" >> $outfilePoolData
    echo "<img src="Figures/WebLogo/Pool_Data/Pool_Data-CA-Probability.png"/><br/></td>" >> $outfilePoolData
  fi
  # C>G
  if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/Pool_Data/Pool_Data-CG-Probability.png" ]]; then
    echo "<td>WARNING: No sequence for C>G </br> </td>" >> $outfilePoolData
  else
    echo "<td><a href="Figures/WebLogo/Pool_Data/Pool_Data-CG.fa">Pool_Data-CG.fa</a><br/>" >> $outfilePoolData
    echo "<img src="Figures/WebLogo/Pool_Data/Pool_Data-CG-Probability.png"/><br/></td>" >> $outfilePoolData
  fi
  # C>T
  if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/Pool_Data/Pool_Data-CT-Probability.png" ]]; then
    echo "<td>WARNING: No sequence for C>T </br> </td>" >> $outfilePoolData
  else
    echo "<td><a href="Figures/WebLogo/Pool_Data/Pool_Data-CT.fa">Pool_Data-CT.fa</a><br/>" >> $outfilePoolData
    echo "<img src="Figures/WebLogo/Pool_Data/Pool_Data-CT-Probability.png"/><br/></td>" >> $outfilePoolData
  fi
  echo "</tr>" >> $outfilePoolData

  # T>A
  echo "<tr>" >> $outfilePoolData
  if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/Pool_Data/Pool_Data-TA-Probability.png" ]]; then
    echo "<td>WARNING: No sequence for T>A </br> </td>" >> $outfilePoolData
  else
    echo "<td><a href="Figures/WebLogo/Pool_Data/Pool_Data-TA.fa">Pool_Data-TA.fa</a><br/>" >> $outfilePoolData
    echo "<img src="Figures/WebLogo/Pool_Data/Pool_Data-TA-Probability.png"/><br/></td>" >> $outfilePoolData
  fi
  # T>C
  if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/Pool_Data/Pool_Data-TC-Probability.png" ]]; then
    echo "<td>WARNING: No sequence for T>C </br> </td>" >> $outfilePoolData
  else
    echo "<td><a href="Figures/WebLogo/Pool_Data/Pool_Data-TC.fa">Pool_Data-TC.fa</a><br/>" >> $outfilePoolData
    echo "<img src="Figures/WebLogo/Pool_Data/Pool_Data-TC-Probability.png"/><br/></td>" >> $outfilePoolData
  fi
  # T>G
  if [[ ! -e "$output_dir/Mutational_Analysis/Figures/WebLogo/Pool_Data/Pool_Data-TG-Probability.png" ]]; then
    echo "<td>WARNING: No sequence for T>G </br> </td>" >> $outfilePoolData
  else
    echo "<td><a href="Figures/WebLogo/Pool_Data/Pool_Data-TG.fa">Pool_Data-TG.fa</a><br/>" >> $outfilePoolData
    echo "<img src="Figures/WebLogo/Pool_Data/Pool_Data-TG-Probability.png"/><br/></td>" >> $outfilePoolData
  fi
  echo "</tr>" >> $outfilePoolData
  echo "</table>" >> $outfilePoolData

  echo "</body></html>" >> $outfilePoolData

fi # End if --poolData

fi # End if [ -d $output_dir/Mutational_Analysis/Figures ]

echo "</body></html>" >> $html

exit 0

