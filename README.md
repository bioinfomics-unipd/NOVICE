# NOVICE
## Analysis of Variability in Coverage

#### Notice:
This code was developed as part of a project carried out during the course of Microbial Metagenomics (Molecular Biology master degree)
at the University of Padova. The project was supervised by Prof. Stefano Campanaro, Dr. Maria Silvia Morlino, Dr. Edoardo Bizzotto,
and Dr. Gabriele Ghiotto.

The project authors are Nicola Gambardella (scriptwriter), Yeganeh Alkhasi (data analysis and debugging) , Silvia Frigo (data analysis and debugging) 
and Konstantina Cheshmedzhieva (GitHub page, manual and report writer).


## What is the purpose of NOVICE?

This script is targetted at analysing the coverage profile of genes predicted in a metagenome-assembled genome (MAG) and defining the coverage variance 
compared to a generated centroid. Additionally, it identifies the most distant genes as a step for their further studying. 

### 1. Input Requirements

The input file must be a csv, tsv or ffn file containing the values of the gene coverage estimated as the average coverage per base inside the gene.

### 2. Command Line

The command line usage of this script is very straightforward and allows the user to choose among several options regarding their exact needs. The command
line would look like this:

`python NOVICE_Build_0_6.py input.csv all mean euclidean prefix`

Let us break down the command line and explore its options:

#### 2.1. Rank, PCA or All

As visible from the sample command line provided above, the user has the options to choose their second argument depending on their needs. If you would 
like to process your data and know only the ranking of the genes, you need to simply put `rank` after the input file. If you would like to extract the 
data only in the form of a PCA, the argument would become `PCA`. Use `all` if you would like to receive information on both the rank and a visualization 
thorugh PCA.

#### 2.2. Mean or Median

After you have chosen what data you need to receive in the end, now it is time to give instructions to the script on which value you need it to use. If 
you would like to have your centroid calculated on the basis of the mean values, set this argument to `mean`. If instead, you would like it to use the 
meadian value to calculate it, just instruct it with `median`.

#### 2.3. Euclidean, Bray-Curtis or Correlation

This argument gives you the freedom to choose between the type of distance you would like to receive your result as. Euclidean distance is the ultimate
go-to distance which represent the space between two dots. The Bray-Curtis distance is often used in natural sciences and gives you a result between 0 
and 1. The correlation distance is a more abstract option as the script will calculate the distance as between two random variables with finite variances.

Depending on your targeted results, you can throw `euclidean`, `braycurtis` or `correlation` as the final argument in your command line.

####  2.4. Prefix

As you may be handling a lot of datasets and would like to receive separate output files with different labels, this argument is your way to put in the 
label you need. For example, you are tackling sample labelled as SAMPLE123. To save the ranking and PCA labelled as SAMPLE123 you just need to insert `SAMPLE123` as the argument called `prefix` in the example command line above. This allows you to process multiple datasets one after another and not over-write the created PCAplots but save them separately.


#### Now you are ready to go and explore the results!



