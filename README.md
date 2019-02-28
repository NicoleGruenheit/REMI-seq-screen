### Analysis tool for screens conducted using the REMI-seq method

This tool can be used for the analysis of fastq files, which were generated by the REMI-seq method (see https://remi-seq.org/ for details on the method). In short, the tool extracts those reads, that contain the end of the integration vector (), the GATC (DpnII) or CATG (NlaIII) site, and a short sequence of genomic DNA. The genomic DNA is then compared to a lookup list of possible insertion points in the *D. discoideum* genome, i.e. sequence tags next to DpnII or NlaIII sites. The number of tags present in a dataset is correlated to the number of respective mutants in the analysed sample. In the simplest case, this tool can just be used to analyse a pool of mutants. However, as is the case with single cell sequencing, there can be drop-outs due to technical reasons. Therefore, it is recommended to sequence at least two better three to four technical replicates of the starting library to determine the technical drop-out rate.

To conduct a full analysis of a starting library and then compare the counts to pools of mutants after a selection process, five steps are necessary:

#### Step 1: Load the experimental design

If you are only analysing one pool, this step can be skipped. Elsewise, the tool needs to know how many different samples there are including the starting library. So if you are just analysing the starting library and the pool at the end of the screen, this number should be set to 2. Next, the full experimental design with filenames for the different samples and replicates is needed. There is an example file (exp_design_small.txt) in the help_files folder:

| Sample | Replicates |     Filename     |  Name  |
| :----: | :--------: | :--------------: | :----: |
| start  |     2      | start1.fastq.zip | start1 |
| start  |     2      | start2.fastq.zip | start2 |
|  test  |     2      | test1.fastq.zip  | test1  |
|  test  |     2      | test2.fastq.zip  | test2  |

Files can be either gzipped or not but it usually makes more sense to not unzip them. Example read files can be found in the folder raw_reads.

Once the experimental design is loaded, you should be able to select the reference sample from the sample names you have provided in the file. This should be the name of the starting library (start in this case). The table at the bottom of the page will display the experimental design so you can check whether it has been loaded correctly. You can now proceed to step 2.

#### Step 2: Extract tags

In this step, you can either extract tags and tag counts from raw fastq files or load precomputed files. If you are analysing files for the first time, you definitely have to you the first option, which can take some time if the fastq files are big. The files needed for the latter option are written at the end of the full analysis and are loaded instantaneously. this step also conducts the lookup step so you have to provide the lookup table (help_files/annotated_positions_03022016) as well. If you are only interested in specific positions, you can filter the lookup table to only contain the positions you are interested in.

The REMI-seq method also can use several different vectors (pGWDI-G1, -G2, -G3, -C4, -G5, -C6, -C7 and -C8 ), which are all based on pLPBLP (Faix et al., 2004) and are obtainable from Dictybase. These vectors all differ by two 6 bp indices (one at each end of the vector), that is contained in the extracted tag. Therefore, this index can be used to distinguish mutants generated in different transformations and also from which side of the vector the tag was extracted. These mutants can be seen as independent and can therefore be used to validate the behaviour of a specific mutant. The indices for the different vectors have to be loaded as well (help_files/indices). Again, if you are only interested in some vectors, you can provided a shortened list.

Finally, you have to provide the path to the folder that contains the all fastq files contained in the experimental design. 

To start the analysis of the fastq files, click the button 'Start Analysis'. If you have precomputed files that you just want to re-analyse, click the button 'Load precomputed file'. The first option will take some time depending on the number and size of the fastq files (give estimate here; maybe add a progress bar). The analysis is finished as soon as a table appears below the options that contains basic stats to each file (e.g. the number of reads and the number of reads containing tags). 

The second option should immediately result in the display of the stats table. If it doesn't, check whether the precomputed files are in the same directory as the tool.

#### Step 3: Generate count tables

Counts for different positions can be summarised in different ways. Step 2 only provides counts for each tag separately, i.e. each side of the vector is also counted separately. To conduct the full analysis, the user can now choose to apply a read count cutoff (tags supported by less reads will be ignored) and to combine the counts in different ways. A cutoff of at least 5 reads resulted in good results in our experience but in some cases it might be better to only look at well supported insertion points so this cutoff can be set to e.g. 100.

Counts for each option are normalised to the sequencing depth apart for the last option. Tables can be downloaded by clicking on the download button below the displayed table. The different ways to combine the counts are:

1. Counts for each insertion point

   ​	This options combines different sides of the vector and different vectors.

2. Counts for different vectors at each insertion point

   ​	This options combines different sides of the vector but provides counts for mutants that were generated using 	different vectors.

3. Counts for each gene

   ​	In most cases, a screen is conducted to identify genes that convey an advantage or disadvantage in a specific 	condition. Therefore, it is prudent to combine insertion sites that are located in the same gene. Intergenic 		insertion sites are ignored in this case.

4. Counts for each gene including 500 bp upstream

   ​	Especially in *D. discoideum* the promoter region of a gene is located within the 500 bp upstream of of the 	gene. Therefore, in most cases an insertion within 500 bp of the start of the gene is also a knock-out. This 	option adds counts of all insertion points within these 500 bp to the counts of a gene.

5. Raw read counts for each insertion point

   ​	This option provides raw counts for each insertion point.

#### Step 4: Check correlation of replicates

One of the main errors in selection experiments is if two biological replicates are too different. This will lead to a very low detection level of over- and underrepresented mutants. Therefore, step 4 is a quality check of all samples. Especially the samples of the starting library (technical replicates) should be highly correlated. If this is not the case, it is usually an indicator for some errors in the experiment. These errors entail e.g. too much or too little selection pressure on one of the replicates. It is prudent to exclude replicates that show any aberrant behaviour. For this analysis, the experimental design file has to be loaded in step 1.

There are also options to exclude mutants with very low read counts, which can decrease the correlation. Also, mutants that show very different read counts in the replicates can be excluded. The plots can be downloaded byt clicking on the download button below the plots.

If the displayed scatter plot show sufficient correlation (usually > 0.6), you can proceed to step 5. 

#### Step 5: Determine over- and underrepresented mutants

