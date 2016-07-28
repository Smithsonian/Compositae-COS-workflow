#Targeted next-gen sequencing for plant phylogenomics: the bioinformatics
##Botany 2016 Workshop

###Authors: 
**Jennifer R. Mandel** (jmandel@memphis.edu)  
**Rebecca B. Dikow** (DikowR@si.edu)  
**Vicki A. Funk** (FunkV@si.edu)  

###Objectives:  
In the course of this workshop, we detail the software and commands you will need to generate and analyze target enrichment data for phylogeny. At each point in the workflow, there are multiple options for software and/or commands, so while we try to present these options, we can't guarantee that we have been exhaustive. The tools are quickly evolving. Other resources attached to this workshop include UNIX cheat sheets, powerpoint about wet-lab procedures and expedition protocols.

###References:
* Mandel *et al.*, 2014
* Mandel *et al.*, 2015

###Outline: 
1. Accessing code via GitHub 
2. Probe design
3. Data wrangling and quality trimming  
4. Assembly  
5. PHYLUCE/orthology detection 
6. Alignment and phylogenetic analysis  
7. Gene tree/species tree methods
8. Visualization   

###1. Accessing code via GitHub 
* Github is the place to go for code, scripts, and software. 
* This project has a GitHub site where there are instructions and the wrapper scripts used for our workflow: [Compositae COS workflow](https://github.com/Smithsonian/Compositae-COS-workflow)
* We recommend installing [GitHub Desktop](https://desktop.github.com) for its ease of use.

###2. Probe design
* The first step when you undertake or even think about undertaking a target enrichment study is to identify all the genomic resources that exist for your group: genomes, transcriptomes, ESTs and to think about the phylogenetic breadth you would like to condiser.

* **Orthology Detection:**  
	+ [OMA](http://omabrowser.org/standalone/): with OMA, you can find *de novo* orthologs from transcriptome data with or without genomes.
	+ [AGALMA](https://bitbucket.org/caseywdunn/agalma): AGALMA is a pipeline that generates orthologs and species trees from transcriptomes. It incorporates [PhyloTreePruner](https://sourceforge.net/projects/phylotreepruner/), which distinguished between in- and out- paralogs.
	+ [Orthofinder](http://www.stevekellylab.com/software/orthofinder): Orthofinder is a program for identifying orthologous protein sequence families. It is written in python and runs as a single command that takes as input a directory of FASTA files, one per species. It outputs a file containing the orthologous groups of genes from these species.
	+ [OrthoMCL](http://orthomcl.org): database that includes ortholog groups for a curated set of taxa. 
	+ [HaMStR](https://sourceforge.net/projects/hamstr/): extends a set of orthologs (e.g. from OrthoMCL) to additional taxa using hidden markov models.

* **Post Orthology Detection:**
	+ Check for approximate copy number and intron boundaries by mapping orthologs to a reference genome. Tools: BWA, TopHat, GSNAP.
	+ Choose single or low copy-number orthologs.
	+ Align orthologs across the taxonomic divergence you want to capture.
	+ [MYcroarray](http://www.mycroarray.com) suggests around 10% divergence to guarantee capture.

###3. Data wrangling and quality trimming 
* Here are some tools that can be used for quality assessement and trimming. You will 
* [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/): FASTQC is a program that allows you to visualize the quality of your raw and/or trimmed data. Very easy to run, even on a laptop or desktop.
* [TrimGalore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/): I like TrimGalore! because it auto-detects adapters.
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic): Commonly used but uses Java, which can be problematic on HPC clusters.
* [Prinseq](http://prinseq.sourceforge.net) Used in Mandel (2014; 2015). 
	+ **Commands:** ```perl prinseq-lite.pl```
	+ **Arguments:**  
	+ ```-fastq <READFILE.fq>```  
	+  ```-out_format 1```
	+  ```-out_good <READFILE_trimmed.fq>``` 
	+  ```–log <READFILE.log>``` 
	+  ```-min_len 40``` 
	+  ```-noniupac```
	+  ```- min_qual_mean 15``` 
	+  ```-lc_method entropy```
	+  ```-lc_threshold 60```
	+  ```-trim_ns_right 10``` 
	+  ```- ns_max_p 20```
* The project gitub site has custom perl script wrappers that encompass a few steps: 
	+ Trimming and BLAST to ESTs (or other resources from which probes were designed), including making a custom BLAST database from your probes and pulling out top hits to the database from your reads using ```Top_Hits.pl```
	+ Wrapper is called: ```Optional_Blast_Wrapper.pl``` and there is a README: ```README for Optional_Blast_ Wrapper.pdf```

 
###4. Assembly 
* The second perl wrapper for this workflow is called ```Optional_Velvet_Prep_Wrapper.pl``` and runs Velvet Optimiser on the BLASTed and filtered reads. 
	+ [Velvet](https://github.com/dzerbino/velvet/tree/master)/[Velvet Optimiser](https://github.com/tseemann/VelvetOptimiser)
	+ Velvet Optimiser runs Velvet with a variety of different kmer values and picks the best assembly.
* Other options for assembly:
	+ [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki): built for transcriptomes but can also work well for UCEs/exons (*i.e.* not contiguous nuclear genomes).
	+ [SPAdes](http://bioinf.spbau.ru/spades): Nice assembler that is great for small genomes and can also work for UCEs/exons.
* It is worth trying more than one assembler on your data to see what works best. The only "danger" of diverging from the perl wrappers is that the output files will look slightly different and may require to adjust scripts downstream. But this is nothing to be afraid of!
* Chloroplast bycatch
	+ [MITObim](https://github.com/chrishah/MITObim): Purpose built to assemble mitochondrial genomes, but the idea works for chloroplasts as well. You give it a reference chloroplast sequence, and iteratively map reads to the reference and *de novo* assemble with [MIRA](http://mira-assembler.sourceforge.net/docs/DefinitiveGuideToMIRA.html).
	+ [Geneious](http://www.geneious.com): not free. Geneious has a number of tools you can use in a GUI interface for mapping and assembly. Handy if you are not comfortable using command line tools. 

###5. PHYLUCE  
###6. Alignment and phylogenetic analysis  
###7. Gene tree/species tree methods
###8. Visualization