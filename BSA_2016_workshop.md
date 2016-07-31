#Targeted next-gen sequencing for plant phylogenomics: the bioinformatics
##Botany 2016 Workshop

###Authors: 
**Jennifer R. Mandel** (jmandel@memphis.edu)  
**Rebecca B. Dikow** (DikowR@si.edu)  
**Vicki A. Funk** (FunkV@si.edu)  
Please contact us with bug reports or other suggested improvements.

###Objectives:  
In the course of this workshop, we detail the software and commands you will need to generate and analyze target enrichment data for phylogeny. At each point in the workflow, there are multiple options for software and/or commands, so while we try to present these options, we can't guarantee that we have been exhaustive. The tools are quickly evolving. Other resources attached to this workshop include powerpoints about wet-lab procedures and expedition protocols.

###References:
* Mandel *et al.*, 2014
* Mandel *et al.*, 2015
* [Compositae COS workflow](https://github.com/Smithsonian/Compositae-COS-workflow)
* Other software links included below

###Outline: 
1. Accessing code via GitHub 
2. Probe design
3. Quality trimming & BLAST filtering  
4. Assembly  
5. PHYLUCE/orthology detection 
6. Alignment  
7. Phylogenetic analysis
8. Visualization   

###1. Accessing code via GitHub 
* Github is one of the best places to go for code, scripts, and software. 
* This project has a GitHub site where there are instructions and the wrapper scripts used for our workflow: [Compositae COS workflow](https://github.com/Smithsonian/Compositae-COS-workflow)
* We recommend installing [GitHub Desktop](https://desktop.github.com) for its ease of use.

###2. Probe design
* The first step when you undertake or even think about undertaking a target (exon) enrichment study is to identify all the genomic resources that exist for your group: genomes, transcriptomes, ESTs and to think about the phylogenetic breadth you would like to consider. You will then need to identify orthologs across the genomes/transcriptomes and pick the best candidates for probe design.

* **Orthology Detection:**  
	+ [OMA](http://omabrowser.org/standalone/): with OMA, you can find *de novo* orthologs from transcriptome data with or without genomes.
	+ [AGALMA](https://bitbucket.org/caseywdunn/agalma): AGALMA is a pipeline that generates orthologs and species trees from transcriptomes. It incorporates [PhyloTreePruner](https://sourceforge.net/projects/phylotreepruner/), which distinguishes between in- and out- paralogs.
	+ [Orthofinder](http://www.stevekellylab.com/software/orthofinder): Orthofinder is a program for identifying orthologous protein sequence families. Input is a directory of FASTA files, one per species. Output is a file containing the orthologous groups of genes.
	+ [OrthoMCL](http://orthomcl.org): database that includes ortholog groups for a curated set of taxa. 
	+ [HaMStR](https://sourceforge.net/projects/hamstr/): extends a set of orthologs (e.g. from OrthoMCL) to additional taxa using hidden markov models.

* **Post Orthology Detection:**
	+ Check for approximate copy number and intron boundaries by mapping orthologs to a reference genome. Tools: [BWA](http://bio-bwa.sourceforge.net), [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml), [GSNAP](http://research-pub.gene.com/gmap/).
	+ Choose single or low copy-number orthologs.
	+ Align orthologs across the taxonomic divergence you want to capture. [MYcroarray](http://www.mycroarray.com) suggests around 10% divergence to guarantee capture. Send alignments to MYcroarray and order probes.

###3. Quality trimming & BLAST filtering 
* The project gitub site has custom perl script wrappers that group some of the steps into single commands. 
	+ ```Optional_Blast_Wrapper.pl``` includes trimming in Prinseq as detailed above and BLASTs the trimmed reads to ESTs or other resources from which probes were designed (our probes were designed from ESTs because it was a long time ago!) and finally pulling out the BLAST top hits. It also includes a README: ```README for Optional_Blast_ Wrapper.pdf```

####Wrapper syntax
+ **Commands:** ```perl Optional_Blast_Wrapper.pl``` 
+ **Arguments:**
+ ```-name <names_file>```
+ ```-evalue <#>	``` Default = 0.00001
+ ```-database <BLAST_database>``` Default, full sequence from which probes were designed(```COS_sunf_lett_saff_all.fasta```)
+ ```-c <T/F>	```	Specify if files are compressed. Default = F/False 

**While you may choose to use the wrapper that includes these steps, it's good to understand what it's doing and what parameters you might think about changing. Software tools are cropping up all the time, and flexibility is key in order to get the most robust and current results.**
	
####Trimming/quality tools 
* There are many tools that can be used for quality assessement and read trimming. Depending on the sequencing platform, you may need to trim adapters as well as low quality bases.
* [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/): FASTQC is a program that allows you to visualize the quality of your raw and/or trimmed data. Very easy to run, even on a laptop or desktop.
* [TrimGalore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/): I like TrimGalore! because it auto-detects adapters.
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic): Commonly used (and built into the PHYLUCE pipeline) but uses Java, which can be problematic on HPC clusters.
* [Prinseq](http://prinseq.sourceforge.net) We used Prinseq in Mandel (2014; 2015) and here are the commands: 
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

####BLAST filtering
* This was originally written to use NCBI BLAST which has been superseded by [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), so I will give both sets of commands.  

	####NCBI BLAST
	+ **Commands:** ```blastall```
	+ **Arguments:** 
	+ ```-i <cleaned_seq_file.fa>``` 
	+ ```-d <Database>``` (e.g. ```COS_sunf_lett_saff_all.fasta```)
	+ ```-o <Outfile.blast>```
	+ ``` -p blastn``` 
	+ ```-e 0.00001```
	+ ```-v 1``` 
	+ ```-b 1``` 
	+ ```-m 8```
	
	####BLAST+ : here is the analogous usage for BLAST+
	+ **Commands:** ```blastn```
	+ **Arguments:** 
	+  ```-query <cleaned_seq_file.fa>``` 
	+  ```-db <Database>``` 
	+  ```-evalue 0.00001``` 
	+  ```-num_descriptions 1``` 
	+  ```-num_alignments 1``` 
	+  ```-out <Outfile.blast>``` 
	+  ```-outfmt 8```
	+  Also check [here](http://www.ncbi.nlm.nih.gov/books/NBK279675/) for all command line options.

####After BLASTing reads to the resources from which probes were designed, we want to pull out only those reads that match best.  
 + **Commands:** ```perl top_hits.pl```  
	+ **Arguments:** 
	+ ```<Outfile.blast>```  
	+ ```<out_file>``` 
	+ ```<# of hits wanted>``` (e.g. 1)
	+ ```<cleaned_seq_file.fa>```
 
###4. Assembly 
* The second perl wrapper for this workflow is called ```Optional_Velvet_Prep_Wrapper.pl``` and prepares BLASTed and filtered reads for Velvet Optimiser. 
+ **Commands:** ```perl Optional_Velvet_Prep_Wrapper.pl``` 
+ **Arguments:** <names_file>
+ The names file should look like this: ```TaxonNameR1=TaxonNameR2=TaxonName```
+ Example on the GitHub site: ```Optional_Velvet_Wrapper_Names_File_Example```
+ The wrapper also runs the two perl scripts ```Pairfq.pl``` and ```Shufflesequences.pl``` which pair and shuffle sequences, important for Velvet. These steps would not be necessary for other assembly programs.
	
####Running VelvetOptimiser
+ [Velvet](https://github.com/dzerbino/velvet/tree/master)/[Velvet Optimiser](https://github.com/tseemann/VelvetOptimiser)
	+ Velvet Optimiser runs Velvet with a variety of different kmer values and picks the best assembly.
	+ **Commands:** ```perl VelvetOptimiser.pl``` 
	+ **Arguments:** 
	+ ```-s 51 -e 71``` (kmer range)
	+  ```-o'-ins_length 350'``` (insert size)
	+  ``` -f'-fasta -shortPaired <ReadName>_paired.fasta -short <ReadName>_singletons.fasta'``` (paired and single files)

* Other options for assembly:
	+ [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki): built for transcriptomes but can also work well for UCEs/exons (*i.e.* not contiguous nuclear genomes).
	+ [SPAdes](http://bioinf.spbau.ru/spades): Nice assembler that is great for small genomes and can also work for UCEs/exons.
* It is worth trying more than one assembler on your data to see what works best. The only "danger" of diverging from the perl wrappers is that the output files will look slightly different and may require to adjust scripts downstream. But this is nothing to be afraid of!

####Chloroplast bycatch  
+ Depending on the coverage depth you achieve, you may be able to assemble whole or mostly complete chloroplast genomes from the off target reads. Here are a few tools:
+ [MITObim](https://github.com/chrishah/MITObim): Purpose built to assemble mitochondrial genomes, but the idea works for chloroplasts as well. You give it a reference chloroplast sequence, and iteratively map reads to the reference and *de novo* assemble with [MIRA](http://mira-assembler.sourceforge.net/docs/DefinitiveGuideToMIRA.html).
	+ [Geneious](http://www.geneious.com): not free. Geneious has a number of tools you can use in a GUI interface for mapping and assembly. Handy if you are not comfortable using command line tools. 

###5. PHYLUCE/orthology detection  
* [PHYLUCE](http://phyluce.readthedocs.org/en/latest/tutorial-one.html) was developed by Brant Faircloth for the analysis of UCE data, but we used parts of this pipeline: LASTZ matching of assemblies to the probe set and selecting only those loci with single matches. The entire PHYLUCE pipeline, however, encompasses everything from trimming to assembly to alignment and matrix generation. 
	+ When we first used PHYLUCE for Mandel *et al.*, 2014, it used Velvet for assembly. Now the new version of PHYLUCE supports Trinity for UCE locus assembly. If you want to use PHYLUCE out of the box on target enrichment data with no modification, you probably should use Trinity because it requires the Trinity headers to be there in the assembly FASTA files. It is not hard to modify, however, so don't let what's hard coded in the software deter you from changing things up.
	+ To do the LASTZ portion of PHYLUCE, you will need to have your probes in a file such as this example on GitHub: ```COS_probes_phyluce.fasta``` 
	+ It's just a FASTA file with the probes listed with the headers labeled for each locus such as:
	+ ```>uce-1_p1```  
	+ ```GCGAAGGGGACGACAAAATCATA```  

+ Then you will run the PHYLUCE ```phyluce_assembly_match_contigs_to_probes.py``` program specifying the directory to your assemblies ```--contigs```, the probe file mentioned above ```--probes```, and an output directory ```--output```. This program in part runs [LASTZ](http://www.bx.psu.edu/~rsharris/lastz/) to locally align contigs to probe sequence. 

* We've done things two ways from here, first, as in Mandel *et al.*, 2014, we considered only those contigs with single matches to the probe set. If a two contigs matched a single probe, then we did not include either contig. What this lead to was a large amount of missing data, because as everyone in this room knows, there are very few truly single copy loci in plants, especially when you are considering a fairly deep divergence.
	+ To continue with this conservative approach, you will run the PHYLUCE ```get_match_counts.py``` and ```get_fastas_from_match_counts.py``` programs and then go onto multiple sequence alignment of your loci.
* For the second approach, we also consider multi-copy loci. What we mean by multi-copy loci is those loci that have multiple probe matches per taxon. For example, for COS-locus-1, most taxa might have only a single matching contig, while a few species may have two, three, or more. We might be interested in including all of these copies rather than none, which is what will happen in PHYLUCE. Obviously this complicates the phylogenetic analyses, but it is a way to start considering copy number and multi-copy genes. 
	+ We take the LASTZ output from the ```phyluce_assembly_match_contigs_to_probes.py``` program and run the ```lastz_uce.pl``` script from the Compositae-COS GitHub.
	+ **Commands:** ```perl lastz_uce.pl```
	+ **Arguments:** 
	+ ``` <names file>```
	+ ```<total # UCEs>```
	+ The output will be ```Final_Uce.fa```, which contains all loci for all taxa.
* We then run [USEARCH](http://www.drive5.com/usearch/) to cluster the matches, first sorting the sequences by length and then clustering at low stringency (65% similarity), then aligning clusters, which may have multiple sequences per species. In USEARCH, first we sort the sequences by length and then cluster to save RAM.
	+ **Commands:** ```usearch -sortbylength <query.fasta> -fastaout <seqs_sorted.fasta>```
	+ **Commands:** ```usearch -cluster_smallmem <seqs_sorted.fasta> -id 0.65 -msaout <output.fasta>```

###6. Alignment
* There are many ways to perform multiple sequence alignment but we prefer [MAFFT](http://mafft.cbrc.jp/alignment/software/). We think it does better with gaps than other options.
	+ To run MAFFT on its own, here are the commands:
	+ ```mafft in.fasta > out.fasta```
	+ MAFFT automatically determines what alignment algorithm to use for your data, and we think it does a good job so there should be no reason to adjust parameters further.
* PHYLUCE has the option to align your loci after orthology detection using ```seqcap_align_2.py```, and during alignment gives the option of trimming alignment ends that only include few taxa, or to trim poorly aligning regions in the middle of loci with [GBLOCKS](http://molevol.cmima.csic.es/castresana/Gblocks.html).
* Deciding whether you want to trim your alignments really depends on your data and your taxa. Since PHYLUCE was developed for UCEs rather than exons (UCEs are designed outside of exon boundaries meaning that the conserved portion of the locus is in the center and divergence increases the further you go out into the flanking region), it makes sense to trim ends. For exon capture, you can try it both ways and see what works best for you.


###7. Phylogenetic Analysis  
* We recommend maximum likelihood analysis with [RAxML](http://sco.h-its.org/exelixis/software.html) or [GARLI](https://code.google.com/archive/p/garli/), for both concatenated matrices and individual loci. 
* For Bayesian trees, [PhyloBayes](http://www.phylobayes.org) seems best for big matrices.

* For species tree analysis, there are many options, some of which we detail here:
	+ [ASTRAL](https://github.com/smirarab/ASTRAL): ASTRAL takes a set of "best tree" gene trees and creates a species tree. It is esay to use because it allows loci to have missing taxa, which is common for large target enrichment datasets (particularly in plants where coy number is a big problem). It scales very well to large datasets. 
	+ [*BEAST](http://beast.bio.ed.ac.uk): *BEAST simultaneously estimates gene and species trees but does not scale well to many loci and does not allow missing taxa.
	+ [BUCKy](http://www.stat.wisc.edu/~ane/bucky/) BUCKy is a method to generate species trees from a distribtion of MrBayes trees for each gene, not just the consensus tree. It scales well to large datasets.


###8. Visualization
* Network diagrams with [SplitsTree](http://www.splitstree.org): We used SplitsTree to generate networks for the multi-copy loci datasets using the ```SuperNetwork``` option.
* Heatmaps with [Gitools](http://www.gitools.org): We used Gitools to create heatmaps showing the distribution of copy number per taxon across all our loci. Input is a CSV file.
* [Densitree](https://www.cs.auckland.ac.nz/~remco/DensiTree/) is another way to visualize gene tree discordance. With a pretty low level of discordance or few different signals, the diagrams can be very informative. If discordance is high, and there are many different signals, they can look a bit chaotic. It all depends on your data!