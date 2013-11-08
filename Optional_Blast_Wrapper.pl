###############
# Code written by: Rishi Masalia, Burke Lab, Univ. of Georgia, August 2013
# Corresponding email: jennifer.r.mandel@gmail.com or masalia@uga.edu 
#
#Goal: This wrapper will take raw sequence reads for a species or individual (fastq). Clean them, both reads 1 and 2 with prinseq-lite.pl (http://prinseq.sourceforge.net/manual.html). BlAST these reads against a database of your choice. Capture the top hit and create a file with the best hit and sequence, which will be assembled farther along in the workflow
#
#Input: perl Optional_Blast_Wrapper.pl [parameters] 
#
#Parameters:
#	-name <file>		This is to specific the names file. Critical file.
#	-evalue <#>		To specific the blast evalue. Default = 0.00001
#	-database <file>	Specify which database to run. Default, Full EST sequence that the probes were designed from. "COS_sunf_lett_saff_all.fasta"
#	-c <T/F>		Specify if files are compressed. Default = F|False 
#
#For instructions on things to do prior to running the script please refer to the README for this script (Optional_Blast_Wrapper.pl)
###############
use Getopt::Long;

system ("mkdir Output/");
system ("mkdir Archived_Files/");

GetOptions( \%params, 'name=s','evalue=i', 'database=s', 'c=s');

if($params{name}){
        $namefile = $params{name};
        $evalue = $params{evalue};
        $zipped = $params{c};
        $database = $params{database};
	open (NAMES, $namefile);
	while (<NAMES>){
        	chomp();
        	$name = $_;

        ##### Cleaning the Sequences ######

       if ($zipped eq "" | $zipped eq "F" | $zipped eq "f"){
                system ("perl ./programs/prinseq-lite.pl -fastq $name.fq -out_format 1 -out_good ./$name.trimmed -out_bad null -log $name.dirty.log -min_len 40 -noniupac -min_qual_mean 15 -lc_method entropy -lc_threshold 60 -trim_ns_right 10 -ns_max_p 20");
        } elsif ($zipped eq "T"| $zipped eq "t"){
                system ("gzip -dc $name.fastq.gz | perl ./programs/prinseq-lite.pl -fastq stdin -out_format 1 -out_good ./$name.trimmed -out_bad null -log $name.dirty.log -min_len 40 -noniupac -min_qual_mean 15 -lc_method entropy -lc_threshold 60 -trim_ns_right 10 -ns_max_p 20");
        } else {
                print "Error in -c parameter\n";
        }

        #### Blast to dataset ###
        if ($evalue eq ""){
		$evalue = "0.00001"; #if evalue is specified use that, if not use default.
	}        
	if ($database eq ""){
                $database = "COS_sunf_lett_saff_all.fasta"; #if database is specified use that, if not use default
        } 
        

	system ("./programs/blastall -i $name.trimmed.fasta -d ./databases/$database -o $name.blast -p blastn -e $evalue -v1 -b1 -m8");

        
	### Get Best Hit ###
        system ("perl ./programs/Top_Hits.pl $name.blast $name.blasted_reads.fasta 1 $name.trimmed.fasta");

        ### Clean up the Folders ####
        system ("mv $name.blast Archived_Files/");
        system ("mv $name.trimmed.fasta Archived_Files/");
        system ("mv $name.dirty.log Archived_Files");
        system ("mv $name.blasted_reads.fasta Output/");


	}
}else {
        print ("*******ERROR: There is no names file present*******************\n");
exit;
}

