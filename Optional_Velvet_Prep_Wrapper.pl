###############
# Written by: Rishi Masalia, Burke Lab, Univ. of Georgia, August 2013
# Corresponding email: jennifer.r.mandel@gmail.com or masalia@uga.edu
#
# Please see README_Optional_Velvet_Prep_Wrapper.pl file for more detailed instructions.
#
# Goal: The purpose of this optional wrapper is to take the output from Top_Hits.pl (end of Optional_Blast_Wrapper, and prepare both Taxon reads (R1 and R2) for Velvet Optimizer.
#
# Usage:
#       perl Optional_Velvet_Prep_Wrapper.pl <names_file>
#
# Names_File: This is a separate file, with just a list of Taxon R1 and R2, as well as the Taxon Name. Each of these elements should be separated by a '=' sign. Example: TaxonNameR1=TaxonNameR2=TaxonName
###############
$file1 = $ARGV[0]; #names file
open (NAME, $file1) or die "\n\n**********\nNo Names File Found\n**********\n\n";
while (<NAME>){
        chomp();
        @pair = split (/=/, $_, 3);
        $forward = $pair[0];
        $reverse = $pair[1];
        $together = $pair[2];
        if ($forward eq "" || $reverse eq "" | $together eq ""){
                print "\n\n**********\n$file1 format incorrect. Correct Format is: R1=R2=Name\n**********\n\n";
                exit;
        }


#### Pairing up Forward and Reverse Sequences, and singletons ###

        system ("perl ./programs/pairfq.pl -f ./Output/$forward.blasted_reads.fasta -r ./Output/$reverse.blasted_reads.fasta -fp $forward.paired.fasta -rp $reverse.paired.fasta -fs $forward.singletons.fasta -rs $reverse.singletons.fasta");

        system ("./programs/shuffleSequences_fasta.pl $forward.paired.fasta $reverse.paired.fasta $together.together_paired.fasta");


### Concatenate Singleton files together ###

        system ("cat $forward.singletons.fasta $reverse.singletons.fasta > $together.together_singletons.fasta");


### Clean Up ###

        system ("mv $together.together_singletons.fasta Output/");
        system ("mv $together.together_paired.fasta Output/");
        system ("mv $forward.paired.fasta Archived_Files/");
        system ("mv $reverse.paired.fasta Archived_Files/");
         system ("mv $forward.singletons.fasta Archived_Files/");
         system ("mv $reverse.singletons.fasta Archived_Files/");

}
