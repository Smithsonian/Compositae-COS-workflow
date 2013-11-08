#!/usr/bin/perl
###################
#Code Written by:  Rishi R. Masalia, Burke Lab, Dept. Plant Biology, Univ. of Georgia
#Corresponding Email: jennifer.r.mandel@gmail.com or  masalia@uga.edu
#
#
#Purpose: To take a blast output and parse the information for a best hit(s) determined by evalue and percent identity, print the output. For this, once the blast output is parsed, only the name taken and compared to an origianl source file to extract name and sequence information from a cleaned Illumina sequence file. 
#
#Commands: perl top_hits.pl <in_file> <out_file> <# of hits wanted> <Cleaned Sequence File>
# The cleaned sequence file is the cleaned fasta file generated after the Prinseq-lite.pl cleaning stage.
#
# This is a necessary script for the Optional_Blast_Wrapper.pl however, if you wish to perform this stage on your own, simply use your favorite Blast Output Parser and match the lines to the Cleaned the Sequence Fasta File. If you wish to use the Optional_Velvet_Wrapper.pl, please make sure the file names are corresponding. 
#
####################

my $in = $ARGV[0];
my $out = $ARGV[1];
my $hits = $ARGV[2];
open (OUT, ">$out"); 

### Open the blast file and get the corresponding number of hits, deteremined by evalue and percent identity ###

my %blasthash = ();
open (INFILE, $in);
while (<INFILE>){
        chomp();
        my $line = $_;
        my @row = split(/\t/,$line);
        my $name = $row[0];
        my $identity = $row[3];
        my $evalue = $row[10];
        $blasthash{$name}{$evalue}{$identity} = $line;
}
close IN;

### Get the # of top hits, by sorting name, evalue and identity ### 
%OrigSeqhash = ();
foreach my $name (sort keys %blasthash){
        my $hits_count = 0;
        foreach my $evalue (sort {$a<=>$b} keys %{$blasthash{$name}}){
                foreach my $identity (sort {$b<=>$a} keys %{$blasthash{$name}{$evalue}}){
                        if ($hits_count < $hits){
                        $OrigSeqhash{$name} = 1; # puts the name into a hash, to be filtered through later
                        }
                my $hits_count++;
                }
        }
}

### Match the name hash with original trimmed file to pull out name and sequence ###

my $cleanseq = $ARGV[3];
open (TRIMMED, $cleanseq);
while(<TRIMMED>){
        chomp();
        $line = $_;
        $part1 = $line;
        $part1 =~ s/\>//g;
        $part1 =~ s/\s.*//g;
        $n = 1 if ($OrigSeqhash{$part1});
                if ($n >=1 and $n <=3){
                        print OUT "$line\n";
                        $n++;
                }




}
exit;
