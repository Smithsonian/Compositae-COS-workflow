################
# Code Written by: Rishi R. Masalia, Burke Lab, Dept. Plant Biology, Univ. of Georgia
# Corresponding Email: jennifer.r.mandel@gmail.com or masalia@uga.edu
#
# Purpose: Finding all copies of COS loci for multigene analysis. To take lastz and contig.fasta output files, and create a large fasta file for Phyluce. This will be comprised of individual "uce" (up to user specified) for all species involved, with corresponding target enrichment sites and corresponding sequence
#
# Commands: perl lastz_use.pl <names file> <total # UCEs>
# Necessary: You must create a "names" file, which has the names of the speices for input into phyluce; eg: phoeb-l or gerb-e -- each of these names should be on separte lines. Additionally, ALL lastz files should end in $name.contigs.lastz, while all contig files should end in, $name.congits.fa and all of these files per species should be in the main directory.
#
# For convience this script archives all species specific information (ie: UCEs per species) and places that information in the "Species_Folder"
#
# Output: This script outputs one catcatenated file, "Final_Uce.fa", which has all UCE/Sequence information for all species involved.
######################################

$file = $ARGV[0];
$UCENUM = $ARGV[1];
open (NAME, $file) || die "Can't open a Names file: $!";
system ("mkdir Species_Folder");
while(<NAME>){
	chomp();
	$name = $_;
	print "********* SPECIES\t $name\n";
	%seen = ();
	system ("mkdir $name");

###### Open lastz file ####

	for ($count = 1; $count <= $UCENUM; $count++){
	
	print "Count is:\t $count\n";
	open (TMP, ">uce-$count.tmp");
	open (LASTZ, "$name.contigs.lastz");
	while(<LASTZ>){
		chomp();
		@line = split(/\s+/, $_);
		$uce = $line[6];
		$node = $line[1];
		$uce =~ s/\>uce-//g;
		$uce =~ s/_.*//g;
		
		if ($uce eq $count){
		push (@array, $node) unless $seen{$node}++;
		}	
	}
	print TMP @array;
	@array = ();
	%seen = ();
	close TMP;

### Create tmp files to house uce-count information ###
	open (TMP1, "uce-$count.tmp");
	open (OUT1, ">uce-$count.outtmp");
	while (<TMP1>){
		chomp();
		$_ =~ s/\>/\n/g;
		print OUT1 "$_";
	}
	close TMP1;
	close OUT1;
	open (TMP2, "uce-$count.outtmp");
	while(<TMP2>){
		chomp();
		$node1 = $_;
		$hash{$node1} = 1;
		
	}

	close TMP1;
	open (CONTIG, "$name.contigs.fa");
	open (FINAL, ">uce-$count.$name.txt");
	
	local $/ = '>';{ 
	while (<CONTIG>){
		$row = $_;
		@line1 = split (/\n/, $row, 2);
		$seq = $line1[1];
		$seq =~ s/\s+//g;
		$seq =~ s/\>//g;
		$header = $line1[0];

		if ($hash{$header}){
			$header =~ s/NODE.*cov/cov/g;
			$newname = $name;
			$newname =~ s/\-/\_/g;
			print FINAL "\>uce-$count\_$newname\|uce-$count\|$header\n$seq\n";
		}


	}
	}
	%hash = ();

	}
system ("rm *.tmp");
system ("rm *.outtmp");
system ("mv *.txt $name/");
system ("find ./$name/ -size 0 -delete");

system ("cat ./$name/*.txt > $name.cat");
system ("tar -zcf $name.tar.gz ./$name/"); 
system ("rm -r $name/");
}
system ("cat *.cat > Final_Uce.fa");
system ("mv *.cat ./Species_Folder/");
system ("mv *.tar.gz ./Species_Folder/");






