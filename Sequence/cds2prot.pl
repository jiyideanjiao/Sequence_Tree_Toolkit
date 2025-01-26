#!/usr/bin/perl -w
use strict;

unless ( $#ARGV >= 0 ){
	print "Usage: perl cds2prot.pl cds.fa > prot.fa\n";
	exit;
}	
my ($id,$protein) = ('','');
my ($i,$j)=(0,0);
my @id_a = ();
my %seq = ();

#reading Seq
foreach my $file (@ARGV){
	open SEQ, "<","$file" or die "Can't open file \"$file\" : $!\n";
	while (<SEQ>){
		chomp;
		if (/^>/){
			$id = $_;
			push @id_a,$_;
		} else {
			$_ = uc($_);
			$_ =~ s/U/T/g;
			$seq{$id} .= $_;
		}
	}
	close SEQ;
}

#translate
foreach my $key (@id_a){
	$j++;

#The 3 line below use to find the sequence id with bad codon.
#$id = (split /\s+/,$key)[0];
#$id =~ s/^>//;
#print STDERR "$j. $id\n";
	for(my $i=0;$i<(length($seq{$key})-2);$i+=3){
		my $codon=substr($seq{$key},$i,3);
		$protein.=&codon2aa($codon);
  	}
	$protein =~ s/(.{80})/$1\n/g;
  	print "$key\n$protein\n";
    	$protein = "";
}
print STDERR "Completed.\n\n";

#codon to amino acid
sub codon2aa{
	my ($codon)=@_;
my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G','AAN'=>'X','GAN'=>'X','CAN'=>'X','TAN'=>'X','NAA'=>'X','NAG'=>'X','NAC'=>'X','NAT'=>'X','NAN'=>'X','AGN'=>'X','GGN'=>'X','CGN'=>'X','TGN'=>'X','NGA'=>'X','NGG'=>'X','NGC'=>'X','NGT'=>'X','NGN'=>'X','ACN'=>'X','GCN'=>'X','CCN'=>'X','TCN'=>'X','NCA'=>'X','NCG'=>'X','NCC'=>'X','NCT'=>'X','NCN'=>'X','ATN'=>'X','GTN'=>'X','CTN'=>'X','TTN'=>'X','NTA'=>'X','NTG'=>'X','NTC'=>'X','NTT'=>'X','NTN'=>'X','ANA'=>'X','ANG'=>'X','ANC'=>'X','ANT'=>'X','ANN'=>'X','GNA'=>'X','GNG'=>'X','GNC'=>'X','GNT'=>'X','GNN'=>'X','CNA'=>'X','CNG'=>'X','CNC'=>'X','CNT'=>'X','CNN'=>'X','TNA'=>'X','TNG'=>'X','TNC'=>'X','TNT'=>'X','TNN'=>'X','NNA'=>'X','NNG'=>'X','NNC'=>'X','NNT'=>'X','NNN'=>'X');
	if(exists $g{$codon}){
		return $g{$codon};
	}else{
		print STDERR "\nBad codon \"$codon\"!!\n";
		print STDERR "Not Completed!\n";
		exit;
	}
}
