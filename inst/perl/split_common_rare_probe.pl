#!/usr/bin/perl

my ($array_file,$cnv_probes) = @ARGV;
#my ($cnv_probes) = "/nfs/ddd0/Tom/CODE/files/031220_031221_wtccc_and_sequenom_probe_ids.txt";

my %common; my %lines;
open IN, "$cnv_probes" || die; while (<IN>) { chomp $_; $common{$_}[0] = "CNV"; } close IN;

my $cName = substr($array_file,0, length($array_file)-4) . "commonCNVprobes.txt"; my $rName = substr($array_file,0, length($array_file)-4) . "rareCNVprobes.txt";
open IN, "$array_file" || die; open OUT1, ">$cName" || die;  open OUT2, ">$rName" || die; 
while (<IN>) { 
chomp $_; my @ar = split("\t", $_); 
	if(defined $common{$ar[7]}) {
		print OUT1 "$_\n";
	} else {
		print OUT2 "$_\n";
	}
} close IN; close OUT1; close OUT2;