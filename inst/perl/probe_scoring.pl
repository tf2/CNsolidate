#/usr/bin/perl

use strict;

# Insert probe scores into the general "bed" file format for the DDD arrays.
# NB. if probe score doesnt exist, true for ~300 probes on each array - score => 0.5

my ($dfile, $ploc) = @ARGV;

my $score_20 =  $ploc . "/31220_scores.txt";
my $score_21 =  $ploc . "/31221_scores.txt";
my $sfile = "";
if($dfile =~m/.*?2531220.*?/) {
	$sfile = $score_20;
} elsif($dfile =~m/.*?2531221.*?/ ) {
	$sfile = $score_21;
} else {
	die "There is a problem with the filenames - please try again!";
}
	my %hash;
	open IN, "$sfile" || die "can't open $sfile $!";
	while (<IN>) {
		chomp $_; my @a=split("\t", $_);
		$hash{$a[0]}[0]=$a[1];
	} close IN;
	open IN, "$dfile" || die  "can't open $dfile $!";
	my @data = <IN>;
	close IN;
	open OUT , ">$dfile" || die  "can't open $dfile $!";
	for my $i (0 .. $#data) {
		chomp $data[$i]; my @a=split("\t", $data[$i]);
		if(defined($hash{$a[9]})) { print OUT "$data[$i]\t$hash{$a[9]}[0]\n"; } else { print OUT "$data[$i]\t0.5\n"; }
	} close OUT;