#/usr/bin/perl

use strict;

my ($dfile, $rfile) = @ARGV;

	open IN, "$dfile" || die  "can't open $dfile $!";
	my $name = substr($dfile,0, length($dfile)-4); open OUT , ">" . $name . ".gff"  || die  "can't open $dfile $!";
	while(<IN>) {
		chomp $_; my @a=split("\t", $_);
		if ($#a < 8) { die "file format is not expected!"; }
		my $color = ""; my $score = $a[8];
       	if ($score > 0.8) {
        	$color = "; color 000000";    
    	} elsif ($score <= 0.8 & $score > 0.6) {
            $color = "; color 0000FF";
    	} elsif ($score <= 0.6 & $score > 0.4) {
       		$color = "; color 00FF00";
       	} elsif ($score <= 0.4 & $score > 0.2) {
            $color = "; color FFFF00";
        } elsif ($score <= 0.2) {
        	$color = "; color FF0000";
        }	
		print OUT "$a[0]\tDNAcopy\t$dfile\t$a[1]\t$a[2]\t$a[3]\t" . "." . "\t" . "." . "\t" . "$color\n";
		
	} close IN;
	open IN, "$rfile" || die  "can't open $rfile $!";
	while(<IN>) {
		chomp $_; my @a=split("\t", $_);
		if ($#a < 8) { die "file format is not expected!"; }
		my $color = "; color FF0000"; my $score = $a[3];
      	# if ($score >=0) {
       	# 	$color = "; color 0000FF";    
    	#}
    	print OUT "$a[0]\tDNAcopy\t$dfile\t$a[1]\t$a[2]\t$a[3]\t" . "." . "\t" . "." . "\t" . "$color\n";	
    	print OUT "$a[0]\tDNAcopy\t$dfile" . "LOOK" . "\t$a[1]\t$a[2]\t$a[8]\t" . "." . "\t" . "." . "\t" . "; color FF0000\n";
    	print OUT "$a[0]\tDNAcopy\t$dfile" . "LOOK" . "\t$a[1]\t$a[2]\t$a[9]\t" . "." . "\t" . "." . "\t" . "; color 0000FF\n";
	} close IN;
	close OUT;
	