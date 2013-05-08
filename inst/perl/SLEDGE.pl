#!/usr/bin/perl

#use strict;
use Getopt::Long;
use threads;
use JSON qw(to_json from_json);
use Data::Dumper;

my ($config, $type, $idir, $odir, $lsf, $loc, $chuck, $project, $self);

GetOptions('conf=s'=> \$config, 'input=s'=> \$idir, 'output=s'=> \$odir, 'lsf=s'=> \$lsf, 'project=s' => \$project, 'rlib=s'=> \$loc, 'type=s'=> \$type, 'self=s' => \$self);

if(!defined $lsf) { $lsf = 0; }
if(!defined $type) { $type = "bed"; }
if(!defined $project) { $project = "team70"; }
if(!defined $self) { $self = "NaN"; }
if(!defined $loc) { $loc = "."; }
if(!defined $idir) { $idir = "."; }
if(!defined $odir) { $odir = $idir; }
if(!defined $idir && !defined $odir) { $idir = "."; $odir = "."; }

print($type);

if (!defined $config) {
	$config= DefaultConfig();
	writeR($config);
	my $len = setProcessDirectories($idir, $odir);
	execute($config, $len);
} else {
	my $configRef = readConfig($config);	
	my $len = setProcessDirectories($configRef->{Process}->{inputdir}, $configRef->{Process}->{outputdir});
	writeR($configRef);
	execute($configRef, $len);
}

sub readConfig {
	my $file = shift;
	open IN, "$file";
	my $js;
	while(<IN>) {
		$js .=$_;
	}
	close IN;
	return(from_json($js));	
}

sub buildAlgos {
	my ($a) = shift;
	my %algo = %{$a}; my @ss;	
	foreach my $k (keys (%algo) ) {
		if($algo{$k} != 0)  {
			push(@ss, $k);
		}
	}
	return(@ss);
}

sub buildOrder {
	my ($aa) = shift;
	my %odi = %{$aa}; my @s;
	foreach $key (sort {$a<=>$b} (keys(%odi))) {
		push(@s, $odi{$key});
	}
	return(@s);
}

sub execute {
	
	my $config = shift;
	my $len = shift;
	
	my $command = "";	
	if($config->{Process}->{lsf} == 1) {
		$command = 'bsub -J' . '"SLEDGE[1-' . $len . ']%' . $config->{Process}->{chuck} . '"' . ' -P"' . $config->{Process}->{project} . '" -q ' . $config->{Process}->{que}
		. '  -R"select[mem>3500] rusage[mem=3500]" -M3500000 -o SLEDGE-%J.%I.out -e SLEDGE-%J.%I.err ' 
		. '"R CMD BATCH --slave ' .  $config->{Process}->{runnerloc} . '"';
		system($command);	
	} else {
		$command = "R CMD BATCH --slave $config->{Process}->{runnerloc}";
		system($command);
	}
	return($command);
}

sub setProcessDirectories {

	my $idir = shift;
	my $odir = shift;
	my $len=0;
	my @com = ();
	opendir DIR, $idir || die "Cant open input directory!!$!";
	my @dircont = readdir(DIR);
	closedir DIR;
	
		system("[ -d $odir ] && rm -r $odir");
		system("mkdir $odir");
		
		if("$config->{Process}->{formatted}" == "gff" || $type eq "gff") {
			for my $i (0 .. $#dircont) {	
				if($dircont[$i] =~ m/.*?gff/) {
					my $samname = $idir . "/" . substr($dircont[$i], 0, length($dircont[$i]));
					my $dirname = $odir . "/" . substr($dircont[$i], 0, length($dircont[$i])-4);
					push(@com, $dirname);
					system("mkdir $dirname");
					system("chmod 777 $dirname");
					system("cp $samname $dirname"); 
					$len++;
				}
			}
		}
		if("$config->{Process}->{formatted}" == "bed" || $type eq "bed" ) {
			for my $i (0 .. $#dircont) {	
				if($dircont[$i] =~ m/.*?bed/) {
					my $samname = $idir . "/" . substr($dircont[$i], 0, length($dircont[$i]));
					my $dirname = $odir . "/" . substr($dircont[$i], 0, length($dircont[$i])-4);
					push(@com, $dirname);
					system("mkdir $dirname");
					system("chmod 777 $dirname");
					system("cp $samname $dirname"); 
					$len++;
				}
			}
		}
		if("$config->{Process}->{formatted}" == "fe" || $type eq "fe" ) {
			for my $i (0 .. $#dircont) {	
				if($dircont[$i] =~ m/.*?txt/) {
					my $samname = $idir . "/" . substr($dircont[$i], 0, length($dircont[$i]));
					my $dirname = $odir . "/" . substr($dircont[$i], 0, length($dircont[$i])-4);
					push(@com, $dirname);
					system("mkdir $dirname");
					system("chmod 777 $dirname");
					system("cp $samname $dirname"); 
					$len++;
				}
			}
		}

	return($len);
}

sub DefaultConfig {

	my $config = {
	
		Details => {
					"Pipeline" => "Combinatorial change point detection --- SLEDGE",
					"Version" => "Version 1.0",
					"Settings" => "Default Configuration File",
					"Author" => "Tomas William Fitzgerald",
					"Email" => 'tf2@sanger.ac.uk',
					"URL" => "mySite.html",
	   	},
		
		Process => { "formatted" => "$type", "configloc" => "$idir/config.json", "runnerloc" => "$idir/runner.R", "libLoc" => "$loc", "project" => "$project", "lsf" => "$lsf", "chuck" => "100", "inputdir" => "$idir", "outputdir" => "$odir", "que" => "long", "wave" => 1, "segMet" => 1, gold=>"../extdata/42Mcalls_all_feb09.txt", gapfile=>"../extdata/gaps_hg19.txt", pipe=>"DDD"  },	
   		Norm => {"wFac" => 1.2, "sFact" => 4.5, "sthes" => 0.68, "sknot" => 1000, "sIt" => 5, "med_norm" => "global", "self" => "$self"},
   		Algorithms => { "ADM2" => 1, "ADM3" => 1, "GADA" => 1, "SMAP" => 1, "SMUG" => 1, "dnacopy" => 1, "fsegment" => 0, "fast" => 0, "safe" => 0, "walk" => 0, "bcp" => 0, "CNCP" => 0 },
   		Sep => { "ADM2" => 0, "ADM3" => 0, "GADA" => 0, "SMAP" => 1, "SMUG" => 1, "dnacopy" => 1, "fsegment" => 1, "fast" => 1, "safe" => 1, "walk" => 0, "bcp" => 1, "CNCP" => 1 }, 
   		Combine => {"Consenues" => 2, "absRatio" => 0.4, "minProbe" => 3, "lRatio" => 0.4, "gRatio" => 0.5, "lest" => 250, "gest" => 200, "tfac" => 2},
   		#Order => { "1" => "spline", "2" => "self", "3" => "wave", "4" => "sep", "5" => "call", "6" => "FeObs", "7" => "mapBreak", "8" => "gap", "9" => "overlap", "10" => "weight", "11"=> "localp", "12"=> "image" },
   		#Order => { "1" => "spline", "2" => "wave", "3" => "self", "4" => "sep", "5" => "call", "6" => "FeObs", "7" => "mapBreak", "8" => "gap", "9" => "overlap", "10" => "weight", "11"=> "localp", "12"=> "image" },
   		Order => { "1" => "spline", "2" => "wave", "3" => "self", "4" => "sep", "5" => "call", "6" => "FeObs", "7" => "mapBreak", "9" => "overlap", "10" => "weight", "11"=> "localp", "12"=> "image" },
   		
   		
   		Settings => {	
   				safe => { "mL" => 3, "mR" => 0.35 },   				
   				fast => { "mL" => 3, "mR" => 0.35, "pT" => 1 },   				
   				GADA => { "alpha" => 0.1, "T" => 3, "mL" => 3, "mR" => 0.35 },
   				ADM2 => {"t" => 5, "np" => 3, "mr" => "0.35", "fact" => 3 },
   				ADM3 => {"t" => 5, "np" => 3, "mr" => "0.35", "nd" => 3 },
   				fsegment => { "sL" => 1000, "cpP" => 2, "mL" => 3, "mR" => 0.35 },   				
   				dnacopy => { "smR" => 2, "oSD" => 4, "sSD" => 2, "trim" => 0.025, "mL" => 3, "mR" => 0.35 },  		
   				SMAP => { "resolution" => 1000, "probC" => 0.05, "p1" => 0.95, "p2" => 0.01, "mL" => 3, "mR" => 0.35, "pT" => 3  },  			
   				SMUG => { "resolution" => 1000, "probC" => 0.05, "p1" => 0.95, "p2" => 0.01, "mL" => 3, "mR" => 0.35, "pT" => 3  },      				
   				bcp => { "resolution" => 1000, "probC" => 0.05, "p1" => 0.95, "p2" => 0.01, "mL" => 3, "mR" => 0.35, "pT" => 3  },    				
   				walk => { "tr1" => 2, "tr2" => 2, "l1" => 25, "l2" => 5, "w1" => 100, "w2" => 500, "w3" => 10, "s" => 500, "mR" => 0.35 },	
   				CNCP => { "fac" => 3, "a" => 3, "n" => 2, "absc" => 0.35, "prob" => 0.4, "dist" => 5000 }, 			
		}	
	};

	my $js = to_json($config);my $json = JSON->new->allow_nonref;	
	my $perl_scalar = $json->decode( $js ); my $pretty_printed = $json->pretty->encode( $perl_scalar );	
	open OUT, ">" . $config->{Process}->{configloc} || die;
	print OUT $pretty_printed;
	close OUT;	
	writeR($config);
	
	return($config);
}

sub writeR {

	my $config = shift;						
	my $c .= makeHeader($config);
	
	my $Command .= ""; $Command .="$c\n
	format = " . "'" . $config->{Process}->{formatted} . "'\n";
	
	my @order = buildOrder($config->{Order});
	for my $i (0 .. $#order) {
	
		if($order[$i] eq "spline") {
			$Command .= "\tjspline(file, idir, format, $config->{Norm}->{sthes}, $config->{Norm}->{sFact}, $config->{Norm}->{sknot}, $config->{Norm}->{sIt})\nformat='bed'\n"; 
			$Command .= "\tdata = inData(file, format)\n";
		}
		if($order[$i] eq "self") {
			$Command .= "\tself(data, file, " . "'" . $config->{Norm}->{self} . "'" .  ")\n"; 
			$Command .= "\tdata = read.table(file)\n";
		}
		if($order[$i] eq "wave") {
			$Command .= "\twave(file, idir, odir, $config->{Norm}->{wFac})\n"; 
			$Command .= "\tdata = read.table(file)\n";
		}
		if($order[$i] eq "call") {
			$Command .= makeParamters($config);
		}
		my @algo = buildAlgos($config->{Algorithms});
		my $s = "c("; for my $i(0 .. $#algo) { $s .= "'" . $algo[$i] . "'" . ","; } chop($s); $s .= ")";
		if($order[$i] eq "FeObs") {	
			$Command .= "\tmakeFObs(file, idir, odir, $config->{Combine}->{Consenues}, $s, $config->{Combine}->{absRatio}," . "'" . $config->{Process}->{gold} . "'" . ")\n";
		}
		if($order[$i] eq "gap") {	
			$Command .= "\tgap(file, idir, odir,$config->{Combine}->{absRatio}, $config->{Combine}->{minProbe}," . "'" .  $config->{Process}->{gapfile} . "'" . ")\n";
		}
		if($order[$i] eq "mapBreak") {	
			$Command .= "\tmapBreakc(file, odir, $config->{Combine}->{absRatio}, $config->{Combine}->{minProbe})\n";
		}
		if($order[$i] eq "overlap") {	
			$Command .= "\toverlap(file, idir, odir," . "'" .  $config->{Process}->{gold} . "'" . ")\n";
		}
		if($order[$i] eq "weight") {	
			$Command .= "\tweight.algorithmsP(file, odir)\n";
		}
		if($order[$i] eq "localp") {	
			$Command .= "\tlocalp(file, odir, $config->{Combine}->{tfac})\n";
		}
		if($order[$i] eq "image") {	
			$Command .= "\tcombin.image(file, odir)\n";
		}
	}
	writeScript($config, $Command);
}


sub makeHeader {

	my $config = shift;
	my $Header = "";
	if($config->{Process}->{lsf}) {	
	$Header .= "
	args<-commandArgs(TRUE)
	x = as.numeric(Sys.getenv('LSB_JOBINDEX'))
	fs = dir('" . $config->{Process}->{inputdir} ."')
	library(SLEDGE,lib.loc= '" . $config->{Process}->{libLoc} . "'" . ")
	files = dir('" . $config->{Process}->{outputdir} ."')
	idir = paste('" . $config->{Process}->{outputdir} . "', '/', files[x], '/', sep='')
	odir =idir
	file = paste(idir, '/', dir(idir), sep='')";
	} else {
	$Header .= "
	library(SLEDGE,lib.loc= '" . $config->{Process}->{libLoc} . "'" . ")
	files = dir('" . $config->{Process}->{outputdir} ."')
	for(y in 1:length(files)) {
	idir = paste('" . $config->{Process}->{outputdir} . "', '/', files[y], '/', sep='')
	odir = idir
	file = paste(idir, '/', dir(idir), sep='')";
	}
	return($Header);
}

sub makeParamters {

	my $config = shift;
	my %para1 = %{$config->{Algorithms}};
	my %para2 = %{$config->{Sep}};
	my @s; my @n; my $command="";
	
	foreach my $k (keys (%para1) ) {
		if($para1{$k} == 1 && $para2{$k} == 1) {
			push(@s, $k);
		} elsif ($para1{$k} == 1 && $para2{$k} == 0) {
			push(@n, $k);
		}
	}

	for my $i (0 .. $#n) {
		if($n[$i] eq "GADA") {
			$command .= "\tGADA_s(data, file, idir, odir, $config->{Settings}->{GADA}->{alpha}, $config->{Settings}->{GADA}->{T}, $config->{Settings}->{GADA}->{mL}, $config->{Settings}->{GADA}->{mR})\n";
		}
		if($n[$i] eq "ADM2") {
			$command .= "\tADM2(data, file, odir, $config->{Settings}->{ADM2}->{t}, $config->{Settings}->{ADM2}->{np}, $config->{Settings}->{ADM2}->{mr}, $config->{Settings}->{ADM2}->{fact})\n";
		}
		if($n[$i] eq "ADM3") {
			$command .= "\tADM3(data, file, odir, $config->{Settings}->{ADM3}->{t}, $config->{Settings}->{ADM3}->{nd}, $config->{Settings}->{ADM3}->{mr}, $config->{Settings}->{ADM3}->{np})\n";
		}
		if($n[$i] eq "SMUG") {
			$command .= "\tSMUG(data, file, idir, odir, $config->{Settings}->{SMUG}->{resolution}, $config->{Settings}->{SMUG}->{probC}, $config->{Settings}->{SMUG}->{p1}, $config->{Settings}->{SMUG}->{p2}, $config->{Settings}->{SMUG}->{mL}, $config->{Settings}->{SMUG}->{mR}, $config->{Settings}->{SMUG}->{pT}," . "'" . $config->{Process}->{gold} . "')" . " " . "\n";
		}
		if($n[$i] eq "SMAP") {
			$command .= "\tSMAP(data, file, idir, odir, $config->{Settings}->{SMAP}->{resolution}, $config->{Settings}->{SMAP}->{probC}, $config->{Settings}->{SMAP}->{p1}, $config->{Settings}->{SMAP}->{p2}, $config->{Settings}->{SMAP}->{mL}, $config->{Settings}->{SMAP}->{mR}, $config->{Settings}->{SMAP}->{pT}," . "'" . $config->{Process}->{gold} . "')" . " " . "\n";
		}
		if($n[$i] eq "fsegment") {
			$command .= "\tFsegment(data, file, idir, odir, $config->{Settings}->{fsegment}->{sL}, $config->{Settings}->{fsegment}->{cpP}, $config->{Settings}->{fsegment}->{mL}, $config->{Settings}->{fsegment}->{mR})\n";
		}
		if($n[$i] eq "dnacopy") {
			$command .= "\tdna.copy(data, file, idir, odir, $config->{Settings}->{dnacopy}->{mL}, $config->{Settings}->{dnacopy}->{mR})\n";
		}
		if($n[$i] eq "safe") {
			$command .= "\tsafe.dec(file, idir, odir, $config->{Settings}->{safe}->{mL}, $config->{Settings}->{safe}->{mR})\n";
		}
		if($n[$i] eq "walk") {
			$command .= "\tjv.walk2(file, idir, odir, $config->{Settings}->{SMUG}->{resolution}, $config->{Settings}->{SMUG}->{probC}, $config->{Settings}->{SMUG}->{p1}, $config->{Settings}->{SMUG}->{p2}, $config->{Settings}->{SMUG}->{mL}, $config->{Settings}->{SMUG}->{mR}, $config->{Settings}->{SMUG}->{pT}," . "'" . $config->{Process}->{gold} . "')" . " " . "\n";
		}
		if($n[$i] eq "CNCP") {
			$command .= "\trcncp(data, file, idir, odir, $config->{Settings}->{CNCP}->{fac}, $config->{Settings}->{CNCP}->{a}, $config->{Settings}->{CNCP}->{n}, $config->{Settings}->{CNCP}->{absc}, $config->{Settings}->{CNCP}->{prob}, $config->{Settings}->{CNCP}->{dist})\n";
		}
		if($n[$i] eq "fast") {
			$command .= "\trFASTCALL(data, file, idir, odir)\n";
		}
	}
	
	$command .= "\tnames= pn.signs(data, file, idir, odir)\n\tfor(y in 1:length(names) ) {\n";	
	$command .= "\t\tdata= read.table(names[y])\n";	
	for my $i (0 .. $#s) {
		
		if($s[$i] eq "GADA") {
			$command .= "\t\tGADA_s(data, names[y], idir, odir, $config->{Settings}->{GADA}->{alpha}, $config->{Settings}->{GADA}->{T}, $config->{Settings}->{GADA}->{mL}, $config->{Settings}->{GADA}->{mR})\n";
		}
		if($s[$i] eq "ADM2") {
			$command .= "\t\tADM2(data, names[y], odir, $config->{Settings}->{ADM2}->{t}, $config->{Settings}->{ADM2}->{np}, $config->{Settings}->{ADM2}->{mr}, $config->{Settings}->{ADM2}->{fact})\n";
		}
		if($s[$i] eq "ADM3") {
			$command .= "\t\tADM3(data, file, odir, $config->{Settings}->{ADM3}->{t}, $config->{Settings}->{ADM3}->{nd}, $config->{Settings}->{ADM3}->{mr}, $config->{Settings}->{ADM3}->{np})\n";
		}
		if($s[$i] eq "SMUG") {
			$command .= "\t\tSMUG(data, names[y], idir, odir, $config->{Settings}->{SMUG}->{resolution}, $config->{Settings}->{SMUG}->{probC}, $config->{Settings}->{SMUG}->{p1}, $config->{Settings}->{SMUG}->{p2}, $config->{Settings}->{SMUG}->{mL}, $config->{Settings}->{SMUG}->{mR}, $config->{Settings}->{SMUG}->{pT}," . "'" . $config->{Process}->{gold} . "')" . " " . "\n";
		}
		if($s[$i] eq "SMAP") {
			$command .= "\t\tSMAP(data, names[y], idir, odir, $config->{Settings}->{SMAP}->{resolution}, $config->{Settings}->{SMAP}->{probC}, $config->{Settings}->{SMAP}->{p1}, $config->{Settings}->{SMAP}->{p2}, $config->{Settings}->{SMAP}->{mL}, $config->{Settings}->{SMAP}->{mR}, $config->{Settings}->{SMAP}->{pT}," . "'" . $config->{Process}->{gold} . "')" . " " . "\n";
		}
		if($s[$i] eq "fsegment") {
			$command .= "\t\tFsegment(data, names[y], idir, odir, $config->{Settings}->{fsegment}->{sL}, $config->{Settings}->{fsegment}->{cpP}, $config->{Settings}->{fsegment}->{mL}, $config->{Settings}->{fsegment}->{mR})\n";
		}
		if($s[$i] eq "dnacopy") {
			$command .= "\t\tdna.copy(data, names[y], idir, odir, $config->{Settings}->{dnacopy}->{mL}, $config->{Settings}->{dnacopy}->{mR})\n";
		}
		if($s[$i] eq "safe") {
			$command .= "\t\tsafe.dec(names[y], idir, odir, $config->{Settings}->{safe}->{mL}, $config->{Settings}->{safe}->{mR})\n";
		}
		if($s[$i] eq "walk") {
			$command .= "\t\tjv.walk2(names[y], idir, odir, $config->{Settings}->{SMUG}->{resolution}, $config->{Settings}->{SMUG}->{probC}, $config->{Settings}->{SMUG}->{p1}, $config->{Settings}->{SMUG}->{p2}, $config->{Settings}->{SMUG}->{mL}, $config->{Settings}->{SMUG}->{mR}, $config->{Settings}->{SMUG}->{pT}," . "'" . $config->{Process}->{gold} . "')" . " " . "\n";
		}
		if($s[$i] eq "CNCP") {
			$command .= "\t\trcncp(data, names[y], idir, odir, $config->{Settings}->{CNCP}->{fac}, $config->{Settings}->{CNCP}->{a}, $config->{Settings}->{CNCP}->{n}, $config->{Settings}->{CNCP}->{absc}, $config->{Settings}->{CNCP}->{prob}, $config->{Settings}->{CNCP}->{dist})\n";
		}
		if($s[$i] eq "fast") {
			$command .= "\t\trFASTCALL(data, names[y], idir, odir)\n";
		}
	}	
	$command .="\t}\n";
	
	return($command);
}

sub writeScript {
	my $config = shift;
	my $command = shift;
	if(!$config->{Process}->{lsf}) { $command .= "\n}\n"; }
	open OUT, ">" . $config->{Process}->{runnerloc} || die;
	print OUT $command;
	close OUT;	
}

