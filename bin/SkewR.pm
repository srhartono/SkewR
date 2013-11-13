package StochHMMToBed;

#################################################################
# StochHMMToBed							#
# This script parse SkewR posterior output file (from STDOUT)	#
# into a concatenated bed file. Bed file coordinate 		#
# will be zero-based [start,end] of StochHMM output.		#
# (Use perldoc for more explanation)				#
#								#
# By Stella Hartono, UC Davis, 2012				#
# This script is provided "as is" without warranty of any kind.	#
#################################################################

use strict; use warnings; 

$SIG{INT} = \&interrupt;

sub main {

	# Check for input sanity and get states from stochhmm output
	my ($input, $colors, $length, $states, $projName, $threshold) = @_;
	my @state = @{$states};
	my %color = %{$colors};
	
	# Open input and output, and print track header into output file
	open (my $in , "<", $input)  or die "Cannot read from $input: $!\n";
	my ($output_name) = getFilename($input);
	open (my $out, ">", "$projName/$output_name.bed") or die "Cannot write to $projName/$output_name: $!\n";
	#print $out "track name=\"$output_name\" description=\"$output_name\" itemRgb=\"On\"\nvisibility=1\n";
	
	# Set up global variables which temporarily store "current" values when processing input file
	my ($current_chr, $current_state);
	my ($start_coor, $previous_coor) = ("-1", "-1"); 
	
	my %data;	# hash that store values for each coordinate
	my @state_pos;	# array to store state position
	
	my $linecount   	= 0; 	   # Variable to store line count
	my $store_data_check	= 0; 	   # Boolean to store data or not (0 is don't store)
	
	# Process the StochHMM posterior probability file
	while (my $line = <$in>) {
		chomp($line);
		$linecount++; # Count line number
		
		# Set up global variables to process the desired StochHMM posterior probability file
		# Line with "Posterior Probabilities Table" means it's the end of the previous data
		# so we stop storing data
		if ($line =~ /Posterior Probabilities Table/) {
			$store_data_check = 0;
		}
	
		# Line with "Sequence" contain new chromosome info.
		elsif ($line =~ /Sequence.+>/) {
		
			# If current state is not initial state (undefined) then we print out the previously stored data
			# and erase the previous data so the hash won't get too big.
			# If we at the end of file, this command won't be invoked. Therefore we will print the
			# last stored data at the end of this loop.
			if (defined($current_state)) {
				# Print each state's peaks coordinate which length is above threshold length
				# Then erase data so the hash won't get too big.
				print_data(\%data, \%color, $length, $current_chr, $out);
				%data = (); 
			}
	
			# Reset global variables to be used in processing the data of the new chromosome.
			undef($current_state);
			$start_coor  = "-1";
		
			# Get the new chromosome stored in global variable
			($current_chr)   = $line =~ /Sequence.+>(\w+)/;
			print "Processing chromosome $current_chr\n";
			die "Fatal error: Undefined chromosome at line $linecount: $line\n" if not defined($current_chr);
		}
	
		# Line with "Position" contains StochHMM state position information we need to store.
		# The subsequent line below "Position" contain data of coordinate and each state's posterior probability
		# Posterior probability below user-inputted StochHMM threshold (not this script's threshold) will be blank
		elsif ($line =~ /^Position/) {
			my ($name, @state_names) = split("\t", $line); # Name is "Position", which we don't use.
	
			# Get the position of states stored in global array state_pos
			# (position means where each user-inputted state is located in the StochHMM output file)
			for (my $i = 0; $i < @state_names; $i++) {
				for (my $j = 0; $j < @state; $j++) {
					$state_pos[$i] = $state[$j] if $state_names[$i] eq $state[$j];
				}
			}
			$store_data_check = 1; # Start storing data at next line
		}
	
		# If boolean store_data_check is 1, we start storing adjacent coordinates and 
		# as a peak if the probability is above this script's user-inputted threshold
		# Btw, posterior probability below user-inputted StochHMM threshold 
		# (not this script's threshold) will be shown blank by StochHMM
		elsif ($store_data_check == 1 and $line =~ /^\d+/) {
	
			die "Fatal error: Undefined chrosome at line $linecount: $line\n" if not defined($current_chr);
			my ($current_coor, @vals) = split("\t", $line);
	
			if (defined($current_coor) and $current_coor =~ /^\d+$/) {
				
				my ($value, $state); # Variable to temporarily store value of the state.
	
				# We assume StochHMM threshold is above 0.5, therefore only 1 state will have a numeric probability.
				for (my $i = 0; $i < @vals; $i++) {
					next if $vals[$i] !~ /^\d+\.?\d*$/; # Skip the state if value is blank.
					$state = $state_pos[$i];
					$value = $vals[$i];
				}
				die "Fatal error: Undefined state at line $linecount: chrom $current_chr $line\n" unless defined($state);
	
				next if $value < $threshold;
				# If current state is initial (not defined), or state changes, or previous coordinate is not adjacent
				# of current coordinate, then we store the temporarily-recorded peak
				# and record the current state, start coordinate of that state, and end coordinate of that state
				if (not defined($current_state) or $current_state ne $state or $current_coor -1 != $previous_coor) {
					delete($data{$start_coor}) if defined($current_state) and $data{$start_coor}{end} - $start_coor + 1 < $length;
					$current_state  = $state;
					$start_coor     = $current_coor;
					$previous_coor  = $current_coor;
					
					# Explanation about how the hash below is stored:
					# A peak at state $current_state, chromosome $current_chr, 
					# starting at coordinate $start_coor, currently ends at $current_coor
					$data{$start_coor}{state} = $current_state;
					$data{$start_coor}{end}   = $current_coor;
				}
	
				# If the current coordinate is adjacent to previous coordinate, then we define
				# new end coordinate for that peak as $current_coor 
				elsif ($previous_coor == $current_coor - 1) {
					$data{$start_coor}{end} = $current_coor;
					$previous_coor		= $current_coor;
				}
				else {
					print "Fatal error: Unknown reason at $linecount: $line\nPlease report this bug to srhartono\@ucdavis.edu";
					$store_data_check = 0;
				}
			}
			else {
				print "Fatal error: Unknown reason at $linecount: $line\nPlease report this bug to srhartono\@ucdavis.edu";
				$store_data_check = 0;
			}
		}
	}
	close $in;
	
	# Print out peaks stored in last data hash
	print_data(\%data, \%color, $length, $current_chr, $out);
	return("$projName/$output_name.bed");
}

sub check_sanity {

        # Check if stochhmm is installed
	my $stochHMMCheck = `which stochhmm`;
	$stochHMMCheck = `which StochHMM` if $stochHMMCheck =~ /^$/;
	die "Please install StochHMM 0.35 or newer (https://github.com/KorfLab/StochHMM) and put it in your \$PATH directoy (e.g. /usr/local/bin or /usr/bin/)\n" if $stochHMMCheck =~ /^$/;
	
	# Check if bedtools is installed
	my $bedtoolsCheck = `which bedtools`;
	die "Please install bedtools version 2 or newer (https://code.google.com/p/bedtools/) before running and put it in your \$PATH directoy (e.g. /usr/local/bin or /usr/bin/)\n" if $bedtoolsCheck =~ /^$/;
	
	# Check bedtools version, and warn user to install correct version if it's not defined or less than 2
	my $bedtoolsVer = `bedtools --version`;
	warn "Please install bedtools version 2 or newer (https://code.google.com/p/bedtools/) before running and put it in your \$PATH directoy (e.g. /usr/local/bin or /usr/bin/)\n" if $bedtoolsVer =~ /^$/;
	if (defined($bedtoolsVer)) {
		($bedtoolsVer) = $bedtoolsVer =~ /bedtools v(\d+)/;
		warn "Please install bedtools version 2 or newer (https://code.google.com/p/bedtools/) before running and put it in your \$PATH directoy (e.g. /usr/local/bin or /usr/bin/)\n" if ($bedtoolsVer) < 2;
	}

	# Help
	print_usage() and die "\n"	if defined($main::opt_h);

	# Help if nothing is input
	print_usage() and die "Fatal Error: Missing input file\n\n" unless defined($main::opt_s) and defined($main::opt_m);
	

	# Variables
	my $seqFile   = $main::opt_s;
	my $modelFile = $main::opt_m;
	my $threshold = defined($main::opt_t) ? $main::opt_t : 0.95;
	my $projName  = defined($main::opt_o) ? $main::opt_o : (print_usage() and die "Fatal Error: Missing Output Directory (-o)\n\n");
	my $threads   = defined($main::opt_z) ? $main::opt_z : 1;

	my $length    = defined($main::opt_l) ? $main::opt_l : 300;
	my $geneFile  = defined($main::opt_g) ? $main::opt_g : (print_usage() and die "Fatal Error: Missing gene file (-g)\n\n");
	my $cpgFile   = defined($main::opt_b) ? $main::opt_b : (print_usage() and die "Fatal Error: Missing CpG file (-b)\n\n");

	# Offset
	die "Please intersect only with either TSS or TTS\n" if defined($main::opt_x) and defined($main::opt_y);
	my $offset = defined($main::opt_x) ? $main::opt_x : defined($main::opt_y) ? $main::opt_y : "0,0";
	if ($offset ne 1) {
		die "-x or -y must be two comma separated integers (e.g. -x -200,500): $offset\n" if $offset !~ /^\d+\,\d+$/;
	}

	# Sequence
	open (my $seqIn, "<", $seqFile) or die "Cannot read from $seqFile: $!\n";
	while (my $line = <$seqIn>) {
		die "Fatal Error: Sequence file does not seem to be Fasta format!\n" if ($line !~ /^>/);
		last;
	}
	close $seqIn;

	# Model
	open (my $modelIn, "<", $modelFile) or die "Cannot read from $modelFile: $!\n";
	while (my $line = <$modelIn>) {
		# Die if first line is not STOCHHMM MODEL FILE
		die "Fatal Error: Model file does not seem to be StochHMM format (did you use version 0.35 or latest?)!\n" if $line !~ /#STOCHHMM MODEL FILE/;
		last; 
	}
	close $modelIn;

	# Threshold
	print_usage() and die "Fatal Error: Threshold must be between number between (including) 0 and 1\n\n" if ($threshold !~ /^\d+\.?\d*$/ or $threshold > 1 or $threshold < 0);

	# Threads
	print_usage() and die "Fatal Error: Threads must be a positive integer\n\n" if ($threads !~ /^\d+$/ or $threads < 0);
	
	# Output directory
	if (not -d $projName) {
		system("mkdir $projName") == 0 or die "Failed to create directory for project name $projName: $!\n";
	}
	else {
		print "Warning: Directory exists, all files will be overwritten. Proceed? (Enter to proceed or Ctrl+C to cancel)";
		<STDIN>;
	}

	# CpG File
	open (my $cpgFileIn, "<", $cpgFile) or die "Cannot read from $cpgFile: $!\n\n";
	while (my $line = <$cpgFileIn>) {
		chomp($line);
		next if $line =~ /track/;
		next if $line =~ /\#/;
		my ($chr, $start, $end) = split("\t", $line);
		die "CpG file does not seem to be a BED 3 format\n" if (not defined($start) or $start !~ /^\d+$/);
		die "CpG file does not seem to be a BED 3 format\n" if (not defined($end) or $end !~ /^\d+$/);
		last;
	}
	close $cpgFileIn;

	# Gene file
	open (my $geneFileIn, "<", $geneFile) or die "Cannot read from $geneFile: $!\n\n";
	while (my $line = <$geneFileIn>) {
		chomp($line);
		next if $line =~ /track/;
		next if $line =~ /\#/;
		my ($chr, $start, $end, $name, $value, $strand) = split("\t", $line);
		die "Gene file does not seem to be a BED 6 format (start coor error)\n" if (not defined($start) or $start !~ /^\d+$/);
		die "Gene file does not seem to be a BED 6 format (end coor error)\n" if (not defined($end) or $end !~ /^\d+$/);
		die "Gene file does not seem to be a BED 6 format (strand error)\n" if (not defined($strand) or ($strand !~ /^\+$/ and $strand !~ /^\-$/));
		last;
	}
	close $geneFileIn;

	# Length
	print_usage() and die "Fatal Error: Length must be a positive integer\n\n" if ($length !~ /^\d+$/ or $length < 0);


	print "
Bedtools folder: $bedtoolsCheck

Input parameters:
Input file	: $seqFile
Model file	: $modelFile
Threshold	: $threshold
Threads		: $threads
Output Dir	: $projName
Min Length	: $length
Gene file	: $geneFile
CpG file	: $cpgFile

";
	return($seqFile, $modelFile, $threshold, $threads, $projName, $length, $geneFile, $cpgFile);
}

sub print_data {
	my ($data, $color, $length, $chr, $out) = @_;
	my %data  = %{$data};
	my %color = %{$color};
	foreach my $coor (sort {$a <=> $b} keys %data) {
		my $state = $data{$coor}{state};
		my $end   = $data{$coor}{end}  ;
		my $diff  = $end - $coor + 1;

		# Get color information for the state
		my $current_color = defined($color{$state}) ? $color{$state} : "0,0,0";
		next if $end - $coor + 1 < $length;
		print $out "$chr\t$coor\t$end\t$state\t$diff\t+\t$coor\t$end\t$current_color\n";
	}
}

sub interrupt {
	print "$0 cancelled\n";
	exit;

}
sub print_usage {
        print "
usage: $0 [options] -s <Sequence> -m <Model> -g <Gene File> -b <CpG file> -o <Output Directory>

Arguments:
-s: Sequence file (Fasta format)
-m: Model file (StochHMM HMM format)
-g: Gene file (BED 6+ format)
-b: CpG file (BED 3+ format)
-o: Output Directory (string)

Options:

1. StochHMM Parameters:
-t: Minimum posterior probability threshold (float [0-1]). Default: 0.9
    Probability threshold is using \"equal or more than\" format (only take probability >= -t)
-z: Thread number (integer)

2. Skew Classes Parameters:
-o: Project name, will be used as output directory (string)
-x: Intersect with N,N bp of Transcription Start Site (default: -500 +1500 of TSS)
    Example: -x -500,1500 is -500 +1500 of TSS
-y: Intersect with Transcription Termintation Site (default: -1500 +500 of TTS)
    Example: -y -1500,500
-l: Minimum length of a SkewR peak to be recorded (integer). Default: 300
Length is using \"equal or more than\" format (only take peak length >= -l)

\n";
}

sub getFilename {
	my ($fh, $type) = @_;
	my (@splitname) = split("\/", $fh);
	my $name = $splitname[@splitname-1];
	pop(@splitname);
	my $folder = join("\/", @splitname);
	@splitname = split(/\./, $name);
	$name = $splitname[0];
	return($name) if not defined($type);
	return($folder, $name) if $type eq "folder";
}


=head1 NAME

package SkewR

=head1 USAGE

usage: $0 [options] -s <Sequence> -m <Model> -g <Gene File> -b <CpG file> -o <Output Directory>

 Arguments:
 -s: Sequence file (Fasta format)
 -m: Model file (StochHMM HMM format)
 -g: Gene file (BED 6+ format)
 -b: CpG file (BED 3+ format)
 -o: Output Directory (string)


Options:

1. StochHMM Parameters:
 -t: Minimum posterior probability threshold (float [0-1]). Default: 0.9
 Probability threshold is using \"equal or more than\" format (only take probability >= -t)
 -z: Thread number (integer). Default: 1

2. Skew Classes Parameters:
 -x: Intersect with N,N bp of Transcription Start Site (default: -500 +1500 of TSS)
     Example: -x -500,1500 is -500 +1500 of TSS
 -y: Intersect with Transcription Termintation Site (default: -1500 +500 of TTS)
     Example: -y -1500,500
 -l: Minimum length of a SkewR peak to be recorded (integer). Default: 300
     Length is using \"equal or more than\" format (only take peak length >= -l)

=head1 SYNOPSIS

Running SkewR on human genome 19 hg19.fa:

Perl bin/RunGC-Skew.pl -s hg19.fa -m model/GC_SKEW_7600.hmm -g hg19_gene.bed -b hg19_cpg.bed -o MyResult

=head1 AUTHOR

Stella R. Hartono (srhartono@ucdavis.edu)

=head1 COPYRIGHT

Copyright 2012 Stella Hartono.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.
A copy of the license is included in the section entitled "GNU
Free Documentation License".

=head1 DISCLAIMER

This script is provided "as is" without warranty of any kind.

=cut

package Intersect;

sub main {

	my ($projName, $splitFastaName, $geneFile, $cpgFile) = @_;

	system("mkdir $projName/temp") if not -d "$projName/temp/";

	my @splitFastaName = @{$splitFastaName};
	my $probBedFile       = "$projName/$projName.bed";
	for (my $i = 0; $i < @splitFastaName; $i++) {
		system("cat $splitFastaName[$i].bed >> $probBedFile") == 0 or die "Failed at merging probability files\n";
	}

	if (-s $probBedFile == 0) {
		print "There are no SkewR peak found!\n";
		exit;
	}

	my $promoterfile = "$projName/temp/promoter.bed";
	my $TTSfile      = "$projName/temp/TTS.bed";
	my $maingeneFile = $promoterfile;
	my $promoterTag  = "-a";
	my ($offset_min, $offset_pos) = (-500, 1500);
	($offset_min, $offset_pos) = split(",", $main::opt_x) if defined($main::opt_x) and $main::opt_x ne 1;
	if (defined($main::opt_y)) {
		if ($main::opt_y eq 1) {
			($offset_min, $offset_pos) = (-1500, 500);
		}
		else {
			($offset_min, $offset_pos) = split(",", $main::opt_y);
		}
		$maingeneFile = $TTSfile;
		$promoterTag  = "-b";
	}

	my ($binDir) = $0 =~ /^(.+)RunGC-SKEW.pl$/;
	$binDir = "./" if not defined($binDir);
	system("$binDir\/bedtools_bed_change.pl $promoterTag -x $offset_min -y $offset_pos -o $maingeneFile -i $geneFile") == 0 or die "Failed to run Bedtools!\n";

	intersect($projName, $maingeneFile, $geneFile, $probBedFile, $cpgFile);
}

sub intersect_Old {
        my ($projName, $maingeneFile, $geneFile, $probBedFile, $cpgFile) = @_;

        print "
Project Name = $projName
Main Gene file = $geneFile
probBedFile = $probBedFile
cpg file = $cpgFile
";

        my $gene_posfile = "$projName/gene_pos.txt";
        my $gene_negfile = "$projName/gene_neg.txt";
        my $bed_GSKEW    = "$projName/bed_GSKEW.txt";
        my $bed_CSKEW    = "$projName/bed_CSKEW.txt";
        my $strongfile   = "$projName/strong.txt";
        my $reversefile  = "$projName/reverse.txt";
        my $weakfile     = "$projName/weak.txt";
        my $bidirectfile = "$projName/bidirect.txt";
        my $noskewfile   = "$projName/noskew.txt";
        my $tempfile     = "$projName/temp.txt";
        my $tempfile2     = "$projName/temp2.txt";

        system("grep '+' $geneFile > $gene_posfile");
        system("grep '-' $geneFile > $gene_negfile");

        system("grep 'G_SKEW' $probBedFile > $bed_GSKEW");
        system("grep 'C_SKEW' $probBedFile > $bed_CSKEW");

        # Strong Skew Class     
        system("bedtools intersect -u -a $gene_posfile -b $bed_GSKEW >  $strongfile");
        system("bedtools intersect -u -a $gene_negfile -b $bed_CSKEW >> $strongfile");
        print "There are no Strong Skew Genek found!\n" if (-s $strongfile == 0);

        # Reverse Skew Class
        system("bedtools intersect -u -a $gene_negfile -b $bed_GSKEW >  $reversefile");
        system("bedtools intersect -u -a $gene_posfile -b $bed_CSKEW >> $reversefile");
        print "There are no Reverse Skew Gene found!\n" if (-s $reversefile == 0);

        # Bidirectional Skew Class
        # Overlap Strong and Reverse Skew Genes to get Bidirectional Skew Genes
        # Which is then kept in Strong Skew Genes and removed from Reverse Skew Genes
        if (-s $strongfile != 0 and -s $reversefile != 0) {
                system("bedtools intersect -u -a $strongfile -b $reversefile > $bidirectfile");
                print "There are no Bidirectional Skew Gene found!\n" if -s $bidirectfile == 0;
                if (-s $bidirectfile != 0) {
                        system("bedtools intersect -v -a $reversefile -b $bidirectfile > $tempfile && mv $tempfile $reversefile");
                }
        }

        system("cat $strongfile $reversefile > $tempfile");

        # No Skew Class
        system("bedtools intersect -v -a $geneFile -b $tempfile > $noskewfile");
        print "There are no No Skew Gene found!\n" if (-s $noskewfile == 0);

        # Weak Skew Class
        if (-s $noskewfile != 0) {
                system("bedtools intersect -u -a $noskewfile -b $cpgFile > $weakfile");
                print "There are no Weak Skew Gene found!\n" if (-s $weakfile == 0);
                if (-s $weakfile != 0) {
                        system("bedtools intersect -v -a $noskewfile -b $weakfile > $tempfile");
                        system("mv $tempfile $noskewfile");
                        print "There are no No Skew Gene found!\n" if (-s $noskewfile == 0);
                }
        }
}

1;
