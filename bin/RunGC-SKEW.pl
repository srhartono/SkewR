#!/usr/bin/perl -w
use strict; use warnings;

BEGIN {
        # Get current directory and push into Perl library
        my $dir;
        $dir = `pwd` . "/bin/";
        $dir =~ s/\n//;
        $dir =~ s/bin\/bin/bin/;
        push (@INC, $dir);

	# Check if SkewR.pm is succesfully pushed
	my $libSkewRCheck = 0;
	for (my $i = 0; $i < @INC; $i++) {
		$libSkewRCheck = 1 if -e "$INC[$i]\/SkewR.pm";
	}

	# Otherwise, ask user to put SkewR.pm to their perl library themselves
	if ($libSkewRCheck == 0) {
		print "Failed to put SkewR.pm into Perl library\nPlease manually copy SkewR.pm yourself to one of these Perl library folder:\n";
		for (my $i = 1; $i <= @INC; $i++) {
			print "$i. $INC[$i]\n";
		}
		die "\n";
	}
}

use Cwd;
use Thread;
use Thread::Queue;
use Getopt::Std;
use SkewR;
use vars qw($opt_s $opt_m $opt_t $opt_v $opt_h $opt_n $opt_z $opt_l $opt_g $opt_b);
getopts("s:m:t:vhn:z:l:g:b:");

my $version = "\nVERSION: 11/1/2013\n\n";
die $version if ($opt_v);

my ($seqFile, $modelFile, $threshold, $threads, $projName, $length, $geneFile, $cpgFile) = StochHMMToBed::check_sanity();

print "Importing Fasta\n";
my @splitFastaName = 0;#@{import_fasta($seqFile, $projName)};
print "Running StochHMM\n";
stochHMM_decode(\@splitFastaName, $modelFile, $threads, $threshold);
print "Deleting Fasta\n";
#delete_fastas(@splitFastaName);
print "Converting Report\n";
convert_report($projName, $length, $geneFile, $cpgFile, $threshold, @splitFastaName);
print "Intersecting Skew Peaks with Gene Info\n";
Intersect::main($projName, \@splitFastaName, $geneFile, $cpgFile);
print "Deleting Reports\n";
#delete_report($projName);

######################################################
#	Run GC-Skew model on whole genome
######################################################
sub get_filename {
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

sub import_fasta {
	my ($fasta, $projName) = @_;
	mkdir $projName if (not -d $projName);

	my $sequence;
	my @splitFastaName;

	my $begin = 1;

	# Split fasta by their chromosomes
	my ($chr, $seq);
	open (my $in, "<", $fasta) or die "Cannot read from $fasta: $!\n";
	while (my $line = <$in>) {
		chomp($line);

		# If current line is fasta header, print chr and fasta
		# Unless it's initial (not defined $chr)
		if ($line =~ /^>/) {
			if (defined($chr)) {
				my $splitFastaName = output_fasta($chr, $seq, $projName);
				push(@splitFastaName, $splitFastaName);
			}
			my @chr = split(" ", $line);
			($chr) = $chr[0] =~ /^>(.+)$/;
			print "\tImporting $chr\n";
			undef($seq);
		}
		else {
			$seq .= uc($line);
		}
	}
	close $in;
	my $splitFastaName = output_fasta($chr, $seq, $projName);
	push(@splitFastaName, $splitFastaName);

	return (\@splitFastaName);
}

sub output_fasta {
	my ($chr, $seq, $projName) = @_;

	my $outputFile = "$projName/$chr";
	open (my $out, ">", "$outputFile.fa") or die "Cannot write to $outputFile.fa: $!\n";
	print $out ">$chr\n$seq\n";
	close $out;

	return ($outputFile);
}

sub stochHMM_decode {
	my ($splitFastaName, $modelName, $threads, $threshold) = @_;
	
	my %arg;
	$arg{SEQ}   	= $splitFastaName;
	$arg{MODEL} 	= $modelName;
	$arg{THREADS}   = $threads;
	$arg{THRESHOLD} = $threshold;
	

	foreach my $arg (sort keys %arg) {
		if ($arg eq "SEQ") {
			for (my $i = 0; $i < @{$arg{$arg}}; $i++) {
				print "SEQ $i\t$arg{$arg}[$i]\n";
			}
		}
		else {
			print "$arg\t$arg{$arg}\n";
		}
	}
	run_stochhmm(\%arg);

	return;
}

sub run_stochhmm{

	my %arg = %{$_[0]};
	$arg{THREADS} = $arg{THREADS} >= @{$arg{SEQ}} ? @{$arg{SEQ}} : $arg{THREADS};

	my $command = "stochhmm -seq FILENAME.fa -model $arg{MODEL} -posterior -threshold $arg{THRESHOLD} > FILENAME.prob";
	my $Q = new Thread::Queue;
	foreach my $seq (@{$arg{SEQ}}) {
		my ($comm) = $command;
		$comm =~ s/FILENAME/$seq/g;
		print "$comm\n";
		$Q->enqueue($comm);
	}
	$Q->end();
	my @threads;

	for (my $i = 0; $i < $arg{THREADS}; $i++){
		$threads[$i] = threads->create(\&worker, $i, $Q);
	}
		
	for (my $i = 0; $i < $arg{THREADS}; $i++){
		$threads[$i]->join();
	}
	
	return;
}

#worker subroutine for run_stochhmm
sub worker {
	my ($thread, $queue) = @_;
	my $tid = threads->tid;
	
	while ($queue->pending) {
		my $command = $queue->dequeue;
		print "processing file $command in thread $tid\n";
		`$command`;
	}
    
    return
}

sub delete_fastas {
	my @splitFastaName = @_;

	foreach my $splitFastaName (@splitFastaName) {
		unlink "$splitFastaName.fa";
	}

	return;
}

#################################################
#	Convert Reports to BED files		#
#################################################

sub convert_report {
	my ($projName, $length, $geneFile, $cpgFile, $threshold, @splitFastaName) = @_;
	my @states = ("C_SKEW", "G_SKEW");
	my %colors = (C_SKEW => "0,0,255", G_SKEW => "255,0,0");
	
	for (my $i = 0; $i < @splitFastaName; $i++) {
		StochHMMToBed::main("$splitFastaName[$i].prob", \%colors, $length, \@states, $projName, $threshold);
	}
}

sub delete_report {
	system('rm -r -f *.prob');
	return;
}

1;
