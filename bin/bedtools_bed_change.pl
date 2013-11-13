#!/usr/bin/perl

use strict; use warnings; use Getopt::Std;

use vars qw($opt_x $opt_y $opt_o $opt_a $opt_b $opt_i $opt_s $opt_f);
getopts("i:x:o:y:absf:");

die "
Usage: $0 [option] -i bed file

options:
-a: Get start of gene (strand specific) and offset accordingly
-b: Get end of gene (strand specific) and offset accordingly
-s: Disable strand specific (default: on)
-x: start offset from current (default: 0)
-y: end offset from current (default: 0)
-o: output (name)
-f: Filter genes less than ths length (default: take all)

E.g. you want all bed file
chr1	5000	6000	name	0	-
to become
chr1	4000	10000	name	0	-
then:

$0 -i foo.bed -x -1000 -y 4000 -o bar.bed 

E.g. +/- 1kb region of TSS (strand specific)
$0 -i foo.bed -x -1000 -y 1000 -o bar.bed

" unless defined($opt_i);

my $input  = $opt_i;
my $x_off  = defined($opt_x) ? $opt_x : 0;
my $y_off  = defined($opt_y) ? $opt_y : 0;
my $output = defined($opt_o) ? $opt_o : "$input.bed";
my $filter = defined($opt_f) ? $opt_f : 0;
die "Filter must be positive integer!\n" unless $filter =~ /^\d+$/;

print "Output = $output\n";
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", $output) or die "Cannot write to $output: $!\n";
while (my $line = <$in>) {
	chomp($line);
	my ($chr, $start, $end, $name, $val, $strand, @others) = split("\t", $line);
	next if $end - $start + 1 < $filter;

	if ($strand ne "-" and $strand ne "+") {
		if (defined($others[0]) and ($others[0] eq "+" or $others[0] eq "-")) {
			my $temp = $strand;
			$strand = $others[0];
			$others[0] = $temp;
		}
		else {
			die "Strand information is incorrect (strand = $strand) at line $line\n";
		}

	}
	my $others = join("\t", @others) if defined($others[0]);
	my ($newstart, $newend);
	if (not defined($opt_a) and not defined($opt_b)) {
		if (not defined($opt_s)) {
			if ($strand eq "+" or $strand eq "1" or $strand eq "F") {
				$newstart = $start + $x_off;
				$newend   = $end   + $y_off;
			}
			elsif ($strand eq "-" or $strand eq "-1" or $strand eq "R") {
				$newstart = $start - $y_off;
				$newend   = $end   - $x_off;
			}
		}
		else {
			$newstart = $start + $x_off;
			$newend   = $end   + $y_off;
		}
	}
	elsif (defined($opt_a)) {
		my $pos = $strand eq "+" ? $start : $end;
		$newstart = $strand eq "+" ? $pos + $x_off : $pos - $y_off;
		$newend   = $strand eq "+" ? $pos + $y_off : $pos - $x_off;
	}
	elsif (defined($opt_b)) {
		my $pos = $strand eq "+" ? $end : $start;
		$newstart = $strand eq "+" ? $pos + $x_off : $pos - $y_off;
		$newend   = $strand eq "+" ? $pos + $y_off : $pos - $x_off;
	}

	if ($newstart < 1 or $newend < 1) {
		print "Skipped [$chr $start $end] because start and/or end is less than 1\n";
		next;
	}
	print $out "$chr\t$newstart\t$newend\t$name\t$val\t$strand\t$others\n" if defined($others);
	print $out "$chr\t$newstart\t$newend\t$name\t$val\t$strand\n" if not defined($others);
}
close $in;
close $out;
