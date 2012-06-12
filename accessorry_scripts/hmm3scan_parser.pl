#!/usr/bin/perl
# file: pfam_parser.pl

use Getopt::Std; 	# read in data at command line using -i type options

$progname=`basename $0`; chomp($progname);
$info.="$progname\nScripted:\tBill Wickstead, 2010\n";
$info.="\tParses domain position and hit quality information from hmmscan search - update for Hmmer3.0";
$usage="Usage:\t\t$progname -i \<infile\> (hmmscan --noali .out file) -o \<outfile\>\n";
$usage.="Options:\t-m \<max number of domains\>\t[default=all]\n";
$usage.="\t\t-s\tmin score [default=-100]\n";
$usage.="\t\t-e\tmax evalue [$default=100]\n";
$usage.="\t\t-G\tdon't squash overlapping domains\n";
$usage.="\t\t-E\tuse envelop estimates of domain position (rather than alignment position)\n";
$usage.="\t\t-l <listfile> ONLY output information on ids from list [default: all]\n";
#$usage.="\t\t-S\tSort output by name\n";
#$usage.="\t\t-N\tNames only (no score and e-val)\n";

sub usage_error {
    print "\n$info\n$usage\n"; exit 1;
}
sub parameter_default {
# returns either $parameters{$_[0]} (if set) or $_[1]
	if (exists $parameters{$_[0]} ) { return $parameters{$_[0]}; }
	else { return $_[1]; }
}
## PERL MATHS SUBS ##
sub log10 { return log($_[0])/log(10); }

## SETUP ##
getopts('i:m:A:SNs:e:f:o:Gl:E',\%parameters);
unless (exists $parameters{'i'} && exists $parameters{'o'}) {  usage_error(); }
$maxhits=parameter_default( 'm', 10000);
$minscore=parameter_default( 's', -100);
$maxeval=parameter_default( 'e', 100);

## PROGRAM STARTS ##
open (IN, $parameters{'i'}) or die "Can't open file: $parameters{'i'}\n";
@data=<IN>;
close (IN);
chomp(@data);

if (exists $parameters{'l'}) {
	open (IN, $parameters{'l'}) or die "Can't open file: $parameters{'l'}\n";
	@list=<IN>;
	close (IN);
	chomp(@list);
	
	for $id (@list) {
		$listname{$id}=1;
	}
}

open (OUT, ">$parameters{'o'}");

for ($n=0; $n<=$#data; $n++) {
	$l=$data[$n];
	if ($l=~m/^Query:\s+(.+)\s+\[L=(\d+)\]$/) {
		$query=$1;
		$len=$2;
		$query=~s/\s//g;
		$length{$query}=$len;
	} elsif ($l=~m/Domain annotation for each model:/) {
		$n++;
		if ($data[$n] eq "") {
			$continue=0;
		} else {
			$continue=1;
		}
		
		if (exists $parameters{'l'} && !exists $listname{$query}) {
			$continue=0;
		}
		
		$m=0; $Nblank=0;
		%model=(); %score=(); %eval=(); %start=(); %stop=(); %key=();
		while ($continue) {
			if ($data[$n] eq "") {
				$Nblank++;
				if ($Nblank>1) {
					$continue=0;
					last;
				}
			} elsif ($data[$n]=~m/^\>\>\s(.+)\s\s(.+)$/) {
				$current=$1;
				$n+=2;
				$Nblank=0;
			} else {
				@tmp=split(/\s+/, $data[$n]);
				if ($tmp[3]>=$minscore && $tmp[6]<=$maxeval) {
					$score{$m}=$tmp[3];
					$eval{$m}=$tmp[6];
					if (exists $parameters{'E'}) {
						$start{$m}=$tmp[13];
						$stop{$m}=$tmp[14];
						$key{$m}=$tmp[15];
					} else {
						$start{$m}=$tmp[10];
						$stop{$m}=$tmp[11];
						$key{$m}=$tmp[12];					
					}

					$model{$m}=$current;
					$m++;
				}
			}
			$n++;
		}
		
		## SORT DOMAINS
		$m=0;
		@sorted_model=(); @sorted_start=(); @sorted_stop=(); @sorted_key=();
		if (scalar(keys %model)>0) {
			foreach $k (sort {$start{$a} <=> $start{$b}} keys %start) {
				$sorted_model[$m]=$model{$k};
				$sorted_start[$m]=$start{$k};
				$sorted_stop[$m]=$stop{$k};
				$sorted_key[$m]=$key{$k};

				if ($m>0 && !exists $parameters{'G'}) {
					if ($sorted_stop[$m-1]>$sorted_start[$m]) { 		## overlap!!
						if ($sorted_stop[$m]-$sorted_start[$m] > $sorted_stop[$m-1]-$sorted_start[$m-1]) {		## keep new model
							$sorted_model[$m-1]=$sorted_model[$m];
							$sorted_start[$m-1]=$sorted_start[$m];
							$sorted_stop[$m-1]=$sorted_stop[$m];
							$sorted_key[$m-1]=$sorted_key[$m];						
						} 		# else allow overwriting of new
					} else { 
						$m++; 
					}
				} else {
					$m++;
				}
			}
		}
		##OUTPUT##
		unless (exists $parameters{'l'} && !exists $listname{$query}) {
			print OUT "$query\t$length{$query}\t";
			if ($#sorted_model>=0) { print OUT $sorted_model[0]; }
			for (my $m=1; $m<=$#sorted_model; $m++) { print OUT "~$sorted_model[$m]"; }
			for (my $m=0; $m<=$#sorted_model; $m++) {
				print OUT "\t$sorted_start[$m]\t$sorted_stop[$m]\t$sorted_key[$m]";
			}
			print OUT "\n";
		}
	}
}

close OUT;
exit;

