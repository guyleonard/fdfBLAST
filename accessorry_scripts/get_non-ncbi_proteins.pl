#!/usr/bin/perl
use Cwd;
use File::Basename;

# Assumes fdfBLAST data-dir structure...

$working_directory = getcwd;
#$bin_directory    = "$working_directory/blast/bin";
####
$genome_set = "11fungi";
$genome_directory = "/home/cs02gl/Desktop/genomes/$genome_set";
print "Genome Directory = $genome_directory\n";
#
$run_ID = "11fungi_1e-10";
$run_dir = "$working_directory";
print "Run Directory = $run_dir\n";
####

&get_genomes;
&get_sequences;

sub get_genomes {
	print "Reading Genomes:\n";
	@genome_names = <$genome_directory/*.fas>;
	foreach $file (@genome_names) {
		( $file, $dir, $ext ) = fileparse( $file, '\..*' );
		#$file = fileparse($file);
		print "\t$file\n";
	}
	return @genome_names;
}


## Use fastacmd to extract sequences from formatdb -o F FASTA databases...
sub get_sequences {
	
	$compo = "all_11fungi_list_ordered_uniq.csv";
	#$split = "plant_split_list_o_uniq.csv";
	
	for ($i = 0; $i <= $#genome_names; $i++ ) {
	#for ($i = 0; $i <= 2; $i++ ) {
	
		$results =  `fastacmd -d $genome_directory\/$genome_names[$i]\.fas -i $run_dir\/$compo > $run_dir\/$genome_names[$i]_fusions\.fas`;

		#$results2 =  `$bin_directory/fastacmd -d $genome_directory\/$genome_names[$i]\.fas -i $run_dir\/$split > $run_dir\/$genome_names[$i]_split_seq\.fas`;
	
	}


}
