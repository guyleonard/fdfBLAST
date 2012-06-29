#!/usr/bin/perl
###########################################
# fdfBLAST June 2012                      #
# Full public release			  #
# Internal version based on 1.9.9.2	  #
$VERSION = "1.2012.beta";    #
###########################################

# Import Modules
use Bio::SearchIO;           # Bioperl for input/output/persing of BLAST etc
use Cwd;                     # Gets pathname of current working directory
use File::Basename;          # Remove path information and extract 8.3 filenames
use GD;                      # Creates PNG images

# use GD::SVG;	       # Creates SVG images - uncomment if you want SVG output
use Math::BigFloat
  ;    # Arbitrary size floating point math package (handles e-values)
use Switch
  ;    # A switch statement for Perl - deprecated, consider switch to given/when
use Time::Local;    # For time elapsed when running different stages

# Global Directory Information
our $WORKING_DIR = getcwd;
our $BLAST_BIN_DIR =
  "/usr/bin";       # Normal install is in /usr/bin/ - change to your dir...

# Other Variables
our $EMPTY = q{};

# Run these subroutines first...
#&run_blast_check;                            # Check for BLAST - You really should just have this installed.
our $GENOME_DIR = &set_genome_dir;    # Set where to look for genome directories
my $core_num, $cores =
  &detect_multi_core;    # Check for single or multi-core machine for BLAST
&which_genomes;          # Check which genome folder to use!
&initial_menu;           # User input for run number or next
&run_ID;                 # Sets run number if none set, sets directory paths
&menu;                   # Main menu

# Set where the genome directories can be found.
sub set_genome_dir {

    ######
    ## USER: You may add locations to this array...
    my @genome_directories = ("$WORKING_DIR/genomes");
    ######
    my $genome_directory = $EMPTY;
    print "Genomes Directory Menu\n**********************\n";
    for ( my $i = 0 ; $i <= $#genome_directories ; $i++ ) {
        print "$i) $genome_directories[$i]\n";
    }
    print "O) Other?\nChoose Genome Directory?\n>:";
    chomp( my $menu_choice = <STDIN> );
    if ( $menu_choice =~ m/O/is ) {
        print "Please enter the location of your genome directory\n>:";
        chomp( $genome_directory = <STDIN> );
        if ( -e $genome_directory && -d $genome_directory ) {

            # Do nothing
        }
        else {
            &set_genome_dir;
        }
    }
    else {
        $genome_directory = "$genome_directories[$menu_choice]";
    }
    print "\nYou chose $genome_directory\n";
    return $genome_directory;
}

# Attempt to detect status of CPU - useful for threads command in BLAST...
sub detect_multi_core {

    my $os_name = $^O;
    my $pprocs  = 1;     # Default
    my $lprocs  = 1;     # Default

    if ( $os_name eq 'linux' ) {

       # Only tested with Linux (Ubuntu), we can fall back to user entry though.

        $pprocs =
          `grep -i "physical id" /proc/cpuinfo | sort -u | wc -l | tr -d '\n'`;
        $lprocs =
          `grep -i "processor" /proc/cpuinfo | sort -u | wc -l | tr -d '\n'`;

    }
    elsif ( $os_name eq 'darwin' ) {

        # Only tested with OSX 10.4+ - will not work with < OSX 10.2
        $pprocs =
`system_profiler | grep -i "Number Of Processors:" | tr -d '[aA-zZ]:[[:punct:]]:\n'`;
        $lprocs =
`system_profiler | grep -i "Total Number Of Cores:" | tr -d '[aA-zZ]:[[:punct:]]:\n'`;

    }
    elsif ( $os_name eq 'MSWin32' ) {

        print "MS Windows is not currently supported.\n";

    }
    else {

        print "Please enter how many physical processors you have\n>:";

        chomp( $pprocs = <STDIN> );
        print
"Please enter how many logical processors (dual core = 2 etc) you have\n>:";
        chomp( $lprocs = <STDIN> );
    }

    if ( $lprocs >= 2 ) {

        $core_num = "$lprocs";
        $cores    = "multi";
    }
    else {
        $core_num = "1";
        $cores    = "single";
    }

    return $core_num, $cores;
}

sub run_blast_check {

    # Very simple check - only if BLAST directory exists.
    # Ideally we would run blast and check version.
    if ( -e $BLAST_BIN_DIR && -d $BLAST_BIN_DIR ) {
        print "Blast Detected\n";
        print `clear`, "\n";
    }
    else {
        print
"\n****\n\nYou do not appear to have LEGACY BLAST installed or it is installed in the wrong location.\nPlease download the LEGACY version from http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download and extract the contents to ../fdf/blast/ or see the Readme file.\n\n****\n";
        print "\nContinue without BLAST? (y/N)\n";
        print ">:";
        chomp( $blast_choice = <STDIN> );

        switch ($blast_choice) {

            case /Y/i {
                print "Blast Detection Overide\n";
                print `clear`, "\n";
                $BLAST_BIN_DIR = "/usr/bin/";
            }
            case /N/i {
                print `clear`, "\n";
                print
"Thanks for using fdfBLAST!\nPlease cite: \"Leonard, G. & Richards, T.A. 2012. Patterns of gene fusion and fission across the fungi. In Preparation.\"\n";
                exit;
            }
            else {
                &run_blast_check;
            }
        }
    }
}

# Retreive list of genome groups from your folder
# Assumes /your_dir/genomes/genome_set_a /your_dir_genome/genome_set_b etc...
sub which_genomes {

    my @all_files   = glob "$GENOME_DIR/*";
    my @folders     = $EMPTY;
    my $menu_choice = $EMPTY;

    foreach my $item (@all_files) {
        if ( -d $item ) {    #Put all folders into array
            push( @folders, $item );
        }
    }

    print "Genome Set Menu\n";
    print "***************\n";
    for ( $i = 0 ; $i <= $#folders ; $i++ ) {
        print "$i) $folders[$i]\n";
    }

    print
"Please enter the number for the directory where your genomes are located.\n>:";

    chomp( $menu_choice = <STDIN> );
    if (   $menu_choice > $#folders
        || $menu_choice < '0'
        || $menu_choice == m/[aA-zZ]/ )
    {
        print `clear`, "\n";
        print "\nIncorrect Menu Choice!\n\n";
        &which_genomes;

    }
    else {
        $GENOME_DIR = "$folders[$menu_choice]";
    }
}

sub initial_menu {
    print `clear`, "\n";
    print "  __     _   __  ___  _       _    ___  _____ \n";
    print " / _| __| | / _|| _ )| |     /_\\  / __||_   _|\n";
    print "|  _|/ _` ||  _|| _ \\| |__  / _ \\ \\__ \\  | |  \n";
    print "|_|  \\__,_||_|  |___/|____|/_/ \\_\\|___/  |_|  \n";
    print "                                        v$VERSION\n";
    print "***********************************************\n";

    print
"Please indicate your previous 'run' number or enter a new number.\nIf no 'run' ID is entered, a number will be generated for you.\n";
}

sub run_ID {
    print ">:";
    chomp( $run_ID = <STDIN> );
    if ( $run_ID eq "" ) {
        $run_ID = time();
        print "Your ID is now: $run_ID\n";
    }
    $run_directory = "$WORKING_DIR/run/$run_ID";

    if ( -e $run_directory && -d $run_directory ) {
        $g2gc_directory = "$GENOME_DIR/g2gc";
        $run_directory  = "$WORKING_DIR/run/$run_ID";

        #
        $gene_hits_directory = "$run_directory/gene_hits";
        $gene_hits_initial   = "$gene_hits_directory/initial";
        $gene_hits_lookup    = "$gene_hits_directory/lookup";

        #
        $gene_hits_results       = "$run_directory/results";
        $gene_hits_differentials = "$gene_hits_results/differentials";
        $gene_hits_duplications  = "$gene_hits_results/duplications";

        #
    }
    else {
        mkdir( $run_directory, 0755 )
          || die "Cannot be the same name as existing directory!";
        $g2gc_directory = "$GENOME_DIR/g2gc";
        $run_directory  = "$WORKING_DIR/run/$run_ID";

        #
        $gene_hits_directory = "$run_directory/gene_hits";
        $gene_hits_initial   = "$gene_hits_directory/initial";
        $gene_hits_lookup    = "$gene_hits_directory/lookup";

        #
        $gene_hits_results       = "$run_directory/results";
        $gene_hits_differentials = "$gene_hits_results/differentials";
        $gene_hits_duplications  = "$gene_hits_results/duplications";
    }
}

# fdfBLAST General Menu
sub menu {
    print `clear`, "\n";

    # Program title
    print "  __     _   __  ___  _       _    ___  _____ \n";
    print " / _| __| | / _|| _ )| |     /_\\  / __||_   _|\n";
    print "|  _|/ _` ||  _|| _ \\| |__  / _ \\ \\__ \\  | |  \n";
    print "|_|  \\__,_||_|  |___/|____|/_/ \\_\\|___/  |_|  \n";
    print "                                        v$VERSION\n";
    print "***********************************************\n";

    # Options for BLAST
    print "Step 1: BLAST Analysis\n---------------------\n";
    print "B) Automated (Options 1+2)\n";
    print "1) FORMATDB Only\t\t2) BLASTALL Only\n\n";

    # Options for fdfBLAST
    print "Steps 2-5: fdfBLAST Options\n---------------------------\n";
    print "3) Comparative Hit Counts\t4) Reciprocal Hit Matching\n";
    print "5) Ranking and Sorting (Identify Putative Fusions)\n\n";

    # Other Options
    print "Other Options\n-------------\n";
    print "8) Change Run ID\t\tQ) Quit\n\n";

    #
    print "$run_ID>:";

    chomp( $menu_choice = <STDIN> );

    switch ($menu_choice) {

        case "1" {
            &menu_one;
            &return_or_quit;
        }
        case "2" {
            &menu_two;
            &return_or_quit;
        }
        case /B/i {
            &menu_one;
            &menu_two;
            &return_or_quit;
        }
        case "3" {
            &menu_three;
            &return_or_quit;
        }
        case "4" {
            &menu_four;
            &return_or_quit;
        }

        # Hidden setting for running Lookup tables only
        # Part of Diagram Step 2. Automatic under option 3)
        case /L/i {
            &menu_l;
            &return_or_quit;
        }
        case "5" {
            &menu_five;
            &return_or_quit;
        }

        # Hidden setting for Duplications - Exeperimental
        case "6" {
            &menu_six;
            &return_or_quit;
        }
        case /A/i {
            &menu_a;
            &return_or_quit;
        }
        case /F/i {
            &menu_f;
            &return_or_quit;
        }
        case "8" {
            print `clear`, "\n";
            print
"Please indicate your desired run number or enter a new number, 'enter' will generate a ID for you\n";
            &run_ID;
            &menu;
        }
        case /Q/i {
            print `clear`, "\n";
            print
"Thanks for using fdfBLAST!\nPlease cite: \"Leonard, G. & Richards, T.A. 2012. Patterns of gene fusion and fission across the fungi. In Preparation.\"\n";
            last;
        }
        else {
            print `clear`, "\n";
            print "*! Wrong Menu Choice - Please Try Again !*\n";
            &menu;
        }
    }
}

sub menu_one {
    $menu_choice = "FormatDB";
    &print_time("start");
    &get_genomes;
    &formatdb;
    &print_time("end");
}

sub menu_two {
    $menu_choice = "BlastAll";
    &print_time("start");
    &get_genomes;
    &blastall;
    &print_time("end");
}

sub menu_three {
    $menu_choice = "Extract Gene Hit Lists";
    &print_time("start");
    &get_genomes;
    &get_g2gc_files;
    &run_gene_hits_helper;
    &run_gene_hits;
    &generate_lookup_tables;
    &print_time("end");
}

sub menu_l {
    $menu_choice = "Lookup Tables";
    &print_time("start");
    &get_genomes;
    &generate_lookup_tables;
    &print_time("end");
}

sub menu_four {
    $menu_choice = "Parse Tables";
    print
      "Limit by Hit Number? (default: 250, set to 500 for testing dataset)\n";
    $hit_limit = "250";
    print ">:";
    chomp( $hit_limit = <STDIN> );
    &print_time("start");
    &get_genomes;
    &differential_new;
    &print_time("end");
}

sub menu_five {
    $menu_choice = "Identify Fusions";
    &print_time("start");
    &parse_lookup;
    &read_domain_file;
    &fusion_scan;
    &print_time("end");
}

# Currently Hidden Menu - Experimental.
sub menu_six {
    $menu_choice = "Identify Duplications";
    &print_time("start");
    &duplication_scan;
    &print_time("end");
}

sub return_or_quit {
    print "\nReturn to (M)enu or (Q)uit\n";
    print ">:";
    chomp( $menu_choice = <STDIN> );
    switch ($menu_choice) {
        case /Q/i {
            print `clear`, "\n";
            print "Thanks for using fdfBLAST, goodbye!ah right\n";
            print
"Thanks for using fdfBLAST!\nPlease cite: \"Leonard, G. & Richards, T.A. 2012. Patterns of gene fusion and fission across the fungi. In Preparation.\"\n";
            last;
        }
        case /M/i {
            print `clear`, "\n";
            &menu;
        }
        else {
            &return_or_quit;
        }
    }
}

sub print_time {
    my $time = shift;
    if ( $time eq "start" ) {
        open my $time_fh, ">>$run_directory/time.txt";
        print $time_fh "Menu Choice = $menu_choice\n";
        $start_time = scalar localtime(time);
        print $time_fh " Start Time = $start_time\n";
        close($time_fh);
    }
    elsif ( $time eq "end" ) {
        open $time_fh, ">>$run_directory/time.txt";
        $end_time = scalar localtime(time);
        print $time_fh "   End Time = $end_time\n";
        &dhms;
        printf $time_fh (
            " Total Time = %4d days%4d hours%4d minutes%4d seconds\n\n",
            @parts[ 7, 2, 1, 0 ]
        );
        close($time_fh);
    }
}

sub dhms {

    my %months = (
        Jan => 1,
        Feb => 2,
        Mar => 3,
        Apr => 4,
        May => 5,
        Jun => 6,
        Jul => 7,
        Aug => 8,
        Sep => 9,
        Oct => 10,
        Nov => 11,
        Dec => 12
    );

    my $beginning = $start_time;
    my $end       = $end_time;

    my @b = split( /\s+|:/, $beginning );
    my @e = split( /\s+|:/, $end );

    my $b =
      timelocal( $b[5], $b[4], $b[3], $b[2], $months{ $b[1] } - 1, $b[-1] );
    my $e =
      timelocal( $e[5], $e[4], $e[3], $e[2], $months{ $e[1] } - 1, $e[-1] );

    my $elapsed = $e - $b;

    @parts = gmtime($elapsed);

    return @parts;

}

sub print_log {
    my $message = shift;
    open my $log_fh, ">>$run_directory/log.txt";
    print $log_fh "$message";
    close($log_fh);
}

## Get all genome names from the genome directory and put into an array
# Files must have .fas as their extension
sub get_genomes {

# Assign all files with extension .fas (a small assumption) to the array @file_names
    @file_names = glob "$GENOME_DIR/*.fas";

    # Loop for all file names in @file_names array
    foreach my $file (@file_names) {
        $file = fileparse($file);
    }

    # This is the number of .fas files in the genome directory
    $genome_num = @file_names;
}

## A method to call the formatdb program from the command line
# and run it for each .fas FASTA formatted genome in a directory
# This will probably need updating to include the new BLAST+ programs
sub formatdb {

    # User output, number of genomes to format

    print "Number of Genomes = $genome_num\n";
    print "FORTMATDB - Formating Databases...\n";

    # For loop, $i up to number of genomes do

    for ( my $i = 0 ; $i < $genome_num ; $i++ ) {
        print "$file_names[$i] \x3E\x3E\x3E";

        # This line calls the formatdb program within a previously set directory
        # it then outputs the databse genomes files to another directory
        $results =
`$BLAST_BIN_DIR/formatdb -t $GENOME_DIR/$file_names[$i] -i $GENOME_DIR/$file_names[$i] -p T -o T`;
        print " Completed\n";
    }
    print "FORMATDB - Finished...\n\n";
}

# This will probably need updating to include the new BLAST+ programs
sub blastall {

# Create the g2gc output folder, if not then error, quit. This is a little harsh...
    mkdir( $g2gc_directory, 0777 );
    #####      || die "This procedure has already been run, please refer to you previous run, $run_ID.\n";
    print
"Number of\: Genomes = $genome_num\tOutput Files = $genome_num x $genome_num = "
      . $genome_num * $genome_num . "\n";
    print "Performing BLASTP\n";

    for ( my $i = 0 ; $i < $genome_num ; $i++ ) {
        for ( my $j = 0 ; $j < $genome_num ; $j++ ) {
            ( $file_i, $dir, $ext ) = fileparse( $file_names[$i], '\..*' );
            ( $file_j, $dir, $ext ) = fileparse( $file_names[$j], '\..*' );
            print "$file_i to $file_j \x3E\x3E\x3E ";

            if ( -e "$g2gc_directory/$file_j\_$file_i.bpo" ) {
                print "already exists...Skipping\n";
            }
            else {
                $results =
`$BLAST_BIN_DIR/blastall -p blastp -d $GENOME_DIR/$file_names[$i] -i $GENOME_DIR/$file_names[$j] -m 0 -a $core_num -F F -o $g2gc_directory/$file_j\_$file_i.bpo`;
                print "Completed\n";
            }
        }
    }
    print "\nBLASTALL Successful\n";
}

# Retrieve list of Blast Parse Output files from directory
sub get_g2gc_files {

    @g2gc_file_names = glob "$g2gc_directory/*.bpo";
    $g2gc_length     = @g2gc_file_names;
}

# Change user submission from 1e-10 to 0.01 etc
sub evaluate {
    $value = shift;

# Perl method bstr will only change '1e' not 'e' to decimal, therefore prefix value with '1'
    if ( $value =~ m/^e/ ) {
        $value = "1" . $value;
        $value = Math::BigFloat->bstr($value);
    }
    if ( $value =~ m/e/ ) {
        $value = Math::BigFloat->bstr($value);
    }
    return $value;
}

sub run_gene_hits_helper {

    # Create the gene_hits output folder, if unable then on error, quit.
    # This could probably be handled better - overwrite? sub directory? etc
    mkdir( $gene_hits_directory, 0755 )
      || die
"This procedure has already been run, please refer to you previous run , $run_ID.\n";

    mkdir( $gene_hits_initial, 0755 );
    print "E Value Upper Limit - e.g 1e-10\n>:";
    chomp( $upper_limit = <STDIN> );
    if ( $upper_limit == "" ) {
        $upper_limit = "1e-10";
    }
    &print_log("Gene Hit List E-Value Range\n");
    &print_log("Upper Value = $upper_limit\n");
    &evaluate($upper_limit);
    $upper_limit = $value;
    print "E Value Lower Limit - e.g 0\n>:";
    chomp( $lower_limit = <STDIN> );
    if ( $lower_limit eq "" ) {
        $lower_limit = "0";
    }
    &print_log("Lower Value = $lower_limit\n");
    &evaluate($lower_limit);
    $lower_limit = $value;
    &print_log("Upper Value = $upper_limit\n");
    &print_log("Lower Value = $lower_limit\n\n");
    print "\nWorking\n";
}

sub run_gene_hits {

    for ( $i = 0 ; $i < $g2gc_length ; $i++ ) {
        $current = $i + 1;
        print "\n$current of $g2gc_length";

        # Assign file name in iterative loop
        $in_file = "$g2gc_file_names[$i]";

        ( $in_filename, $dir, $ext ) = fileparse( $in_file, '\..*' );

        #$in_filename = fileparse($in_file);
        print "    $in_filename";

        # Internal iterators for 2D array
        # X = Blast sequence number, Y = Hits against query
        $x_pos = 0;

        #@comparison_length = ();
        undef(@comparison_length);

        my $search = new Bio::SearchIO(
            '-format' => 'blast',
            '-file'   => $in_file
        );
        while ( my $result = $search->next_result ) {
            $y_pos                   = 0;
            $query_accession         = $result->query_accession;
            $query_length            = $result->query_length;
            $query_genome[$x_pos][0] = "$query_accession";
            $query_genome[$x_pos][1] = "$query_length";
            while ( my $hit = $result->next_hit ) {
                $hit_accession    = $hit->accession;
                $hit_length       = $hit->length;
                $hit_significance = $hit->significance;

                while ( my $hsp = $hit->next_hsp ) {
                    $percent_identity = $hsp->percent_identity;

                    # Range of query that hits the subject gene
                    @query_range = $hsp->range('query');
                    #####@hit_range   = $hsp->range('hit');

                    # We don't really need 13 decimal places, round to none..
                    $percent_identity = sprintf( "%.0f", $percent_identity );
                }

                $comparison_genome[$x_pos][$y_pos] =
"$hit_accession,$hit_length,$hit_significance,$percent_identity,$query_range[0],$query_range[1],";
                $y_pos++;
            }

            push( @comparison_length, $y_pos );
            $x_pos++;
        }
        open my $out_fh, ">$gene_hits_initial/$in_filename.csv";
        ##### Surely we don't need this twice? I don't know Guy, what were you trying to do?
        $query_genome_length = $x_pos;

        #Append hit numbers to query_genome array
        for ( $x = 0 ; $x < $query_genome_length ; $x++ ) {

            $comparison_genome_length = $comparison_length[$x];
            $count                    = 0;
            for ( $d = 0 ; $d < $comparison_genome_length ; $d++ ) {

                $evalue = "$comparison_genome[$x][$d]";
                $evalue =~ m/(.*?)\,(\d*)\,(.*?)\,(.*?)\,/;
                $evalue = $3;
                &evaluate( $value = "$evalue" );
                $evalue = $value;
                if (   ( $evalue <= $upper_limit )
                    && ( $evalue >= $lower_limit ) )
                {
                    $count++;
                }
            }
            $query_genome[$x][2] = $count;
        }

        for ( $x = 0 ; $x < $query_genome_length ; $x++ ) {

            print $out_fh
              "$query_genome[$x][2],$query_genome[$x][0],$query_genome[$x][1],";

            $comparison_genome_length2 = $query_genome[$x][2];

            for ( $y = 0 ; $y < $comparison_genome_length2 ; $y++ ) {

                $evalue = "$comparison_genome[$x][$y]";
                $evalue =~ m/(.*?)\,(\d*)\,(.*?)\,(.*?)\,/;
                $evalue = $3;
                &evaluate( $value = "$evalue" );
                $evalue = $value;

                ## 0 is the lower limit as the lower the E-value,
                ## or the closer it is to zero, the more "significant" the match is
                ## therefore
                if (   ( $evalue <= $upper_limit )
                    && ( $evalue >= $lower_limit ) )
                {
                    print $out_fh "$comparison_genome[$x][$y]";
                }
            }
            print $out_fh "\n";
        }
        close($out_fh);

        undef(@comparison_genome);
        undef(@query_genome);
    }
    print "\nFinished\n";

}

# A method to create tables of the recorded hits for each genome group
sub generate_lookup_tables {
    print "Generating Lookup Tables";

    # Check to see if the initial hit lists exist
    if ( -e $gene_hits_initial && -d $gene_hits_initial ) {

        # Number of gene hit files
        $gene_hits_length = @gene_hits_file_names;

        # Edited gene hit lists get their own directory
        mkdir( $gene_hits_lookup, 0755 );

        # Internal counter
        $run = 0;

        # As we have called get_genomes in menu_l we can use
        # genome_num and file_names - which I need to redo for strictures!!!
        # Two for loops to iterate through each gene hits initial file
        # and add the values to a 2d array for output to file
        for ( $i = 0 ; $i < $genome_num ; $i++ ) {
            for ( $j = 0 ; $j < $genome_num ; $j++ ) {

                ( $file_i, $dir, $ext ) = fileparse( $file_names[$i], '\..*' );
                ( $file_j, $dir, $ext ) = fileparse( $file_names[$j], '\..*' );
                open my $in_fh, "<$gene_hits_initial/$file_i\_$file_j\.csv";

              # We're using while with an iterator instead of foreach to
              # read in the file line by line adding each column to the 2d array
                $x = 0;
                while (<$in_fh>) {
                    my ($line) = $_;
                    chomp($line);
                    @temp = split( /,/, $line );
                    $AoA[$x][$j] = $temp[0];
                    $x++;
                }

                # Then we add up the values in each column for each file
                # using a for loop
                $total_in_array_column = 0;
                for my $i ( 0 .. $#AoA ) {
                    $total_in_array_column += $AoA[$i][$j];
                }

                # Then those values are added plainly to an array,
                # the sum of all values is calculated and if they are
                # equal to 0 then there are no hits, else we continue...
                push( @check_for_zero, $total_in_array_column );
                $total_of_array_columns += $_ for @check_for_zero;
                if ( $total_of_array_columns == 0 ) {
                    print
"\nThere are no hits at all, please try some different e-values...\n";
                    last;
                }
                ###
            }    # End second (internal) for loop
            close($in_fh);
            open my $out_fh, ">$gene_hits_lookup/$file_i\.csv";

            # Here we use two for loops to iterate through the 2D array
            # for $i from 0 to length of @AoA
            for my $i ( 0 .. $#AoA ) {

                # A reference to the row at position $i
                $aref = $AoA[$i];

                # Length of row
                $n          = @$aref - 1;
                $line_total = 0;

                # Second for, $j from 0 to length of row
                for my $j ( 0 .. $n ) {

                    # Sum the individual numbers of each row
                    $line_total = $line_total + $AoA[$i][$j];

                    # Print out each element of the row
                    print $out_fh "$AoA[$i][$j],";
                }
                print $out_fh "\n";
            }

            # Increment run
            $run = $run + 1;

            # Reset array, so length is not kept at largest
            # @AoA = ();
            undef(@AoA);
            print ".";
        }
        close($out_fh);
    }
    print "\n";
}

sub parse_lookup {

    print "Lower Ratio Value (default: 0.1)\n>:";
    chomp( $lower_ratio = <STDIN> );
    if ( $lower_ratio == "" ) {
        $lower_ratio = "0.1";
    }
    print "Higher Ratio Value (default: 1.0)\n>:";
    chomp( $higher_ratio = <STDIN> );
    if ( $higher_ratio == "" ) {
        $higher_ratio = "1.0";
    }
}

sub differential_new {

    # Output directory
    mkdir( $gene_hits_results,       0755 );
    mkdir( $gene_hits_differentials, 0755 );
    mkdir( $gene_hits_duplications,  0755 );

    print "\nDifferential Extraction\n";
    $verbose          = 0;
    @variables        = @_;
    $limit            = $variables[0];
    @lookup_filenames = glob "$gene_hits_lookup/*.csv";
    $lookup_number    = @lookup_filenames;

    for ( $k = 0 ; $k < $lookup_number ; $k++ ) {
        #####
        print "K$k\n" if $verbose == 1;
        #####
        for ( $i = 0 ; $i < $lookup_number ; $i++ ) {

# An option to include self-genome comparison could go here, e.g. Genome A to Genome A
# currently they are excluded...not useful for synapomorphy analysis!
            if ( $k != $i ) {

                # This is where the main comparisons occur
                # e.g. a-b, a-c & a-d
                # uses $k and $i
                #print "\n$lookup_filenames[$k]\t$lookup_filenames[$i]\n";
                $file_k = fileparse( $lookup_filenames[$k] );
                $file_i = fileparse( $lookup_filenames[$i] );

                #####
                print "Processing\t$file_k\tand\t$file_i\n";
                #####
                &comparison( $lookup_filenames[$k], $lookup_filenames[$i], $k,
                    $i, $lookup_filenames[$k] );
            }
        }
    }
    &remove_single_recips;
}

sub comparison {

    @variables     = @_;
    $ifile         = $variables[0];
    $jfile         = $variables[1];
    $ivar          = $variables[2];
    $jvar          = $variables[3];
    $kfile         = $variables[4];
    $kfile         = fileparse($kfile);
    $ifile_compare = fileparse($ifile);
    $verbose       = 0;

    open my $query_fh, "<$ifile";
    #####
    print "\tOpening: $ifile\n" if $verbose == 1;
    #####
    @query_array = ();
    while (<$query_fh>) {
        my ($line) = $_;
        chomp($line);
        my ($line_number) = $.;
        @temp = split( /,/, $line );
        push( @query_array, "$temp[$jvar]" );
    }

    $query_array_length = @query_array;
    for ( $l = 0 ; $l < $query_array_length ; $l++ ) {
        if ( $query_array[$l] == "0" && $ifile_compare eq $kfile ) {
            $line_number = $l + 1;
            &get_accession( $ifile, $jfile, $line_number );
            $missing_genome = fileparse($jfile);
            open my $nohits_fh, ">>$gene_hits_results/no_hits_$missing_genome";
            print $nohits_fh "$accession,\n";
        }
        elsif ( $query_array[$l] == "1" ) {

            $line_number = $l + 1;
            &get_accession( $ifile, $jfile, $line_number );
            #####
            print "Q$accession_array[0],$accession_array[1]" if $verbose == 1;
            #####
            &get_hit_info( $ifile, $jfile, $line_number, "0" );
            #####
            print
" hits S\_$return_array[0],$return_array[1],$return_array[2],$return_array[3]\n"
              if $verbose == 1;
            #####
            &scan_subject_initial( $accession_array[0], $ifile, $jfile,
                $return_array[0], "0" );
            #####
            print "NStatus = $scan_array[0]\n\n" if $verbose == 1;
            #####

            if ( $scan_array[0] eq "recip" ) {

                open my $non_diff_fh,
                  ">>$gene_hits_duplications/duplications_$kfile";

                &get_line_number_start(
                    "$gene_hits_duplications/duplications_$kfile");

                print $non_diff_fh
"$last_line_number;$accession_array[0];$accession_array[1];$return_array[0];$return_array[1];$return_array[2];$return_array[3]\n";
            }
        }
        elsif ( $query_array[$l] >= "2" && $query_array[$l] <= $hit_limit ) {
            $verbose = 0;
            #####
            print "\# Query Hits = $query_array[$l]\n" if $verbose == 2;
            #####
            for ( $n = 0 ; $n < $query_array[$l] ; $n++ ) {

                $line_number = $l + 1;
                &get_accession( $ifile, $jfile, $line_number );
                $query_number = $n + 1;
                #####
                print "Q$query_number\_$accession_array[0],$accession_array[1]"
                  if $verbose == 2;
                #####
                &get_hit_info( $ifile, $jfile, $line_number, $n );
                #####
                print
" hits S\_$return_array[0],$return_array[1],$return_array[2],$return_array[3]\n"
                  if $verbose == 1;
                #####
                &scan_subject_initial( $accession_array[0], $ifile, $jfile,
                    $return_array[0], $n );

                #####
                print "Status = $scan_array[0]\n\n" if $verbose == 2;
                #####

                if ( $scan_array[0] eq "recip" ) {

                    open my $differentials_fh,
                      ">>$gene_hits_differentials/differentials_$kfile";

                    &get_line_number_start(
                        "$gene_hits_differentials/differentials_$kfile");

                    if ( $last_accession eq $accession_array[0] ) {
                        $last_line_number = $last_line_number + 1;
                    }
                    elsif ( $last_accession ne $accession_array[0] ) {
                        $last_line_number = 0;
                    }

                    print $differentials_fh
"$last_line_number;$accession_array[0];$accession_array[1];$return_array[0];$return_array[1];$return_array[2];$return_array[3];$return_array[4];$return_array[5]\n";
                }    #endif
            }    #endfor
        }    #endelsif

    }    #endfor
    #print $non_diff_fh "\n";
}

sub remove_single_recips {

    @diff_file_names = glob "$gene_hits_differentials/*.csv";
    $num_files       = @diff_file_names;

    for ( $q = 0 ; $q < $num_files ; $q++ ) {
        print "\nRemoving singles from $diff_file_names[$q]\n";
        open my $diff_fh, "<$diff_file_names[$q]";

        ( $file, $dir, $ext ) = fileparse( $diff_file_names[$q], '\..*' );
        @temp_array = ();
        while (<$diff_fh>) {
            my ($line) = $_;
            chomp($line);
            my ($line_number) = $.;
            push( @temp_array, "$line" );
        }

        $temp_array_length = @temp_array;

        for ( $r = 0 ; $r < $temp_array_length ; $r++ ) {

            @line_n  = split( /;/, $temp_array[$r] );
            @line_n1 = split( /;/, $temp_array[ $r + 1 ] );

            $ln_n  = $line_n[0];
            $ln_n1 = $line_n1[0];

            if ( $ln_n eq 0 && $ln_n1 eq 0 ) {

                # Do Nothing
            }
            else {
                open my $output_fh, ">>$dir$file\_fixed$ext";
                foreach my $element (@line_n) {
                    print $output_fh $element . "\;";
                }

                print $output_fh "\n";
            }
        }
        print "End\n";
        close($output_fh);
    }
}

sub fusion_scan {

    $verbose = 0;

    print "\nPerforming Fusion Scans\n";

    @fusion_filenames = glob "$gene_hits_differentials/*fixed.csv";
    $fusion_number    = @fusion_filenames;

    print "Number of files = $fusion_number\n";

    for ( my $i = 0 ; $i < $fusion_number ; $i++ ) {

        undef(@fusion);
        open my $fusion_fh, "<$fusion_filenames[$i]";
        $input_file = fileparse( $fusion_filenames[$i] );
        print "Opening $fusion_filenames[$i]\n";

        while (<$fusion_fh>) {
            my ($line) = $_;
            chomp($line);
            push( @fusion, "$line" );
        }

        $fusion_array_length = @fusion;

        for ( my $j = 0 ; $j < $fusion_array_length ; $j++ ) {

            @array_line = split( /;/, $fusion[$j] );
            $line_number = $array_line[0];

            if ( $line_number eq 0 ) {

                $query_accession = $array_line[1];
                $query_length    = $array_line[2];

                $subject_accession       = $array_line[3];
                $subject_length          = $array_line[4];
                $subject_evalue          = $array_line[5];
                $subject_hit_range_start = $array_line[7];
                $subject_hit_range_end   = $array_line[8];
                $match_length =
                  $subject_hit_range_end - $subject_hit_range_start;
                push( @subject_hits,
"$subject_accession,$subject_hit_range_start,$subject_hit_range_end,$subject_length"
                ) if $match_length > 50 && $subject_length > 50;
                #####
                print "\nQ $query_accession -> $subject_accession\n"
                  if $verbose == 1;
                #####

                for ( my $k = $j + 1 ; $k < $fusion_array_length ; $k++ ) {

                    @next_array_line = split( /;/, $fusion[$k] );
                    $next_line_number = $next_array_line[0];

                    if ( $next_line_number ne 0 ) {

                        $next_query_accession = $next_array_line[1];

                        $next_subject_accession       = $next_array_line[3];
                        $next_subject_length          = $next_array_line[4];
                        $next_subject_evalue          = $next_array_line[5];
                        $next_subject_hit_range_start = $next_array_line[7];
                        $next_subject_hit_range_end   = $next_array_line[8];

             # > 50 should be user selectable limit for now it is hard-coded
             # to exclude all sequences with <=50 bases.
             # Another limit is Matched Length (end - start) this should also be
             # more than 50.
                        $next_match_length =
                          $next_subject_hit_range_end -
                          $next_subject_hit_range_start;

                        if (   $next_line_number ne 0
                            && $query_accession eq $next_query_accession
                            && $next_subject_length > 50
                            && $next_match_length > 50 )
                        {
                            print
"\n$next_line_number ne 0 && $query_accession eq $next_query_accession && $next_subject_length > 50 && $next_match_length > 50\n"
                              if $verbose == 1;
                            push( @subject_hits,
"$next_subject_accession,$next_subject_hit_range_start,$next_subject_hit_range_end,$next_subject_length,$next_subject_evalue"
                            );
                        }
                    }
                    else {
                        last;
                    }
                }
                &rank_sort;
                undef(@subject_hits);
            }
        }
        close($fusion_fh);
    }
}

sub rank_sort {
    $subject_hit_length = @subject_hits;

    # Get fusion ORF length and half
    $fusionORF = $query_length;
    $half_fusion_length = sprintf( "%.0f", $query_length / 2 );
    print "Fusion Length - $fusionORF / 2 = $half_fusion_length\n"
      if $verbose == 1;

    for ( $q = 0 ; $q < $subject_hit_length ; $q++ ) {
        @unfused = split( /,/, $subject_hits[$q] );

        #print "Q = $q\n";
        if ( $unfused[1] <= $half_fusion_length ) {
            $pos1 = "L";
        }
        elsif ( $unfused[1] >= $half_fusion_length ) {
            $pos1 = "R";
        }
        if ( $unfused[2] <= $half_fusion_length ) {
            $pos2 = "L";
        }
        elsif ( $unfused[2] >= $half_fusion_length ) {
            $pos2 = "R";
        }

        #print "$pos1 - $unfused[1] , $unfused[2] - $pos2\n";
        if ( "$pos1$pos2" eq "LL" ) {

            #print "$pos1$pos2\t$unfused[0] ORF is Leftmost\n";
            push( @leftmost, "$subject_hits[$q]" );
        }
        elsif ( "$pos1$pos2" eq "RR" ) {

            #print "$pos1$pos2\t$unfused[0] ORF is Rightmost\n";
            push( @rightmost, "$subject_hits[$q]" );
        }
        elsif ( "$pos1$pos2" eq "LR" ) {

            #print "$pos1$pos2\t$unfused[0] ORF is more left\n";
            push( @middles, "$subject_hits[$q]" );
        }
        elsif ( "$pos1$pos2" eq "RL" ) {

            #print "$pos1$pos2\t$unfused[0] ORF is more right\n";
            push( @middles, "$subject_hits[$q]" );
        }
    }

    $left_num   = @leftmost;
    $right_num  = @rightmost;
    $middle_num = @middles;
    print "L:$left_num\tR:$right_num\tM:$middle_num\n" if $verbose == 1;

    if ( $left_num == 0 && $right_num == 0 && $middle_num == 0 ) {

        # Discard
    }
    elsif ($left_num == 0 && $right_num == 0
        or $right_num == 0 && $middle_num == 0
        or $left_num == 0  && $middle_num == 0 )
    {

        # Discard
    }
    elsif ( $left_num == 0 && $right_num >= 1 and $middle_num >= 1 ) {

        &ignore_orthologues;
        print "\t\tXX $continue XX\n" if $verbose == 1;
        if ( $continue eq "yes" ) {
            &middle_and_right;
        }
    }
    elsif ( $right_num == 0 && $left_num >= 1 and $middle_num >= 1 ) {

        &ignore_orthologues;
        print "\t\tXX $continue XX\n" if $verbose == 1;
        if ( $continue eq "yes" ) {
            &middle_and_left;
        }

    }
    elsif ( $middle_num == 0 && $left_num >= 1 and $right_num >= 1 ) {

        print "011\n" if $verbose == 1;
        &left_and_right;

    }
    elsif ( $left_num >= 1 && $right_num >= 1 and $middle_num >= 1 ) {

        #print "\nXX - All - XX\n";
        &ignore_orthologues;
        print "\t\tXX $continue XX\n" if $verbose == 1;
        if ( $continue eq "yes" ) {
            print "111\n" if $verbose == 1;
            &left_and_right;
            &middle_and_left;
            &middle_and_right;
        }
    }
    undef(@middles);
    undef(@leftmost);
    undef(@rightmost);
}

sub left_and_right {
    print "\nLR\n" if $verbose == 1;
    for ( my $a = 0 ; $a < $left_num ; $a++ ) {
        print "A:$a\n" if $verbose == 1;
        @unfused_one = split( /,/, $leftmost[$a] );
        my $one_end = $unfused_one[2];

        for ( my $b = 0 ; $b < $right_num ; $b++ ) {
            print "\tB:$b\n" if $verbose == 1;
            @unfused_two = split( /,/, $rightmost[$b] );
            my $two_start = $unfused_two[1];

            print "\t$query_accession -> $unfused_one[0] + $unfused_two[0]\n"
              if $verbose == 1;

            $ratio = $one_end / $two_start;
            $ratio = sprintf( "%.2f", $ratio );
            print "\tRlr: $one_end / $two_start = $ratio\n" if $verbose == 1;

            if ( $ratio >= $lower_ratio && $ratio <= $higher_ratio ) {
                undef(@subject_hits);
                push( @subject_hits,
"$unfused_one[0],$unfused_one[1],$unfused_one[2],$unfused_one[3]"
                );
                push( @subject_hits,
"$unfused_two[0],$unfused_two[1],$unfused_two[2],$unfused_two[3]"
                );
                $subject_hit_length = @subject_hits;
                &generate_image;
            }
        }
    }
}

sub ignore_orthologues {

    print "Ignoring Potential Orthologues\n" if $verbose == 1;
    $continue = "yes";

    for ( my $a = 0 ; $a < $middle_num ; $a++ ) {
        @unfused_middles = split( /,/, $middles[$a] );
        my $middle_length = $unfused_middles[3];
        print "\t$query_accession - $middle_evalue\n" if $verbose == 1;
        $length_ratio = $middle_length / $query_length;

        if ( $length_ratio <= "0.95" ) {

            #$continue = "yes";
            push( @cont, "yes" );
            print
"\t\t$unfused_middles[0] - $middle_length / $query_length = $length_ratio - yes\n"
              if $verbose == 1;
        }
        elsif ( $length_ratio > "0.95" ) {

            push( @cont, "no" );
            print
"\t\t$unfused_middles[0] - $middle_length / $query_length = $length_ratio - no\n"
              if $verbose == 1;
        }

    }
    foreach my $item (@cont) {
        print $item . "," if $verbose == 1;
    }
    print "\n" if $verbose == 1;

    $continue = "no" if ( grep /^no$/, @cont );

    undef(@cont);
    return $continue;
}

sub middle_and_left {
    print "\nML\n" if $verbose == 1;
    for ( my $a = 0 ; $a < $left_num ; $a++ ) {
        print "A:$a\n" if $verbose == 1;
        @unfused_one = split( /,/, $leftmost[$a] );

        my $one_end = $unfused_one[2];

        for ( my $b = 0 ; $b < $middle_num ; $b++ ) {
            print "\tB:$b\n" if $verbose == 1;
            @unfused_two = split( /,/, $middles[$b] );

            my $two_start = $unfused_two[1];
            $two_match_length = $unfused_two[2] - $unfused_two[1];
            print "$query_accession -> $unfused_one[0] + $unfused_two[0]\n"
              if $verbose == 1;

            $ratio = $one_end / $two_start;
            $ratio = sprintf( "%.2f", $ratio );

            if ( $two_start <= $one_end ) {
                print "\t$two_start <= $one_end Overlap!\n" if $verbose == 1;
            }
            else {
                print "\tRml: $one_end / $two_start = $ratio\n"
                  if $verbose == 1;
                if ( $ratio >= $lower_ratio && $ratio <= $higher_ratio ) {
                    undef(@subject_hits);
                    push( @subject_hits,
"$unfused_one[0],$unfused_one[1],$unfused_one[2],$unfused_one[3]"
                    );
                    push( @subject_hits,
"$unfused_two[0],$unfused_two[1],$unfused_two[2],$unfused_two[3]"
                    );
                    $subject_hit_length = @subject_hits;
                    &generate_image;
                }
            }
        }
    }
}

sub middle_and_right {
    print "\nMR\n" if $verbose == 1;
    for ( my $a = 0 ; $a < $middle_num ; $a++ ) {
        print "A:$a\n" if $verbose == 1;
        @unfused_one = split( /,/, $middles[$a] );

        my $one_end = $unfused_one[2];

        for ( my $b = 0 ; $b < $right_num ; $b++ ) {
            print "\tB:$a\n" if $verbose == 1;
            @unfused_two = split( /,/, $rightmost[$b] );

            my $two_start = $unfused_two[1];
            $two_match_length = $unfused_two[2] - $unfused_two[1];
            print "$query_accession -> $unfused_one[0] + $unfused_two[0]\n"
              if $verbose == 1;

            $ratio = $one_end / $two_start;
            $ratio = sprintf( "%.2f", $ratio );

            if ( $one_end >= $two_start ) {
                print "\t$one_end >= $two_start Overlap!\n" if $verbose == 1;
            }
            else {
                print "\tRmr: $one_end / $two_start = $ratio\n"
                  if $verbose == 1;
                if ( $ratio >= $lower_ratio && $ratio <= $higher_ratio ) {
                    undef(@subject_hits);
                    push( @subject_hits,
"$unfused_one[0],$unfused_one[1],$unfused_one[2],$unfused_one[3]"
                    );
                    push( @subject_hits,
"$unfused_two[0],$unfused_two[1],$unfused_two[2],$unfused_two[3]"
                    );
                    $subject_hit_length = @subject_hits;
                    &generate_image;
                }
            }
        }
    }
}

## These could be call to print_log with a filename variable...
sub print_comp_list {
    my $message = shift;
    open my $log_fh, ">>$run_directory/composite_list.csv";
    print $log_fh "$message";
    close($log_fh);
}

sub print_split_list {
    my $message = shift;
    open my $log_fh, ">>$run_directory/split_list.csv";
    print $log_fh "$message";
    close($log_fh);
}

sub read_domain_file {
    print "Name of all-domains file?\n>:";
    chomp( $domain_in = <STDIN> );
    open my $domains_fh, "$run_directory\/$domain_in";
    while (<$domains_fh>) {
        $line = $_;
        push( @domains, $line );
    }
    close($domains_fh);
}

sub generate_image {

    $query_acc_size = ( length $query_accession );
    push( @names_length, $query_acc_size );

    $number_of_hits = $subject_hit_length;

    $unfused_one_length = ( length $unfused_one[0] );
    push( @names_length, $unfused_one_length );
    $unfused_two_length = ( length $unfused_two[0] );
    push( @names_length, $unfused_two_length );

    # Get longest gene accession and * 6 for gene accession length in 'pixels'
    my $highest;
    for (@names_length) {
        $highest = $_ if $_ > $highest;
    }
    $name_length = $highest * 6;

    # Image specific
    $padding_left  = 100;
    $padding_right = 100;
    $padding       = $padding_left + $padding_right;

    # place to start drawing sequences from
    $left_pos = $padding_left + $name_length;

    # Image Width(length?)
    $image_length = $query_length + $padding + $name_length;
    if ( $image_length < 500 ) {
        $image_length += 250;
    }

    # number of hits * 50 + 50 for query + 20 for padding
    $image_height = ( $number_of_hits * 50 ) + 50 + 20;

    # set positions
    $font_vpos      = 2;
    $bar_vpos       = 20;
    $accession_vpos = 15;
    $start          = 0;
    $end            = 0;

    # create a new image
    $im = new GD::Image( $image_length, $image_height, 1 )
      ;    # Not paletted - allows for alpha transparency of domains

    # allocate some specific colors
    $white = $im->colorAllocate( 255, 255, 255 );
    $black = $im->colorAllocate( 0,   0,   0 );
    $red   = $im->colorAllocate( 255, 0,   0 );
    $blue  = $im->colorAllocate( 0,   21,  181 );

    # make the background non-transparent and interlaced
    $im->transparent(-1);    # no transparency
    $im->interlaced('true');
    $im->alphaBlending(1);

    # White background
    $im->filledRectangle( 0, 0, $image_length, $image_height, $white );

    # Put a black frame around the picture
    $im->rectangle( 0, 0, $image_length - 1, $image_height - 1, $black );

    $actual_length = "";
    &draw_query;
    &draw_subjects;
    &watermark;

    @unfused_one = split( /,/, $subject_hits[0] );
    @unfused_two = split( /,/, $subject_hits[1] );

    # We only really need 1dp for the folders, 2dp ratio is printed in image...
    $ratio = sprintf( "%.1f", $ratio );

    $query_accession_dir = "$gene_hits_differentials/$query_accession";

    &print_comp_list("$query_accession,\n");
    &print_split_list("$unfused_one[0],\n$unfused_two[0],\n");

    if ( -e $query_accession_dir && -d $query_accession_dir ) {

        $ratio_dir = "$query_accession_dir/$ratio";
        if ( -e $ratio_dir && -d $ratio_dir ) {

            open my $out_fh, '>',
"$ratio_dir/$query_accession\_\_$unfused_one[0]\_\_$unfused_two[0].png";
            binmode $out_fh;
            print $out_fh $im->png(9);
            close($out_fh);
        }
        else {
            mkdir( $ratio_dir, 0755 );
            open $out_fh, '>',
"$ratio_dir/$query_accession\_\_$unfused_one[0]\_\_$unfused_two[0].png";
            binmode $out_fh;
            print $out_fh $im->png(9);
            close($out_fh);
        }
    }
    else {

        mkdir( $query_accession_dir, 0755 );
        $ratio_dir = "$query_accession_dir/$ratio";
        if ( -e $ratio_dir && -d $ratio_dir ) {

            open $out_fh, '>',
"$ratio_dir/$query_accession\_\_$unfused_one[0]\_\_$unfused_two[0].png";
            binmode $out_fh;
            print $out_fh $im->png(9);
            close($out_fh);
        }
        else {
            mkdir( $ratio_dir, 0755 );
            open $out_fh, '>',
"$ratio_dir/$query_accession\_\_$unfused_one[0]\_\_$unfused_two[0].png";
            binmode $out_fh;
            print $out_fh $im->png(9);
            close($out_fh);
        }
    }
}

sub draw_query {
    $bar_length = $query_length;
    $accession  = $query_accession;
    &scale_bars;

    #
    $end       = $query_length;
    $end_scale = $query_length;
    &length_bars;

    #
    &bar($blue);
    &accession_name($blue);
    &draw_domain;
}

sub draw_subjects {

    for ( $c = 0 ; $c < $number_of_hits ; $c++ ) {

        @temp          = split( /,/, $subject_hits[$c] );
        $left_pos      = $left_pos + $temp[1];
        $bar_vpos      = $bar_vpos + 50;
        $bar_length    = ( $temp[2] - $temp[1] );
        $actual_length = $temp[3];
        $font_vpos     = $font_vpos + 50;

        #
        $colour_ratio = $bar_length / $actual_length;
        $colour_ratio = sprintf( "%.2f", $colour_ratio );

        &get_colour;

        #
        &bar($ratio_colour);

        #
        $start     = $temp[1];
        $end       = $bar_length;
        $end_scale = $end + $start;
        &scale_bars;

        #
        &length_bars;

        #
        $accession = $temp[0];
        $left_pos  = $left_pos - $temp[1];
        ##$left_pos = $padding_left + $name_length;
        $accession_vpos = $accession_vpos + 50;
        &accession_name($ratio_colour);
        ###############################
        # I think this section fixes the second domain out of sync bug!
        if ( $c eq 0 ) {

            #print "***C*** = $c\n";
            &draw_domain;
        }
        elsif ( $c ge 1 ) {

            #$left_pos = $padding_left + $name_length + 400;
            $left_pos = $padding_left + $temp[1] + $name_length;

            #print "---C--- = $c\n";
            &draw_domain;
        }
        ###############################
    }
}

sub draw_domain {

    for ( my $i = 0 ; $i <= $#domains ; $i++ ) {

        my @temp = split( /\t/, $domains[$i] );

        if ( $accession eq $temp[0] ) {
            print "$accession = $temp[0]\n";

            my @temp2 = split( /\~/, $temp[2] );
            for ( my $j = 0 ; $j <= $#temp2 ; $j++ ) {

                $var = $j * 3 + 3;

                $domain       = $temp2[$j];
                $domain_start = $temp[$var];
                $domain_end   = $temp[ $var + 1 ];
                $domain_type  = $temp[ $var + 2 ];
                chomp($domain_type)
                  ;    # Remove end of line, causing last domain to be ignored.

                print
"\t$#temp2 - $j - $domain\t$domain_start\t$domain_end\t$domain_type\n";

                &domain_colour($domain);
                if ( $domain_type eq ".." ) {
                    &domain_start_open;
                    &domain_middle;
                    &domain_end_open;
                }
                elsif ( $domain_type eq "[." ) {
                    &domain_start_closed;
                    &domain_middle;
                    &domain_end_open;
                }
                elsif ( $domain_type eq ".]" ) {
                    &domain_start_open;
                    &domain_middle;
                    &domain_end_closed;
                }
                elsif ( $domain_type eq "[]" ) {
                    &domain_start_closed;
                    &domain_middle;
                    &domain_end_closed;
                }
            }
        }
    }
}

sub domain_end_open {

    # make a polygon
    $hpos = $left_pos + $domain_end - 12;

    $poly = new GD::Polygon;
    $poly->addPt( $hpos,     $bar_vpos );
    $poly->addPt( $hpos + 4, $bar_vpos - 4 );
    $poly->addPt( $hpos,     $bar_vpos - 4 );
    $poly->addPt( $hpos,     $bar_vpos );
    $poly->addPt( $hpos,     $bar_vpos + 8 );
    $poly->addPt( $hpos + 4, $bar_vpos + 4 );
    $poly->addPt( $hpos,     $bar_vpos );

    # draw the polygon, filling it with a color
    $im->filledPolygon( $poly, $dom_colour );
}

sub domain_start_open {

    # make a polygon
    $hpos = $left_pos + $domain_start + 9;

    $poly = new GD::Polygon;
    $poly->addPt( $hpos + 4, $bar_vpos - 4 );    #40
    $poly->addPt( $hpos,     $bar_vpos - 4 );    #00
    $poly->addPt( $hpos + 4, $bar_vpos - 4 );    #44
    $poly->addPt( $hpos,     $bar_vpos );        #08
    $poly->addPt( $hpos + 4, $bar_vpos + 4 );    #412
    $poly->addPt( $hpos,     $bar_vpos + 8 );    #016
    $poly->addPt( $hpos + 4, $bar_vpos + 8 );    #416

    # draw the polygon, filling it with a color
    $im->filledPolygon( $poly, $dom_colour );

}

sub domain_end_closed {

    $im->setAntiAliased($dom_colour);
    $im->filledArc(
        $left_pos + $domain_end - 13,
        $bar_vpos + 2,
        13, 13, 270, 90, $dom_colour, gdArc
    );
}

sub domain_start_closed {

    $im->filledArc(
        $left_pos + $domain_start + 14,
        $bar_vpos + 2,
        13, 13, 90, 270, $dom_colour, gdArc
    );
    $im->setAntiAliased($dom_colour);
}

sub domain_middle {
    $im->filledRectangle(
        $left_pos + $domain_start + 14,
        $bar_vpos - 4,
        $left_pos + $domain_end - 13,
        $bar_vpos + 8, $dom_colour
    );    # + 14/-13 for caps
    $domain_name_length = length($domain) * 3;
    $dark_grey = $im->colorAllocate( 83, 83, 83 );
    if ( $domain_name_length < ( $domain_end - $domain_start ) ) {
        $im->string( gdSmallFont,
            $left_pos +
              $domain_start +
              ( ( ( $domain_end - $domain_start ) / 2 ) - $domain_name_length )
              - 1,
            $bar_vpos - 4 - 1,
            "$domain",
            $dark_grey
        );    # shadow
        $im->string( gdSmallFont,
            $left_pos +
              $domain_start +
              ( ( ( $domain_end - $domain_start ) / 2 ) - $domain_name_length ),
            $bar_vpos - 4,
            "$domain",
            $white
        );
    }
}

sub domain_colour {
    my $temp = shift;

    # Add the ascii values of the domain text to generate
    # colours specific to each name...
    my @domain = unpack( "C*", "$temp" );
    $sum = eval { join '+', @domain };
    $sum = $sum / 2;

    $r = int( ( ( $sum * 10 ) - 40 ) - 240 );
    $g = int( ( ( $sum * 10 ) - 40 ) - 120 );
    $b = int( ( ( $sum * 10 ) - 40 ) - 80 );
    $a = "45";

    #$dom_colour = $im->colorAllocate( $r, $g, $b );
    $dom_colour = $im->colorAllocateAlpha( $r, $g, $b, $a );

    return $dom_colour;
}

sub scale_bars {
    $scale_quotient = int( $bar_length / 100 );

    #Large Scale Bars
    for ( $a = 0 ; $a <= $scale_quotient ; $a++ ) {
        $scale = $a * 100;
        $im->string( gdSmallFont, $scale + $left_pos,
            $font_vpos, "$scale", $black );
        $im->rectangle(
            $scale + $left_pos,
            $bar_vpos - 5,
            $scale + $left_pos + 1,
            $bar_vpos - 1, $black
        );
    }
}

sub length_bars {

    $end_padded = $end + $left_pos;

    # Start
    $im->string( gdSmallFont, $left_pos, $font_vpos + 30, "$start", $red );
    $im->rectangle(
        $left_pos,
        $bar_vpos + 6,
        $left_pos + 1,
        $bar_vpos + 10, $red
    );

    #End
    $im->string( gdSmallFont, $end_padded, $font_vpos + 30, "$end_scale",
        $red );
    $im->rectangle(
        $end_padded - 1,
        $bar_vpos + 6,
        $end_padded, $bar_vpos + 10, $red
    );

#Length End
#$im->string( gdSmallFont, $end_padded, $font_vpos, "$bar_length", $black );
#$im->rectangle( $end_padded - 1, $bar_vpos - 5, $end_padded, $bar_vpos - 1, $black );

}

sub bar {
    $colour = shift;
    $im->filledRectangle(
        $left_pos, $bar_vpos,
        $bar_length + $left_pos,
        $bar_vpos + 5, $colour
    );
}

sub accession_name {
    my $colour = shift;
    if ( $actual_length == "" ) {
        $im->string( gdSmallFont, $padding_left - 50,
            $accession_vpos, "$accession", $colour );
        $im->string( gdTinyFont,
            $padding_left - 50,
            $accession_vpos + 12,
            "L=$bar_length", $colour
        );
    }
    else {
        $im->string( gdSmallFont, $padding_left - 50,
            $accession_vpos, "$accession", $colour );
        $im->string( gdTinyFont,
            $padding_left - 50,
            $accession_vpos + 20,
            "ML=$bar_length", $colour
        );
        $im->string( gdTinyFont,
            $padding_left - 50,
            $accession_vpos + 12,
            "L =$actual_length", $colour
        );
    }
}

sub watermark {

    $dark_grey  = $im->colorAllocate( 83,  83,  83 );
    $light_grey = $im->colorAllocate( 192, 192, 192 );
    $dark_red   = $im->colorAllocate( 146, 84,  83 );
    $dark_green = $im->colorAllocate( 129, 158, 107 );

    $darker_grey  = $im->colorAllocate( 52,  50,  51 );
    $darker_red   = $im->colorAllocate( 123, 59,  59 );
    $purple       = $im->colorAllocate( 96,  84,  112 );
    $dark_blue    = $im->colorAllocate( 33,  135, 204 );
    $darker_green = $im->colorAllocate( 43,  166, 22 );
    $pink         = $im->colorAllocate( 255, 0,   255 );

    $logo_left_x = $image_length - $padding_right;
    $logo_left_y = $image_height - 20;

    $fdfBLAST = "fdfBLAST v$VERSION";

    $im->string( gdSmallFont,
        $image_length - 100,
        $logo_left_y + 7,
        "$fdfBLAST", $light_grey
    );

    $im->string( gdSmallFont,
        $image_length - 180,
        $logo_left_y + 7,
        "Ratio = $ratio", $light_grey
    );

    $im->string( gdTinyFont,
        $image_length - 425,
        $logo_left_y + 10,
        "Length Match %", $light_grey
    );

    $im->filledRectangle(
        $image_length - 350,
        $logo_left_y + 10,
        $image_length - 320,
        $logo_left_y + 17,
        $darker_grey
    );
    $im->filledRectangle(
        $image_length - 320,
        $logo_left_y + 10,
        $image_length - 290,
        $logo_left_y + 17,
        $darker_red
    );
    $im->filledRectangle(
        $image_length - 290,
        $logo_left_y + 10,
        $image_length - 260,
        $logo_left_y + 17,
        $purple
    );
    $im->filledRectangle(
        $image_length - 260,
        $logo_left_y + 10,
        $image_length - 230,
        $logo_left_y + 17,
        $dark_blue
    );
    $im->filledRectangle(
        $image_length - 230,
        $logo_left_y + 10,
        $image_length - 200,
        $logo_left_y + 17,
        $darker_green
    );
    $im->filledRectangle(
        $image_length - 200,
        $logo_left_y + 10,
        $image_length - 190,
        $logo_left_y + 17, $pink
    );

    $im->string( gdTinyFont,
        $image_length - 343,
        $logo_left_y + 10,
        "<40", $light_grey
    );
    $im->string( gdTinyFont,
        $image_length - 317,
        $logo_left_y + 10,
        "40-60", $light_grey
    );
    $im->string( gdTinyFont,
        $image_length - 287,
        $logo_left_y + 10,
        "60-70", $light_grey
    );
    $im->string( gdTinyFont,
        $image_length - 257,
        $logo_left_y + 10,
        "70-80", $light_grey
    );
    $im->string( gdTinyFont,
        $image_length - 229,
        $logo_left_y + 10,
        "80-100", $light_grey
    );
    $im->string( gdTinyFont,
        $image_length - 197,
        $logo_left_y + 10,
        "!", $light_grey
    );
}

sub get_colour {

    if ( $colour_ratio le "0.4" ) {

        # Black
        $r = int( 42 - ( 40 - ( 10 * ( 10 * 1 ) ) ) );
        $g = int( 40 - ( 40 - ( 10 * ( 10 * 1 ) ) ) );
        $b = int( 41 - ( 40 - ( 10 * ( 10 * 1 ) ) ) );

        #print "Black $r, $g, $b\n";
        $ratio_colour = $im->colorAllocate( $r, $g, $b );
    }
    elsif ( $colour_ratio gt "0.4" && $colour_ratio le "0.6" ) {

        # Red
        $r = int( 123 - ( 60 - ( 10 * ( 10 * 1 ) ) ) );
        $g = int( 59 -  ( 60 - ( 10 * ( 10 * 1 ) ) ) );
        $b = int( 59 -  ( 60 - ( 10 * ( 10 * 1 ) ) ) );

        #print "Red $r, $g, $b\n";
        $ratio_colour = $im->colorAllocate( $r, $g, $b );
    }
    elsif ( $colour_ratio gt "0.6" && $colour_ratio le "0.7" ) {

        # Purple
        $r = int( 96 -  ( 70 - ( 10 * ( 10 * 1 ) ) ) );
        $g = int( 84 -  ( 70 - ( 10 * ( 10 * 1 ) ) ) );
        $b = int( 112 - ( 70 - ( 10 * ( 10 * 1 ) ) ) );

        #print "Purple $r, $g, $b\n";
        $ratio_colour = $im->colorAllocate( $r, $g, $b );
    }
    elsif ( $colour_ratio gt "0.7" && $colour_ratio le "0.8" ) {

        # Blue
        $r = int( 33 -  ( 80 - ( 10 * ( 10 * 1 ) ) ) );
        $g = int( 135 - ( 80 - ( 10 * ( 10 * 1 ) ) ) );
        $b = int( 204 - ( 80 - ( 10 * ( 10 * 1 ) ) ) );

        #print "Blue $r, $g, $b\n";
        $ratio_colour = $im->colorAllocate( $r, $g, $b );
    }
    elsif ( $colour_ratio gt "0.8" && $colour_ratio le "1.0" ) {

        # Green
        $r = int( 43 -  ( 100 - ( 10 * ( 10 * 1 ) ) ) );
        $g = int( 166 - ( 100 - ( 10 * ( 10 * 1 ) ) ) );
        $b = int( 22 -  ( 100 - ( 10 * ( 10 * 1 ) ) ) );

        #print "Green $r, $g, $b\n";
        $ratio_colour = $im->colorAllocate( $r, $g, $b );
    }
    else {

        # Bright Pink
        $ratio_colour = $im->colorAllocate( 255, 0, 255 );
    }
    return $ratio_colour;
}

sub duplication_scan {

    $verbose = 0;
    #####
    print "\nPerforming Duplication Scans\n";
    #####

    @duplication_filenames = glob "$gene_hits_duplications/*fixed.csv";
    $duplication_number    = @duplication_filenames;

    for ( $i = 0 ; $i < $duplication_number ; $i++ ) {
        open my $duplication_fh, "<$duplication_filenames[$i]";
        #####
        #print "Opening $duplication_filenames[$i]...\n" if $verbose == 1;
        #####
        $input_file = fileparse( $duplication_filenames[$i] );
        print "$input_file";
        $pass = 0;
        while (<$duplication_fh>) {
            my ($line)        = $_;
            my ($line_number) = $.;

            #print "\nLN\:$line_number\n";
            chomp($line);

            #$pass = 0;
            if ( $line ne "" ) {

                push( @duplication_array, "$line" );
            }
            elsif ( $line eq "" && $pass == 0 ) {
                $duplicated_genes = $line_number;
                $pass++;
            }
        }
        close($duplication_fh);
        $duplication_array_length = @duplication_array;

        for ( $k = 0 ; $k < $duplicated_genes - 1 ; $k++ ) {
            $count        = 0;
            @temp         = split( /\;/, $duplication_array[$k] );
            $query_gene   = $temp[1];
            $query_length = $temp[2];
            #####
            print "$k\t$query_gene\t" if $verbose == 1;
            #####
            for ( $j = 0 ; $j < $duplication_array_length ; $j++ ) {
                @temp2          = split( /\;/, $duplication_array[$j] );
                $subject_gene   = $temp2[3];
                $subject_length = $temp2[4];
                $subject_evalue = $temp[5];
                $subject_ident  = $temp[6];
                $query_count    = $temp2[1];
                #####
                print "$subject_gene\n" if $verbose == 1;
                #####
                if ( $query_gene eq $query_count ) {
                    $count++;
                    push( @subject,
"$subject_gene,$subject_length,$subject_evalue,$subject_ident,"
                    );
                }
            }
            $output_file = fileparse( $duplication_filenames[$i] );
            #####
#print "Output to $gene_hits_duplications/$count\_hits_$output_file\n" if $verbose == 1;
            #####
            $hits = $count + 1;
            open my $output_fh,
              ">>$gene_hits_duplications/$hits\_hits_$output_file";
            print $output_fh "$query_gene,$query_length,";
            foreach my $value (@subject) {
                print $output_fh "$value";
            }
            print $output_fh "\n";
            #####
            print "Hits = $count\n" if $verbose == 1;
            #####
            @subject = ();
            print ".";
        }

        @duplication_array = ();
        print "\n";
        close($output_fh);
    }
}

sub scan_subject_initial {
    $verbose       = 0;
    @scan_array    = ();
    @variables     = @_;
    $query_hits    = $variables[0];
    $init_file     = $variables[1];
    $init_file_two = $variables[2];

    #####$init_file =~ m/(.*?\/lookup\/)(.*?\..*?)(\.csv)/;
    #####$init_file = $2;
    ( $init_file, $dir, $ext ) = fileparse( $init_file, '\..*' );
    #####$init_file_two =~ m/(.*?\/lookup\/)(.*?\..*?)(\.csv)/;
    #####$init_file_two     = $2;
    ( $init_file_two, $dir, $ext ) = fileparse( $init_file_two, '\..*' );
    $subject_accession = $variables[3];
    $loop_position     = $variables[4];

    #####
    # print "Opening $gene_hits_initial/$init_file_two\_$init_file\.txt\n";
    #####

    open my $in_fh, "<$gene_hits_initial/$init_file_two\_$init_file\.csv";

    while (<$in_fh>) {
        my ($line) = $_;
        chomp($line);
        @temp = split( /,/, $line );
        #####
        print "if $temp[1] == $subject_accession\n" if $verbose == 1;
        #####
        if ( $temp[1] eq $subject_accession ) {
            #####
            # print "$query_hits != $temp[0]\n";
            print "\# Reciprocal Hits = $temp[0]\n" if $verbose == 1;
            ######
            for ( $r = 0 ; $r < $temp[0] ; $r++ ) {
                $loop_position = $r + 1;
                $location      = ( 6 * $loop_position ) - 3;
                #####
                print
                  "Loc = $location\tQH = $query_hits \<\-\> $temp[$location]\n"
                  if $verbose == 1;
                #####
                if ( $query_hits eq $temp[$location] ) {

                    push( @scan_array, "recip" );
                    push( @scan_array,
                        "$temp[(6 * $loop_position) - 3]",
                        "$temp[(6 * $loop_position) - 2]",
                        "$temp[(6 * $loop_position) - 1]",
                        "$temp[(6 * $loop_position)]",
                        "$temp[(6 * $loop_position) + 1]",
                        "$temp[(6 * $loop_position) + 2]" );
                }
                else {
                    push( @scan_array, "norecip" );

                }
            }
        }
    }
    return @scan_array;
}

sub get_line_number_start {

    # This method opens a file and gets the last
    # line number and returns the value
    # This will probably be exponentially slow. Hrmmm.
    $last_line_number = 0;

    @variables   = @_;
    $file_handle = $variables[0];

    open my $file, "<$file_handle";

    while (<$file>) {

        #$last_line_number = $. if eof;
        #####$last_line        = $_ if eof;
        $last_accession = $_ if eof;
        @temp             = split( /;/, $last_accession );
        $last_accession   = $temp[1];
        $last_line_number = $temp[0];
    }
    close($file);
    #####chomp($last_line);
    ######$last_line_number = $last_line_number + 1;
    return "$last_line_number";
    return "$last_accession";
}

sub get_accession {
    $verbose = 0;

    # First we need to open the corresponding file in the initial dir
    # To do that we need the file name (no ext or dir) of the lookup file
    @variables          = @_;
    $init_file          = $variables[0];
    $init_file_two      = $variables[1];
    $lookup_line_number = $variables[2];

    #####$init_file =~ m/(.*?\/lookup\/)(.*?\..*?)(\.csv)/;
    #####$init_file = $2;
    ( $init_file, $dir, $ext ) = fileparse( $init_file, '\..*' );

    #####$init_file_two =~ m/(.*?\/lookup\/)(.*?\..*?)(\.csv)/;
    #####$init_file_two = $2;
    ( $init_file_two, $dir, $ext ) = fileparse( $init_file_two, '\..*' );

    # Set initial directory and open file

    open my $access_fh, "<$gene_hits_initial\/$init_file\_$init_file_two\.csv";
    #####
# print "Opening $gene_hits_initial\/$init_file\_$init_file_two\.txt\n" if $verbose == 1;
    #####

    # While file, read each line and compare line numbers
    #
    while (<$access_fh>) {
        my ($line) = $_;
        chomp($line);
        my ($initial_line_number) = $.;

        # When line numbers match, assign accession and then return
        if ( $lookup_line_number == $initial_line_number ) {

            # Split the line on , and retrieve accession (pos 1)
            @temp             = split( /,/, $line );
            $accession        = $temp[1];
            $accession_length = $temp[2];
            @accession_array  = ( "$accession", "$accession_length" );
        }

        # Else do nothing
    }
    close($access_fh);
    return $accession_array;
}

sub get_hit_info {

    $verbose      = 0;
    @return_array = ();

    # First we need to open the corresponding file in the initial dir
    # To do that we need the file name (no ext or dir) of the lookup file
    @variables          = @_;
    $init_file          = $variables[0];
    $init_file_two      = $variables[1];
    $lookup_line_number = $variables[2];
    $location           = $variables[3];
    $location           = $location + 1;

    ( $init_file, $dir, $ext ) = fileparse( $init_file, '\..*' );

    ( $init_file_two, $dir, $ext ) = fileparse( $init_file_two, '\..*' );

    # Set initial directory and open file

    open my $access_fh, "<$gene_hits_initial\/$init_file\_$init_file_two\.csv";

    # While file, read each line and compare line numbers
    #
    while (<$access_fh>) {
        my ($line) = $_;
        chomp($line);
        my ($initial_line_number) = $.;

        # When line numbers match, assign accession and then return
        if ( $lookup_line_number == $initial_line_number ) {

            if ( $location == 1 ) {

                # Split the line on ',' and retrieve accession (pos 1)
                @temp                  = split( /,/, $line );
                $hit_accession         = $temp[3];
                $length                = $temp[4];
                $evalue                = $temp[5];
                $percent_identity      = $temp[6];
                $query_hit_range_start = $temp[7];
                $query_hit_range_end   = $temp[8];
                @return_array          = (
                    "$hit_accession",         "$length",
                    "$evalue",                "$percent_identity",
                    "$query_hit_range_start", "$query_hit_range_end"
                );
            }
            else {

                $location              = ( $location * 6 );
                @temp                  = split( /,/, $line );
                $hit_accession         = $temp[ $location - 3 ];
                $length                = $temp[ $location - 2 ];
                $evalue                = $temp[ $location - 1 ];
                $percent_identity      = $temp[$location];
                $query_hit_range_start = $temp[ $location + 1 ];
                $query_hit_range_end   = $temp[ $location + 2 ];
                @return_array          = (
                    "$hit_accession",         "$length",
                    "$evalue",                "$percent_identity",
                    "$query_hit_range_start", "$query_hit_range_end"
                );
            }
        }
    }
    close($access_fh);
    return @return_array;
}
