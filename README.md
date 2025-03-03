![DOI](https://zenodo.org/badge/3666/guyleonard/fdfBLAST.png)

1. Prerequisites / Dependencies
2. Initial Setup
3. Running fdFBLAST
4. PFAM Domains
5. Installing Dependencies
6. Example Genomes
7. Other Notes
8. Accessory Scripts

Q. Why are there several versions of the script?

A.fdfBLAST.pl - current working version
fdfBLAST_1.20120702.pl - current working version with version number
fdfBLAST_1.20120704.pl - experimental update, consider buggy
fdfBLAST_old.pl - This is the version used for manuscript preparation, deposited here for reproducibility/openness

Prerequisites
--------------
BLAST Legacy Executables

PERL Modules
  Math::BigFloat
  BioPerl
  GD

HMMER3 and PFAM-A

Optional
  GD::SVG

Initial Setup
-------------
Don't panic. There is generally very little for you to do, everything should be explained on screen and you have only a few options onscreen in menus.

1) Download the script from https://github.com/guyleonard/fdfBLAST and create a folder called 'fdf' or whatever you prefer and put the script there.

2) Download your genomes and place them in a folder under a directory callled 'genomes'.
	This folder can be anywhere, e.g. fdf/genomes/your_genomes, or you can enter the directory manually when running fdf or edit the script under the Sub set_genome_dir.
2a) Make sure the files all have a .fas extension.
	It helps to go with a consistent naming scheme so for example 'Homo sapiens' I would rename to 'Hsapiens.fas' or similar.
2b) fdfBLAST is currently unable to deal with ncbi style or other complex fasta definition lines.
	Therefore, you need to either use the tool REFGEN to alter them or manually edit them with 'sed' or a regex enabled text editor, so that you have a unique (across all genomes) accession number and some way of identifying which genome it came from e.g. ">12345Hs"

3) Install BLAST - Skip if you already have this installed.
	The executables can be global or in a sub folder of the fdf directory called 'blast'

4) Create a directory called 'run' under your fdf directory.


Running fdFBLAST
----------------

fdfBLAST will present you with a series of options on the command line via a series of menus. Follow them as they appear.

i) $ perl fdfBLAST.pl
ii) Answer all the menu options from the options given

1) Step 1 and 2 of the program deal with running formatdb and blastp.
	This will only need to be run once, per set of genomes.
	It checks for previous blast comparisons, so you cann add a new genome to an analysis and only have to perform the new set of analyses. 

2) Steps 3 and 4 are the main parts of fdfBLAST and should run without any interaction
	There are a couple of options to change if you wish but I would leave them at the defaults for the first run.

3) Step 5 is the stage that outputs fusion candidates.
	Initially it will output alignment diagrams (in the differentials folder).
	In order to map PFAM domains you will need to run it again, you can leave the file entry blank initially.
	Go to PFAM Domains, Step 1. Then re-run this step.

PFAM Domains
------------

You can do these operations anywhere, as long as your put the final resulting file back into your particluar run directory.
It involves using a few command line tools, to process the data.

1) Use the program hmmpress on the database you downloaded to get it ready for hmmscan.
	$ hmmpress Pfam-A.hmm

2) In your 'run' folder there will be two files: composite_list.csv and split_list.csv
	These represent the gene accessions for the total of the fusion output
	Composite = potential fused domains, split = individual domains in potential fusion

2a) You need to use the command line program "sort" to sort each file.
	$ sort composite_list.csv >composite_list_sorted.csv

2b) Then use the program "uniq -u" to remove any duplicates.
	$ uniq -u composite_list_sorted.csv >composite_list_sorted_uniq.csv

2c) Then use "cat" to put them together.
	$ cat composite_list_sorted_uniq.csv split_list_sorted_uniq.csv >all_domains.csv

3) Use the accessory script "get_non_ncbi_seqs.pl" to extract all the sequencess based on the accession in the file you created in step 2c
	You should put this script and the file from 2c in your current run directory
	You will need to make some ammendments to this file in a text editor.
	$genome_set = the name of your genome folder, just the name no directory information
	$genome_directory = the directory path goes here, no trailing slash
	$run_ID = your current run ID

3a) The script will produce a set of .fas sequence files named after each of the genomes.
	Use "cat" to put them all together and then you will need to run PFAM or CDD on them.
	e.g. $ cat genomeA.fas genomeB.fas genomeC.fas >all_domain_seqs.fas

4) Use hmmscan to identify the domains associated with each sequence.
	$ hmmscan --cpu 6 --noali all_domain_seqs.fas >output.out

5) Use the other accessory script hmm3scan_parser.pl to create a modified version of the above output.
	$ perl hmm3scan_parser.pl -i <infile> -o <outfile>

6) Go to fdfBLAST 3) above with the output file.

NB - You can also use the CDD database via RPSBLAST. However, I have not included any scripts for this to work. As long as you can get the output from RPSBLAST and CDD into the same tab delimited format that is generated by hmm3scan_parser.pl then fdfBLAST will be able to read the file.

e.g. the file should look like this
Accession	Length	Domain1~Domain2~etc	Start	End		Start	End
201166		777	Pkinase~Pkinase		73	242	..	73	242	[.
201490		300	Caveolin~FA_hydroxylase	92	141	..	151	248	[]

Installing Dependencies
-----------------------
Legacy BLAST can be found here http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

Perl Libraries
You can install them on the command line thus:
    sudo perl -MCPAN -e shell
followed by 
    install Math::BigFloat
You can search for other packages using 
    i /package::name/
There may be dependencies I have forgotten about related to GD that you need to install as part of debian/linux distro.

PFAM can be found here ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam26.0/Pfam-A.hmm.gz

HMMER can be found here http://hmmer.janelia.org/

Example Genomes
---------------

Other Notes
-----------
The fdfBLAST program was developed under a free and open-source Ubuntu Linux environment (http://www.ubuntu.com), although it should be capable of running on any UNIX based operating system (including Apple OSX) where the following programs are available for your system. Perl v5.10.0 (http://www.perl.org) was used to write the main portions of the script. The legacy NCBI BLAST executables (version 2.2.21) are used to perform the BLASTp and RPS-BLAST {Marchler-Bauer, 2005} searches and to prepare the proteome databases ready for use.

A new version of the local BLAST executables is now available, called BLAST+ executables, but they were not released at the inception of the program – there is no reason why fdfBLAST could not be adapted to utilise these in a future release – although there is no pressing need to change the script currently since no particular advantage is conferred. Later portions of the script employ the use of the program HMMER3 (http://hmmer.org) and a copy of the Sanger Institute’s PFAM database (Pfam 24.0 - October 2009; Finn, 2010). and the NCBI Conserved Domain Database (CDD version 2.2.2).

Accessorry Scripts
------------------
hmm3scan_parser.pl -- Copyright Bill Wickstead.
