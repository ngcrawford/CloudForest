#!/usr/bin/perl -w
#-------------------------------------------------------------------------#
# This is
 $program_name =    "mraic_mod.pl";
 $version =         "1.4.4";
 $lastChanges =	    "05/12/2009 01:16:31 AM CEST"; # JN
# made by
 $author =          "Johan Nylander";
#
# Email:            jnylander at users dot sorceforge dot net
#
# Description:	    Script for calculating Akaike weights, AIC, AIc,
#                   and BIC values for nucleotide substitution models.
#                   Likelihood scores are calculated using PHYML
#                   (Guindon and Gascuel, 2003).
#
# Usage:            Interactivly or by passing a filename: 'mraic.pl filename'
#                   Input data is DNA characters in PHYML format.
#
# Requirements:     PHYML v. 3.0. (http://www.atgc-montpellier.fr/phyml/) needs
#                   to be installed!
 $PHYML = "phyml";	#  Optionally, specify the full PATH to phyml here
#
# Version history:  Version 1.0 Oct 2004.
#                   Version 1.1 Nov 2004: Branch lengths are now
#                   included in total number of parameters in calculations.
#                   AICc values are only printed in Nparams/Nchars < 40.
#                   Version 1.2 Jan 2005: The script can now be run on
#                   Windows machines, at least on XP/NT (used unlink
#                   instead of rm). [Thanks to Johannes Lundberg].
#                   The script now allows interleaved format as well
#                   as sequential input. [J.L]
#                   The script does not sleep the system by default.
#                   Uncomment the line containing "sleep" if you want
#                   the progress to be printed.
#                   Version 1.3 Jan 2005: Bug fix. The commands
#                   for two of the models invoked by the -modeltest flag
#                   were incorrect. [Thanks to Pär Larsson]
#                   Version 1.3.1 Feb 2005: Bug fix. MbBlock for K2PIG
#                   now printed correctly. [Thanks to Torsten Eriksson]
#                   Version 1.3.2 Mar 2005: Bug fix. MbBlock for SYMI
#                   now printed correctly.
#                   Version 1.3.3 Mar 2005: Minor change. Made the code
#                   slightly more generous to differences in the input
#                   format. [Thanks to Cymon Cox]
#                   Version 1.4.0 Mar 2005: Bug fixes. AICc calculations
#                   failed (giving negative AICc) when there were more
#                   parameters than data. For those cases, branch lengths
#                   are now not included in number of parameters and warnings
#                   are given.
#                   Output tables now sorted correctly. [Thanks to Cymon Cox]
#                   Version 1.4.1 Feb 2005: Added some error checking for
#                   handling failure in PHYML optimization.
#                   Version 1.4.2 June 2006: Minor changes. Added some code
#                   for handling the case where previous phyml output is in
#                   the working directory. Thanks to Per and Ivar Erixon.
#                   Version 1.4.3 April 2007: The output from PHYML (>v.2.4.4)
#                   changed. The likelihood is now read from the stat file
#                   instead of the lk file. Fix by Olaf Bininda-Emonds.
#                   Version 1.4.4 May 2009. The output from PHYML (>v.3.)
#                   changed. Commands adjusted to work with the new PHYML
#                   version. Thanks to Liliana Davalos.
#                   
#--------------------------------------------------------------------------#

use Getopt::Long;
use File::Basename;
use File::Spec;

# Globals
$interactive              = NO;                # No interactive by default
$debugLevel               = 0;                 # Set value equal to '1' for debug information (written to file)
$useModeltest             = FALSE;             # Do not compare 56 model by default  
$nmodels                  = 24;                # Default is to compare the 24 MrBayes models
$useBrlens                = TRUE;              # Include branch lengths in total number of parameters
$BAratio                  = 40;                # Burnham Anderson ratio nchar/nparams
$changedNparams           = FALSE;             # If the BA ratio is hit, the script changes the number of
                                               # parameters to not include branch lengths (and sets
                                               # $changedNparams = TRUE)
$mraic_outfile_ending     = '.MrAIC.txt';      # File ending for the output file
$phyml_tree_file_ending   = '_phyml_tree.txt'; # Note: file endings might change with different phyml version!
$phyml_stat_file_ending   = '_phyml_stats.txt';
$command_file_ending      = '_phyml_command_file';



# Start

# Check if debug is on
if ($debugLevel > 0) {
	$debugFile = "debugfile.MrAIC.txt";
	open DEBUG, '>', "$debugFile";
}

my ($help, $infilename, $output, $mt);
 
#-- prints usage if no command line parameters are passed or there is an unknown
#   parameter or help option is passed
usage() if ( @ARGV < 1 or ! GetOptions('help|h' => \$help, 'infile=s' => \$infilename, 'output_dir=s' => \$output, 'modeltest' => \$mt) or defined $help );
         
if ($mt) {
    $useModeltest = TRUE;
    $nmodels = 56;
}

if ($help) {
    helpInformation();
    exit(0);
}

if (! $output) {
    $output = '';
}
 
sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "usage: modeltest [--infile INFILE] [--output_dir DIR] [--modeltest] [--help|-h]\n";
  exit;
}

# Handle arguments or get necessary info before running
if ($infilename) {

    # Check filenames to prevent overwriting
    check_pers_files($infilename);

    # Get nchar and ntax
    getDimensions();

    # Guess input format (interleaved/sequential)
    $sequential = guess_input_format($infilename, $ntax, $nchar);
    
    # Test if number of parameters is larger than the data size
    if ($maxNumParameters > $nchar) {
        $useBrlens = FALSE;
        $maxNumParametersBrLens = $maxNumParameters;
        $maxNumParameters = 10; # 10: Number of params in GTRIG
        $changedNparams = TRUE;
    }
}
else {

    # If no file name on command line. Perform a more interactive run

    $interactive = YES;
    printf "\n[%s version %s by %s.]\n",$program_name, $version, $author;
    print "[At any time, type h for help of q to quit!]\n\n";
    if ($useModeltest) {
        do {
            print "Enter name of infile: ";
            chomp ($infilename = <STDIN>);
        }
        until ($infilename ne '');
        if ($infilename eq 'h' or $infilename eq 'H') {
            helpInformation();
            exit(0);
        }
        elsif ($infilename eq 'q' or $infilename eq 'Q') {
            exit(0);
        }
        elsif (check_pers_files($infilename)) {
            # check if old files exists
        }
        else {
            die "\a\nERROR!\n\nCan't find file '$infilename'.\n\n" unless -e $infilename;
        }
    }

	# Get model set
	do {
		print "Use only MrBayes models? (y/n) ";
		chomp ($_ = <STDIN>);
	}
	until ((/^[\s]*y/i) or (/^[\s]*n/i) or (/^[\s]*q/i) or (/^[\s]*h/i));
	if (/^[\s]*q/i) {
		exit(0);
	}
	elsif (/^[\s]*h/i) {
		helpInformation();
		exit(0);
	}
	elsif (/^[\s]*n/i) {
		$useModeltest = TRUE;
		$nmodels = 56;
	}

	# Get ntax and nchar
	getDimensions();

    # Guess input format (sequential/interleaved)
    $sequential = guess_input_format($infilename, $ntax, $nchar);

	# Confirm input
	print  "\nConfirming input\n";
	print  "\tInfile: $infilename\n";
    print  "\tSequences seems to be in";
	if ($sequential) {
        print " sequential format.\n";
    }
    else {
        print " interleaved format.\n";
    }
	print  "\tNumber of taxa: $ntax\n";
	print  "\tNumber of characters: $nchar\n";
	print  "\tMax. num. of parameters (incl. branch lengths): $maxNumParameters\n";
    if ($debugLevel > 0) { # DEBUG 
		printf DEBUG "\nin run() ntax: $ntax\n";
		printf DEBUG "\nin run() nchar: $nchar\n";
		printf DEBUG "\nin run() maxNumParameters: $maxNumParameters\n";
		$ratio_n_over_K = ($nchar/$maxNumParameters);
		printf DEBUG "\nin run() ratio_n_over_K: $ratio_n_over_K\n";
	}
	if ($maxNumParameters > $nchar) {
		$maxNumParametersBrLens = $maxNumParameters;
		$maxNumParameters = 10;
		$useBrlens = FALSE;
		$changedNparams = TRUE;
		$estMinBASize = ($BAratio * $maxNumParametersBrLens);
		$estMinSize = ($maxNumParametersBrLens);
		print  "\n\tWARNING:\n";
		print  "\tNumber of parameters exceeds the number of columns\n";
		print  "\tin the matrix! You need more sequence data!\n";
		printf "\t(For %g terminal branches you probably need more\n", $ntax;
		printf "\tthan %g characters, and an approx. data size for\n", $estMinSize;
		printf "\thaving a nchar/nparam ratio over %g is about %g bp.)\n", $BAratio, $estMinBASize;
		print  "\n\tWill run the script anyway but will not include\n";
		print  "\tbranch lengths as number of parameters.\n";
		print  "\n\tNumber of parameters is now set to $maxNumParameters\n";
		printf "\t(the number of free parameters in the GTR+I+G model).\n";
		print  "\n\tContinuing fitting substitution models might\n";
		print  "\tnot be a good idea.\n\n";
		if ($debugLevel > 0) { # DEBUG
			printf DEBUG "\nin run() revised maxNumParameters: $maxNumParameters\n";
			$revisedRatio_n_over_K = ($nchar/$maxNumParameters);
			printf DEBUG "\nin run() revised ratio_n_over_K: $revisedRatio_n_over_K\n";
		}
	}
	printf "\tComparing %g models.\n", $nmodels;
	print "\nHit return to start! ";
	chomp ($_ = <STDIN>);
	if (/^[\s]*q/i) {
		exit(0);
	}
	elsif (/^[\s]*h/i) {
		helpInformation();
		exit(0);
	}
	print "Run!\n\n";
}

# Run the script
run();

# Print the result
MrAICprint();


# Sub routines

#---------------- guess_input_format() -------------------------#
# Try to see if file is interleaved or not
#
sub guess_input_format {

    my ($infile_name, $ntax, $nchar) = @_;
    my $max_line_length = 0;
    my $nlines = 0;

    open my $FI, "<", $infile_name or die "\a\nERROR!\n\nCould not open infile for reading: $!\n";
    while (<$FI>) {
        chomp();
        my $line = $_;
        next if (/^\s+$/);
        $nlines++;
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        $line =~ s/\s+/ /g;
        if(length($line) > $max_line_length) {
            $max_line_length = length($line);
        }
    }
    close($FI) or warn "Could not close file handle: $! \n";

    if ($max_line_length < $nchar) { # File is probably interleaved if
        return(0);                   # no line is as long or longer than
    }                                # number of characters.
    else {                           # Else it's probably sequential
        return(1);
    }

} # end of guess_input_format


#---------------- check_pers_files() ---------------------------#
# Check if there are files from previous runs to avoid overwriting
#
sub check_pers_files {

    my ($infile_name) = @_;
    my $tree_file     = "$infile_name" . "$phyml_tree_file_ending";
    my $stat_file     = "$infile_name" . "$phyml_stat_file_ending";
    my $mraic_file    = "$infile_name" . "$mraic_outfile_ending";
    my @check_files   = ($tree_file, $stat_file, $mraic_file);
    my @pers_files    = ();

    foreach (@check_files) {
        if (-e $_) {
            push @pers_files, $_;
        }
    }

    if (scalar @pers_files) {
        print "\a\nWarning!\n";
        print "There is already phyml (and/or MrAIC) output for '$infilename' in the working directory:\n";
        print "(@pers_files)\n";
        print "Phyml/MrAIC will overwrite unless you move/rename the file(s).\n\n";
        print "Quitting.\n\n";
        exit(0);
    }
    else {
        return 0;
    }
}


#---------------- run() ---------------------------#
# Run PHYML
#
sub run {

    # File names
    $phyml_command_file  = "$infilename" . "_mraic" . "$command_file_ending";
    $phyml_tree_txt_file = "$infilename" .  "$phyml_tree_file_ending";
    $phyml_stat_txt_file = "$infilename" .  "$phyml_stat_file_ending";

	# Initialize the hashes
    if ($useModeltest eq TRUE) {
        SetModeltestMODELnr();
        SetModeltestMODELcmd();
        SetModeltestMODELdf();
    }
    else {
        SetMbMODELnr();
        SetMbMODELcmd();
        SetMbMODELdf();
        SetMbMODELmbcmd();
    }

    # Run PHYML and collect the output for nmodels
    $k = 0;	# Variable for progress printing
    for ($i=0; $i<$nmodels; $i++) {

        # Create a command file for PHYML
        open COMMANDS, '>', $phyml_command_file or die "\a\nERROR!\n\nCan't open command file: $!\n";
            print COMMANDS "$infilename\n$MODELcmd{$i}";
        close COMMANDS;

        # Print progress. Not optimal, needs a 'sleep' to be visable, see below.
        $k = ($i + 1);
        printf STDOUT "\nAnalysis %g (of %g). Model %s\n", $k, $nmodels, $MODELnr{$i};
        print  STDOUT "--------------------------------------------------------------\n";
        if ($debugLevel > 0) { # DEBUG
            print  DEBUG "--------------------------------------------------------------\n";
            printf DEBUG "\nAnalysis %g (of %g). Model %s\n", $k, $nmodels, $MODELnr{$i};
            print  DEBUG "--------------------------------------------------------------\n";
        }
        #system "sleep 0,5";	# DEBUG: Needed to prevent PHYML owerwriting the display. Uncomment if needed.

        # Run PHYML through the shell
        system "$PHYML < $phyml_command_file";

        # Get the lnL in the *stats.txt file
        if (($statFileName) = glob($phyml_stat_txt_file)) {
            # check indirectly if PHYML is installed
        }
        else {
            unlink ("$phyml_command_file");
            die "\a\nERROR!\n\nCannot find output from PHYML.\nIs PHYML installed (and named \"phyml\")?\n\n";
        }

        open (STATFILE, '<', $statFileName) or die "\a\nERROR!\n\ncannot open $a for reading: $!";
            $tmpLnL = 999.999;
            while (<STATFILE>) {
                chomp;
                $tmpLnL = $1 if ($_ =~ m/Log-likelihood:\s+(.+)/i); # O.B-E.
            }
        close(STATFILE);
        if ($tmpLnL == 999.999) {
            die "\a\nERROR!\n\nCould not read the likelihood from PHYML.\n\n";
        }
        else {
            $MODELlnL{$i} = $tmpLnL;
            if ($debugLevel == 1) { # DEBUG
                printf DEBUG "\nin run() lnL: $MODELlnL{$i}\n";
            }
        }

        # Get the tree in the *tree.txt file
        ($treeFileName) = glob($phyml_tree_txt_file);
        open (TREEFILE, '<', $treeFileName) or die "\a\nERROR!\n\nCannot open $a for reading: $!";
        $tmpTree = "";
            while (<TREEFILE>) {
                chomp($tmpTree = $_);
            }
        close(TREEFILE);
        if ($tmpTree eq "") {
            die "\a\nERROR!\n\nCould not read the tree from PHYML.\n\n";
        }
        else {
            $MODELtree{$i} = $tmpTree;
        }

		# Remove phyml_*.txt and phyml_commands to avoid overwrite warnings by PHYML
        unlink ($phyml_command_file);
        unlink ($phyml_tree_txt_file);
        unlink ($phyml_stat_txt_file);

		# Calculate AIC
		AIC($i);

		# Calculate AICc
		AICc($i);

		# Calculate BIC
		BIC($i);
	}

	# Get the minAIC, model name, tree, and mbblock
	$AIC_min = $MODELAIC{0};
	$j=0;
	for ($i=0; $i<$nmodels; $i++) {
		if ($MODELAIC{$i} <= $AIC_min) {

			# Get the min AIC value
			$AIC_min = $MODELAIC{$i};
			$j = $i;

			# Get the name of min AIC model
			$minimumAICmodel = $MODELnr{$j};

			# Get the tree from the min AIC model
			$minimumAICtree = $MODELtree{$j};

			# Get the mbcmd for the min AIC model
			$minimumAICmbcmd = $MODELmbcmd{$j};
		}
		if ($debugLevel == 1) { # DEBUG
			printf DEBUG "\nin run minAICmodel: $minimumAICmodel minAIC: $AIC_min\n";
		}
	}

	# Get the minAICc, model name, tree, and mbblock
	$AICc_min = $MODELAICc{0};
	$j=0;
	for ($i=0; $i<$nmodels; $i++) {
		if ($MODELAICc{$i} <= $AICc_min) {

			# Get the min AICc value
			$AICc_min = $MODELAICc{$i};
			$j = $i;

			# Get the name of min AICc model
			$minimumAICcmodel = $MODELnr{$j};

			# Get the tree from the min AICc model
			$minimumAICctree = $MODELtree{$j};

			# Get the mbcmd for the min AICc model
			$minimumAICcmbcmd = $MODELmbcmd{$j};
		}
		if ($debugLevel == 1) { # DEBUG
			printf DEBUG "\nin run minAICcmodel: $minimumAICcmodel minAICc: $AICc_min\n";
		}
	}

	# Get the minBIC, model name, tree, and mbblock
	$BIC_min = $MODELBIC{0};
	$j=0;
	for ($i=0; $i<$nmodels; $i++) {
		if ($MODELBIC{$i} <= $BIC_min) {

			# Get the min BIC value
			$BIC_min = $MODELBIC{$i};
			$j = $i;

			# Get the name of min BIC model
			$minimumBICmodel = $MODELnr{$j};

			# Get the tree from the min BIC model
			$minimumBICtree = $MODELtree{$j};

			# Get the mbcmd for the min BIC model
			$minimumBICmbcmd = $MODELmbcmd{$j};
		}
		if ($debugLevel == 1) { # DEBUG
			printf DEBUG "\nin run minBICmodel: $minimumBICmodel minBIC: $BIC_min\n";
		}
	}

	# Calculate AIC weights
	$sumAICExp = 0.0;
	for ($i=0; $i<$nmodels; $i++) {
		$deltaAIC{$i} = $MODELAIC{$i} - $AIC_min;
		$sumAICExp += exp((-0.5) * ($deltaAIC{$i}));
	}
	if ($debugLevel == 1) { # DEBUG
		printf DEBUG "\nin run sumAICExp: $sumAICExp\n";
	}
	for ($i=0; $i<$nmodels; $i++) {
		$AkaikeAICWeight{$i} = exp((-0.5) * ($deltaAIC{$i})) / $sumAICExp;
	}
	# Calculate AICc weights
	$sumAICcExp = 0.0;
	for ($i=0; $i<$nmodels; $i++) {
		$deltaAICc{$i} = $MODELAICc{$i} - $AICc_min;
		$sumAICcExp += exp((-0.5) * ($deltaAICc{$i}));
	}
	if ($debugLevel == 1) { # DEBUG
		printf DEBUG "\nin run sumAICcExp: $sumAICcExp\n";
	}
	for ($i=0; $i<$nmodels; $i++) {
		$AkaikeAICcWeight{$i} = exp((-0.5) * ($deltaAICc{$i})) / $sumAICcExp;
	}
	# Calculate BIC weights
	$sumBICExp = 0.0;
	for ($i=0; $i<$nmodels; $i++) {
		$deltaBIC{$i} = $MODELBIC{$i} - $BIC_min;
		$sumBICExp += exp((-0.5) * ($deltaBIC{$i}));
	}
	if ($debugLevel == 1) { # DEBUG
		printf DEBUG "\nin run sumBICExp: $sumBICExp\n";
	}
	for ($i=0; $i<$nmodels; $i++) {
		$AkaikeBICWeight{$i} = exp((-0.5) * ($deltaBIC{$i})) / $sumBICExp;
	}

	# Print final notes
	print  "\n\n\n\n";
	print  "------------------------------------------------------------------\n";
	printf "Normal end of %s\n",$program_name;
	print  "------------------------------------------------------------------\n";
	printf "Results saved to file: %s.MrAIC.txt\n\n",$infilename;
	printf "Tree from min AIC model (%s) saved to file: %s.AIC-%s.tre\n",$minimumAICmodel, $infilename, $minimumAICmodel;
	printf "Tree from min AICc model (%s) saved to file: %s.AICc-%s.tre\n",$minimumAICcmodel, $infilename, $minimumAICcmodel;
	printf "Tree from min BIC model (%s) saved to file: %s.BIC-%s.tre\n",$minimumBICmodel, $infilename, $minimumBICmodel;
	print  "------------------------------------------------------------------\n";
	print  "\n";
	if ($debugLevel > 0) {
		close(DEBUG);
	}
}


#---------------- MrAICprint() ---------------------------#
# Print the results
#
sub MrAICprint {

    # Set output filename
    $filenameonly = basename($infilename);
    $MrAIC_output = File::Spec->join($output, "$filenameonly" . "$mraic_outfile_ending");
    
    #Get the date
    $datestring = prettydate();

    # Print to the outputfile
	open (OUT, '>', $MrAIC_output) or die "\a\nERROR!\n\nCould not create output file $!\n";
		printf OUT "\nOutput from %s version %s by %s\n",$program_name, $version, $author;
		print  OUT "-------------------------------------------------------------\n";
        printf OUT "$datestring\n";        
		printf OUT "\nInput data from file \"%s\" (%s taxa, %s characters)\n", $infilename, $ntax, $nchar;
		printf OUT "\nMinimum AIC  model: %s", $minimumAICmodel;
		printf OUT "\nMinimum AICc model: %s", $minimumAICcmodel;
		printf OUT "\nMinimum BIC  model: %s", $minimumBICmodel;
		print  OUT "\n";

		# See if the number of parameters are small compared to sample size
		if($changedNparams eq TRUE) {
			$estMinBASize = ($BAratio * $maxNumParametersBrLens);
			$estMinSize = ($maxNumParametersBrLens);
			printf  OUT "\n\nWARNING:\n";
			printf  OUT "Number of parameters is larger than the sample size!\n";
			printf  OUT "This causes problems in AICc calculations and you need to get\n";
			printf  OUT "more sequence data!\n";
			printf  OUT "\t(For %g terminal branches you probably need more\n", $ntax;
			printf  OUT "\tthan %g characters, and an approx. data size for\n", $estMinSize;
			printf  OUT "\thaving a nchar/nparam ratio over %g is about %g bp.)\n", $BAratio, $estMinBASize;
			printf  OUT "Calculations were done anyway while not considering branch\n";
			printf  OUT "lengths in the total number of parameters. Number of\n";
			printf  OUT "parameters was set to $maxNumParameters (number of free parameters\n";
			printf  OUT "in the GTR+I+G model).\n\n";
		}

		# Output sorted by AICc
		elsif (($nchar/$maxNumParameters) < $BAratio) { # Esentially testing if nchar < 400
			printf  OUT "\n\nNote: number of parameters compared to sample size\n";
			printf  OUT "is low (Nchar/Nparams < %g). Consider the use of AICc.\n", $BAratio;

			# Output sorted by AICc
			$cumAICcWeight = 0.0;
			printf  OUT "\n\n(Output sorted by AICc. w: Akaike weight, cw: cumulative w.)";
			printf  OUT "\n\nModel\tdf\tlnL\t\tAICc\t\twAICc\tcwAICc";
			foreach $i (sort {$MODELAICc{$a} <=> $MODELAICc{$b}} keys %MODELAICc) {
				$cumAICcWeight += $AkaikeAICcWeight{$i};
				printf OUT "\n%s\t%g\t%.4f\t%.4f\t%.4f\t%.4f", $MODELnr{$i}, $MODELdf{$i}, $MODELlnL{$i}, $MODELAICc{$i}, $AkaikeAICcWeight{$i}, $cumAICcWeight;
			}
			print  OUT "\n\n";
		}

		# Output sorted by AIC
		$cumAICWeight = 0.0;
		printf  OUT "\n\n(Output sorted by AIC. w: Akaike weight, cw: cumulative w.)";
		printf  OUT "\n\nModel\tdf\tlnL\t\tAIC\t\twAIC\tcwAIC";
		foreach $i (sort {$MODELAIC{$a} <=> $MODELAIC{$b}} keys %MODELAIC) {
			$cumAICWeight += $AkaikeAICWeight{$i};
			printf OUT "\n%s\t%g\t%.4f\t%.4f\t%.4f\t%.4f", $MODELnr{$i}, $MODELdf{$i}, $MODELlnL{$i}, $MODELAIC{$i}, $AkaikeAICWeight{$i}, $cumAICWeight;
		}
		print  OUT "\n\n";

		# Output sorted by BIC
		$cumBICWeight = 0.0;
		printf  OUT "\n\n(Output sorted by BIC. w: Akaike weight, cw: cumulative w.)";
		printf  OUT "\n\nModel\tdf\tlnL\t\tBIC\t\twBIC\tcwBIC";
		foreach $i (sort {$MODELBIC{$a} <=> $MODELBIC{$b}} keys %MODELBIC) {
			$cumBICWeight += $AkaikeBICWeight{$i};
			printf OUT "\n%s\t%g\t%.4f\t%.4f\t%.4f\t%.4f", $MODELnr{$i}, $MODELdf{$i}, $MODELlnL{$i}, $MODELBIC{$i}, $AkaikeBICWeight{$i}, $cumBICWeight;
		}
		print  OUT "\n\n";

		# Print MrBayes block
		if ($useModeltest eq FALSE) {
			my $dataSetName = "$infilename";
			my $partitionName = "Dummy";
			print  OUT "\n\n";
			print  OUT "To specify default values for best AIC, AICc, or BIC models in MrBayes:\n";
			print  OUT "\n";
			printf OUT "[Mrbayes block for the best AIC model (%s)]\n", $minimumAICmodel;
			print  OUT "\nBEGIN MRBAYES;\n";
			printf OUT " Charset %s = 1 - %s;\n",$dataSetName, $nchar;
			printf OUT " Partition %s = 1:%s;\n",$partitionName, $dataSetName;
			printf OUT " Set partition = %s;\n",$partitionName;
			printf OUT "%s\n", $minimumAICmbcmd;
			print  OUT "END;\n";
			print  OUT "\n\n";
			printf OUT "[Mrbayes block for the best AICc model (%s)]\n", $minimumAICcmodel;
			print  OUT "\nBEGIN MRBAYES;\n";
			printf OUT " Charset %s = 1 - %s;\n",$dataSetName, $nchar;
			printf OUT " Partition %s = 1:%s;\n",$partitionName, $dataSetName;
			printf OUT " Set partition = %s;\n",$partitionName;
			printf OUT "%s\n", $minimumAICcmbcmd;
			print  OUT "END;\n";
			print  OUT "\n\n";
			printf OUT "[Mrbayes block for the best BIC model (%s)]\n", $minimumBICmodel;
			print  OUT "\nBEGIN MRBAYES;\n";
			printf OUT " Charset %s = 1 - %s;\n",$dataSetName, $nchar;
			printf OUT " Partition %s = 1:%s;\n",$partitionName, $dataSetName;
			printf OUT " Set partition = %s;\n",$partitionName;
			printf OUT "%s\n", $minimumBICmbcmd;
			print  OUT "END;\n";
		}
		print  OUT "\n\n";
		print  OUT "-------------------------------------------------------------\n";
		print  OUT "End of Output\n\n";
	close(OUT);
	
	$infilename = File::Spec->join($output, "$filenameonly");
	
	# Print trees to treefiles
	open (AICTREE, ">$infilename.AIC-$minimumAICmodel.tre") or die "\a\nERROR!\n\nCould not create tree file $!\n";
		print  AICTREE "$minimumAICtree\n";
	close(AICTREE);
	open (AICcTREE, ">$infilename.AICc-$minimumAICcmodel.tre") or die "\a\nERROR!\n\nCould not create tree file $!\n";
		print  AICcTREE "$minimumAICctree\n";
	close(AICcTREE);
	open (BICTREE, ">$infilename.BIC-$minimumBICmodel.tre") or die "\a\nERROR!\n\nCould not create tree file $!\n";
		print  BICTREE "$minimumBICtree\n";
	close(BICTREE);
}


#------------- prettydate() --------------------#
# Print time stamp
sub prettydate {
   @_ = localtime(shift || time);
   return(sprintf("%02d:%02d %02d/%02d/%04d", @_[2,1], $_[4]+1, $_[3], $_[5]+1900));
} 


#------------- AIC($i) -------------------------#
# Calculate AIC
# AIC = -2lnL + 2K
# L: max likelihood, K: number of free parameters
#
sub AIC {
	my ($model_i) = @_;
	$MODELAIC{$model_i} = (((-2.0)*($MODELlnL{$model_i})) + (2.0*($MODELdf{$model_i})));
	if ($debugLevel == 1) { # DEBUG
		print DEBUG "\nin AIC() AIC: $MODELAIC{$model_i}\n";
	}
}


#------------- AICc($i) -------------------------#
# Calculate AICc
# AICc = -2lnL + 2K + 2K(K+1)/n-K-1
# L: max likelihood, K: number of free parameters, n: sample size
#
sub AICc {
	my ($model_i) = @_;
	$MODELAICc{$model_i} = (((-2.0)*($MODELlnL{$model_i})) + (2.0*($MODELdf{$model_i})) + ((2.0*($MODELdf{$model_i})*(($MODELdf{$model_i})+1.0))/($nchar-($MODELdf{$model_i})-1.0)));
	if ($debugLevel == 1) { # DEBUG
		print DEBUG "\nin AICc() AICc: $MODELAICc{$model_i}\n";
	}
}


#------------- BIC($i) -------------------------#
# Calculate BIC
# BIC = -2lnL + Kln(n)
# L: max likelihood, K: number of free parameters, n: sample size
#
sub BIC {
	my ($model_i) = @_;
	$MODELBIC{$model_i} = (((-2.0)*($MODELlnL{$model_i})) + ((log($nchar)) * ($MODELdf{$model_i})));
	if ($debugLevel == 1) { # DEBUG
		print DEBUG "\nin BIC() BIC: $MODELBIC{$model_i}\n";
	}
}


#------------- getDimensions() -------------------------#
# Get the number of taxa and characters from the infile
#
sub getDimensions {
	open (INFILE,"< $infilename") or die "\a\nERROR!\n\nCannot open $a for reading: $!";
	chomp (my $firstline = <INFILE>);
	$_ = $firstline;
	close(INFILE);
	if (/\s*(\d+)\s+(\d+)/) {
		$ntax = $1;
		$nchar = $2;
		$nbranch = ((2*$ntax)-3);
		$maxNumParameters = $nbranch + 10;	# Number of branches plus parameters in GTR+I+G
	}
	else {
		die "\a\nERROR!\n\nCould not read number of taxa and/or characters\n\n"
	}
}


#------------- helpInformation() ----------------#
# Print help information
#
sub helpInformation {
	print  "\n\n";
	printf "This is %s version %s by %s\n", $program_name, $version, $author;
	printf "Last changes: %s\n\n", $lastChanges;
	print  "Usage:  Interactively or by passing arguments as below\n\n";
	printf "   %s infile\n\n", $program_name;
	printf "   %s -modeltest infile\n", $program_name;
	print  "   \n";
	print  "Output:  A text file \"<infile>.MrAIC.txt\" with AIC, AICc,\n";
	print  "   BIC, Akaike weights, and MrBayes blocks for specifying\n";
	print  "   the best models. Three files \"*.tre\" with trees (in\n";
	print  "   Newick format) obtained under the minimum AIC, AICc, and\n";
	print  "   BIC models, respectively.\n";
	print  "   \n";
	printf "Description:  %s is a Perl script for calculating AIC,\n", $program_name;
	print  "   AICc, BIC, and  Akaike weights for nucleotide substitution\n";
	print  "   models. Likelihood scores under different models are\n";
	print  "   estimated using PHYML (Guindon and Gasquel, 2003).\n";
	print  "   Input is DNA data in PHYML format (see README.html).\n";
	print  "   If the argument '-modeltest' is parsed, 56\n";
	print  "   models (the ones tested in Modeltest) are evaluated in PHYML.\n";
	print  "   Default is to test the 24 models that can be specified in\n";
	print  "   MrBayes v3. These are JC, F81, K2P (aka K80), HKY, SYM\n";
	print  "   and GTR, each combined with Propinv (I) and/or Gamma (G).\n";
	print  "   \n";
	print  "Requirements:  PHYML must be installed on your system (named\n";
	print  "   \"phyml\") and be in the PATH. Alternatively, user might edit\n";
	print  "    mraic.pl to specify the full path to the PHYML binary.\n";
	print  "   \n";
	print  "Notes: In this script, sample size (n) used in AICc and BIC is\n";
	print  "   assumed to be the number of characters in the data matrix.\n";
	print  "   This is probably NOT correct when it comes to phylogenetic\n";
	print  "   analyses (Nylander, 2004), but serve as an approximation to\n";
	print  "   the true n.\n";
	print  "   Branch lengths are included in the total number of parameters\n";
	print  "   for each model when calculating the AIC, AICc, and BIC.\n";
	print  "   When the number of characters is low and/or when the number\n";
	print  "   of parameters is low, or when the ratio nchar/nparams is close\n";
	print  "   to 1, this causes problems in the AICc calculations.\n";
	print  "   AICc is not designed to handle the latter case and if\n";
	print  "   encountered, users are strongly encouraged to get more\n";
	print  "   sequence data. MrAIC.pl will run these data anyway, but\n";
	print  "   ignoring the number of branches as parameters.\n";
	print  "   Furthermore, all calculations of Akaike weights, AIC, AICc,\n";
	print  "   and BIC are all dependent on PHYML's ability to find the\n";
	print  "   maximum of the likelihood under each model. More elaborate\n";
	print  "   searches might be necessary to get more correct assessment\n";
	print  "   of the ML for some data sets!\n";
	print  "\n\n";
}


#-------------- SetMbMODELmbcmd() -------------------------------------------------#
# Set the MrBayes command for model. Key=model number, value=command.
#
sub SetMbMODELmbcmd {
	my $equal = "rates=equal";
	my $propinv = "rates=propinv";
	my $gamma = "rates=gamma";
	my $invgamma = "rates=invgamma";
	my $Prset = "Prset applyto=(1)";
	my $Lset = "Lset applyto=(1)";
	my $nst1 = "nst=1";
	my $nst2 = "nst=2";
	my $nst6 = "nst=6";
	my $shapepr = "shapepr=Uniform(0.1,50.0)";
	my $pinvarpr = "pinvarpr=Uniform(0.0,1.0)";
	my $eqfreqpr = "statefreqpr=Fixed(Equal)";
	my $uneqfreqpr = "statefreqpr=Dirichlet(1.0,1.0,1.0,1.0)";
	my $revmatpr = "revmatpr=Dirichlet(1.0,1.0,1.0,1.0,1.0,1.0)";
	my $tratiopr = "tratiopr=Beta(1.0,1.0)";
	$MODELmbcmd{0}	=	" $Lset $nst1 $equal;\n $Prset $eqfreqpr;";	#JC69
	$MODELmbcmd{1}	=	" $Lset $nst1 $propinv;\n $Prset $eqfreqpr $pinvarpr";	#JC69I
	$MODELmbcmd{2}	=	" $Lset $nst1 $gamma;\n $Prset $eqfreqpr $shapepr;";	#JC69G
	$MODELmbcmd{3}	=	" $Lset $nst1 $invgamma;\n $Prset $eqfreqpr $shapepr $pinvarpr;";	#JC69IG
	$MODELmbcmd{4}	=	" $Lset $nst1 $equal;\n $Prset $uneqfreqpr;";	#F81
	$MODELmbcmd{5}	=	" $Lset $nst1 $propinv;\n $Prset $uneqfreqpr $pinvarpr;";	#F81I
	$MODELmbcmd{6}	=	" $Lset $nst1 $gamma;\n $Prset $uneqfreqpr $shapepr;";	#F81G
	$MODELmbcmd{7}	=	" $Lset $nst1 $invgamma;\n $Prset $uneqfreqpr $shapepr $pinvarpr;";	#F81IG
	$MODELmbcmd{8}	=	" $Lset $nst2 $equal;\n $Prset $eqfreqpr $tratiopr;";	#K2P
	$MODELmbcmd{9}	=	" $Lset $nst2 $propinv;\n $Prset $eqfreqpr $tratiopr $pinvarpr;";	#K2PI
	$MODELmbcmd{10}	=	" $Lset $nst2 $gamma;\n $Prset $eqfreqpr $tratiopr $shapepr;";	#K2PG
	$MODELmbcmd{11}	=	" $Lset $nst2 $invgamma;\n $Prset $eqfreqpr $tratiopr $shapepr $pinvarpr;";	#K2PIG
	$MODELmbcmd{12}	=	" $Lset $nst2 $equal;\n $Prset $uneqfreqpr $tratiopr;";	#HKY
	$MODELmbcmd{13}	=	" $Lset $nst2 $propinv;\n $Prset $uneqfreqpr $tratiopr $pinvarpr;";	#HKYI
	$MODELmbcmd{14}	=	" $Lset $nst2 $gamma;\n $Prset $uneqfreqpr $tratiopr $shapepr;";	#HKYG
	$MODELmbcmd{15}	=	" $Lset $nst2 $invgamma;\n $Prset $uneqfreqpr $tratiopr $shapepr $pinvarpr;";	#HKYIG
	$MODELmbcmd{16}	=	" $Lset $nst6 $equal;\n $Prset $revmatpr $eqfreqpr;";	#SYM
	$MODELmbcmd{17}	=	" $Lset $nst6 $propinv;\n $Prset $revmatpr $eqfreqpr $pinvarpr;";	#SYMI
	$MODELmbcmd{18}	=	" $Lset $nst6 $gamma;\n $Prset $revmatpr $eqfreqpr;";	#SYMG
	$MODELmbcmd{19}	=	" $Lset $nst6 $invgamma;\n $Prset $revmatpr $eqfreqpr $shapepr $pinvarpr;";	#SYMIG
	$MODELmbcmd{20}	=	" $Lset $nst6 $equal;\n $Prset $revmatpr $uneqfreqpr;";	#GTR
	$MODELmbcmd{21}	=	" $Lset $nst6 $propinv;\n $Prset $revmatpr $uneqfreqpr $pinvarpr;";	#GTRI
	$MODELmbcmd{22}	=	" $Lset $nst6 $gamma;\n $Prset $revmatpr $uneqfreqpr $shapepr;";	#GTRG
	$MODELmbcmd{23}	=	" $Lset $nst6 $invgamma;\n $Prset $revmatpr $uneqfreqpr $shapepr $pinvarpr;";	#GTRIG
}


#-------------- SetMbMODELnr() -------------------------------------------------#
# Set the model number. Key=model number, value=model name.
#
sub SetMbMODELnr {
	$MODELnr{0}	=	"JC69";
	$MODELnr{1}	=	"JC69I";
	$MODELnr{2}	=	"JC69G";
	$MODELnr{3}	=	"JC69IG";
	$MODELnr{4}	=	"F81";
	$MODELnr{5}	=	"F81I";
	$MODELnr{6}	=	"F81G";
	$MODELnr{7}	=	"F81IG";
	$MODELnr{8}	=	"K2P";	# aka K80
	$MODELnr{9}	=	"K2PI";
	$MODELnr{10}	=	"K2PG";
	$MODELnr{11}	=	"K2PIG";
	$MODELnr{12}	=	"HKY";
	$MODELnr{13}	=	"HKYI";
	$MODELnr{14}	=	"HKYG";
	$MODELnr{15}	=	"HKYIG";
	$MODELnr{16}	=	"SYM";
	$MODELnr{17}	=	"SYMI";
	$MODELnr{18}	=	"SYMG";
	$MODELnr{19}	=	"SYMIG";
	$MODELnr{20}	=	"GTR";
	$MODELnr{21}	=	"GTRI";
	$MODELnr{22}	=	"GTRG";
	$MODELnr{23}	=	"GTRIG";
}


#---------- SetMbMODELcmd() -------------------------------------------#
# Set model commands. Key=model name, value=command string to PHYML.
#
sub SetMbMODELcmd {
	#my $I="V\nY\n";
	#my $Eq="F\n0.25\n0.25\n0.25\n0.25\n";
    my $I = q{};
    if ($sequential) {
        $I = "I\n";
    }
    else {
        $I = "";
    }
	$MODELcmd{0}	=	"$I+\nM\nM\nM\nM\nM\nR\nY\n";	# JC69
	$MODELcmd{1}	=	"$I+\nM\nM\nM\nM\nM\nV\nY\nR\nY\n";	# JC69I
	$MODELcmd{2}	=	"$I+\nM\nM\nM\nM\nM\nY\n";	# JC69G
	$MODELcmd{3}	=	"$I+\nM\nM\nM\nM\nM\nV\nY\nY\n";	# JC69IG
	$MODELcmd{4}	=	"$I+\nM\nM\nM\nM\nM\nM\nM\nF\nR\nY\n";	# F81
	$MODELcmd{5}	=	"$I+\nM\nM\nM\nM\nM\nM\nM\nF\nV\nY\nR\nY\n";	# F81I
	$MODELcmd{6}	=	"$I+\nM\nM\nM\nM\nM\nM\nM\nF\nY\n";	# F81G
	$MODELcmd{7}	=	"$I+\nM\nM\nM\nM\nM\nM\nM\nF\nV\nY\nY\n";	# F81IG
	$MODELcmd{8}	=	"$I+\nM\nM\nM\nM\nM\nM\nT\nY\nR\nY\n";	# K2P
	$MODELcmd{9}	=	"$I+\nM\nM\nM\nM\nM\nM\nT\nY\nR\nV\nY\nY\n";	# K2PI
	$MODELcmd{10}	=	"$I+\nM\nM\nM\nM\nM\nM\nT\nY\nY\n";	# K2PG
	$MODELcmd{11}	=	"$I+\nM\nM\nM\nM\nM\nM\nT\nY\nV\nY\nY\n";	# K2PIG
	$MODELcmd{12}	=	"$I+\nF\nT\nY\nR\nY\n";	# HKY
	$MODELcmd{13}	=	"$I+\nF\nT\nY\nR\nV\nY\nY\n";	# HKYI
	$MODELcmd{14}	=	"$I+\nF\nT\nY\nY\n";	# HKYG
	$MODELcmd{15}	=	"$I+\nF\nT\nY\nV\nY\nY\n";	# HKYIG
	$MODELcmd{16}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nY\n";	# SYM
	$MODELcmd{17}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n";	# SYMI
	$MODELcmd{18}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nY\n";	# SYMG
	$MODELcmd{19}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n";	# SYMIG
	$MODELcmd{20}	=	"$I+\nM\nM\nM\nF\nR\nY\n";	# GTR
	$MODELcmd{21}	=	"$I+\nM\nM\nM\nF\nR\nV\nY\nY\n";	# GTRI
	$MODELcmd{22}	=	"$I+\nM\nM\nM\nF\nY\n";	# GTRG
	$MODELcmd{23}	=	"$I+\nM\nM\nM\nF\nV\nY\nY\n";	# GTRIG
}


#-------------- SetMbMODELdf() -------------------------------------------------#
# Set degrees of freedom for models. Key=model name, value=degrees of freedom
# Branch length parameters are included in number of parameters!
#
sub SetMbMODELdf {
	if ($useBrlens eq TRUE) {
		$MODELdf{0}	=	0 + $nbranch;	# JC69
		$MODELdf{1}	=	1 + $nbranch;	# JC69I
		$MODELdf{2}	=	1 + $nbranch;	# JC69G
		$MODELdf{3}	=	2 + $nbranch;	# JC69IG
		$MODELdf{4}	=	3 + $nbranch;	# F81
		$MODELdf{5}	=	4 + $nbranch;	# F81I
		$MODELdf{6}	=	4 + $nbranch;	# F81G
		$MODELdf{7}	=	5 + $nbranch;	# F81IG
		$MODELdf{8}	=	1 + $nbranch;	# K2P
		$MODELdf{9}	=	2 + $nbranch;	# K2PI
		$MODELdf{10}	=	2 + $nbranch;	# K2PG
		$MODELdf{11}	=	3 + $nbranch;	# K2PIG
		$MODELdf{12}	=	4 + $nbranch;	# HKY
		$MODELdf{13}	=	5 + $nbranch;	# HKYI
		$MODELdf{14}	=	5 + $nbranch;	# HKYG
		$MODELdf{15}	=	6 + $nbranch;	# HKYIG
		$MODELdf{16}	=	5 + $nbranch;	# SYM
		$MODELdf{17}	=	6 + $nbranch;	# SYMI
		$MODELdf{18}	=	6 + $nbranch;	# SYMG
		$MODELdf{19}	=	7 + $nbranch;	# SYMIG
		$MODELdf{20}	=	8 + $nbranch;	# GTR
		$MODELdf{21}	=	9 + $nbranch;	# GTRI
		$MODELdf{22}	=	9 + $nbranch;	# GTRG
		$MODELdf{23}	=	10 + $nbranch;	# GTRIG
	}
	else {
		$MODELdf{0}	=	0;	# JC69
		$MODELdf{1}	=	1;	# JC69I
		$MODELdf{2}	=	1;	# JC69G
		$MODELdf{3}	=	2;	# JC69IG
		$MODELdf{4}	=	3;	# F81
		$MODELdf{5}	=	4;	# F81I
		$MODELdf{6}	=	4;	# F81G
		$MODELdf{7}	=	5;	# F81IG
		$MODELdf{8}	=	1;	# K2P
		$MODELdf{9}	=	2;	# K2PI
		$MODELdf{10}	=	2;	# K2PG
		$MODELdf{11}	=	3;	# K2PIG
		$MODELdf{12}	=	4;	# HKY
		$MODELdf{13}	=	5;	# HKYI
		$MODELdf{14}	=	5;	# HKYG
		$MODELdf{15}	=	6;	# HKYIG
		$MODELdf{16}	=	5;	# SYM
		$MODELdf{17}	=	6;	# SYMI
		$MODELdf{18}	=	6;	# SYMG
		$MODELdf{19}	=	7;	# SYMIG
		$MODELdf{20}	=	8;	# GTR
		$MODELdf{21}	=	9;	# GTRI
		$MODELdf{22}	=	9;	# GTRG
		$MODELdf{23}	=	10;	# GTRIG
	}
}


#-------------- SetModeltestMODELnr() -------------------------------------------------#
# Set the model number. Key=model number, value=model name.
#
sub SetModeltestMODELnr {
	$MODELnr{0}	=	"JC69";
	$MODELnr{1}	=	"JC69I";
	$MODELnr{2}	=	"JC69G";
	$MODELnr{3}	=	"JC69IG";
	$MODELnr{4}	=	"F81";
	$MODELnr{5}	=	"F81I";
	$MODELnr{6}	=	"F81G";
	$MODELnr{7}	=	"F81IG";
	$MODELnr{8}	=	"K80";
	$MODELnr{9}	=	"K80I";
	$MODELnr{10}	=	"K80G";
	$MODELnr{11}	=	"K80IG";
	$MODELnr{12}	=	"HKY";
	$MODELnr{13}	=	"HKYI";
	$MODELnr{14}	=	"HKYG";
	$MODELnr{15}	=	"HKYIG";
	$MODELnr{16}	=	"TrNef";
	$MODELnr{17}	=	"TrNefI";
	$MODELnr{18}	=	"TrNefG";
	$MODELnr{19}	=	"TrNefIG";
	$MODELnr{20}	=	"TrN";	# aka TN93
	$MODELnr{21}	=	"TrNI";
	$MODELnr{22}	=	"TrNG";
	$MODELnr{23}	=	"TrNIG";
	$MODELnr{24}	=	"K3P";	# aka K81
	$MODELnr{25}	=	"K3PI";
	$MODELnr{26}	=	"K3PG";
	$MODELnr{27}	=	"K3PIG";
	$MODELnr{28}	=	"K3Puf";
	$MODELnr{29}	=	"K3PufI";
	$MODELnr{30}	=	"K3PufG";
	$MODELnr{31}	=	"K3PufIG";
	$MODELnr{32}	=	"TIMef";
	$MODELnr{33}	=	"TIMefI";
	$MODELnr{34}	=	"TIMefG";
	$MODELnr{35}	=	"TIMefIG";
	$MODELnr{36}	=	"TIM";
	$MODELnr{37}	=	"TIMI";
	$MODELnr{38}	=	"TIMG";
	$MODELnr{39}	=	"TIMIG";
	$MODELnr{40}	=	"TVMef";
	$MODELnr{41}	=	"TVMefI";
	$MODELnr{42}	=	"TVMefG";
	$MODELnr{43}	=	"TVMefIG";
	$MODELnr{44}	=	"TVM";
	$MODELnr{45}	=	"TVMI";
	$MODELnr{46}	=	"TVMG";
	$MODELnr{47}	=	"TVMIG";
	$MODELnr{48}	=	"SYM";
	$MODELnr{49}	=	"SYMI";
	$MODELnr{50}	=	"SYMG";
	$MODELnr{51}	=	"SYMIG";
	$MODELnr{52}	=	"GTR";
	$MODELnr{53}	=	"GTRI";
	$MODELnr{54}	=	"GTRG";
	$MODELnr{55}	=	"GTRIG";
}


#---------- SetModeltestMODELcmd() -------------------------------------------#
# Set model commands. Key=model name, value=command string to PHYML.
#
sub SetModeltestMODELcmd {
	#my $I="V\nY\n";
	#my $Eq="E\n0.25\n0.25\n0.25\n0.25\n";
    my $I = q{};
    if ($sequential) {
        $I = "I\n";
    }
    else {
        $I = "";
    }
	$MODELcmd{0}	=	"$I+\nM\nM\nM\nM\nM\nR\nY\n";	# JC69
	$MODELcmd{1}	=	"$I+\nM\nM\nM\nM\nM\nV\nY\nR\nY\n";	# JC69I
	$MODELcmd{2}	=	"$I+\nM\nM\nM\nM\nM\nY\n";	# JC69G
	$MODELcmd{3}	=	"$I+\nM\nM\nM\nM\nM\nV\nY\nY\n";	# JC69IG
	$MODELcmd{4}	=	"$I+\nM\nM\nM\nM\nM\nM\nM\nF\nR\nY\n";	# F81
	$MODELcmd{5}	=	"$I+\nM\nM\nM\nM\nM\nM\nM\nF\nV\nY\nR\nY\n";	# F81I
	$MODELcmd{6}	=	"$I+\nM\nM\nM\nM\nM\nM\nM\nF\nY\n";	# F81G
	$MODELcmd{7}	=	"$I+\nM\nM\nM\nM\nM\nM\nM\nF\nV\nY\nY\n";	# F81IG
	$MODELcmd{8}	=	"$I+\nM\nM\nM\nM\nM\nM\nT\nY\nR\nY\n";	# K2P
	$MODELcmd{9}	=	"$I+\nM\nM\nM\nM\nM\nM\nT\nY\nR\nV\nY\nY\n";	# K2PI
	$MODELcmd{10}	=	"$I+\nM\nM\nM\nM\nM\nM\nT\nY\nY\n";	# K2PG
	$MODELcmd{11}	=	"$I+\nM\nM\nM\nM\nM\nM\nT\nY\nV\nY\nY\n";	# K2PIG
	$MODELcmd{12}	=	"$I+\nF\nT\nY\nR\nY\n";	# HKY
	$MODELcmd{13}	=	"$I+\nF\nT\nY\nR\nV\nY\nY\n";	# HKYI
	$MODELcmd{14}	=	"$I+\nF\nT\nY\nY\n";	# HKYG
	$MODELcmd{15}	=	"$I+\nF\nT\nY\nV\nY\nY\n";	# HKYIG
	$MODELcmd{16}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n010020\n1.00\n1.00\n1.00\nR\nY\n";	# TrNef
	$MODELcmd{17}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n010020\n1.00\n1.00\n1.00\nR\nV\nY\nY\n";	# TrNefI
	$MODELcmd{18}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n010020\n1.00\n1.00\n1.00\nY\n";	# TrNefG
	$MODELcmd{19}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n010020\n1.00\n1.00\n1.00\nV\nY\nY\n";	# TrNefIG
	$MODELcmd{20}	=	"$I+\nM\nM\nF\nT\nY\nR\nY\n";	# TrN
	$MODELcmd{21}	=	"$I+\nM\nM\nF\nT\nY\nR\nV\nY\nY\n";	# TrNI
	$MODELcmd{22}	=	"$I+\nM\nM\nF\nT\nY\nY\n";	# TrNG
	$MODELcmd{23}	=	"$I+\nM\nM\nF\nT\nY\nV\nY\nY\n";	# TrNIG
	$MODELcmd{24}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012210\n1.00\n1.00\n1.00\nR\nY\n";	# K3P
	$MODELcmd{25}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012210\n1.00\n1.00\n1.00\nR\nV\nY\nY\n";	# K3PI
	$MODELcmd{26}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012210\n1.00\n1.00\n1.00\nY\n";	# K3PG
	$MODELcmd{27}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012210\n1.00\n1.00\n1.00\nV\nY\nY\n";	# K3PIG
	$MODELcmd{28}	=	"$I+\nM\nM\nM\nM\nF\nK\n012210\n1.00\n1.00\n1.00\nR\nY\n";	# K3Puf
	$MODELcmd{29}	=	"$I+\nM\nM\nM\nM\nF\nK\n012210\n1.00\n1.00\n1.00\nR\nV\nY\nY\n";	# K3PufI
	$MODELcmd{30}	=	"$I+\nM\nM\nM\nM\nF\nK\n012210\n1.00\n1.00\n1.00\nY\n";	# K3PufG
	$MODELcmd{31}	=	"$I+\nM\nM\nM\nM\nF\nK\n012210\n1.00\n1.00\n1.00\nV\nY\nY\n";	# K3PufIG
	$MODELcmd{32}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012230\n1.00\n1.00\n1.00\n1.00\nR\nY\n";	# TIMef
	$MODELcmd{33}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012230\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n";	# TIMefI
	$MODELcmd{34}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012230\n1.00\n1.00\n1.00\n1.00\nY\n";	# TIMefG
	$MODELcmd{35}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012230\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n";	# TIMefIG
	$MODELcmd{36}	=	"$I+\nM\nM\nM\nM\nF\nK\n012230\n1.00\n1.00\n1.00\n1.00\nR\nY\n";	# TIM
	$MODELcmd{37}	=	"$I+\nM\nM\nM\nM\nF\nK\n012230\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n";	# TIMI
	$MODELcmd{38}	=	"$I+\nM\nM\nM\nM\nF\nK\n012230\n1.00\n1.00\n1.00\n1.00\nY\n";	# TIMG
	$MODELcmd{39}	=	"$I+\nM\nM\nM\nM\nF\nK\n012230\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n";	# TIMIG
	$MODELcmd{40}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nY\n";	# TVMef
	$MODELcmd{41}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n";	# TVMefI
	$MODELcmd{42}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nY\n";	# TVMefG
	$MODELcmd{43}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n";	# TVMefIG
	$MODELcmd{44}	=	"$I+\nM\nM\nM\nM\nF\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nY\n";	# TVM
	$MODELcmd{45}	=	"$I+\nM\nM\nM\nM\nF\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n";	# TVMI
	$MODELcmd{46}	=	"$I+\nM\nM\nM\nM\nF\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nY\n";	# TVMG
	$MODELcmd{47}	=	"$I+\nM\nM\nM\nM\nF\nK\n012314\n1.00\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n";	# TVMIG
	$MODELcmd{48}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nY\n";	# SYM
	$MODELcmd{49}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n";	# SYMI
	$MODELcmd{50}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nY\n";	# SYMG
	$MODELcmd{51}	=	"$I+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n";	# SYMIG
	$MODELcmd{52}	=	"$I+\nM\nM\nM\nF\nR\nY\n";	# GTR
	$MODELcmd{53}	=	"$I+\nM\nM\nM\nF\nR\nV\nY\nY\n";	# GTRI
	$MODELcmd{54}	=	"$I+\nM\nM\nM\nF\nY\n";	# GTRG
	$MODELcmd{55}	=	"$I+\nM\nM\nM\nF\nV\nY\nY\n";	# GTRIG

}


#-------------- SetModeltestMODELdf() -------------------------------------------------#
# Set degrees of freedom for models. Key=model name, value=degrees of freedom
# Branch length parameters are included in number of parameters!
#
sub SetModeltestMODELdf {
	if ($useBrlens eq TRUE) {
		$MODELdf{0}	=	0 + $nbranch;	# JC69
		$MODELdf{1}	=	1 + $nbranch;	# JC69I
		$MODELdf{2}	=	1 + $nbranch;	# JC69G
		$MODELdf{3}	=	2 + $nbranch;	# JC69IG
		$MODELdf{4}	=	3 + $nbranch;	# F81
		$MODELdf{5}	=	4 + $nbranch;	# F81I
		$MODELdf{6}	=	4 + $nbranch;	# F81G
		$MODELdf{7}	=	5 + $nbranch;	# F81IG
		$MODELdf{8}	=	1 + $nbranch;	# K80
		$MODELdf{9}	=	2 + $nbranch;	# K80I
		$MODELdf{10}	=	2 + $nbranch;	# K80G
		$MODELdf{11}	=	3 + $nbranch;	# K80IG
		$MODELdf{12}	=	4 + $nbranch;	# HKY
		$MODELdf{13}	=	5 + $nbranch;	# HKYI
		$MODELdf{14}	=	5 + $nbranch;	# HKYG
		$MODELdf{15}	=	6 + $nbranch;	# HKYIG
		$MODELdf{16}	=	2 + $nbranch;	# TrNef
		$MODELdf{17}	=	3 + $nbranch;	# TrNefI
		$MODELdf{18}	=	3 + $nbranch;	# TrNefG
		$MODELdf{19}	=	4 + $nbranch;	# TrNefIG
		$MODELdf{20}	=	5 + $nbranch;	# TrN
		$MODELdf{21}	=	6 + $nbranch;	# TrNI
		$MODELdf{22}	=	6 + $nbranch;	# TrNG
		$MODELdf{23}	=	7 + $nbranch;	# TrNIG
		$MODELdf{24}	=	2 + $nbranch;	# K3P
		$MODELdf{25}	=	3 + $nbranch;	# K3PI
		$MODELdf{26}	=	3 + $nbranch;	# K3PG
		$MODELdf{27}	=	4 + $nbranch;	# K3PIG
		$MODELdf{28}	=	5 + $nbranch;	# K3Puf
		$MODELdf{29}	=	6 + $nbranch;	# K3PufI
		$MODELdf{30}	=	6 + $nbranch;	# K3PufG
		$MODELdf{31}	=	7 + $nbranch;	# K3PufIG
		$MODELdf{32}	=	3 + $nbranch;	# TIMef
		$MODELdf{33}	=	4 + $nbranch;	# TIMefI
		$MODELdf{34}	=	4 + $nbranch;	# TIMefG
		$MODELdf{35}	=	5 + $nbranch;	# TIMefIG
		$MODELdf{36}	=	6 + $nbranch;	# TIM
		$MODELdf{37}	=	7 + $nbranch;	# TIMI
		$MODELdf{38}	=	7 + $nbranch;	# TIMG
		$MODELdf{39}	=	8 + $nbranch;	# TIMIG
		$MODELdf{40}	=	4 + $nbranch;	# TVMef
		$MODELdf{41}	=	5 + $nbranch;	# TVMefI
		$MODELdf{42}	=	5 + $nbranch;	# TVMefG
		$MODELdf{43}	=	6 + $nbranch;	# TVMefIG
		$MODELdf{44}	=	7 + $nbranch;	# TVM
		$MODELdf{45}	=	8 + $nbranch;	# TVMI
		$MODELdf{46}	=	8 + $nbranch;	# TVMG
		$MODELdf{47}	=	9 + $nbranch;	# TVMIG
		$MODELdf{48}	=	5 + $nbranch;	# SYM
		$MODELdf{49}	=	6 + $nbranch;	# SYMI
		$MODELdf{50}	=	6 + $nbranch;	# SYMG
		$MODELdf{51}	=	7 + $nbranch;	# SYMIG
		$MODELdf{52}	=	8 + $nbranch;	# GTR
		$MODELdf{53}	=	9 + $nbranch;	# GTRI
		$MODELdf{54}	=	9 + $nbranch;	# GTRG
		$MODELdf{55}	=	10 + $nbranch;	# GTRIG
	}
	else {
		$MODELdf{0}	=	0;	# JC69
		$MODELdf{1}	=	1;	# JC69I
		$MODELdf{2}	=	1;	# JC69G
		$MODELdf{3}	=	2;	# JC69IG
		$MODELdf{4}	=	3;	# F81
		$MODELdf{5}	=	4;	# F81I
		$MODELdf{6}	=	4;	# F81G
		$MODELdf{7}	=	5;	# F81IG
		$MODELdf{8}	=	1;	# K80
		$MODELdf{9}	=	2;	# K80I
		$MODELdf{10}	=	2;	# K80G
		$MODELdf{11}	=	3;	# K80IG
		$MODELdf{12}	=	4;	# HKY
		$MODELdf{13}	=	5;	# HKYI
		$MODELdf{14}	=	5;	# HKYG
		$MODELdf{15}	=	6;	# HKYIG
		$MODELdf{16}	=	2;	# TrNef
		$MODELdf{17}	=	3;	# TrNefI
		$MODELdf{18}	=	3;	# TrNefG
		$MODELdf{19}	=	4;	# TrNefIG
		$MODELdf{20}	=	5;	# TrN
		$MODELdf{21}	=	6;	# TrNI
		$MODELdf{22}	=	6;	# TrNG
		$MODELdf{23}	=	7;	# TrNIG
		$MODELdf{24}	=	2;	# K3P
		$MODELdf{25}	=	3;	# K3PI
		$MODELdf{26}	=	3;	# K3PG
		$MODELdf{27}	=	4;	# K3PIG
		$MODELdf{28}	=	5;	# K3Puf
		$MODELdf{29}	=	6;	# K3PufI
		$MODELdf{30}	=	6;	# K3PufG
		$MODELdf{31}	=	7;	# K3PufIG
		$MODELdf{32}	=	3;	# TIMef
		$MODELdf{33}	=	4;	# TIMefI
		$MODELdf{34}	=	4;	# TIMefG
		$MODELdf{35}	=	5;	# TIMefIG
		$MODELdf{36}	=	6;	# TIM
		$MODELdf{37}	=	7;	# TIMI
		$MODELdf{38}	=	7;	# TIMG
		$MODELdf{39}	=	8;	# TIMIG
		$MODELdf{40}	=	4;	# TVMef
		$MODELdf{41}	=	5;	# TVMefI
		$MODELdf{42}	=	5;	# TVMefG
		$MODELdf{43}	=	6;	# TVMefIG
		$MODELdf{44}	=	7;	# TVM
		$MODELdf{45}	=	8;	# TVMI
		$MODELdf{46}	=	8;	# TVMG
		$MODELdf{47}	=	9;	# TVMIG
		$MODELdf{48}	=	5;	# SYM
		$MODELdf{49}	=	6;	# SYMI
		$MODELdf{50}	=	6;	# SYMG
		$MODELdf{51}	=	7;	# SYMIG
		$MODELdf{52}	=	8;	# GTR
		$MODELdf{53}	=	9;	# GTRI
		$MODELdf{54}	=	9;	# GTRG
		$MODELdf{55}	=	10;	# GTRIG
	}
}