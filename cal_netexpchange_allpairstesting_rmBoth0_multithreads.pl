#=============================================================================#
#Author: Jun-Yeong Lee
#Date: 03162020
#Usage: perl cal_netexpchange_allpairstesting_rmBoth0_multithreads.pl [#threads] [table]
#Examp: perl cal_netexpchange_allpairstesting_rmBoth0_multithreads.pl 16 normcnt_deseq2
#Description: Calculate net expression change of each Y vs each O.
#             Each sample should have "y" or "o" as the first word. (e.g. y.6.M.GSM12345; y=group; 6=age; M=sex)
#             transpose.pl is needed to run this script. A input table will be transposed.
#             Calculate after removing if gene exp in both are 0.
#=============================================================================#

use strict;
use List::Util qw/ min max /;
use Parallel::ForkManager;
my $thread = $ARGV[0];		# you can change thread number.

#############################################
###### CHECK this part to fit your data.
my $sSet = 0;
my @aFile = glob("*".$ARGV[1]."*");	# input files
#############################################

my $ForkManager = Parallel::ForkManager->new($thread);
foreach my $sIn (@aFile) {
	$ForkManager->start and next;

	my $sName = substr $sIn, 0, -4;

	unless (-f $sName.".txt_transpose.txt") {
		system ("perl transpose.pl ".$sName.".txt;\n");
	}


	open h_out, ">netexpchange_allpairstesting_rmBoth0_".$sName.".txt";

	print h_out "sampleC\tsampleD\tC_AUC\tD_AUC\tdAUC\(D-C\)\n";




	# read cnt table
	open h_in, $sName.".txt_transpose.txt" or die;
	my $sHead = <h_in>;
	chomp ($sHead); $sHead =~ s/\n//g; $sHead =~ s/\r//g;

	#############################################
	###### CHECK this part to fit your data.
	my @aHead;
	if ($sSet == 0) {
		@aHead = split /\t/, $sHead;

	} elsif ($sSet == 1) {
		my @aHead1 = split /\t/, $sHead;
		foreach my $sCC (@aHead1) {
			my @aTT = split /\|/, $sCC;		# change this to find a gene symbol
			push @aHead, $aTT[2];			# change this
		}
	}
	#############################################


	open hTempY, ">".$sName.".txt_transpose.txt_TempY";
	open hTempO, ">".$sName.".txt_transpose.txt_TempO";

	print hTempY $aHead[0];
	print hTempO $aHead[0];

	for (my $TT = 1; $TT < scalar @aHead; $TT++) {
		print hTempY "\t".$aHead[$TT];
		print hTempO "\t".$aHead[$TT];
	}
	print hTempY "\n";
	print hTempO "\n";

	while (my $sLine = <h_in>) {
		chomp ($sLine); $sLine =~ s/\n//g; $sLine =~ s/\r//g;
		my @aLine = split /\t/, $sLine;
		my @aGroup = split /\./, $aLine[0];
		if ($aGroup[0] eq "c" || $aGroup[0] eq "y" ) {
			print hTempY $sLine."\n";
		} elsif ($aGroup[0] eq "d" || $aGroup[0] eq "o" ) {
			print hTempO $sLine."\n";
		}
	}

	close hTempY;
	close hTempO;

	close h_in;



	open hinTempY, $sName.".txt_transpose.txt_TempY";
	my $sHeadTY = <hinTempY>;


	while (my $sLine1 = <hinTempY>) {
		chomp ($sLine1); $sLine1 =~ s/\n//g; $sLine1 =~ s/\r//g;
		my @aLine1 = split /\t/, $sLine1;
		my @aAll;
		for (my $i = 1; $i <= $#aHead; $i++) {
			push @aAll, $aLine1[$i];
		}


		open hinTempO, $sName.".txt_transpose.txt_TempO";
		my $sHeadTO = <hinTempO>;

		while (my $sLine = <hinTempO>) {
			chomp ($sLine); $sLine =~ s/\n//g; $sLine =~ s/\r//g;
			my @aLine = split /\t/, $sLine;
			my @aTempAll;
			for (my $i = 1; $i <= $#aHead; $i++) {
					push @aTempAll, $aLine[$i];
			}

		################ dAUC ################

			my @aConAll;
			my @aTarAll;

			for (my $i = 0; $i <= $#aAll; $i++) {
		#		print $aAll[$i]."\t".$aTempAll[$i]."\n"; <STDIN>;
				unless ($aAll[$i] == 0 && $aTempAll[$i] == 0) {
					push @aConAll, $aAll[$i];
					push @aTarAll, $aTempAll[$i];
		#			print $aAll[$i]."\t".$aTempAll[$i]."\n"; <STDIN>;
				}
			}

			my @aMinAll;
			my @aTAll = (@aConAll, @aTarAll);

			foreach my $sA (@aTAll) {
				if ($sA > 0) { push @aMinAll, $sA; }
			}

			my $sMinAll = min (@aMinAll);
			my $sMaxAll = max (@aMinAll);
			my $sAucAll;
			my $sConAucAll;

	#		print $aLine1[0]."\t".$aLine[0]."\t".$sMaxAll."\t".$sMinAll."\t".$s1."\t".$s2."\t".$s3."\n";

			for (my $i = 0; $i < scalar @aConAll; $i++) {
				if ($aConAll[$i] > 0) {
					$sConAucAll = $sConAucAll + log_(10,$aConAll[$i]) - log_(10,$sMinAll);
				}
				if ($aTarAll[$i] > 0) {
					$sAucAll = $sAucAll + log_(10,$aTarAll[$i]) - log_(10,$sMinAll);
				}
			}

	#		print $sMaxAll."\t".@aConAll."\n";

			my $sLogMaxMinAll = log_(10,$sMaxAll) - log_(10,$sMinAll);

			$sAucAll = ($sAucAll / $sLogMaxMinAll) / scalar @aConAll;
			$sConAucAll = ($sConAucAll / $sLogMaxMinAll) / scalar @aConAll;

			my $sdAucAll = $sAucAll - $sConAucAll;


			print h_out $aLine1[0]."\t".$aLine[0]."\t";
			print h_out $sConAucAll."\t".$sAucAll."\t".$sdAucAll."\n";



		}
		close h_inTempO;
	}
	close hinTempY;

	close h_out;


	#system ("rm ".$sIn."_transpose.txt;\n");

	system ("rm ".$sName.".txt_transpose.txt_TempO;\n");
	system ("rm ".$sName.".txt_transpose.txt_TempY;\n");
	$ForkManager->finish;
}
$ForkManager-> wait_all_children;



##########################
sub log_ {
	return log($_[1]) / log($_[0]);
}

















