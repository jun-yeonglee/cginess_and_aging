#=============================================================================#
#Author: Jun-Yeong Lee
#Date: 03012020
#Usage: perl cal_WCX_rmBoth0_nbynYO.pl [0 or 1] [table]
#Examp: perl cal_WCX_rmBoth0_nbynYO.pl 0 normcnt_deseq2_GSE121539_mouse_infoAdded.txt
#Description: Perform Wilcoxon rank sum test and print out p-value between each Y vs each O.
#             Each sample should have "y" or "o" as the first word. (e.g. y.6.M.GSM12345; y=group; 6=age; M=sex)
#             transpose.pl is needed to run this script. A input table will be transposed.
#             Calculate after removing if gene exp in both are 0.
################################################################################################################
############# USE "1" as ARGV[0] when gene name should be separated (e.g. NM_XXXXXX|0|SYMBOL). see below note to fit.
############# if separation is NOT needed, use "0" as ARGV[0];
################################################################################################################
#=============================================================================#

use strict;
use Statistics::Test::WilcoxonRankSum;
use List::Util qw/ min max /;
use Parallel::ForkManager;
my $thread = 20;		# you can change thread number.

my $sSet = 0;
my @aFile = glob("GSE*.read.merged.deseq2_norm_e*_infoAdded_sex*_cgi[pnb].txt");	# input files

my $ForkManager = Parallel::ForkManager->new($thread);
foreach my $sIn (@aFile) {
	$ForkManager->start and next;
	my $sName = substr $sIn, 0, -4;


	unless (-f $sName.".txt_transpose.txt") {
		system ("perl transpose.pl ".$sName.".txt;\n");
	}

	open h_outWCX, ">nbynCD_wcx.rs_rmBoth0_".$sName.".txt";

	print h_outWCX "sampleC\tsampleD\tWCX\n";



	# read cnt table
	open h_in, $sIn."_transpose.txt" or die;
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
			push @aHead, $aTT[3];			# change this
			print $aTT[3]."\n";
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

		################ WCX ################

			my @aMedAll;
			my @aTarAll;

			for (my $i = 0; $i <= $#aAll; $i++) {
		#		print $aAll[$i]."\t".$aTempAll[$i]."\n"; <STDIN>;
				if ($aAll[$i] != 0 and $aTempAll[$i] != 0) {
					push @aMedAll, $aAll[$i];
					push @aTarAll, $aTempAll[$i];
		#			print $aAll[$i]."\t".$aTempAll[$i]."\n"; <STDIN>;
				}
			}

			my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
			$wilcox_test->load_data(\@aMedAll, \@aTarAll);
			my $sProb = $wilcox_test->probability();
			print h_outWCX $aLine1[0]."\t".$aLine[0]."\t".$sProb."\n";

		}
		close h_inTempO;
	}
	close hinTempY;

	close h_outWCX;



	system ("rm ".$sName.".txt_transpose.txt_TempO;\n");
	system ("rm ".$sName.".txt_transpose.txt_TempY;\n");
	$ForkManager->finish;
}
$ForkManager-> wait_all_children;

#system ("rm ".$sIn."_transpose.txt;\n");


















