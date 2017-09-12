use strict;
use warnings;
use POSIX qw(tmpnam);
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use List::MoreUtils qw(uniq);
use Time::localtime;
use List::Util qw( min max );

## ======================================
## Usage: see -h
## ======================================

sub usage{
  warn <<END;
  Usage:
  Run by typing: perl Primer_Design_For_Certain_Seq.pl -cfg [Parameter config file] -input [Seq fasta file] -out [Primer output] 
    Required params:
	-c|cfg							[s]	Parameter config file
	-i|input						[s]	Flanking Seq fasta file
	-o|out							[s]	Primer output
    Example: perl Primer_Design_For_Certain_Seq.pl -cfg [Parameter config file] -input [Seq fasta file] -out [Primer output] 
END
  exit;
}
## ======================================
## Get options
## ======================================

my %opt;
%opt = (
	'help'				=> undef,
	'debug'				=> undef,
	'cfg'				=> undef,
	'input'				=> undef,
	'out'				=> undef
);

die usage() if @ARGV == 0;
GetOptions (
  'h|help'				=> \$opt{help},
  'debug'				=> \$opt{debug},
  'c|cfg=s'				=> \$opt{cfg},
  'i|input=s'			=> \$opt{input},
  'o|out=s'				=> \$opt{out}
) or die usage();

#check input paramaters
die usage() if $opt{help};
die usage() unless ( $opt{cfg} );
die usage() unless ( $opt{input} );
die usage() unless ( $opt{out} );

########
#Main Function
########
my %complements;
$complements{"A"} = "T";
$complements{"C"} = "G";
$complements{"G"} = "C";
$complements{"T"} = "A";
$complements{"N"} = "N";
$complements{"a"} = "t";
$complements{"c"} = "g";
$complements{"g"} = "c";
$complements{"t"} = "a";
$complements{"n"} = "n"; 

if (-e "Primer3.Log.txt"){
		system "rm Primer3.Log.txt";
}
if (-e "Blast.log.txt"){
		system "rm Blast.log.txt";
}
if (-e $opt{out}){
		system "rm $opt{out}";
}
if (-e "tmp.dir"){
	system "rm -r tmp.dir";
}
system "mkdir tmp.dir";
my %params = &CFG_Parser($opt{cfg});
print "$opt{cfg} has been read! \n";
#print $params{'max_dist_for_sequencing'}."\n";
my %seq = &Seq_Parser($opt{input});
print "$opt{input} has been parsed! \n";


foreach my $snpid (sort keys %seq){
	my @h = split /\_/,$snpid;
	my $chr = $h[1];
	my @locations = split /-/,$h[2];
	my $start = $locations[0];
	print "                      \n";
	print "                      \n";
	print "########################################################\n";
	print "           Designing Primer for $snpid                  \n";
	print "########################################################\n";
	my $Seq = $seq{$snpid};
	my $input1 = $Seq; 
	(my $hout1_forward,my $hout1_reverse) = &Primer3_Run($snpid,$input1);
	my @forwards = uniq(@{$hout1_forward});
	my @reverses = uniq(@{$hout1_reverse});
	#print "Step1: \n$forwards[0]\n$reverses[0]\n";
	die "No Forward Primers can be designed by Primer3! \n" if (!@forwards);
	die "No Reverse Primers can be designed by Primer3! \n" if (!@reverses);
	print "Forward and Reverse Primers Candidates has been designed by Primer3! \n";
	my $hblastouts_forward; my $hblastouts_reverse;
	$hblastouts_forward = &qualified_blast_run(\@forwards);
	$hblastouts_reverse = &qualified_blast_run(\@reverses);
	#my $hcompetitive_reverses_blast_outs;my $hcommon_forwards_blast_outs;my $hcommon_reverses_blast_outs;
	print "No Forward Primers Blast outs are qualified! \n" if (!defined($hblastouts_forward));
	print "No Reverse Primers Blast outs are qualified! \n" if (!defined($hblastouts_reverse));
	if (!defined($hblastouts_forward ) and !defined($hblastouts_reverse) ){
		print "No Primers Blast outs are qualified! \n";
		next;
	} else {
		print "Primers Blast outs are qualified! \n";
	}
	###Determine Unique Primer Pairs###
	my %Unique_Primer_Pairs;
	if ($hblastouts_forward and $hblastouts_reverse){
		if ( (%{$hblastouts_forward} and %{$hblastouts_reverse})){
			%Unique_Primer_Pairs = &Unique_Primers($hblastouts_forward, $hblastouts_reverse,$chr); ###Filter and output unique Primer Pair 
			print " Primer Pairs have been uniqued! \n";
		} 
	}else {
		print "No input for common_forwards_blast and competitive_reverses_blast! \n";
	}
	print "No Unique Primers Pairs are found! \n" if (!%Unique_Primer_Pairs);
	open(OUT,">>$opt{out}");
	print "Outputing .... >> $opt{out} \n";
	print OUT "##The PCR Primers of $snpid \n";
	print OUT "Forward_Primer\tStart\tEnd\tGC\tTm\tReverse_Primer\tStart\tEnd\tGC\tTm\tProduct_Size\n";
	if (%Unique_Primer_Pairs){
		my $i=0;
		foreach my $id1 (keys %Unique_Primer_Pairs){
			my @infos1 = &Extract_Print_Info($id1,$start);
			foreach my $id2 (@{$Unique_Primer_Pairs{$id1}}){
				my @infos2 = &Extract_Print_Info($id2,$start);
				foreach my $info1 (@infos1){
					foreach my $info2 (@infos2){
						my @g1 = split /\t/,$info1;
						my @g2 = split /\t/,$info2;
						my $product_size = $g2[2] - $g1[1] + 1;
						print OUT $info1."\t".$info2."\t".$product_size."\n";
						$i++;
					}
				}
			}
		}
		my $num1 = $i/2;
		print "$num1 Primer Pairs with Forward Competitive Primer has been designed! \n";
	}
	close OUT;
}









########
#Subfunction
########


sub CFG_Parser{
	#Parse the configure file and obtain parameters for Primer3 and Blastn.
	my $cfgfile = shift;
	open(IN,$cfgfile);
	my %params;
	while(<IN>){
		chomp;s/\r//;
		next unless /\=/;
		my @g = split /\s*\=\s*/,$_;
		$params{$g[0]} = $g[1];
	}
	return %params;
}

sub Seq_Parser{
	#Parse the input sequence fasta file and obtain the Left, Right flanking sequene and SNP alleles.
	my $fa = shift;
	my %seq;
	open(FA,$fa);
	my $id;
	while(<FA>){
		chomp;s/\r//;
		if (/^>/){
			$id = $_;
			$id =~ s/^>//;
		} else {
			$_ =~ s/\s//mg;
			$seq{$id} .= $_;
		}
	}
	close FA;
	return %seq;
}


sub Primer3_Run {
	#Output Primer3 input file and Run Primer3, and then read output as array.
	my $snpid = shift;
	my $seq = shift;
	die "The length of Right Flanking sequence is smaller than $params{'PRIMER_MIN_SIZE'}! " if (length($seq) < $params{"PRIMER_MIN_SIZE"});
	if (-e "Primer3.tmp.input"){
		system "rm Primer3.tmp.input";
	}
	my $tmp_input = "Primer3.tmp.input";
	my $len = length($seq);
	open (OUT,">$tmp_input");
	print OUT "SEQUENCE_ID"."=".$snpid."\n";
	print OUT "SEQUENCE_TEMPLATE"."=".$seq."\n";
	#print OUT "SEQUENCE_TARGET"."="."0\,0"."\n";
	print OUT "PRIMER_TASK"."="."generic"."\n";
	print OUT "PRIMER_PICK_LEFT_PRIMER"."="."1"."\n";
	print OUT "PRIMER_PICK_INTERNAL_OLIGO"."="."0"."\n";
	print OUT "PRIMER_PICK_RIGHT_PRIMER"."="."1"."\n";
	print OUT "PRIMER_OPT_SIZE"."=".$params{"PRIMER_OPT_SIZE"}."\n";
	print OUT "PRIMER_MIN_SIZE"."=".$params{"PRIMER_MIN_SIZE"}."\n";
	print OUT "PRIMER_MAX_SIZE"."=".$params{"PRIMER_MAX_SIZE"}."\n";
	print OUT "PRIMER_MAX_TM"."=".$params{"PRIMER_MAX_TM"}."\n";
	print OUT "PRIMER_MIN_TM"."=".$params{"PRIMER_MIN_TM"}."\n";
	print OUT "PRIMER_OPT_TM"."=".$params{"PRIMER_OPT_TM"}."\n";
	print OUT "PRIMER_MAX_GC"."=".$params{"PRIMER_MAX_GC"}."\n";
	print OUT "PRIMER_MIN_GC"."=".$params{"PRIMER_MIN_GC"}."\n";
	print OUT "PRIMER_MAX_NS_ACCEPTED"."="."1"."\n";
	print OUT "PRIMER_PRODUCT_SIZE_RANGE"."=".$params{"min_pcr_product_size"}."-".$params{"max_pcr_product_size"}."\n";
	print OUT "P3_FILE_FLAG"."="."1"."\n";
	print OUT "PRIMER_EXPLAIN_FLAG"."="."1"."\n";
	print OUT "PRIMER_THERMODYNAMIC_PARAMETERS_PATH"."=".$params{"PRIMER_THERMODYNAMIC_PARAMETERS_PATH"}."\n";
	print OUT "PRIMER_MAX_SELF_ANY"."=".$params{"PRIMER_MAX_SELF_ANY"}."\n";
	print OUT "PRIMER_PAIR_MAX_COMPL_ANY"."=".$params{"PRIMER_PAIR_MAX_COMPL_ANY"}."\n";
	#print OUT "PRIMER_PAIR_MAX_DIFF_TM"."=".$params{"PRIMER_PAIR_MAX_DIFF_TM"}."\n";;
	print OUT "=\n";
	close OUT;
	print "The Primer3 input file 'Primer3.tmp.input' for Primer has been generated! \n";
	if (-e "$snpid.for"){
		system "rm $snpid.for";
	}
	if (-e "$snpid.rev"){
		system "rm $snpid.rev";
	}
	system "$params{'primer3'} Primer3.tmp.input > Primer3.Log.txt 2>&1 ";
	print "Primer3 is done! \n";
	my $forward = "$snpid.for";
	my @forward_out;
	my $i=0;
	open(IN,$forward);
	while(<IN>){
		chomp;s/\r//;
		$i++;
		if ($i==1){
			next;
		}
		next if /#/;
		$_ =~ s/\s+/\t/g;
		$_ =~ s/^\t//g;
		my @g = split /\t/,$_;
		shift @g;
		my $seq = $g[0];
		my $line =join "\t", @g;
		push @forward_out, $line;
	}
	close IN;
	
	my $reverse = "$snpid.rev";
	my @reverse_out;
	my $ii=0;
	open(IN,$reverse);
	while(<IN>){
		chomp;s/\r//;
		$ii++;
		if ($ii==1){
			next;
		}
		next if /#/;
		$_ =~ s/\s+/\t/g;
		$_ =~ s/^\t//g;
		my @g = split /\t/,$_;
		shift @g;
		my $start = $g[1]-$g[2]+1;###This is important! For reverse primer the start point is the 5'end in the reverse direction and 3'end in the positive direction. 
		$g[1] = $start;
		my $line =join "\t", @g;
		push @reverse_out, $line;
	}
	close IN;
	if (-e "$snpid.for"){
		system "mv $snpid.for tmp.dir";
	}
	if (-e "$snpid.rev"){
		system "mv $snpid.rev tmp.dir";
	}
	return (\@forward_out,\@reverse_out);
}



sub qualified_blast_run {
	#Make Blast query fasta file and Filter Blast Otput according to mismatch number and SNP overlap.
	my $query = shift;
	my $queryfile = "Blast.input.fa";
	open(IN,">".$queryfile);
	my $out = "tmp.blast.out";
	if (-e "tmp.blast.out"){
		system "rm tmp.blast.out";
	}
	my @primers = @{$query};
	foreach my $primer (@primers){
		my @g = split /\t/,$primer;
		my $seq = $g[0];
		my $id = join "_",@g;
		print IN ">$id\n$seq\n";
	}
	close IN;
	system "$params{'blastn'} -db $params{'reference_genome'} -query $queryfile -out $out -evalue $params{'evalue'} -num_threads $params{'num_cpus'} -outfmt \"6 std gaps nident\" -dust no -gapopen 4 -gapextend 2 -penalty -2 -reward 2 -word_size $params{'word_size'} -max_target_seqs 50 > Blast.log.txt 2>&1";
	#print "$params{'blastn'} -db $params{'reference_genome'} -query $queryfile -out $out -evalue $params{'evalue'} -num_threads $params{'num_cpus'} -outfmt \"6 std gaps nident\" -dust no -gapopen 4 -gapextend 2 -penalty -2 -reward 2 -word_size $params{'word_size'} -max_target_seqs 50 > Blast.log.txt 2>&1\n";
	my %outs;
	open (OUT,"$out");
	while(<OUT>){
		chomp;s/\r//;
		my @g = split /\t/,$_;
		my $start = $g[8];
		my $end = $g[9];
		if ($end < $start) {
			$g[8] = $end;
			$g[9] = $start;
		}
		my @ots = &Filter_Blast(@g);
		if (@ots){
			my $id = shift @g;
			my $line = join "\t",@g;
			push @{$outs{$id}}, $line;
		} else {
			next;
		}
	}
	close OUT;
	return (\%outs);
}


sub Filter_Blast {
	#Filter Blast Otput according to mismatch number and SNP overlap
	my ($id, $chr, $a, $b, $c, $d, $e, $f, $pos1, $pos2, $score1, $score2, $g, $h) = @_;
	my @raw = ($id, $chr, $a, $b, $c, $d, $e, $f, $pos1, $pos2, $score1, $score2, $g, $h);
	my @g = split /\_/,$id;
	my $seq = $g[0];
	my $len = length($seq);
	if (($h>$len-$params{"min_mismatches_close3p"}) or ($f<=$len-$params{"min_dist3p"} and $h>$len-$params{"min_mismatches"})){
		return (@raw);
	} else {
		return (());
	}
}



sub Unique_Primers {
	#Get Unique Primers and test the Tm difference and Cross-Dimer between Forward and Reverse Primers!.
	(my $hforward, my $hreverse, my $target_chr) = @_;
	my %outs;
	my %forwards = %{$hforward}; my %reverses = %{$hreverse};
	foreach my $id1 (keys %forwards){
		my $hforward_blast_outs = \@{$forwards{$id1}};
		foreach my $id2 (keys %reverses){
			my $hreverse_blast_outs = \@{$reverses{$id2}};
			if (&Unique_Blast($hforward_blast_outs,$hreverse_blast_outs,$target_chr)){
				#print "$id1\n$id2\n";
				my @TMs1 = &Extract_TM($id1);
				my @TMs2 = &Extract_TM($id2);
				my @Seqs1 = &Extract_Seq($id1);
				my @Seqs2 = &Extract_Seq($id2);
				my $max_TM_diff = 0;
				my $Dimer = 0;
				foreach my $TM1 (@TMs1){
					foreach my $TM2 (@TMs2){
						if ( abs($TM2-$TM1) > $max_TM_diff ){
							$max_TM_diff = abs($TM2-$TM1);
						}
					}
				}
				foreach my $seq1 (@Seqs1){
					foreach my $seq2 (@Seqs2){
						if ( &crossdimer($seq1,$seq2)==1 ){
							$Dimer = 1;
							#print "Cross-Dimer between Forward and Reverse Primers exists! $seq1 $seq2\n";
						}
					}
				}
				foreach my $seq1 (@Seqs1){
						if ( &dimer($seq1,$seq1)==1 ){
							$Dimer = 1;
							#print "Self-Dimer for Forward Primers exists! $seq1 \n";
						}
				}
					foreach my $seq2 (@Seqs2){
						if ( &dimer($seq2,$seq2)==1 ){
							$Dimer = 1;
							#print "Self-Dimer for Reverse Primers exists! $seq2 \n";
						}
					}
				if (($max_TM_diff < 5) and ($Dimer ==0)){
					push @{$outs{$id1}},$id2;
				} else {
					next;
					#print "TM difference between Left Primer and Right Primer is larger than 5! Or dimer Existed! \n";
					#print "TM difference between Left Primer and Right Primer is larger than 5!\n $id1\n$id2\n";
				}
			}
		}
	}
	return (%outs);
}

sub Unique_Blast {
	#Get the primers pairs with unique relative position in the whole genome.
	(my $hforward_blast_outs, my $hreverse_blast_outs,my $target_chr) = @_;
	my @forward_blast_outs = @{$hforward_blast_outs};
	my @reverse_blast_outs = @{$hreverse_blast_outs};
	my %test_unique;
	my %best_unique;
	foreach my $out1 (@forward_blast_outs){
		my ($chr1, $a1, $b1, $c1, $d1, $e1, $f1, $pos11, $pos12, $score11, $score12, $g1, $h1) = split /\t/,$out1;
		foreach my $out2 (@reverse_blast_outs ){
			my ($chr2, $a2, $b2, $c2, $d2, $e2, $f2, $pos21, $pos22, $score21, $score22, $g2, $h2) = split /\t/,$out2;
			my $max_dis = ($pos22-$pos11+1);
			if (($chr1 eq $chr2) and ($max_dis <= $params{"max_dist_for_sequencing"}) and ($max_dis > 0)){
				$test_unique{$chr2}++;
				if ( ( $max_dis <= $params{'max_pcr_product_size'}) and ($max_dis >= $params{'min_pcr_product_size'}) ){
					$best_unique{$chr2}++;
				}
			}
		}
	}
	if (!%test_unique || !%best_unique){
		return 0;
	}elsif ((scalar(keys %test_unique) == 1) and (scalar(keys %best_unique) == 1) and exists($test_unique{$target_chr}) and exists($best_unique{$target_chr}) and ($test_unique{$target_chr} == 1) and ($best_unique{$target_chr} == 1)){
		return 1;
	} else {
		return 0;
	}
	
}

sub Extract_TM {
	#Extract Tm values from the Primer3 results.
	my $id = shift;
	my @TMs;
	if ($id =~ /\-/){
		my @g1 = split /\-/,$id;
		foreach my $v (@g1){
			my @g2 = split /\_/,$v;
			push @TMs, $g2[5];
		}
	} else {
		my @g2 = split /\_/,$id;
		push @TMs, $g2[5];
	}
	return (@TMs);
}

sub Extract_Seq {
	#Extract Sequence from the Primer3 results.
	my $id = shift;
	my @seqs;
	if ($id =~ /\-/){
		my @g1 = split /\-/,$id;
		foreach my $v (@g1){
			my @g2 = split /\_/,$v;
			push @seqs, $g2[0];
		}
	} else {
		my @g2 = split /\_/,$id;
		push @seqs, $g2[0];
	}
	return (@seqs);
}


sub Extract_Print_Info {
	#Extract Sequence, position, GC content and Tm values from the Primer3 results. It can distinguish the forward and reverse primer.
	(my $id, my $s) = @_;
	my @outs;
	my @g2 = split /\_/,$id;
	my $seq = $g2[0];
	my $start = $s+$g2[1]+1;
	my $end = $start +$g2[2] -1+1;
	my $GC = $g2[4];
	my $TM = $g2[5];
	my $line = join "\t", ($seq,$start,$end,$GC,$TM);
	push @outs, $line;
	return (@outs);
}

sub Complement {
	#Return the complemental sequence. 
	my $seq = shift;
	my @g = split //,$seq;
	my @out;
	foreach my $s (@g){
		push @out, $complements{$s};
	}
	my $line = join "",@out;
	return $line;
}

sub crossdimer {
	#Dimers: A primer self-dimer is formed by intermolecular interactions between two of the same primers, 
	#i.e., forward vs. forward or reverse vs. reverse. Cross-dimers are formed by intermolecular interactions
	#between a forward and a reverse primer. Primer dimers with more than 4 consecutive base pairings are not passed,
	#and will not be considered as potential primer pairs.
	(my $seq1, my $seq2) = @_;
	my $total_len1 = length($seq1);
	my $total_len2 = length($seq2);
	my $index = 0;
	my $reverse_complemennt_3end_seq1 = &Reverse_Complement(substr($seq1,($total_len1-3),3));
	my $end_seq2 = substr($seq2,($total_len2-3),3);
	if ($reverse_complemennt_3end_seq1 eq $end_seq2){
		$index=1;
	}
	if ($index==1){
		return 1;
	} else {
		return 0;
	}
}

sub dimer {
	#Dimers: A primer self-dimer is formed by intermolecular interactions between two of the same primers, 
	#i.e., forward vs. forward or reverse vs. reverse. Cross-dimers are formed by intermolecular interactions
	#between a forward and a reverse primer. Primer dimers with more than 4 consecutive base pairings are not passed,
	#and will not be considered as potential primer pairs.
	(my $seq1, my $seq2) = @_;
	my $total_len = length($seq1);
	my $index = 0;
	for(my $i=0;$i<=$total_len-3;$i++){
		my $reverse_complemennt_seq = &Reverse_Complement(substr($seq1,$i,3));
		if ($seq2=~/$reverse_complemennt_seq/){
			$index = 1;
		}
	}
	if ($index==1){
		return 1;
	} else {
		return 0;
	}
}

sub Reverse_Complement {
	my $seq = shift;
	my @g = split //,$seq;
	my @out;
	foreach my $s (reverse @g){
		push @out, $complements{$s};
	}
	my $line = join "",@out;
	return $line;
}
