
#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
no warnings 'experimental::smartmatch';
use List::Util qw(shuffle);

my $vcf_imp_f = '';
my $vcf_seq_f = '';

my $var_hash = '';

my $extra = '';

GetOptions(
    "dosage-file=s"   => \$vcf_imp_f,
    "geno-file=s"     => \$vcf_seq_f,
    "site-list"       => \$var_hash,
    "extra"           => \$extra
);

sub getline {
    my $infile = shift;

    if ($infile and my $line = <$infile>) {
	chomp $line;
	return split(/\t/, $line);
    }
    else{
	return ();
    }
}

if( $vcf_imp_f eq '' || $vcf_seq_f eq '' ){
	die "specify both --dosage-file and --geno-file";
}

if( not -e $vcf_imp_f ){
	 die "$vcf_imp_f doesn't exist";
}
if( not -e $vcf_seq_f ){
	die "$vcf_seq_f doesn't exist";
}

open(my $vcf_seq, "-|", "zgrep -v '##' $vcf_seq_f | cut -f1-2,4-5,9- ");
open(my $vcf_imp, "-|", "zgrep -v '##' $vcf_imp_f | cut -f1-2,4-5,9- ");

my $nxcol = 5;

my @line_seq = getline($vcf_seq);
my @line_imp = getline($vcf_imp);

my @iid_seq = @line_seq[ $nxcol .. $#line_seq ];
my @iid_imp = @line_imp[ $nxcol .. $#line_imp ];

my %imp_map = map { $_ => 1 } @iid_imp;
my %seq_map = map { $_ => 1 } @iid_seq;

my @idx_imp;
my @idx_seq;

my %index_seq;
@index_seq{ @iid_seq } = (0..$#iid_seq);

for (my $i = 0; $i <= $#iid_imp; $i++ ){
    if( exists $index_seq{ $iid_imp[$i] } ){
        push @idx_imp, $i + $nxcol;
	push @idx_seq, $index_seq{ $iid_imp[$i] } + $nxcol;
    }
}

#Make hash for sequenced file
my %var_hash;

@line_seq = getline($vcf_seq);
while ($#line_seq>0){
   my $seqkey=join(":",@line_seq[1..3]);
   $var_hash{$seqkey}=[@line_seq];
   @line_seq=getline($vcf_seq);
}

#print "@{[%var_hash]}\n";


@line_imp = getline($vcf_imp);
print "CHR\tPOS\tREF\tALT\tN\tAC\tAF\tIF\tRsq\tRsq_MaCH"; #\tRsq_unph";

if($extra){
    print "\tED0\tVD0\tED1\tVD1";
}

print "\n";


while ( $#line_imp > 0) {
    my $impkey=join(":",@line_imp[1..3]);
    if(exists $var_hash{$impkey} ){
        my @line_seq=@{$var_hash{$impkey}};	
	my @format_seq = split(/:/, $line_seq[4]);
	my @format_imp = split(/:/, $line_imp[4]);
	print join("\t",@format_seq);
	print("\n");


	next if not "GT" ~~ @format_seq;
	next if not "DS" ~~ @format_imp;
	
	my ($GT_i) = grep { $format_seq[$_] eq "GT" } 0..$#format_seq;
	my ($DS_i) = grep { $format_imp[$_] eq "DS" } 0..$#format_imp;
	
	my ($SG,$SD,$SDD,$SGG,$SDG,$N) = (0.0)x6;
	for( my $i = 0; $i <= $#idx_imp; $i++){
	    my @geno_seq = split(/:/,$line_seq[ $idx_seq[$i] ]);
	    my @geno_imp = split(/:/,$line_imp[ $idx_imp[$i] ]);
	    my @GT = split(/[|\/]/, $geno_seq[$GT_i] );
	    next if '.' ~~ @GT;
	    my $G = $GT[0] + $GT[1];
	    my $D = $geno_imp[$DS_i];
	    	$SG += $G;
	    	$SD += $D;
	    	$SDD += $D*$D;
	    	$SGG += $G*$G;
	    	$SDG += $D*$G;

	    	$N += 1;
	}
        my ($EG, $VG, $ED, $VD, $CovDG) = (0.0)x5;
	if ($N > 1){
            $EG = $SG/$N;
	    $VG = ($SGG - $N*$EG*$EG)/($N - 1);
	    $ED = $SD/$N;
	    $VD = ($SDD - $N*$ED*$ED)/($N - 1);
	    $CovDG = ($SDG - $N*$ED*$EG)/($N - 1);
	}
	print join("\t", @line_seq[0..3]) . "\t$N\t$SG";
	printf("\t%.3f\t%.3f",0.5*$EG,0.5*$ED);
	
	my ($Rsq, $R2m, $Rsq_unph) = (0.0)x3;
	if( $SG > 0 and $VG*$VD > 0 and $ED > 0 and $EG > 0 ){
	    $Rsq = $CovDG*$CovDG / ( $VD*$VG );
	    $R2m = $VD / ( $ED * ( 1 - 0.5*$ED ) );
	    printf("\t%.3f\t%.3f",$Rsq,$R2m);
	}else{
	    print "\tNA\tNA";
	}
	print "\n";
	@line_imp = getline($vcf_imp);
    }
    else{
	@line_imp = getline($vcf_imp);    	
    }
}
close($vcf_seq);
close($vcf_imp);
