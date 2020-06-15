#!/usr/bin/perl -w
use strict;
my $main = "/storage/chen/home/jw29";
my $var_file = "HGMD/data/HGMD_2014_2016_USH2A_updated";
my $temp = "RP_gene";
my $freq_cut = 1;
my $output = $var_file."_maxAF$freq_cut"."_bin_pt_combi";
my $R = "/storage/chen/Software/R-3.5.0/bin/R";
my $STGD = "HGMD/data/dxdb_stage1_all_sort_freq_clean_NULL-2-15-2017_USH2A_add_screen_organized.gz";
my $total = 107 * 2;
open(OUTPUT,">$output") || die ("$output");
#1       94488103        .       T       C       3.0     0.00396825396825397     het:1|hom:1;ANN=C|intron_variant|MODIFIER|ABCA4|ABCA4|transcript|NM_000350.2|protein_coding|32/49|c.4668-596A>G||||||   FBP_144:1/1:5:0:5:.|RKK_374:0/1:4:1:5:.|        HGMD_None 

#print OUTPUT "chr\tpos\t.\tref\talt\tac\taf\tannotation\tpt_indo\tHGMD\tfisher_p_greater_max\tfisher_p_greater_all\tfisher_p_less_max\tfisher_p_less_all\tmax_freq\tmax_ac\tmax_an\tfreq\tac\tan\thom\tmax_hom\n";

print OUTPUT "Var\tHGMD_ref\tAC_STGD\tAF_STGD\tAC_Max\tAF_Max\tAN_Max\tNFE_AC\tNFE_AF\tNFE_AN\tMax_bino_p_greater\tNfe_bino_p_greater\tMax_bino_p_less\tNfe_bino_p_less\tHom\tMax_bino_af_p_greater\tNFE_bino_af_p_greater\n";

my $tabix  ="tabix";
my $gnomad = "/storage/chen/tmp/disease_gene/gnomad/gnomad_pop_freq_new_sort.gz";

#my $gnomad = "gnomad_pop_freq_new_sort.gz";
#my $gnomad = "$main/gnomad/data/gnomad_exomes_r2.0.1.sites_sort_5-26-2017.gz";
#NO      1       13445   C       A       1       0.000315258511979823    DDX11L1|DDX11L1 nonsynonymous_SNV       RP:solved_RP:SRF_980
#YES     1       1993722 C       T       2.0     0.00537634408602150538  PRKCZ   intron_variant  unsolved_LCA:1275:1/1:6:0:6:.


my $dis_freq = sqrt(1*(0.11)/4000);

open(INPUT,$var_file) || die ("$var_file");
my %hash;
while(my $line = <INPUT>){
if($line =~ /^#/){
next;
}
chomp $line;
my @info = split(/\s+/,$line);
my @info2 = split(/:/,$info[0]);
my $region = "$info2[0]:$info2[1]-$info2[1]";
my $ref = uc($info2[2]);
my $alt = uc($info2[3]);
my $temp_line = "$main/STGD/data/tabix_$temp";
#print "$region\n";
` $tabix $gnomad $region > $temp_line`;

#print "$gnomad\n";
open(INPUT1,$temp_line) || die ("$temp_line");
my $freq = 0;
my $max_freq =0;
my $max_ac =0;
my $max_an =123136*2;
my $ac=0;
my $an=123136*2;
my $hom=0;
my $max_hom=0;
my $flag_wes=0;
while(my $line1 = <INPUT1>){
my @info1 = split(/\s+/,$line1);
#1	879176	G	T	AC:4	AF:1.64101e-05	AN:243752	AC_AFR:0	AC_AMR:0	AC_ASJ:0	AC_EAS:0	AC_FIN:0	AC_NFE:4	AC_OTH:0	AC_SAS:0	AC_Male:2	AC_Female:2	AN_AFR:14948	AN_AMR:33546	AN_ASJ:9790	^CAN_EAS:17204	AN_FIN:22154	AN_NFE:109878	AN_OTH:5454	AN_SAS:30778	AN_Male:134250	AN_Female:109502	AF_AFR:0.00000e+00	AF_AMR:0.00000e+00	AF_ASJ:0.00000e+00	AF_EAS:0.00000e+00	AF_FIN:0.00000e+00	AF_NFE:3.64040e-05	AF_OTH:0.00000e+00	AF_SAS:0.00000e+00	AF_Male:1.48976e-05	AF_Female:1.82645e-05	Hom_AFR:0	Hom_AMR:0	Hom_ASJ:0	Hom_EAS:0	Hom_FIN:0	Hom_NFE:0	Hom_OTH:0	Hom_SAS:0	Hom_Male:0	Hom_Female:0	Hom:0	POPMAX:NFE	AC_POPMAX:4	AN_POPMAX:109878	AF_POPMAX:3.64040e-05
#1       94487251        G       T       AC:7    AF:2.85679e-05  AN:245030       AC_AFR:0        AC_AMR:3        AC_ASJ:0        AC_EAS:0        AC_FIN:0   AC_NFE:3 AC_OTH:1        AC_SAS:0        AC_Male:4       AC_Female:3     AN_AFR:15022    AN_AMR:33512    AN_ASJ:9756     AN_EAS:17202    AN_FIN:22260    AN_NFE:111148       AN_OTH:5464     AN_SAS:30666    AN_Male:134252  AN_Female:110778        AF_AFR:0.00000e+00      AF_AMR:8.95202e-05      AF_ASJ:0.00000e+00 AF_EAS:0.00000e+00       AF_FIN:0.00000e+00      AF_NFE:2.69910e-05      AF_OTH:1.83016e-04      AF_SAS:0.00000e+00      AF_Male:2.97947e-05     AF_Female:2.70812e-05       Hom_AFR:0       Hom_AMR:0       Hom_ASJ:0       Hom_EAS:0       Hom_FIN:0       Hom_NFE:0       Hom_OTH:0       Hom_SAS:0       Hom_Male:0 Hom_Female:0     Hom:0   POPMAX:AMR      AC_POPMAX:3     AN_POPMAX:33512 AF_POPMAX:8.95202e-05
print "$info1[0]\t$info1[1]\t$info1[2]\t$info1[3]\n";

if(($info2[0] eq "$info1[0]") && ($info2[1] == $info1[1])&&(uc($info1[2]) eq $ref)&&(uc($info1[3])eq $alt)){
my @ac = split(/:/,$info1[12]);
my @af = split(/:/,$info1[32]);
my @an = split(/:/,$info1[22]);
my @ac_max = split(/:/,$info1[-3]);
my @an_max = split(/:/,$info1[-2]);
my @af_max = split(/:/,$info1[-1]);
my @hom = split(/:/,$info1[-5]);
#my @maxhom = split(/:/,$info1[-1]);
$ac=$ac[1];
$freq=$af[1];
$an=$an[1];
$max_ac = $ac_max[1];
$max_an = $an_max[1];
$max_freq = $af_max[1];
$hom = $hom[1];
#$max_hom = $maxhom[1];
$flag_wes=1;
last;
}
}
#my $gnomad = "$main/gnomad/data/gnomad_exomes_r2.0.1.sites_sort_5-26-2017.gz";
#NO      1       13445   C       A       1       0.000315258511979823    DDX11L1|DDX11L1 nonsynonymous_SNV       RP:solved_RP:SRF_980
#YES     1       1993722 C       T       2.0     0.00537634408602150538  PRKCZ   intron_variant  unsolved_LCA:1275:1/1:6:0:6:.
`rm $temp_line`;
my $new_region = "chr$region";
` $tabix $STGD $new_region > $temp_line`;
open(INPUT1,$temp_line) || die ("$temp_line") ;
my $ac_STGD=0;
my $af_STGD=0;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\t/,$line1);
$info1[0] =~ s/chr//g;
#print "$info1[0]\t$info1[1]\t$info1[3]\t$info1[4]\n";
if(($info2[0] eq "$info1[0]") && ($info2[1] == $info1[1])&&(uc($info1[3]) eq $ref)&&(uc($info1[4])eq $alt)){
#print "enter_RP\t$line1\n";
$ac_STGD=$info1[5];
$af_STGD=$info1[6];
#print "$ac_STGD\t$af_STGD\n";

last;
}
}
print "$ac_STGD\t$af_STGD\n";
my $pt_freq;
my $bin_p = "NA";
my $direction = "greater";
#print "$line\t$ac_STGD\t$af_STGD\t$max_ac\t$max_freq\t$max_an\t$ac\t$freq\t$an\t$hom\n";
#exit;

if($max_freq eq "\."){
$max_freq=0;
}
#$max_an=123136*2;
#$pt_freq = $max_freq/0.01;

$pt_freq = $max_freq/$dis_freq;
if($max_freq <= $dis_freq){
$bin_p = bin_test($pt_freq,$ac_STGD,$total,$direction);
}
my $bin_p1 = "NA";

if($freq eq "\."){
$freq=0;
}
my $pt_freq1 = $freq/$dis_freq;
if($freq <=$dis_freq){
$bin_p1 = bin_test($pt_freq1,$ac_STGD,$total,$direction);
}

my $STGD_rest = $total - $ac_STGD;
#my $max_rest = $max_an - $max_ac;
#my $fisher_p = fisher_test($ac_STGD,$STGD_rest,$max_ac, $max_rest,"less"); 

#my $rest = $an - $ac;
#my $fisher_p1 = fisher_test($ac_STGD,$STGD_rest,$ac, $rest,"less"); 


my $bin_p_less = "NA";
$direction = "less";

if($max_freq <=$dis_freq){
$bin_p_less = bin_test($pt_freq,$ac_STGD,$total,$direction);
}


my $bin_p_less1 = "NA";
if($freq <= $dis_freq){
$bin_p_less1 = bin_test($pt_freq1,$ac_STGD,$total,$direction);
}

my $bin_p_greater = "NA";
$direction = "greater";

#if($max_freq <= $dis_freq){
$bin_p_greater = bin_test($max_freq,$ac_STGD,$total,$direction);
#}

my $bin_p_greater1 = "NA";

#if($freq <=$dis_freq){
$bin_p_greater1 = bin_test($freq,$ac_STGD,$total,$direction);



print OUTPUT "$line\t$ac_STGD\t$af_STGD\t$max_ac\t$max_freq\t$max_an\t$ac\t$freq\t$an\t$bin_p\t$bin_p1\t$bin_p_less\t$bin_p_less1\t$hom\t$bin_p_greater\t$bin_p_greater1\n";

}

sub bin_test {
my $R_script = "STGD/scripts/R_$ARGV[0].R";
my $R_out = "STGD/scripts/R_$ARGV[0].out";

open (OUTFILE, ">$R_script");
print OUTFILE "pvalues<-binom.test(", "$_[1]",",$_[2]",", p= $_[0]",", alternative=\"$_[3]\" ,conf.level= 0.95)\$p\.value\n";
print OUTFILE "write(pvalues\,\"$R_out\")\n";
close OUTFILE;
`$R < $R_script --no-save --slave`;
open INFILE, "<$R_out";
my $file = <INFILE>;
my @results = split(/\s+/,$file);
`rm $R_script $R_out`;
close INFILE;
return($results[0]);
}

sub chi2_p {
my $R_script = "STGD/scripts/R_$ARGV[0].R";
my $R_out = "STGD/scripts/R_$ARGV[0].out";

open (OUTFILE, ">$R_script");
print OUTFILE "pvalues<-pchisq(", "$_[0]",", df=", "$_[1]",", lower.tail = FALSE)\n";
print OUTFILE "write(pvalues\,\"$R_out\")";
close OUTFILE;
`$R < $R_script --no-save --slave`;
open INFILE, "<$R_out";
my $file = <INFILE>;
my @results = split(/\s+/,$file);
`rm $R_script $R_out`;
close INFILE;
return($results[0]);
}


sub fisher_test {
my $R_script = "STGD/scripts/R_$ARGV[0].R";
my $R_out = "STGD/scripts/R_$ARGV[0].out";
#print "enter\n";
open (OUTFILE, ">$R_script");
print OUTFILE	"Tea = matrix(c($_[0]\,$_[1]\,$_[2]\,$_[3])\,nrow=2\,dimnames =list(Guess = c(\"Milk\"\,\"Tea\")\,Truth = c(\"Milk\"\,\"Tea\")))\n";
print OUTFILE	"pvalues<-fisher.test(Tea\, alternative=\"$_[4]\")\n";
print OUTFILE "write(pvalues\$p\.value,(\"$R_out\"))\n";
close OUTFILE;
`$R < $R_script --no-save --slave`;
open INFILE, "<$R_out";
my $file = <INFILE>;
my @results = split(/\s+/,$file);
`rm $R_script $R_out`;
close INFILE;
return($results[0]);
}
