#!/usr/bin/perl -w
use strict;
my $main = "/home";
#my $pt_freq = $ARGV[0]; # ac/(number of pt *2)
# $pt_freq = $pt_freq *2;
my $R_comd = "/storage/chen/Software/R-3.5.0/bin/R";
#my $dom_dis_freq = 0.001;
my $dom_dis_freq = 0.001;
my $rec_dis_freq = 0.0001;
#my @freq_factor = (0.1, 0.2, 0.5,0.7, 0.9, 1,2,5,10);
my @freq_factor = (0.1, 0.2, 0.5,0.7, 0.9, 1,2,5);
my $label = "conv";
my @size = (50,80,100,120,200,250,300,400,500,600,700,800,900,1000);
#my $normal_freq = 0.00005;
my $normal_freq = 0.0001;


my $output1 = "HGMD/data/power_test_pf_sampling_bias_dom_test_new_ground_truth_new_value_$normal_freq"."_$label";
open(OUTPUT1,">$output1");

print OUTPUT1 "AF\t";

        for my $s (@size){
            print OUTPUT1 "$s\t";

        }
            print OUTPUT1 "\n";

	for my $freq_factor (@freq_factor){  # allele freq
		 print OUTPUT1 "$freq_factor\t";

             for my $s (@size){   # pt cohort aize
                   my $size = $s;
                   
                   my $power = binProb($size,$normal_freq,$dom_dis_freq,"dom_power",$freq_factor);
                   print OUTPUT1 "$power\t";
               }
                   print OUTPUT1 "\n";
         }
close(OUTPUT1);
##############
#@normal_freq = (0.000005,0.00001,0.00002,0.00005, 0.0001);

$output1 = "HGMD/data/spec_test_pf_sampling_bias_dom_test_new_ground_truth_new_value_$normal_freq"."_$label";
open(OUTPUT1,">$output1");
print OUTPUT1 "AF\t";

        for my $s (@size){
            print OUTPUT1 "$s\t";

        }
            print OUTPUT1 "\n";

	for my $freq_factor (@freq_factor){  # allele freq
		 print OUTPUT1 "$freq_factor\t";

             for my $s (@size){   # pt cohort aize
                   my $size = $s;
                   
                   my $power = binProb($size,$normal_freq,$dom_dis_freq,"dom_spec",$freq_factor);
                   print OUTPUT1 "$power\t";
               }
                   print OUTPUT1 "\n";
         }
close(OUTPUT1);

####################
#@normal_freq = (0.0001, 0.0002, 0.0005, 0.001, 0.005);
#$normal_freq = 0.0005;
$normal_freq = 0.001;

$output1 = "HGMD/data/power_test_pf_sampling_bias_rec_test_new_ground_truth_new_value_$normal_freq"."_$label";
open(OUTPUT1,">$output1");

print OUTPUT1 "AF\t";

        for my $s (@size){
            print OUTPUT1 "$s\t";

        }
            print OUTPUT1 "\n";

	for my $freq_factor (@freq_factor){  # allele freq
		 print OUTPUT1 "$freq_factor\t";

             for my $s (@size){   # pt cohort aize
                   my $size = $s * 2;
                   
                   my $power = binProb($size,$normal_freq,$rec_dis_freq,"rec_power",$freq_factor);
                   print OUTPUT1 "$power\t";
               }
                   print OUTPUT1 "\n";
         }
close(OUTPUT1);
################


 #@normal_freq = (0.00005, 0.0001, 0.0002, 0.0005,0.0008, 0.001);
$output1 = "HGMD/data/spec_test_pf_sampling_bias_rec_test_new_ground_truth_new_value_$normal_freq"."_$label";
open(OUTPUT1,">$output1");
print OUTPUT1 "AF\t";

        for my $s (@size){
            print OUTPUT1 "$s\t";

        }
            print OUTPUT1 "\n";

	for my $freq_factor (@freq_factor){  # allele freq
		 print OUTPUT1 "$freq_factor\t";

             for my $s (@size){   # pt cohort aize
                   my $size = $s*2;
                   
                   my $power = binProb($size,$normal_freq,$rec_dis_freq,"rec_spec",$freq_factor);
                   print OUTPUT1 "$power\t";
               }
                   print OUTPUT1 "\n";
         }
close(OUTPUT1);

###############


sub bin_test {
my $R_script = "STGD/scripts/R_pt_AD$ARGV[0].R";
my $R_out = "STGD/scripts/R_pt_AD$ARGV[0].out";

open (OUTFILE, ">$R_script");
print OUTFILE "pvalues<-binom.test(", "$_[1]",",$_[2]",", p= $_[0]",", alternative=\"$_[3]\" ,conf.level= 0.95)\$p\.value\n";
print OUTFILE "write(pvalues\,\"$R_out\")\n";
close OUTFILE;
`$R_comd < $R_script --no-save --slave`;
open INFILE, "<$R_out";
my $file = <INFILE>;
my @results = split(/\s+/,$file);
`rm $R_script $R_out`;
close INFILE;
return($results[0]);
}

sub binProb {
my $R_script = "STGD/scripts/R_pt_AD$ARGV[0].R";
my $R_out = "STGD/scripts/R_pt_AD$ARGV[0].out";

open (OUTFILE, ">$R_script");
print OUTFILE "set.seed(123)\n";
if($_[3] eq "dom_power"){
print OUTFILE "g_array<-dbinom(0:$_[0]\, $_[0]\, ($_[1]\*2\*$_[4]))\n";
}elsif($_[3] eq "rec_power"){
print OUTFILE "g_array<-dbinom(0:$_[0]\, $_[0]\, $_[1]\*$_[4])\n";
}elsif($_[3] eq "dom_spec"){
print OUTFILE "g_array<-dbinom(0:$_[0]\, $_[0]\, $_[1]\*$_[4]/\(1-sqrt\(1-$_[2]\)\)\)\n";
}elsif($_[3] eq "rec_spec"){
print OUTFILE "g_array<-dbinom(0:$_[0]\, $_[0]\, $_[1]\*$_[4]/sqrt\($_[2]\)\)\n";
}
print OUTFILE "b=c(0:$_[0])\n";

####if($_[3] eq "dom_power"){
####print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/\(1-sqrt\(1-$_[2]\)\)\)<=0.05,1,0)\n";
####}elsif($_[3] eq "dom_spec"){
####print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/\(1-sqrt\(1-$_[2]\)\)\)>0.05,1,0)\n";
####}elsif($_[3] eq "rec_power"){
####print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/sqrt\($_[2]\)\)<=0.05,1,0)\n";
####}elsif($_[3] eq "rec_spec"){
####print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/sqrt\($_[2]\)\)>0.05,1,0)\n";
####}
if($_[3] eq "dom_power"){
print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/\(1-sqrt\(1-$_[2]\)\)\)<=0.05,1,0)\n";
}elsif($_[3] eq "dom_spec"){
print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/\(1-sqrt\(1-$_[2]\)\)\)>0.05,1,0)\n";
}elsif($_[3] eq "rec_power"){
print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/sqrt\($_[2]\)\)<=0.05,1,0)\n";
}elsif($_[3] eq "rec_spec"){
print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/sqrt\($_[2]\)\)>0.05,1,0)\n";
}
print OUTFILE "c=c(-1:($_[0]-1))\n";
if($_[3] eq "dom_power"){
print OUTFILE "factor1=ifelse(pbinom(c, size=$_[0], prob=$_[1]*2,lower.tail=F)>0.05,1,0)\n";
}elsif($_[3] eq "dom_spec"){
print OUTFILE "factor1=ifelse(pbinom\(c, size=$_[0], prob=$_[1]*2,lower.tail=F\)<=0.05,1,0)\n";
}elsif($_[3] eq "rec_power"){
print OUTFILE "factor1=ifelse(pbinom\(c, size=$_[0], prob=$_[1],lower.tail=F\)>0.05,1,0)\n";
}elsif($_[3] eq "rec_spec"){
print OUTFILE "factor1=ifelse(pbinom\(c, size=$_[0], prob=$_[1],lower.tail=F\)<=0.05,1,0)\n";
}
#if($_[3] eq "dom_power"){
#print OUTFILE "factor1=ifelse(pbinom\(b, size=$_[0], prob=$_[1]*2\)<0.95,1,0)\n";
#}elsif($_[3] eq "dom_spec"){
#print OUTFILE "factor1=ifelse(pbinom\(b, size=$_[0], prob=$_[1]*2\)>=0.95,1,0)\n";
#}elsif($_[3] eq "rec_power"){
#print OUTFILE "factor1=ifelse(pbinom\(b, size=$_[0], prob=$_[1]\)<0.95,1,0)\n";
#}elsif($_[3] eq "rec_spec"){
#print OUTFILE "factor1=ifelse(pbinom\(b, size=$_[0], prob=$_[1]\)>=0.95,1,0)\n";
#}
if($_[3] =~ /power/){
print OUTFILE "factor2=factor*factor1\n";
}elsif($_[3] =~ /spec/){
print OUTFILE "factor2=ifelse((factor+factor1)>0,1,0)\n";
}

print OUTFILE "val=sum(g_array\*factor2)\n";
print OUTFILE "write(val\,ncolumns=$_[0]\,file=\"$R_out\")\n";
close OUTFILE;
`$R_comd < $R_script --no-save --slave`;

open INFILE, "<$R_out";
my $file = <INFILE>;
my @results = split(/\s+/,$file);
`rm $R_script $R_out`;
close INFILE;
return($results[0]);
}
