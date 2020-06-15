#!/usr/bin/perl -w
use strict;
my $main = "/home";
#my $pt_freq = $ARGV[0]; # ac/(number of pt *2)
# $pt_freq = $pt_freq *2;
my $R_comd = "/storage/chen/Software/R-3.5.0/bin/R";
#my $dom_normal_freq = 0.00005;
#my $rec_normal_freq = 0.0005;
#my $size = 50;
my $size = 100;
my $label = "conv";
my @dom_normal_freq = (1/1000000,2/1000000,3/1000000,4/1000000,5/1000000,6/1000000,7/1000000,8/1000000,9/1000000,0.00001,0.00002,0.00003,0.00004,0.00005,0.00006,0.00007,0.00008,0.00009,0.0001,0.00015,0.0002,0.0003,0.0004,0.0005,0.001);
my @rec_normal_freq = (1/1000000,5/100000,7/100000,1/10000,2/10000,1.5/10000,3/10000,0.0004,0.0005,0.0006,0.0007,0.0008,0.001,0.0015,0.002,0.003,0.005);

my @dis_freq = (1/200,1/1000,1/2000,1/5000,1/10000);

my $output1 = "HGMD/data/power_test_af_cutoff_dom_test_new_ground_truth_new_value_$size"."_$label";
open(OUTPUT1,">$output1");

print OUTPUT1 "AF\t";

        for my $dom_normal_freq (@dom_normal_freq){
            print OUTPUT1 "$dom_normal_freq\t";

        }
            print OUTPUT1 "\n";

	for my $dis_freq (@dis_freq){  # allele freq
		 print OUTPUT1 "$dis_freq\t";

#             for my $s (@size){   # pt cohort aize
	       for my $dom_normal_freq (@dom_normal_freq){ 
                #  my $size = $s;
                   
                   my $power = binProb($size,$dom_normal_freq,$dis_freq,"dom_power");
                   print OUTPUT1 "$power\t";
               }
                   print OUTPUT1 "\n";
         }
close(OUTPUT1);
##############
#@normal_freq = (0.000005,0.00001,0.00002,0.00005, 0.0001);

$output1 = "HGMD/data/spec_test_af_cutoff_dom_test_new_ground_truth_new_value_$size"."_$label";
open(OUTPUT1,">$output1");
print OUTPUT1 "AF\t";

#        for my $s (@size){
 #           print OUTPUT1 "$s\t";

  #      }

 for my $dom_normal_freq (@dom_normal_freq){
            print OUTPUT1 "$dom_normal_freq\t";

        }
            print OUTPUT1 "\n";

	for my $dis_freq (@dis_freq){  # allele freq
		 print OUTPUT1 "$dis_freq\t";

         #    for my $s (@size){   # pt cohort aize
         #          my $size = $s;
            for my $dom_normal_freq (@dom_normal_freq){         
                   my $power = binProb($size,$dom_normal_freq,$dis_freq,"dom_spec");
                   print OUTPUT1 "$power\t";
               }
                   print OUTPUT1 "\n";
         }
close(OUTPUT1);

####################
#@normal_freq = (0.0001, 0.0002, 0.0005, 0.001, 0.005);
# @dis_freq = (0.001, 0.0005, 0.0001, 0.00005,0.00001);
$size = $size*2;
@dis_freq = (1/2000, 1/5000,1/10000, 1/20000, 1/100000);

$output1 = "HGMD/data/power_test_af_cutoff_rec_test_new_ground_truth_new_value_$size"."_$label";
open(OUTPUT1,">$output1");

print OUTPUT1 "AF\t";

       # for my $s (@size){
        #    print OUTPUT1 "$s\t";

        #}
          for my $rec_normal_freq (@rec_normal_freq){
            print OUTPUT1 "$rec_normal_freq\t";

        }

            print OUTPUT1 "\n";

	for my $dis_freq (@dis_freq){  # allele freq
		 print OUTPUT1 "$dis_freq\t";
		for my $rec_normal_freq (@rec_normal_freq){
#             for my $s (@size){   # pt cohort aize
                   #$size = $s * 2;
                   
                   my $power = binProb($size,$rec_normal_freq,$dis_freq,"rec_power");
                   print OUTPUT1 "$power\t";
               }
                   print OUTPUT1 "\n";
         }
close(OUTPUT1);
################


# @normal_freq = (0.00005, 0.0001, 0.0002, 0.0005,0.0008, 0.001);
$output1 = "HGMD/data/spec_test_af_cutoff_rec_test_new_ground_truth_new_value_$size"."_$label";
open(OUTPUT1,">$output1");
print OUTPUT1 "AF\t";

#        for my $s (@size){
 #           print OUTPUT1 "$s\t";

  #      }
 
  for my $rec_normal_freq (@rec_normal_freq){
            print OUTPUT1 "$rec_normal_freq\t";

        }

            print OUTPUT1 "\n";

	for my $dis_freq (@dis_freq){  # allele freq
		 print OUTPUT1 "$dis_freq\t";

#             for my $s (@size){   # pt cohort aize
 #                  my $size = $s*2;
  
for my $rec_normal_freq (@rec_normal_freq){                 
                   my $power = binProb($size,$rec_normal_freq,$dis_freq,"rec_spec");
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
my $R_script = "STGD/scripts/R_pt_AD_dis$ARGV[0].R";
my $R_out = "STGD/scripts/R_pt_AD_dis$ARGV[0].out";

open (OUTFILE, ">$R_script");
print OUTFILE "set.seed(123)\n";
if($_[3] eq "dom_power"){
print OUTFILE "g_array<-dbinom(0:$_[0]\, $_[0]\, ($_[1]\*2))\n";
}elsif($_[3] eq "rec_power"){
print OUTFILE "g_array<-dbinom(0:$_[0]\, $_[0]\, $_[1])\n";
}elsif($_[3] eq "dom_spec"){
print OUTFILE "g_array<-dbinom(0:$_[0]\, $_[0]\, $_[1]/\(1-sqrt\(1-$_[2]\)\)\)\n";
}elsif($_[3] eq "rec_spec"){
print OUTFILE "g_array<-dbinom(0:$_[0]\, $_[0]\, $_[1]/sqrt\($_[2]\)\)\n";
}
print OUTFILE "b=c(0:$_[0])\n";
#####if($_[3] eq "dom_power"){
#####print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/\(1-sqrt\(1-$_[2]\)\)\)<=0.05,1,0)\n";
####}elsif($_[3] eq "dom_spec"){
####print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/\(1-sqrt\(1-$_[2]\)\)\)>0.05,1,0)\n";
####}elsif($_[3] eq "rec_power"){
####print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/sqrt\($_[2]\)\)<=0.05,1,0)\n";
####}elsif($_[3] eq "rec_spec"){
####print OUTFILE "factor=ifelse(pbinom\(b, size=$_[0], prob=$_[1]/sqrt\($_[2]\)\)>0.05,1,0)\n";
####}
#if($_[3] eq "dom_power"){
#print OUTFILE "sd=($_[1]*2 - $_[1]/(1-sqrt(1-$_[2])))/-2.326\n";
#print OUTFILE "cutoff=qnorm(0.05,mean=$_[1]/(1-sqrt(1-$_[2])),sd=sd)\n";
#print OUTFILE "factor=ifelse(pbinom(b, size=$_[0], prob=cutoff   )<=0.05,1,0)\n";
#}elsif($_[3] eq "dom_spec"){
#print OUTFILE "sd=($_[1]*2 - $_[1]/(1-sqrt(1-$_[2])))/-2.326\n";
#print OUTFILE "cutoff=qnorm(0.05,mean=$_[1]/(1-sqrt(1-$_[2])),sd=sd)\n";
#print OUTFILE "factor=ifelse(pbinom(b, size=$_[0], prob=cutoff  )>0.05,1,0)\n";
#}elsif($_[3] eq "rec_power"){
#print OUTFILE "sd=($_[1] - $_[1]/sqrt($_[2]))/-2.326\n";
#print OUTFILE "cutoff=qnorm(0.05,mean=$_[1]/sqrt($_[2]),sd=sd)\n";
#print OUTFILE "factor=ifelse(pbinom(b, size=$_[0], prob=cutoff )<=0.05,1,0)\n";
#}elsif($_[3] eq "rec_spec"){
#print OUTFILE "sd=($_[1] - $_[1]/sqrt($_[2]))/-2.326\n";
#print OUTFILE "cutoff=qnorm(0.05,mean=$_[1]/sqrt($_[2]),sd=sd)\n";
#print OUTFILE "factor=ifelse(pbinom(b, size=$_[0], prob=cutoff   )>0.05,1,0)\n";

#}
#print OUTFILE "val=sum(g_array\*factor)\n";
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
if($_[3] =~ /power/){
print OUTFILE "factor2=factor*factor1\n";
}elsif($_[3] =~ /spec/){
print OUTFILE "factor2=ifelse((factor+factor1)>0,1,0)\n";
}
print OUTFILE "val=sum(g_array\*factor2)\n";

print OUTFILE "write(val\,ncolumns=$_[0]\,file=\"$R_out\")\n";
close OUTFILE;
`$R_comd < $R_script --no-save --slave`;
#exit;
open INFILE, "<$R_out";
my $file = <INFILE>;
my @results = split(/\s+/,$file);
`rm $R_script $R_out`;
close INFILE;
return($results[0]);
}
