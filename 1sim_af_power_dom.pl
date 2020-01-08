#!/usr/bin/perl -w
use strict;
my $main = "/home";
my $pt_freq = $ARGV[0]; # ac/(number of pt *2)
 $pt_freq = $pt_freq *2;
my $output = "HGMD/data/power_test_pf_dom_$ARGV[0]"."_test_new_ground_truth_new";
my $output1 = "HGMD/data/power_test_pf_dom_$ARGV[0]"."_test_new_ground_truth_new_value";
my $R_comd = "R-3.5.0/bin/R";
my $dis_freq = 0.001;
open(OUTPUT,">$output");
open(OUTPUT1,">$output1");
print OUTPUT "#$pt_freq\n";
print OUTPUT1 "#$pt_freq\n";

my @size = (10,20,40,80,100,120,200,250,300,400,500,600,700,800,900,1000);
my @normal_freq = (0.00001,0.00005,0.0001,0.0005);

print OUTPUT "AF\t";
print OUTPUT1 "AF\t";

        for my $s (@size){
            print OUTPUT "$s\t";
            print OUTPUT1 "$s\t";

        }
            print OUTPUT "\n";
            print OUTPUT1 "\n";

	for my $i (@normal_freq){  # allele freq
		 print OUTPUT "$i\t";
		 print OUTPUT1 "$i\t";

             for my $s (@size){   # pt cohort aize
                   my $size = $s;
                   my @allele_x_obs=@{rbinom(50, $size, $pt_freq)};
                   my $n=0;
                   my $normal_freq = $i; 
                   my $pt_freq_exp = $normal_freq / (1-sqrt(1-$dis_freq)); 
                   my $obser_exp = int($pt_freq*$size);
                   my $p_origin=bin_test($pt_freq_exp, $obser_exp, $size, "less");

                   for my $allele_x_obs (@allele_x_obs){
                   my $p = bin_test($pt_freq_exp,$allele_x_obs, $size,"less");
                   if($p < 0.05){
                   $n++;
                   }
                   }
                   my $power = $n/50;
                   print OUTPUT "$power:$p_origin\t";
                   print OUTPUT1 "$power\t";
                   print "$pt_freq_exp\t$obser_exp\t$size\t$p_origin\n"; 

               }
                   print OUTPUT "\n";
                   print OUTPUT1 "\n";
         }

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

sub rbinom {
my $R_script = "STGD/scripts/R_pt_AD$ARGV[0].R";
my $R_out = "STGD/scripts/R_pt_AD$ARGV[0].out";

open (OUTFILE, ">$R_script");
print OUTFILE "g_array<-rbinom($_[0]\,$_[1]\, $_[2])\n";
print OUTFILE "write(g_array\,ncolumns=$_[0]\,file=\"$R_out\")\n";
close OUTFILE;
`$R_comd < $R_script --no-save --slave`;

open INFILE, "<$R_out";
my $file = <INFILE>;
my @results = split(/\s+/,$file);
`rm $R_script $R_out`;
close INFILE;
return(\@results);
}
