#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Find;
use File::Path;
use Cwd;
use File::Basename;
#use diagnostics;

 ########################### WELCOME  ##########################################
#                                                                   			#
#                          sgRNAcas9(Cas12a/cpf1)                               # 
# ---a tool for fast designing CRISPR sgRNA with high specificity   			#
#                                                                        		#
# AUTHOR  : Zhongping Xu                                                        #
# Email   : zhongpingxu\@aliyun.com                                 			#
# Homepage: http://tiramisutes.github.io/2017/01/13/CRISPR-Designer.html        #
#                 												     			#
# Huazhong Agricultural University   											#
# Version: sgRNAcas9_3.0.5                                          			#
# Begin       : 2013.12.9                                           			#
# LAST REVISED: 2020.7.13                                           			#
 ###############################################################################

my $oldtime = time();

my ($Inputfile_Fasta, $truncat, $GC_l, $GC_m, $Genome, $Option, $Type, $Seqmap_vesion, $Num_mismatch, $offset_s, $offset_e, $path);

GetOptions( "i=s" => \$Inputfile_Fasta,        #Input file
            "x=i" => \$truncat,                #Length of sgRNA[23]
            "l=i" => \$GC_l,                   #The minimum value of GC content [35]
			"m=i" => \$GC_m,                   #The maximum value of GC content [65] 
	        "g=s" => \$Genome,                 #The reference genome sequence
			"o=s" => \$Option,                 #Searching CRISPR target sites using DNA strands based option(s/a/b)
			"t=s" => \$Type,                   #Type of gRNA searching mode(s/p) 
			"v=s" => \$Seqmap_vesion,          #Operation system [w, for windows; l, for linux-64; u, for linux-32;  m, for MacOSX-64; a, for MacOSX-32]
			"n=i" => \$Num_mismatch,           #Maximum number of mismatches [5]
			"s=i" => \$offset_s,               #The minimum value of sgRNA offset [-2]
            "e=i" => \$offset_e,               #The maximum value of sgRNA offset [32]
            "p=s" => \$path,                   #Output path, this path must exists
          );

#default
$truncat ||= "23";  
$GC_l ||= "35";                                #35 % < GC% < 65 %
$GC_m ||= "65";
$Seqmap_vesion ||= "l";                        #linux
$Num_mismatch ||="5";                          #number of mismatches: 5
$offset_s ||="-3";                             #sgRNA offset: -2 to 32 bp
$offset_e ||="33";
my $dir_default = getcwd;                      #default output
$path ||= $dir_default;

my $dir =$path;
my $faname = basename($Inputfile_Fasta);

mkdir("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname",0755)||die "Can't create directory: Directory exists at $dir. Please delete, move or rename the exist directory before you run this program.$!" ;

open  (LOG, ">>$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/sgRNAcas9.Log.txt") || die "Can't open sgRNAcas9.Log.txt for writing!" ."\n";
print  LOG "################################# Log ###########################################".        "\n\n";
#print "Writing Log information.                                                                      " ."\n";
print  LOG "#                              sgRNAcas9                                                  " ."\n";
print  LOG "#     ---a tool for fast designing CRISPR sgRNA with high specificity                     " ."\n";          
print  LOG "#                                                                                         " ."\n";
print  LOG "#       contact:  Xie Shengsong, Email: ssxieinfo\@gmail.com                                .\n\n";
######


#mkdir("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/A.Final_report",0755)||die "can't create directory: $!" ;

print "\n\tWelcome to sgRNAcas9\n";
print "\t---a tool for fast designing CRISPR sgRNA with high specificity\n";
print "\t---------------------------------------------------------\n";
print "Version   : 3.0"."\n";
print "Copyright : Free software"."\n";
print "Author    : Shengsong Xie"."\n";
print "Email     : ssxieinfo\@gmail.com"."\n";
print "Homepage  : www.biootools.com"."\n";

my $local_time;
$local_time = localtime();
print "Today     : $local_time\n\n";
print  LOG "# Time, begin at $local_time."."\n";
print  LOG "# Usage: perl $0 -i $Inputfile_Fasta -x $truncat -l $GC_l -m $GC_m -g $Genome -o $Option -t $Type -v $Seqmap_vesion -n $Num_mismatch -s $offset_s -e $offset_e -p $path 2>log.txt\n\n";

################################### format seq ###################################较为可简单，代码完全理解，基础代码，不需要修改
print "Start sgRNAcas9 program........\n";
print "Step1: Format target sequences.\n";
print  LOG "# Start sgRNAcas9 program........\n";
print  LOG "# Step1: Format target sequences.\n";

open (Inseq, $Inputfile_Fasta) || die "Can't open $Inputfile_Fasta for reading!\n";
open (FASTA, ">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/TargetSeq.fa") || die "Can't open TargetSeq.fa for writing!\n";

my $TmpTit="";
my $TmpSeq="";
my $TmpTit_S="";
my $TmpTit_A="";
my $TmpSeq_Acomp;
my $TmpSeq_ARevcomp;

if ($Option eq "b") {              #B=both:sense and anti-sense strand

    while (<Inseq>) {

	  chomp $_; 
	 
      if (/^>(\S+)/){
	   
	    if ($TmpTit && $TmpSeq) {
	
			 $TmpTit_S=$TmpTit."_S";
		     $TmpTit_A=$TmpTit."_A";
		
			 $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
		     $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/;  	
		     $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
				
		     print FASTA ">$TmpTit_S\n$TmpSeq\n";	
	         print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n";
		    }
	
	         $TmpTit=$1;	
	         $TmpSeq="";
	
	   }elsif(/(\w+)/){
	
	         $TmpSeq.=$1;
		
	  }
	 
	}

       if ($TmpTit && $TmpSeq) {

	        my $TmpTit_S=$TmpTit."_S";
            my $TmpTit_A=$TmpTit."_A";

	        my $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
            $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/;  	
            my $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
     
            print FASTA ">$TmpTit_S\n$TmpSeq\n";
            print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n"; 
	
  }
}elsif ($Option eq "a") {           #Non-template: A: anti-sense strand

	while (<Inseq>) {

		chomp $_;

	if (/^>(\S+)/){
	
	    if ($TmpTit && $TmpSeq) {
		  	  
		     $TmpTit_A=$TmpTit."_A";
		
			 $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
		     $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/;  	
		     $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
				
		     print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n";	
		  }
	
	         $TmpTit=$1;	
	         $TmpSeq="";
	
	  }elsif(/(\w+)/){
	
	         $TmpSeq.=$1;
		
	  }
	}

       if ($TmpTit && $TmpSeq) {

            my $TmpTit_A=$TmpTit."_A";

            my $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
            $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/; 
            my $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
   
            print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n";
   
  }
}elsif ($Option eq "s" || $Option eq "") {     #Template: S: sense strand

	while (<Inseq>) {

		chomp $_;

	if (/^>(\S+)/){
	
	   if ($TmpTit && $TmpSeq) {
		  	  
		    $TmpTit_S=$TmpTit."_S";
				
		    print FASTA ">$TmpTit_S\n$TmpSeq\n";	
		  
	    }
	
	        $TmpTit=$1;	
	        $TmpSeq="";
	
	  }elsif(/(\w+)/){
	
	    	$TmpSeq.=$1;
		
	  } 
	}

      if ($TmpTit && $TmpSeq) {
	
		   my $TmpTit_S=$TmpTit."_S";

           print FASTA ">$TmpTit_S\n$TmpSeq\n";
    
  }
}
close(Inseq);
close(FASTA);

########################### Find CRISPR targets-single #################################已经理解
print "Step2: Find CRISPR targets.\n";
print  LOG "# Step2: Find CRISPR targets.\n";

open(IntPut, "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/TargetSeq.fa") || die "Can't open TargetSeq.fa for reading!\n";
#report_protospacer_single.txt为CRISPR.targets_S.txt和CRISPR.targets_A.txt文件合并结果。
open(OutPut, ">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/report_protospacer_single.txt") || die "Can't open report_protospacer_single.txt for writing!\n";
open(OutPut1, ">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_single.fa") || die "Can't open CRISPR.targets_single.fa for writing!\n";
#open(OutPut2, ">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/full_length_sgRNA.txt") || die "Can't open full_length_sgRNA.txt for writing!\n";
open(OutPut3, ">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_S.txt") || die "Can't open CRISPR.targets_single.fa for writing!\n";
open(OutPut4, ">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_A.txt") || die "Can't open CRISPR.targets_single.fa for writing!\n";

#print OutPut "Candidate sgRNA, Pattern: GGX18NGG, GX19NGG, X20NGG, GC% >=$GC %\n\n";
print OutPut "sgRID\t"."Start\t"."End\t"."CRISPR_target_sequence(5'-3')\t"."Length(nt)\t"."GC%\n";

my $ID=""; 
my $seq="";
my $probe_id="";

while(<IntPut>) {
	chomp $_;

	if (/^>(\S+)/){   #\S 表示匹配非空白字符，即fa文件中序列ID号
		analysis(); $ID=$1;

  }else{
		$seq=$_;
		
  }
}
analysis();
close OutPut;
close OutPut1;
close OutPut3;
close OutPut4;
#至此生成report_protospacer_single.txt，CRISPR.targets_S.txt，CRISPR.targets_A.txt和CRISPR.targets_single.fa四个文件。
############################# Find CRISPR targets-pairs #################################
open( PA, "<$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_A.txt" ) || die "can't open CRISPR.targets_A.txt!";
open( PS, "<$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_S.txt" ) || die "can't open CRISPR.targets_S.txt!";
open( PAIRS1, ">>$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/report_protospacer_pairs.xls" ) || die "can't open report_protospacer_pairs.xls!";
open( PAIRS2, ">>$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_pairs.fa" ) || die "can't open CRISPR.targets_pairs.fa!";

my $sgRID_A;
my $Start_A;
my $End_A;
my $target_seq_A ="";
my $Pattern_A;
my $GC_A;

my $sgRID_S;
my $Start_S;
my $End_S;
my $target_seq_S ="";
my $Pattern_S;
my $GC_S;

my $len_target_seq_A=0;
my $len_target_seq_S=0;

my $ID_A;
my $ID_S;

print PAIRS1 "\t\tPaired-gRNA\n";
print PAIRS1 "sgRID_S\ttarget_seq_S\tStart_S\tEnd_S\tGC%_S\t\tsgRID_A\ttarget_seq_A\tStart_A\tEnd_A\tGC%_A\tsgRNA_offset(bp)\n";

while ( <PA> ) {
	chomp $_; 
	(my $sgRID_A, my $Start_A, my $End_A, my $target_seq_A,	my $Pattern_A, my $GC_A, my $emp_A)=split/\t/, $_;

	$len_target_seq_A = length($target_seq_A);   #计算sgRNA靶标碱基长度

	next if $len_target_seq_A ne ($truncat+4);  #长度相同则继续后面程序，不等就退出；cas9和cpf1因PAM长度不同，在此处需要修改后面数字

	$ID_A = $sgRID_A;
	$ID_A=~ s/_A_(\d+)//m;
	
	seek PS, 0, 0;

	while ( <PS> ) {
    chomp $_; 
		(my $sgRID_S, my $Start_S, my $End_S, my $target_seq_S,	my $Pattern_S, my $GC_S, my $emp_S)=split/\t/, $_;

		next if $target_seq_S eq "";
		$len_target_seq_S = length($target_seq_S);
		next if $len_target_seq_S ne ($truncat+4);               

		$ID_S = $sgRID_S;
		$ID_S=~ s/_S_(\d+)//m;   #匹配Gh_D08G0270_S_1中的Gh_D08G0270 ID号
		
		my $offset_value = $Start_S -$End_A;  #是否成对计算公式

		if (($ID_A eq $ID_S) and ($offset_value > "$offset_s" and $offset_value < "$offset_e")) {   #$offset_s和$offset_e为参数指定值：-2 to 32 bp or 5 to 35 bp

			print PAIRS1 "$sgRID_A"."\t"."$target_seq_A"."\t"."$Start_A"."\t"."$End_A"."\t"."$GC_A"."\t<->\t";
			print PAIRS1 "$sgRID_S"."\t"."$target_seq_S"."\t"."$Start_S"."\t"."$End_S"."\t"."$GC_S"."\t"."$offset_value"."\n";
			print PAIRS2 ">"."$sgRID_A"."\n"."$target_seq_A"."\n";
			print PAIRS2 ">"."$sgRID_S"."\n"."$target_seq_S"."\n";

		}
  }
}
close(PA);
close(PS);
close(PAIRS1);
close(PAIRS2);

########################### Unique pairs sgR fasta seq ######################
my %seq;
my $title;
my $infile="$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_pairs.fa";

open (IN,"$infile") || die "can't open $infile!";
while (<IN>){
	$_=~s/\n//;
	$_=~s/\r//;
	if ($_=~/>/){
		$title=$_;
		$title=~s/>//;
	}
	else{
		$seq{$_}=$title;
	}
}
close IN;
#remove the abundant sequences
my @seq=keys (%seq);
my @uniqueseq;
my $find=0;
foreach (@seq){
	$find=0;
	my $seq=uc($_);
	foreach (@uniqueseq){
		if ($seq=~/$_/){
			$_=$seq;#replace with longer seq
			$find=1;
		}
		if ($_=~/$seq/){
			$find=1;
		}
	}
	if ($find==0){
		push @uniqueseq,$seq;
	}
}
#outout the final result
open (OUT,">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/unique_pairs.fa");
foreach (@uniqueseq){
	print OUT ">$seq{$_}\n$_\n";
}
close OUT;

################################# seqmap：Mapping Millions Of Short Sequences To The Genome #########################找到的靶标短序列比对到基因组
#“s”模式用CRISPR.targets_single.fa靶标，"p"模式用unique_pairs.fa靶标
print "Step3: Evaluate CRISPR potential off-target effect.\n";
print "Step3-1: Whole genome mapping.\n\n";
print  LOG "# Step3: Evaluate CRISPR potential off-target effect.\n";
print  LOG "# Step3-1: Whole genome mapping.\n";

if ($Seqmap_vesion eq "l" and $Type eq "s") {       #Linux-64

	my @args = ("./Seqmap/seqmap-1.0.12-linux-64","$Num_mismatch", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_single.fa", "$Genome", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt", "/output_all_matches");
	system(@args);

	if($? == -1) {
		die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "l" and $Type eq "p") {

	my @args = ("./Seqmap/seqmap-1.0.12-linux-64","$Num_mismatch", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/unique_pairs.fa", "$Genome", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt", "/output_all_matches");
	system(@args);

	if($? == -1) {
		die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "w" and $Type eq "s") {   #windows 32, 64-bit

	my @args = ("./Seqmap/seqmap-1.0.12-windows.exe","$Num_mismatch", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_single.fa", "$Genome", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt", "/output_all_matches");
	system(@args);

	if($? == -1) {
		die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "w" and $Type eq "p") {   

	my @args = ("./Seqmap/seqmap-1.0.12-windows.exe","$Num_mismatch", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/unique_pairs.fa", "$Genome", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt", "/output_all_matches");
	system(@args);

	if($? == -1) {
		die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "m" and $Type eq "s") {    #macOSX-64

	my @args = ("./Seqmap/seqmap-1.0.12-mac-64","$Num_mismatch", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_single.fa", "$Genome", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
		die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "m" and $Type eq "p") {   

	my @args = ("./Seqmap/seqmap-1.0.12-mac-64","$Num_mismatch", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/unique_pairs.fa", "$Genome", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
	  die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "a" and $Type eq "s") {    #macOSX-32

	my @args = ("./Seqmap/seqmap-1.0.12-mac","$Num_mismatch", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_single.fa", "$Genome", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
	  die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "a" and $Type eq "p") {   

	my @args = ("./Seqmap/seqmap-1.0.12-mac","$Num_mismatch", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/unique_pairs.fa", "$Genome", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
	  die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "u" and $Type eq "s") {    #Linux-32

	my @args = ("./Seqmap/seqmap-1.0.12-linux","$Num_mismatch", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_single.fa", "$Genome", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
	  die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "u" and $Type eq "p") {   

	my @args = ("./Seqmap/seqmap-1.0.12-linux","$Num_mismatch", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/unique_pairs.fa", "$Genome", "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
	  die "system @args failed: $?";
	}

}

################################## count_sgR_Targetsite #########################从比对结果中挑选NGG和NAG靶标结果输出到report_sgRid.txt文件
print "\nFinished mapping candidate CRISPR target sequences to genome: $Genome.\n";
print "\nStep3-2: Count the No.of potential off-target cleavage sites.\n";
print  LOG "# Finished mapping candidate CRISPR target sequences to genome: $Genome.\n";
print  LOG "# Step3-2: Count the No.of potential off-target cleavage sites.\n";

open(Input, "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt") || die "Can't open seqmap_output.txt for reading!\n";
open(Out, ">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/report_sgRid.txt") || die "Can't open report_sgRid.txt for writing!\n";

print Out "sgRID\t"."Total_No.OT\n";

my %ID;
my $len1=$truncat-23;
my $len2=$truncat-22;
my $len3=$truncat-21;

while (<Input>) {

	chomp $_; 
	
	(my $trans_id, my $trans_coord, my $target_seq, my $probe_id, my $probe_seq, my $num_mismatch, my $strand)=split/\t/, $_; 
	#$target_seq为比对到基因组上的序列，$probe_seq为找到的靶标序列;下面脚本需要修改，作用是只保留NGG和NAG的结果
	
	next if ($trans_id eq "trans_id" or (substr($target_seq,$len1,1) ne "T") or (substr($target_seq,$len2,1) ne "T") or (substr($target_seq,$len3,1) ne "T") ) ;
    #                 匹配到表头                          第一个碱基不是T                          第二个碱基不是T                     第三个碱基不是T

    while ($probe_id=~ /(\w[\w-]*)/g) {    #\w：任意单词字符，[_0-9a-zA-Z]
		$ID{$1}++;
	}
}

foreach (sort {$ID{$a}<=>$ID{$b}}keys %ID) {
	
	next if $_=~ /probe_id/;
    
	print Out "$_\t$ID{$_}\n";

}
close Input;
close Out;

################################## format_seqmap #################################比对seqmap_output.txt结果中分类统计OT和POT脱靶情况
print "Step3-3: Format and sort whole genome mapping result.\n";
print  LOG "# Step3-3: Format and sort whole genome mapping result.\n";

mkdir("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/A.Sort_OT_byID",0755)||die "can't create directory: $!" ;

my $trans_id;
my $trans_coord;
my $target_seq;
my $probe_id1;
my $probe_seq;
my $num_mismatch;
my $New_num_mismatch;
my $strand;

my $sgR_ID;

my $countsgR;

my $i=0;
my $j=0;

my $len4= $truncat-23; #0
my $len5= $truncat-22; #1
my $len6= $truncat-21; #2
my $len7= $truncat-20; #3

open(Seq, "<$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt" ) || die "can't open seqmap_output.txt!";
open(sgR, "<$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/report_sgRid.txt" ) || die "can't open report_sgRid.txt!";

open(Out3, ">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/search_OT.txt") || die "Can't open search_OT.txt for writing!\n";
print Out3 "\n\t\t------------------------------------ off-target analysis(OT:off-target) -----------------------------------\n\n";

while (<Seq>) {
	chomp $_; 

	(my $trans_id, my $trans_coord, my $target_seq, my $probe_id1, my $probe_seq, my $num_mismatch, my $strand)=split/\t/, $_;

	my $New_num_mismatch = $num_mismatch;
   
	seek sgR, 0, 0;
	
	while (<sgR>) {
		chomp $_; 
	
		(my $sgR_ID, my $countsgR)=split/\t/, $_;

		if ($probe_id1=~/$sgR_ID$/m) {
			$j++;
			    
			my @seq1 = split //, $probe_seq;
			my @seq2 = split //, $target_seq;
			
			my @Iden;
			$i++;
			
			for my $n (0..@seq1-1) {
				if   ( ($seq1 [$n] eq 'A' && $seq2 [$n] eq 'A') 
				    || ($seq1 [$n] eq 'T' && $seq2 [$n] eq 'T') 
				    || ($seq1 [$n] eq 'C' && $seq2 [$n] eq 'C') 
				    || ($seq1 [$n] eq 'G' && $seq2 [$n] eq 'G') ) {
				
				push @Iden, '-';
				}
				else {
				
					push @Iden, lc($seq2[$n]);
				
				}
			}
		#next if 为假时只是退出某一次循环 
   		#discard  以下主要是为了去掉TTTN的PAM区错配比对的基因组比对结果
	    next if ($trans_id eq "trans_id" or (substr($target_seq,$len4,1) ne "T") or (substr($target_seq,$len5,1) ne "T") or (substr($target_seq,$len6,1) ne "T") ) ;
		
			open (Out2,">>$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/A.Sort_OT_byID/$sgR_ID.txt");
			
	   #  处理例如TTTA和TTTG这种PAM区的错配                   TTT-[N]                                        
	    if($num_mismatch > 0 and (substr($target_seq,$len7,1) ne substr($probe_seq,$len7,1))){
	
	    	$New_num_mismatch=$num_mismatch-1;  #TTTN中N允许不同，所以这几种情况序列比对上的错配不算作靶标错配，故错配值减去1
	
				print Out2 "$probe_id1\t"."POT$i\t"."@Iden"."\t$New_num_mismatch"."M\t"."$trans_id\t"."$trans_coord\t"."$strand\n";
	
	    	#print Out3  "\t\t"."20           13 12      8 7           1 N G G\n";
	    	print Out3 "$probe_id1\t"."@seq1"."\n"."POT$i\t\t"."@seq2"."\n"."POT$i\t\t"."@Iden"."\t$New_num_mismatch"."M\t"."$trans_id\t"."$trans_coord\t"."$strand\n";
	
			}else {
	
				print Out2 "$probe_id1\t"."POT$i\t"."@Iden"."\t$num_mismatch"."M\t"."$trans_id\t"."$trans_coord\t"."$strand\n";
	
	    	print Out3 "$probe_id1\t"."@seq1"."\n"."POT$i\t\t"."@seq2"."\n"."POT$i\t\t"."@Iden"."\t$num_mismatch"."M\t"."$trans_id\t"."$trans_coord\t"."$strand\n";
	
			}
    }
  }
}
close(Seq);
close(sgR);
close(Out2);
close(Out3);

################################### subroutine ###################################
sub analysis {  
	if($ID eq "" || $seq eq "") {
	  return();
	}
	#my $seq =~ s/\n//g;
	my $SEQRNA  = uc($seq);
	my $SEQDNA  = $SEQRNA;
	$SEQDNA  =~ tr/U/T/;

	my $len = length($seq);

	my $idx_s;
	my $i = 0;

	for($idx_s=-1; $idx_s<length($seq); $idx_s++) {

		my $lentruncat = $truncat+4;  #27
		my $lentruncat2 = $truncat+1; #24

		if (substr($SEQDNA, $idx_s+1, $lentruncat) =~ /TTT[ATCG]{$lentruncat2}/) {    #TTTN-X?  
			my $sgRidx_s = $idx_s+2;                                  #For sense strand,sgRidx_s表达的是所匹配sgRNA序列的起始位置
			my $sgR = $&;  # $& 表示与格式匹配的字符串,此处即表示匹配到的27bp的sgRNA序列
			#print $sgR."\n";
			my $SgR1 = substr($sgR,3,24);   #去掉PAM区TTTN的前三个TTT碱基
			#print $SgR1."\n";
			my $sgRAT = &AT($sgR);          # 调用子函数AT来计算AT含量
			#my $sgRGC = &GC($sgR);
			my $sgRGC = &GC($SgR1);

			#my $SgRsd = substr($sgR,8,12);                           #5-X12-NGG-3
			#my $SgRseed ="$SgRsd"."TTN";
			my $SgRNt = substr($sgR,4,$truncat);                      #truncat序列，TTTN后面的23个碱基
			my $sgRlen = length($sgR);
			my $sgRidx_e = $sgRidx_s+$sgRlen-1;                       #sgRidx_e表示的是所匹配sgRNA序列的终止位置
			$i++;
			
			my $sgRidx_A_s =$len-$sgRidx_e+1;                         #For anti-sense strand
			my $sgRidx_A_e =$len-$sgRidx_s+1;
			
			my $sgrna = $SgRNt;
			$sgrna =~ tr/T/U/;
            #         $sgR 表示匹配到的27bp序列，此处设置GC含量上下限
				if (!($sgR =~ /T{8,15}/g)){       #/g 修饰词声明一个全局匹配——也就是说，在该字串里匹配尽可能多的次数。8个T的GC含量为35%，15个的为65%
					
					if ($ID=~ /._S$/ and ($sgRGC >= $GC_l and $sgRGC < $GC_m)) {    #For sense strand
           
						print OutPut "$ID\_$i\t"."$sgRidx_s\t"."$sgRidx_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";   #输出为CRISPR.targets_S.txt文件形式
						print OutPut3 "$ID\_$i\t"."$sgRidx_s\t"."$sgRidx_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";
						
						print OutPut1 ">$ID\_$i\n"."$sgR\n";                                                    #输出为CRISPR.targets_single.fa文件形式
						#print OutPut2 ">$ID\_$i\n"."$sgRNA"."\n";

					}elsif ($ID=~ /._A$/ and ($sgRGC >= $GC_l and $sgRGC < $GC_m) ) {   #For anti-sense strand
           
						print OutPut "$ID\_$i\t"."$sgRidx_A_s\t"."$sgRidx_A_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";
						print OutPut4 "$ID\_$i\t"."$sgRidx_A_s\t"."$sgRidx_A_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";
						
						print OutPut1 ">$ID\_$i\n"."$sgR\n";
						#print OutPut2 ">$ID\_$i\n"."$sgRNA"."\n";
    			}
				}
		}
	}
}

########count_ot#########
my $Dir3_ot = "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/A.Sort_OT_byID" ;
my $Ext3_ot = "txt" ; #file type

my ($OTname_ot, $OTID_ot, $pattern_ot, $mismatch_ot,$chrosome_ot,$location_ot,$strand_ot);

opendir(DH, "$Dir3_ot") or die "Can't open: $!\n" ;
my @list_ot = grep {/$Ext3_ot$/ && -f "$Dir3_ot/$_" } readdir(DH) ;
closedir(DH) ;
chdir($Dir3_ot) or die "Can't cd dir: $!\n" ;

my $file_ot;

open(OutPut, ">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/report_count_total_ot.txt") || die "Can't open report_count_total_ot.txt for writing!\n";

print OutPut "ID\t0M(on-/off-)\t1M\t2M\t3M\t4M\t5M\tTotal OT number\n";  #cpf1_$faname_sgRNAcas9_report.xls结果中OT统计结果

foreach my $file_ot (@list_ot){
    open(FH, "$file_ot") or die "Can't open: $!\n" ;
	$file_ot =~ s/.txt//m;

my $j=0;
my $e=0;
my $count_0M=0;
my $count_1M=0;
my $count_2M=0;
my $count_3M=0;
my $count_4M=0;
my $count_5M=0;
my $line1_ot;
my $count_error;

while($line1_ot=<FH>){
    $j++;
	($OTname_ot, $OTID_ot, $pattern_ot, $mismatch_ot,$chrosome_ot,$location_ot,$strand_ot)=split(/\t/, $line1_ot);
	#pig_NR1H4_cds_A_4	POT264	- - - - - - - t - - - - - - - t g - - - g - -	3M	10	9370885	+
	if ($mismatch_ot eq '0M'){
	 ++$count_0M;
	} elsif ($mismatch_ot eq '1M'){
	 ++$count_1M;
	} elsif ($mismatch_ot eq '2M'){
	 ++$count_2M;
	} elsif ($mismatch_ot eq '3M'){
	 ++$count_3M;
	} elsif ($mismatch_ot eq '4M'){
	 ++$count_4M;
	} elsif ($mismatch_ot eq '5M'){
	 ++$count_5M;
	} else {
	 print OutPut "Error-Strange mismatch numbers:$mismatch_ot\n";
	 ++$count_error;
	}	
}

$e=$j-1;

    print OutPut $file_ot."\t";  #统计结果写入report_count_total_ot.txt文件
    print OutPut "$count_0M\t";
	print OutPut "$count_1M\t";
	print OutPut "$count_2M\t";
	print OutPut "$count_3M\t";
    print OutPut "$count_4M\t";
	print OutPut "$count_5M\t";

if ($count_0M == 0) {

	print OutPut "$j\n";

}else {

    print OutPut "$e\n";
}
}

close OutPut;

#########
#########
my $File1_ot = "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/report_protospacer_single.txt";
my $File2_ot = "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/report_count_total_ot.txt";

my ($sgRID_ot,	$Start_ot, $End_ot, $CRISPR_target_sequence_ot,$Length_ot,$GC_ot);
my ($ID_ot, $M0_ot, $M1_ot,$M2_ot, $M3_ot, $M4_ot,$M5_ot,$OT_ot);

open( PA, "$File1_ot" ) || die "can't open $File1_ot!";
open( PS, "$File2_ot" ) || die "can't open $File2_ot!";
open( PAIRS1, ">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/single_and_ot.txt" ) || die "can't open single_and_ot.txt!";

while ( <PA> ) {
    chomp();
($sgRID_ot,$Start_ot, $End_ot, $CRISPR_target_sequence_ot,$Length_ot,$GC_ot)=split/\t/, $_;

    seek PS, 0, 0;  #seek 设置文件的当前位置，http://blog.chinaunix.net/uid-20648367-id-3736970.html，http://blog.csdn.net/wzhwho/article/details/6181449

    while (  <PS> ) {
    chomp();
   ($ID_ot, $M0_ot, $M1_ot,$M2_ot, $M3_ot, $M4_ot,$M5_ot,$OT_ot)=split/\t/, $_;

if ($ID_ot eq $sgRID_ot) {

print PAIRS1 $sgRID_ot."\t".$Start_ot."\t".$End_ot."\t".$CRISPR_target_sequence_ot."\t".$Length_ot."\t".$GC_ot."\t".$ID_ot."\t".$M0_ot."\t".$M1_ot."\t".$M2_ot."\t".$M3_ot."\t".$M4_ot."\t".$M5_ot."\t".$OT_ot."\n";

}
}
}
close(PA);
close(PS);
close(PAIRS1);

#########
my $File1 = "$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/single_and_ot.txt";

my ($Risk_evaluation);

open( PA, "$File1" ) || die "can't open $File1!";
open( PAIRS1, ">$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/cpf1_$faname._sgRNAcas9_report.xls" ) || die "can't open cpf1_$faname._sgRNAcas9_report.xls!";

print PAIRS1 "sgRID\t"."Start\t"."End\t"."Protospacer_sequence+PAM(5'-3')\t"."Length(nt)\t"."GC%_of_Protospacer\t"."Protospacer+PAM(OT)\t"."0M(on-/off-)\t"."1M\t"."2M\t"."3M\t"."4M\t"."5M\t"."Total_No.of_OT\t"."Risk_evaluation\n";

while ( <PA> ) {
    chomp();
my ($sgRID_POT, $Start_POT, $End_POT, $CRISPR_target_sequence_POT,$Length_POT,$GC_POT,$ID_POT, $M0_POT, $M1_POT,$M2_POT, $M3_POT, $M4_POT,$M5_POT,$OT_POT);
($sgRID_POT, $Start_POT, $End_POT, $CRISPR_target_sequence_POT,$Length_POT,$GC_POT,$ID_POT, $M0_POT, $M1_POT,$M2_POT, $M3_POT, $M4_POT,$M5_POT,$OT_POT)=split/\t/, $_;


if ($M0_POT == 0 ) {

$Risk_evaluation = "Discard"; 

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t".""."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t".$Risk_evaluation."\n";

}elsif ($M0_POT >1 ) {

$Risk_evaluation = "Repeat_sites_or_bad?";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t".$Risk_evaluation."\n";

}elsif ($M1_POT ==0 and $M2_POT ==0 and $M3_POT <=10 ){   #可优化

$Risk_evaluation = "Best";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t".$Risk_evaluation."\n";

}elsif ($M1_POT ==0 and $M2_POT ==0 ){

$Risk_evaluation = "Low_risk";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t".$Risk_evaluation."\n";

}elsif ($M1_POT>=1 ) {

$Risk_evaluation = "High_risk";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t".$Risk_evaluation."\n";

}elsif ($M1_POT==1 and $M2_POT >=1 ){

$Risk_evaluation = "Moderate_risk";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t".$Risk_evaluation."\n";

}else{

$Risk_evaluation = "Unclassified";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t".$Risk_evaluation."\n";
}
}
close(PA);
close(PS);
close(PAIRS1);

#######subroutine to calculate GC% content
sub GC { 

    my $seq2 = shift @_; 

    $seq2 = uc($seq2);

    my $seq2DNA;
    $seq2DNA = $seq2;
    $seq2DNA =~ tr/U/T/;

    my $seq2length = length($seq2DNA);

    my $count = 0;
    for (my $i = 0; $i < $seq2length; $i++) {
			my $sub = substr($seq2DNA,$i,1);
			if ($sub =~ /G|C/i) {
	    	$count++;
			}
    }
    my $gc = sprintf("%.1f",$count * 100 /$seq2length);
    return $gc;
}


#######subroutine to calculate AT% content
sub AT { 

    my $seq2 = shift @_; 

    $seq2 = uc($seq2);

    my $seq2DNA;
    $seq2DNA = $seq2;
    $seq2DNA =~ tr/U/T/;

    my $seq2length = length($seq2DNA);

    my $count = 0;
    for (my $i = 0; $i < $seq2length; $i++) {
			my $sub = substr($seq2DNA,$i,1);
			if ($sub =~ /A|T/i) {
	    	$count++;
			}
    }
    my $gc = sprintf("%.1f",$count * 100 /$seq2length);
    return $gc;
}

close IntPut;
close OutPut;
#exit;


#delete Temporary files.
#unlink ("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/seqmap_output.txt")||die "Can't delete seqmap_output.txt file";
#unlink ("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/search_OT.txt")||die "Can't delete search_OT.txt file";
#unlink ("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/unique_pairs.fa")||die "Can't delete unique_pairs.fa file";
#unlink ("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_A.txt")||die "Can't delete CRISPR.targets_A.txt file";
#unlink ("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_S.txt")||die "Can't delete CRISPR.targets_S.txt file";
#unlink ("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_pairs.fa")||die "Can't delete CRISPR.targets_pairs.fa file";
#unlink ("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/CRISPR.targets_single.fa")||die "Can't delete CRISPR.targets_single.fa file";
#unlink ("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/report_sgRid.txt")||die "Can't delete report_sgRid.txt file";
#unlink ("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/report_protospacer_single.txt")||die "Can't delete report_protospacer_single.txt";
#unlink ("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/report_count_total_ot.txt")||die "Can't delete report_count_total_ot.txt";
#unlink ("$dir/cpf1.sgRNAcas9.report_$truncat.$Option.$faname/single_and_ot.txt")||die "Can't delete single_and_ot.txt";


######################################## Running Time #################################################
my $oo = time() -  $oldtime;
printf  "\nTotal time consumption is $oo second.\n";

printf LOG "# Total time consumption is $oo seconds." ."\n"; 
print "Your job is done, open cpf1.sgRNAcas9.report_$truncat.$Option directory to check result.". "\n";
print  LOG "# Your job is done, open cpf1.sgRNAcas9.report_$truncat.$Option directory to check result.". "\n";

print  LOG "################################# END ###########################################"."\n\n";
