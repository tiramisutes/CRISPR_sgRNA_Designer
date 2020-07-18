#!/bin/bash

# sysinfo_page - A script to filter BE (ABE and CBE) sgRNA from Cas9 sgRNAs.

# USAGE: ./sgRNAcas9toBE.sh ABE|CBE Ghir_A12G025340.fa.sgRNAcas9_report.xls ABE|CBE_sgRNA.xls

##### Constants
editor="$1"
input="$2"
output="$3"


if [ "${editor}" == "ABE" ]
then
	cat ${input} \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,substr($4,2,6);} else if($2%3==1) {print $0,substr($4,4,6);} else {print $0,substr($4,3,6);}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,substr($24,1,2),"#",substr($24,3,4);} else if($2%3==1) {print $0,substr($24,5,2),"#",substr($24,1,4);} else {print $0,substr($24,1,1),substr($24,6,1),substr($24,2,4);}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {gsub(/A/,"G",$27); print}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,$25$27;} else if($2%3==1) {print $0,$27$25;} else {print $0,$25$27$26;}}' \
	| awk -F"\t" 'BEGIN{OFS="\t";codon["ATA"]="Ile";codon["ATC"]="Ile";codon["ATT"]="Ile";codon["ATG"]="Met";codon["ACA"]="Thr";codon["ACC"]="Thr";codon["ACG"]="Thr";codon["ACT"]="Thr";codon["AAC"]="Asn";codon["AAT"]="Asn";codon["AAA"]="Lys";codon["AAG"]="Lys";codon["AGC"]="Ser";codon["AGT"]="Ser";codon["AGA"]="Arg";codon["AGG"]="Arg";codon["CTA"]="Leu";codon["CTC"]="Leu";codon["CTG"]="Leu";codon["CTT"]="Leu";codon["CCA"]="Pro";codon["CCC"]="Pro";codon["CCG"]="Pro";codon["CCT"]="Pro";codon["CAC"]="His";codon["CAT"]="His";codon["CAA"]="Gln";codon["CAG"]="Gln";codon["CGA"]="Arg";codon["CGC"]="Arg";codon["CGG"]="Arg";codon["CGT"]="Arg";codon["GTA"]="Val";codon["GTC"]="Val";codon["GTG"]="Val";codon["GTT"]="Val";codon["GCA"]="Ala";codon["GCC"]="Ala";codon["GCG"]="Ala";codon["GCT"]="Ala";codon["GAC"]="Asp";codon["GAT"]="Asp";codon["GAA"]="Glu";codon["GAG"]="Glu";codon["GGA"]="Gly";codon["GGC"]="Gly";codon["GGG"]="Gly";codon["GGT"]="Gly";codon["TCA"]="Ser";codon["TCC"]="Ser";codon["TCG"]="Ser";codon["TCT"]="Ser";codon["TTC"]="Phe";codon["TTT"]="Phe";codon["TTA"]="Leu";codon["TTG"]="Leu";codon["TAC"]="Tyr";codon["TAT"]="Tyr";codon["TAA"]="Stop";codon["TAG"]="Stop";codon["TGC"]="Cys";codon["TGT"]="Cys";codon["TGA"]="Stop";codon["TGG"]="Trp";} {print $0,codon[substr($24,1,3)]codon[substr($24,4,3)],codon[substr($28,1,3)]codon[substr($28,4,3)]}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,substr($4,2,9);} else if($2%3==1) {print $0,substr($4,4,6);} else {print $0,substr($4,3,9);}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,substr($31,1,2),substr($31,9,1),substr($31,3,6);} else if($2%3==1) {print $0,"#","#",substr($31,1,6);} else {print $0,substr($31,1,1),substr($31,8,2),substr($31,2,6);}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {gsub(/A/,"G",$34); print}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,$32$34$33;} else if($2%3==1) {print $0,$34;} else {print $0,$32$34$33;}}' \
	| awk -F"\t" 'BEGIN{OFS="\t";codon["ATA"]="Ile";codon["ATC"]="Ile";codon["ATT"]="Ile";codon["ATG"]="Met";codon["ACA"]="Thr";codon["ACC"]="Thr";codon["ACG"]="Thr";codon["ACT"]="Thr";codon["AAC"]="Asn";codon["AAT"]="Asn";codon["AAA"]="Lys";codon["AAG"]="Lys";codon["AGC"]="Ser";codon["AGT"]="Ser";codon["AGA"]="Arg";codon["AGG"]="Arg";codon["CTA"]="Leu";codon["CTC"]="Leu";codon["CTG"]="Leu";codon["CTT"]="Leu";codon["CCA"]="Pro";codon["CCC"]="Pro";codon["CCG"]="Pro";codon["CCT"]="Pro";codon["CAC"]="His";codon["CAT"]="His";codon["CAA"]="Gln";codon["CAG"]="Gln";codon["CGA"]="Arg";codon["CGC"]="Arg";codon["CGG"]="Arg";codon["CGT"]="Arg";codon["GTA"]="Val";codon["GTC"]="Val";codon["GTG"]="Val";codon["GTT"]="Val";codon["GCA"]="Ala";codon["GCC"]="Ala";codon["GCG"]="Ala";codon["GCT"]="Ala";codon["GAC"]="Asp";codon["GAT"]="Asp";codon["GAA"]="Glu";codon["GAG"]="Glu";codon["GGA"]="Gly";codon["GGC"]="Gly";codon["GGG"]="Gly";codon["GGT"]="Gly";codon["TCA"]="Ser";codon["TCC"]="Ser";codon["TCG"]="Ser";codon["TCT"]="Ser";codon["TTC"]="Phe";codon["TTT"]="Phe";codon["TTA"]="Leu";codon["TTG"]="Leu";codon["TAC"]="Tyr";codon["TAT"]="Tyr";codon["TAA"]="Stop";codon["TAG"]="Stop";codon["TGC"]="Cys";codon["TGT"]="Cys";codon["TGA"]="Stop";codon["TGG"]="Trp";} {print $0,codon[substr($31,1,3)]codon[substr($31,4,3)]codon[substr($31,7,3)],codon[substr($35,1,3)]codon[substr($35,4,3)]codon[substr($35,7,3)]}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($29==$30) {print $0,"4to7: No";} else {print $0,"Yes";}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($36==$37) {print $0,"4to9: No";} else {print $0,"Yes";}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($29 !~ /Arg|Lys|His/ && $30 ~ /Arg|Lys|His/ || $29 ~ /Arg|Lys|His/ && $30 !~ /Arg|Lys|His/) {print $0,"PH>7";} else if($29 !~ /Asp|Glu/ && $30 ~ /Asp|Glu/ || $29 ~ /Asp|Glu/ && $30 !~ /Asp|Glu/) {print $0,"PH<7";} else if($29 !~ /Stop/ && $30 ~ /Stop/ || $29 ~ /Stop/ && $30 !~ /Stop/) {print $0,"Find Stop";} else {print $0,"PH=7";}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if(substr($4,5,1) ~ /A/) {print $0,substr($4,5,1)"5th is A";} else {print $0,"No A in 5th";}}' \
	> ${output}
elif [ "${editor}" == "CBE" ]
then
	cat ${input} \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,substr($4,2,6);} else if($2%3==1) {print $0,substr($4,4,6);} else {print $0,substr($4,3,6);}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,substr($24,1,2),"#",substr($24,3,4);} else if($2%3==1) {print $0,substr($24,5,2),"#",substr($24,1,4);} else {print $0,substr($24,1,1),substr($24,6,1),substr($24,2,4);}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {gsub(/C/,"T",$27); print}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,$25$27;} else if($2%3==1) {print $0,$27$25;} else {print $0,$25$27$26;}}' \
	| awk -F"\t" 'BEGIN{OFS="\t";codon["ATA"]="Ile";codon["ATC"]="Ile";codon["ATT"]="Ile";codon["ATG"]="Met";codon["ACA"]="Thr";codon["ACC"]="Thr";codon["ACG"]="Thr";codon["ACT"]="Thr";codon["AAC"]="Asn";codon["AAT"]="Asn";codon["AAA"]="Lys";codon["AAG"]="Lys";codon["AGC"]="Ser";codon["AGT"]="Ser";codon["AGA"]="Arg";codon["AGG"]="Arg";codon["CTA"]="Leu";codon["CTC"]="Leu";codon["CTG"]="Leu";codon["CTT"]="Leu";codon["CCA"]="Pro";codon["CCC"]="Pro";codon["CCG"]="Pro";codon["CCT"]="Pro";codon["CAC"]="His";codon["CAT"]="His";codon["CAA"]="Gln";codon["CAG"]="Gln";codon["CGA"]="Arg";codon["CGC"]="Arg";codon["CGG"]="Arg";codon["CGT"]="Arg";codon["GTA"]="Val";codon["GTC"]="Val";codon["GTG"]="Val";codon["GTT"]="Val";codon["GCA"]="Ala";codon["GCC"]="Ala";codon["GCG"]="Ala";codon["GCT"]="Ala";codon["GAC"]="Asp";codon["GAT"]="Asp";codon["GAA"]="Glu";codon["GAG"]="Glu";codon["GGA"]="Gly";codon["GGC"]="Gly";codon["GGG"]="Gly";codon["GGT"]="Gly";codon["TCA"]="Ser";codon["TCC"]="Ser";codon["TCG"]="Ser";codon["TCT"]="Ser";codon["TTC"]="Phe";codon["TTT"]="Phe";codon["TTA"]="Leu";codon["TTG"]="Leu";codon["TAC"]="Tyr";codon["TAT"]="Tyr";codon["TAA"]="Stop";codon["TAG"]="Stop";codon["TGC"]="Cys";codon["TGT"]="Cys";codon["TGA"]="Stop";codon["TGG"]="Trp";} {print $0,codon[substr($24,1,3)]codon[substr($24,4,3)],codon[substr($28,1,3)]codon[substr($28,4,3)]}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,substr($4,2,9);} else if($2%3==1) {print $0,substr($4,4,6);} else {print $0,substr($4,3,9);}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,substr($31,1,2),substr($31,9,1),substr($31,3,6);} else if($2%3==1) {print $0,"#","#",substr($31,1,6);} else {print $0,substr($31,1,1),substr($31,8,2),substr($31,2,6);}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {gsub(/C/,"T",$34); print}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($2%3==0) {print $0,$32$34$33;} else if($2%3==1) {print $0,$34;} else {print $0,$32$34$33;}}' \
	| awk -F"\t" 'BEGIN{OFS="\t";codon["ATA"]="Ile";codon["ATC"]="Ile";codon["ATT"]="Ile";codon["ATG"]="Met";codon["ACA"]="Thr";codon["ACC"]="Thr";codon["ACG"]="Thr";codon["ACT"]="Thr";codon["AAC"]="Asn";codon["AAT"]="Asn";codon["AAA"]="Lys";codon["AAG"]="Lys";codon["AGC"]="Ser";codon["AGT"]="Ser";codon["AGA"]="Arg";codon["AGG"]="Arg";codon["CTA"]="Leu";codon["CTC"]="Leu";codon["CTG"]="Leu";codon["CTT"]="Leu";codon["CCA"]="Pro";codon["CCC"]="Pro";codon["CCG"]="Pro";codon["CCT"]="Pro";codon["CAC"]="His";codon["CAT"]="His";codon["CAA"]="Gln";codon["CAG"]="Gln";codon["CGA"]="Arg";codon["CGC"]="Arg";codon["CGG"]="Arg";codon["CGT"]="Arg";codon["GTA"]="Val";codon["GTC"]="Val";codon["GTG"]="Val";codon["GTT"]="Val";codon["GCA"]="Ala";codon["GCC"]="Ala";codon["GCG"]="Ala";codon["GCT"]="Ala";codon["GAC"]="Asp";codon["GAT"]="Asp";codon["GAA"]="Glu";codon["GAG"]="Glu";codon["GGA"]="Gly";codon["GGC"]="Gly";codon["GGG"]="Gly";codon["GGT"]="Gly";codon["TCA"]="Ser";codon["TCC"]="Ser";codon["TCG"]="Ser";codon["TCT"]="Ser";codon["TTC"]="Phe";codon["TTT"]="Phe";codon["TTA"]="Leu";codon["TTG"]="Leu";codon["TAC"]="Tyr";codon["TAT"]="Tyr";codon["TAA"]="Stop";codon["TAG"]="Stop";codon["TGC"]="Cys";codon["TGT"]="Cys";codon["TGA"]="Stop";codon["TGG"]="Trp";} {print $0,codon[substr($31,1,3)]codon[substr($31,4,3)]codon[substr($31,7,3)],codon[substr($35,1,3)]codon[substr($35,4,3)]codon[substr($35,7,3)]}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($29==$30) {print $0,"4to7: No";} else {print $0,"Yes";}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($36==$37) {print $0,"4to9: No";} else {print $0,"Yes";}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if($29 !~ /Arg|Lys|His/ && $30 ~ /Arg|Lys|His/ || $29 ~ /Arg|Lys|His/ && $30 !~ /Arg|Lys|His/) {print $0,"PH>7";} else if($29 !~ /Asp|Glu/ && $30 ~ /Asp|Glu/ || $29 ~ /Asp|Glu/ && $30 !~ /Asp|Glu/) {print $0,"PH<7";} else if($29 !~ /Stop/ && $30 ~ /Stop/ || $29 ~ /Stop/ && $30 !~ /Stop/) {print $0,"Find Stop";} else {print $0,"PH=7";}}' \
	| awk -F"\t" 'BEGIN{OFS="\t"} {if(substr($4,5,1) ~ /C/) {print $0,substr($4,5,1)"5th is C";} else {print $0,"No C in 5th";}}' \
	> ${output}
else
	echo -e "##################################################################\nPleast select ABE or CBE\n##################################################################\n"
fi
