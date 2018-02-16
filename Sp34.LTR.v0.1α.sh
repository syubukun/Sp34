#!/bin/bash

cat $1|sed -e 's/name=AP/name=PR/g'|sed -e 's/name=RVP/name=PR/g'|sed -e 's/name=RNaseH/name=RH/g'|sed -e 's/name=rve/name=RH/g'|sed -e 's/name=RVT/name=RT/g' >tmp1.txt

cat tmp1.txt|sed '/^#.*/d'|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="repeat_region")print "##"$0;else print "#"$0}'|tr -d '\n'|sed -e 's/##/\n/g' >tmp2.txt

mkdir sep1
cat tmp2.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0!~"protein_match")print $0}'|tr '#' '\n'|sed '/^$/d' >sep1/LTRharvest_solo.txt
cat tmp2.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"protein_match")print $0}' >LTRharvest_par_and_full_tmp1.txt
cat LTRharvest_par_and_full_tmp1.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"name=RT"&&$0~"name=PR"&&$0~"name=INT"&&$0~"name=RH")print $0}'|sed -e 's/^.*Parent=/Parent=/g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $0"#"}' >full_list.txt

cat full_list.txt|awk 'BEGIN{ORS="\\\|"}{print $0}'|sed -e 's/\\|$//g' >sep1/full.grep.txt
cat LTRharvest_par_and_full_tmp1.txt|grep -f sep1/full.grep.txt|tr '#' '\n' >sep1/LTRharvest_full.txt
cat LTRharvest_par_and_full_tmp1.txt|grep -v -f sep1/full.grep.txt|tr '#' '\n' >sep1/LTRharvest_par.txt

mkdir sep2
cat sep1/LTRharvest_solo.txt|awk 'BEGIN{FS="\t";OFS="\t"}{print $0";status=solo"}' >sep2/LTRharvest_solo2.txt
cat sep1/LTRharvest_par.txt|awk 'BEGIN{FS="\t";OFS="\t"}{print $0";status=par"}' >sep2/LTRharvest_par2.txt
cat sep1/LTRharvest_full.txt|awk 'BEGIN{FS="\t";OFS="\t"}{print $0";status=full"}' >sep2/LTRharvest_full2.txt

mkdir sep3
cat sep2/LTRharvest_full2.txt sep2/LTRharvest_par2.txt sep2/LTRharvest_solo2.txt|awk 'BEGIN{FS="\t"}{if($1=="Sp34_Chr1"&&($3=="LTR_retrotransposon"||$3=="long_terminal_repeat"||$3=="protein_match"||$3=="RR_tract"))print}'|sort -k4,4 -n >sep3/LTRharvest_Chr1.txt
cat sep2/LTRharvest_full2.txt sep2/LTRharvest_par2.txt sep2/LTRharvest_solo2.txt|awk 'BEGIN{FS="\t"}{if($1=="Sp34_Chr2"&&($3=="LTR_retrotransposon"||$3=="long_terminal_repeat"||$3=="protein_match"||$3=="RR_tract"))print}'|sort -k4,4 -n >sep3/LTRharvest_Chr2.txt
cat sep2/LTRharvest_full2.txt sep2/LTRharvest_par2.txt sep2/LTRharvest_solo2.txt|awk 'BEGIN{FS="\t"}{if($1=="Sp34_Chr3"&&($3=="LTR_retrotransposon"||$3=="long_terminal_repeat"||$3=="protein_match"||$3=="RR_tract"))print}'|sort -k4,4 -n >sep3/LTRharvest_Chr3.txt
cat sep2/LTRharvest_full2.txt sep2/LTRharvest_par2.txt sep2/LTRharvest_solo2.txt|awk 'BEGIN{FS="\t"}{if($1=="Sp34_Chr4"&&($3=="LTR_retrotransposon"||$3=="long_terminal_repeat"||$3=="protein_match"||$3=="RR_tract"))print}'|sort -k4,4 -n >sep3/LTRharvest_Chr4.txt
cat sep2/LTRharvest_full2.txt sep2/LTRharvest_par2.txt sep2/LTRharvest_solo2.txt|awk 'BEGIN{FS="\t"}{if($1=="Sp34_Chr5"&&($3=="LTR_retrotransposon"||$3=="long_terminal_repeat"||$3=="protein_match"||$3=="RR_tract"))print}'|sort -k4,4 -n >sep3/LTRharvest_Chr5.txt
cat sep2/LTRharvest_full2.txt sep2/LTRharvest_par2.txt sep2/LTRharvest_solo2.txt|awk 'BEGIN{FS="\t"}{if($1=="Sp34_ChrX"&&($3=="LTR_retrotransposon"||$3=="long_terminal_repeat"||$3=="protein_match"||$3=="RR_tract"))print}'|sort -k4,4 -n >sep3/LTRharvest_ChrX.txt

mkdir MGEScan_compare
cat $2|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"Sp34")print "MGEScan_"$1,$2,$5}'|grep 'Chr1'|sort -k2,2 -n|sed -e 's/Sp34_//g' >MGEScan_compare/MGEScan_Chr1.table.txt
cat $2|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"Sp34")print "MGEScan_"$1,$2,$5}'|grep 'Chr2'|sort -k2,2 -n|sed -e 's/Sp34_//g' >MGEScan_compare/MGEScan_Chr2.table.txt
cat $2|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"Sp34")print "MGEScan_"$1,$2,$5}'|grep 'Chr3'|sort -k2,2 -n|sed -e 's/Sp34_//g' >MGEScan_compare/MGEScan_Chr3.table.txt
cat $2|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"Sp34")print "MGEScan_"$1,$2,$5}'|grep 'Chr4'|sort -k2,2 -n|sed -e 's/Sp34_//g' >MGEScan_compare/MGEScan_Chr4.table.txt
cat $2|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"Sp34")print "MGEScan_"$1,$2,$5}'|grep 'Chr5'|sort -k2,2 -n|sed -e 's/Sp34_//g' >MGEScan_compare/MGEScan_Chr5.table.txt
cat $2|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"Sp34")print "MGEScan_"$1,$2,$5}'|grep 'ChrX'|sort -k2,2 -n|sed -e 's/Sp34_//g' >MGEScan_compare/MGEScan_ChrX.table.txt
cat sep3/LTRharvest_Chr1.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="LTR_retrotransposon")print $4,$5,$9}'|sed -e 's/;Parent=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $3,$1,$2}'|sed -e 's/ID=LTR_retrotransposon/LTRharvest_Chr1_/g' >MGEScan_compare/LTRharvest_Chr1.table.txt
cat sep3/LTRharvest_Chr2.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="LTR_retrotransposon")print $4,$5,$9}'|sed -e 's/;Parent=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $3,$1,$2}'|sed -e 's/ID=LTR_retrotransposon/LTRharvest_Chr2_/g' >MGEScan_compare/LTRharvest_Chr2.table.txt
cat sep3/LTRharvest_Chr3.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="LTR_retrotransposon")print $4,$5,$9}'|sed -e 's/;Parent=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $3,$1,$2}'|sed -e 's/ID=LTR_retrotransposon/LTRharvest_Chr3_/g' >MGEScan_compare/LTRharvest_Chr3.table.txt
cat sep3/LTRharvest_Chr4.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="LTR_retrotransposon")print $4,$5,$9}'|sed -e 's/;Parent=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $3,$1,$2}'|sed -e 's/ID=LTR_retrotransposon/LTRharvest_Chr4_/g' >MGEScan_compare/LTRharvest_Chr4.table.txt
cat sep3/LTRharvest_Chr5.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="LTR_retrotransposon")print $4,$5,$9}'|sed -e 's/;Parent=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $3,$1,$2}'|sed -e 's/ID=LTR_retrotransposon/LTRharvest_Chr5_/g' >MGEScan_compare/LTRharvest_Chr5.table.txt
cat sep3/LTRharvest_ChrX.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="LTR_retrotransposon")print $4,$5,$9}'|sed -e 's/;Parent=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $3,$1,$2}'|sed -e 's/ID=LTR_retrotransposon/LTRharvest_ChrX_/g' >MGEScan_compare/LTRharvest_ChrX.table.txt

cd MGEScan_compare
cat MGEScan_Chr1.table.txt LTRharvest_Chr1.table.txt|sort -k 2,2 -n >Chr1_marged.table.txt
cat MGEScan_Chr2.table.txt LTRharvest_Chr2.table.txt|sort -k 2,2 -n >Chr2_marged.table.txt
cat MGEScan_Chr3.table.txt LTRharvest_Chr3.table.txt|sort -k 2,2 -n >Chr3_marged.table.txt
cat MGEScan_Chr4.table.txt LTRharvest_Chr4.table.txt|sort -k 2,2 -n >Chr4_marged.table.txt
cat MGEScan_Chr5.table.txt LTRharvest_Chr5.table.txt|sort -k 2,2 -n >Chr5_marged.table.txt
cat MGEScan_ChrX.table.txt LTRharvest_ChrX.table.txt|sort -k 2,2 -n >ChrX_marged.table.txt

cat Chr1_marged.table.txt|grep -A1 -B1 'MGE'|grep -v '-'|awk 'BEGIN{FS="\t";OFS="\t"}{a[NR]=$1;b[NR]=$2;c[NR]=$3}{if($2-c[NR-1]<0)print a[NR-1],b[NR-1],c[NR-1],$0}' >Chr1.tmp.txt
cat Chr2_marged.table.txt|grep -A1 -B1 'MGE'|grep -v '-'|awk 'BEGIN{FS="\t";OFS="\t"}{a[NR]=$1;b[NR]=$2;c[NR]=$3}{if($2-c[NR-1]<0)print a[NR-1],b[NR-1],c[NR-1],$0}' >Chr2.tmp.txt
cat Chr3_marged.table.txt|grep -A1 -B1 'MGE'|grep -v '-'|awk 'BEGIN{FS="\t";OFS="\t"}{a[NR]=$1;b[NR]=$2;c[NR]=$3}{if($2-c[NR-1]<0)print a[NR-1],b[NR-1],c[NR-1],$0}' >Chr3.tmp.txt
cat Chr4_marged.table.txt|grep -A1 -B1 'MGE'|grep -v '-'|awk 'BEGIN{FS="\t";OFS="\t"}{a[NR]=$1;b[NR]=$2;c[NR]=$3}{if($2-c[NR-1]<0)print a[NR-1],b[NR-1],c[NR-1],$0}' >Chr4.tmp.txt
cat Chr5_marged.table.txt|grep -A1 -B1 'MGE'|grep -v '-'|awk 'BEGIN{FS="\t";OFS="\t"}{a[NR]=$1;b[NR]=$2;c[NR]=$3}{if($2-c[NR-1]<0)print a[NR-1],b[NR-1],c[NR-1],$0}' >Chr5.tmp.txt
cat ChrX_marged.table.txt|grep -A1 -B1 'MGE'|grep -v '-'|awk 'BEGIN{FS="\t";OFS="\t"}{a[NR]=$1;b[NR]=$2;c[NR]=$3}{if($2-c[NR-1]<0)print a[NR-1],b[NR-1],c[NR-1],$0}' >ChrX.tmp.txt

cat Chr1.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"MGEScan")print $4;else if($4~"MGEScan")print $1}' >Remove_from_LTRharvest.Chr1.list.txt
cat Chr2.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"MGEScan")print $4;else if($4~"MGEScan")print $1}' >Remove_from_LTRharvest.Chr2.list.txt
cat Chr3.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"MGEScan")print $4;else if($4~"MGEScan")print $1}' >Remove_from_LTRharvest.Chr3.list.txt
cat Chr4.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"MGEScan")print $4;else if($4~"MGEScan")print $1}' >Remove_from_LTRharvest.Chr4.list.txt
cat Chr5.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"MGEScan")print $4;else if($4~"MGEScan")print $1}' >Remove_from_LTRharvest.Chr5.list.txt
cat ChrX.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~"MGEScan")print $4;else if($4~"MGEScan")print $1}' >Remove_from_LTRharvest.ChrX.list.txt

cat Remove_from_LTRharvest.Chr1.list.txt|sed -e 's/LTRharvest_Chr1_/LTR_retrotransposon/g'|sed -e 's/$/;/g'|awk 'BEGIN{ORS="\\\|"}{print}'|sed -e 's/\\|$//g' >Remove_from_LTRharvest.Chr1.grep.script.txt
cat Remove_from_LTRharvest.Chr2.list.txt|sed -e 's/LTRharvest_Chr2_/LTR_retrotransposon/g'|sed -e 's/$/;/g'|awk 'BEGIN{ORS="\\\|"}{print}'|sed -e 's/\\|$//g' >Remove_from_LTRharvest.Chr2.grep.script.txt
cat Remove_from_LTRharvest.Chr3.list.txt|sed -e 's/LTRharvest_Chr3_/LTR_retrotransposon/g'|sed -e 's/$/;/g'|awk 'BEGIN{ORS="\\\|"}{print}'|sed -e 's/\\|$//g' >Remove_from_LTRharvest.Chr3.grep.script.txt
cat Remove_from_LTRharvest.Chr4.list.txt|sed -e 's/LTRharvest_Chr4_/LTR_retrotransposon/g'|sed -e 's/$/;/g'|awk 'BEGIN{ORS="\\\|"}{print}'|sed -e 's/\\|$//g' >Remove_from_LTRharvest.Chr4.grep.script.txt
cat Remove_from_LTRharvest.Chr5.list.txt|sed -e 's/LTRharvest_Chr5_/LTR_retrotransposon/g'|sed -e 's/$/;/g'|awk 'BEGIN{ORS="\\\|"}{print}'|sed -e 's/\\|$//g' >Remove_from_LTRharvest.Chr5.grep.script.txt
cat Remove_from_LTRharvest.ChrX.list.txt|sed -e 's/LTRharvest_ChrX_/LTR_retrotransposon/g'|sed -e 's/$/;/g'|awk 'BEGIN{ORS="\\\|"}{print}'|sed -e 's/\\|$//g' >Remove_from_LTRharvest.ChrX.grep.script.txt

cd ..
mkdir sep4
cd sep4
cat ../sep3/LTRharvest_Chr1.txt|grep -v -f ../MGEScan_compare/Remove_from_LTRharvest.Chr1.grep.script.txt >LTRharvest_Chr1.gff.txt
cat ../sep3/LTRharvest_Chr2.txt|grep -v -f ../MGEScan_compare/Remove_from_LTRharvest.Chr2.grep.script.txt >LTRharvest_Chr2.gff.txt
cat ../sep3/LTRharvest_Chr3.txt|grep -v -f ../MGEScan_compare/Remove_from_LTRharvest.Chr3.grep.script.txt >LTRharvest_Chr3.gff.txt
cat ../sep3/LTRharvest_Chr4.txt|grep -v -f ../MGEScan_compare/Remove_from_LTRharvest.Chr4.grep.script.txt >LTRharvest_Chr4.gff.txt
cat ../sep3/LTRharvest_Chr5.txt|grep -v -f ../MGEScan_compare/Remove_from_LTRharvest.Chr5.grep.script.txt >LTRharvest_Chr5.gff.txt
cat ../sep3/LTRharvest_ChrX.txt|grep -v -f ../MGEScan_compare/Remove_from_LTRharvest.ChrX.grep.script.txt >LTRharvest_ChrX.gff.txt

cd ..
mkdir sep5
cd sep5
cat ../sep4/LTRharvest_Chr1.gff.txt|sed -e 's/$/;/g' >LTRharvest.Chr1.gff.txt
cat ../sep4/LTRharvest_Chr2.gff.txt|sed -e 's/$/;/g' >LTRharvest.Chr2.gff.txt
cat ../sep4/LTRharvest_Chr3.gff.txt|sed -e 's/$/;/g' >LTRharvest.Chr3.gff.txt
cat ../sep4/LTRharvest_Chr4.gff.txt|sed -e 's/$/;/g' >LTRharvest.Chr4.gff.txt
cat ../sep4/LTRharvest_Chr5.gff.txt|sed -e 's/$/;/g' >LTRharvest.Chr5.gff.txt
cat ../sep4/LTRharvest_ChrX.gff.txt|sed -e 's/$/;/g' >LTRharvest.ChrX.gff.txt

cd ..
cat OR/Sp34.genome.v7.7_ltrdigest_output.renamed.gff.txt|head -13 >header.gff.txt

mkdir MGEtoLTR
cd MGEtoLTR
cat ../OR/ltr.renamed.out.txt|grep 'Chr1'|sort -n -k2,2|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$3,$4,$5,$1}'|sed -e 's/_[0-9]*$//g'|nl|sed -e 's/^  *//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $6,"LTRharvest","repeat_region",$2,$5,".","?",".","ID=repeat_region"$1"\n"$6,"LTRharvest","LTR_retrotransposon",$2,$5,".","?",".","ID=LTR_retrotransposon"$1";Parent=repeat_region"$1"\n"$6,"LTRharvest","long_terminal_repeat",$2,$3,".","?",".","Parent=LTR_retrotransposon"$1"\n"$6,"LTRharvest","long_terminal_repeat",$4,$5,".","?",".","Parent=LTR_retrotransposon"$1}' >MGEtoLTR.Chr1.tmp1.gff.txt
cat ../OR/ltr.renamed.out.txt|grep 'Chr2'|sort -n -k2,2|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$3,$4,$5,$1}'|sed -e 's/_[0-9]*$//g'|nl|sed -e 's/^  *//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $6,"LTRharvest","repeat_region",$2,$5,".","?",".","ID=repeat_region"$1"\n"$6,"LTRharvest","LTR_retrotransposon",$2,$5,".","?",".","ID=LTR_retrotransposon"$1";Parent=repeat_region"$1"\n"$6,"LTRharvest","long_terminal_repeat",$2,$3,".","?",".","Parent=LTR_retrotransposon"$1"\n"$6,"LTRharvest","long_terminal_repeat",$4,$5,".","?",".","Parent=LTR_retrotransposon"$1}' >MGEtoLTR.Chr2.tmp1.gff.txt
cat ../OR/ltr.renamed.out.txt|grep 'Chr3'|sort -n -k2,2|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$3,$4,$5,$1}'|sed -e 's/_[0-9]*$//g'|nl|sed -e 's/^  *//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $6,"LTRharvest","repeat_region",$2,$5,".","?",".","ID=repeat_region"$1"\n"$6,"LTRharvest","LTR_retrotransposon",$2,$5,".","?",".","ID=LTR_retrotransposon"$1";Parent=repeat_region"$1"\n"$6,"LTRharvest","long_terminal_repeat",$2,$3,".","?",".","Parent=LTR_retrotransposon"$1"\n"$6,"LTRharvest","long_terminal_repeat",$4,$5,".","?",".","Parent=LTR_retrotransposon"$1}' >MGEtoLTR.Chr3.tmp1.gff.txt
cat ../OR/ltr.renamed.out.txt|grep 'Chr4'|sort -n -k2,2|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$3,$4,$5,$1}'|sed -e 's/_[0-9]*$//g'|nl|sed -e 's/^  *//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $6,"LTRharvest","repeat_region",$2,$5,".","?",".","ID=repeat_region"$1"\n"$6,"LTRharvest","LTR_retrotransposon",$2,$5,".","?",".","ID=LTR_retrotransposon"$1";Parent=repeat_region"$1"\n"$6,"LTRharvest","long_terminal_repeat",$2,$3,".","?",".","Parent=LTR_retrotransposon"$1"\n"$6,"LTRharvest","long_terminal_repeat",$4,$5,".","?",".","Parent=LTR_retrotransposon"$1}' >MGEtoLTR.Chr4.tmp1.gff.txt
cat ../OR/ltr.renamed.out.txt|grep 'Chr5'|sort -n -k2,2|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$3,$4,$5,$1}'|sed -e 's/_[0-9]*$//g'|nl|sed -e 's/^  *//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $6,"LTRharvest","repeat_region",$2,$5,".","?",".","ID=repeat_region"$1"\n"$6,"LTRharvest","LTR_retrotransposon",$2,$5,".","?",".","ID=LTR_retrotransposon"$1";Parent=repeat_region"$1"\n"$6,"LTRharvest","long_terminal_repeat",$2,$3,".","?",".","Parent=LTR_retrotransposon"$1"\n"$6,"LTRharvest","long_terminal_repeat",$4,$5,".","?",".","Parent=LTR_retrotransposon"$1}' >MGEtoLTR.Chr5.tmp1.gff.txt
cat ../OR/ltr.renamed.out.txt|grep 'ChrX'|sort -n -k2,2|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$3,$4,$5,$1}'|sed -e 's/_[0-9]*$//g'|nl|sed -e 's/^  *//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $6,"LTRharvest","repeat_region",$2,$5,".","?",".","ID=repeat_region"$1"\n"$6,"LTRharvest","LTR_retrotransposon",$2,$5,".","?",".","ID=LTR_retrotransposon"$1";Parent=repeat_region"$1"\n"$6,"LTRharvest","long_terminal_repeat",$2,$3,".","?",".","Parent=LTR_retrotransposon"$1"\n"$6,"LTRharvest","long_terminal_repeat",$4,$5,".","?",".","Parent=LTR_retrotransposon"$1}' >MGEtoLTR.ChrX.tmp1.gff.txt

cat ../header.gff.txt MGEtoLTR.Chr1.tmp1.gff.txt >MGEtoLTR.Chr1.tmp2.gff.txt
cat ../header.gff.txt MGEtoLTR.Chr2.tmp1.gff.txt >MGEtoLTR.Chr2.tmp2.gff.txt
cat ../header.gff.txt MGEtoLTR.Chr3.tmp1.gff.txt >MGEtoLTR.Chr3.tmp2.gff.txt
cat ../header.gff.txt MGEtoLTR.Chr4.tmp1.gff.txt >MGEtoLTR.Chr4.tmp2.gff.txt
cat ../header.gff.txt MGEtoLTR.Chr5.tmp1.gff.txt >MGEtoLTR.Chr5.tmp2.gff.txt
cat ../header.gff.txt MGEtoLTR.ChrX.tmp1.gff.txt >MGEtoLTR.ChrX.tmp2.gff.txt

cd ..
mkdir GT
cd GT
cp ../OR/Sp34.genome.v7.7.fa .
gt suffixerator -tis -des -dna -ssp -db Sp34.genome.v7.7.fa -indexname Sp34.genome.v7.7.fa
gt extractseq Sp34.genome.v7.7.fa >Sp34.genome.v7.7.log.txt
gt ltrdigest -hmms  ../OR/hmm/*.hmm -outfileprefix MGEScan_Chr1_gt -encseq Sp34.genome.v7.7.fa -matchdescstart <../MGEtoLTR/MGEtoLTR.Chr1.tmp2.gff.txt >MGEtoLTR.Chr1.gff.txt
gt ltrdigest -hmms  ../OR/hmm/*.hmm -outfileprefix MGEScan_Chr2_gt -encseq Sp34.genome.v7.7.fa -matchdescstart <../MGEtoLTR/MGEtoLTR.Chr2.tmp2.gff.txt >MGEtoLTR.Chr2.gff.txt
gt ltrdigest -hmms  ../OR/hmm/*.hmm -outfileprefix MGEScan_Chr3_gt -encseq Sp34.genome.v7.7.fa -matchdescstart <../MGEtoLTR/MGEtoLTR.Chr3.tmp2.gff.txt >MGEtoLTR.Chr3.gff.txt
gt ltrdigest -hmms  ../OR/hmm/*.hmm -outfileprefix MGEScan_Chr4_gt -encseq Sp34.genome.v7.7.fa -matchdescstart <../MGEtoLTR/MGEtoLTR.Chr4.tmp2.gff.txt >MGEtoLTR.Chr4.gff.txt
gt ltrdigest -hmms  ../OR/hmm/*.hmm -outfileprefix MGEScan_Chr5_gt -encseq Sp34.genome.v7.7.fa -matchdescstart <../MGEtoLTR/MGEtoLTR.Chr5.tmp2.gff.txt >MGEtoLTR.Chr5.gff.txt
gt ltrdigest -hmms  ../OR/hmm/*.hmm -outfileprefix MGEScan_ChrX_gt -encseq Sp34.genome.v7.7.fa -matchdescstart <../MGEtoLTR/MGEtoLTR.ChrX.tmp2.gff.txt >MGEtoLTR.ChrX.gff.txt

cat MGEtoLTR.Chr1.gff.txt|grep -v '#'|sed -e 's/LTRharvest/MGEScan_LTR/g'|sed -e 's/LTRdigest/MGEScan_LTR/g'|sed -e 's/name=AP/name=PR/g'|sed -e 's/name=RVP/name=PR/g'|sed -e 's/name=RNaseH/name=RH/g'|sed -e 's/name=rve/name=RH/g'|sed -e 's/name=RVT/name=RT/g' >MGEtoLTR2.Chr1.gff.txt
cat MGEtoLTR.Chr2.gff.txt|grep -v '#'|sed -e 's/LTRharvest/MGEScan_LTR/g'|sed -e 's/LTRdigest/MGEScan_LTR/g'|sed -e 's/name=AP/name=PR/g'|sed -e 's/name=RVP/name=PR/g'|sed -e 's/name=RNaseH/name=RH/g'|sed -e 's/name=rve/name=RH/g'|sed -e 's/name=RVT/name=RT/g' >MGEtoLTR2.Chr2.gff.txt
cat MGEtoLTR.Chr3.gff.txt|grep -v '#'|sed -e 's/LTRharvest/MGEScan_LTR/g'|sed -e 's/LTRdigest/MGEScan_LTR/g'|sed -e 's/name=AP/name=PR/g'|sed -e 's/name=RVP/name=PR/g'|sed -e 's/name=RNaseH/name=RH/g'|sed -e 's/name=rve/name=RH/g'|sed -e 's/name=RVT/name=RT/g' >MGEtoLTR2.Chr3.gff.txt
cat MGEtoLTR.Chr4.gff.txt|grep -v '#'|sed -e 's/LTRharvest/MGEScan_LTR/g'|sed -e 's/LTRdigest/MGEScan_LTR/g'|sed -e 's/name=AP/name=PR/g'|sed -e 's/name=RVP/name=PR/g'|sed -e 's/name=RNaseH/name=RH/g'|sed -e 's/name=rve/name=RH/g'|sed -e 's/name=RVT/name=RT/g' >MGEtoLTR2.Chr4.gff.txt
cat MGEtoLTR.Chr5.gff.txt|grep -v '#'|sed -e 's/LTRharvest/MGEScan_LTR/g'|sed -e 's/LTRdigest/MGEScan_LTR/g'|sed -e 's/name=AP/name=PR/g'|sed -e 's/name=RVP/name=PR/g'|sed -e 's/name=RNaseH/name=RH/g'|sed -e 's/name=rve/name=RH/g'|sed -e 's/name=RVT/name=RT/g' >MGEtoLTR2.Chr5.gff.txt
cat MGEtoLTR.ChrX.gff.txt|grep -v '#'|sed -e 's/LTRharvest/MGEScan_LTR/g'|sed -e 's/LTRdigest/MGEScan_LTR/g'|sed -e 's/name=AP/name=PR/g'|sed -e 's/name=RVP/name=PR/g'|sed -e 's/name=RNaseH/name=RH/g'|sed -e 's/name=RVT/name=RT/g' >MGEtoLTR2.ChrX.gff.txt

cd ..
mkdir MGE
cd MGE
cat ../GT/MGEtoLTR2.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="repeat_region")print "##"$0;else print "#"$0}'|tr -d '\n'|sed -e 's/##/\n/g' >MGEtoLTR.Chr1.oneline.txt
cat ../GT/MGEtoLTR2.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="repeat_region")print "##"$0;else print "#"$0}'|tr -d '\n'|sed -e 's/##/\n/g' >MGEtoLTR.Chr2.oneline.txt
cat ../GT/MGEtoLTR2.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="repeat_region")print "##"$0;else print "#"$0}'|tr -d '\n'|sed -e 's/##/\n/g' >MGEtoLTR.Chr3.oneline.txt
cat ../GT/MGEtoLTR2.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="repeat_region")print "##"$0;else print "#"$0}'|tr -d '\n'|sed -e 's/##/\n/g' >MGEtoLTR.Chr4.oneline.txt
cat ../GT/MGEtoLTR2.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="repeat_region")print "##"$0;else print "#"$0}'|tr -d '\n'|sed -e 's/##/\n/g' >MGEtoLTR.Chr5.oneline.txt
cat ../GT/MGEtoLTR2.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="repeat_region")print "##"$0;else print "#"$0}'|tr -d '\n'|sed -e 's/##/\n/g' >MGEtoLTR.ChrX.oneline.txt

cat MGEtoLTR.Chr1.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0!~"protein_match")print $0}'|tr '#' '\n'|sed '/^$/d' >MGEtoLTR.Chr1.solo.gff.txt
cat MGEtoLTR.Chr2.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0!~"protein_match")print $0}'|tr '#' '\n'|sed '/^$/d' >MGEtoLTR.Chr2.solo.gff.txt
cat MGEtoLTR.Chr3.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0!~"protein_match")print $0}'|tr '#' '\n'|sed '/^$/d' >MGEtoLTR.Chr3.solo.gff.txt
cat MGEtoLTR.Chr4.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0!~"protein_match")print $0}'|tr '#' '\n'|sed '/^$/d' >MGEtoLTR.Chr4.solo.gff.txt
cat MGEtoLTR.Chr5.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0!~"protein_match")print $0}'|tr '#' '\n'|sed '/^$/d' >MGEtoLTR.Chr5.solo.gff.txt
cat MGEtoLTR.ChrX.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0!~"protein_match")print $0}'|tr '#' '\n'|sed '/^$/d' >MGEtoLTR.ChrX.solo.gff.txt

cat MGEtoLTR.Chr1.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"protein_match")print $0}' >MGEtoLTR.Chr1.par_and_full.tmp.txt
cat MGEtoLTR.Chr2.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"protein_match")print $0}' >MGEtoLTR.Chr2.par_and_full.tmp.txt
cat MGEtoLTR.Chr3.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"protein_match")print $0}' >MGEtoLTR.Chr3.par_and_full.tmp.txt
cat MGEtoLTR.Chr4.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"protein_match")print $0}' >MGEtoLTR.Chr4.par_and_full.tmp.txt
cat MGEtoLTR.Chr5.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"protein_match")print $0}' >MGEtoLTR.Chr5.par_and_full.tmp.txt
cat MGEtoLTR.ChrX.oneline.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"protein_match")print $0}' >MGEtoLTR.ChrX.par_and_full.tmp.txt

cat MGEtoLTR.Chr1.par_and_full.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"name=RT"&&$0~"name=PR"&&$0~"name=INT"&&$0~"name=RH")print $0}'|sed -e 's/^.*Parent=/Parent=/g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $0"#"}' >MGEtoLTR.Chr1.list.txt
cat MGEtoLTR.Chr2.par_and_full.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"name=RT"&&$0~"name=PR"&&$0~"name=INT"&&$0~"name=RH")print $0}'|sed -e 's/^.*Parent=/Parent=/g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $0"#"}' >MGEtoLTR.Chr2.list.txt
cat MGEtoLTR.Chr3.par_and_full.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"name=RT"&&$0~"name=PR"&&$0~"name=INT"&&$0~"name=RH")print $0}'|sed -e 's/^.*Parent=/Parent=/g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $0"#"}' >MGEtoLTR.Chr3.list.txt
cat MGEtoLTR.Chr4.par_and_full.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"name=RT"&&$0~"name=PR"&&$0~"name=INT"&&$0~"name=RH")print $0}'|sed -e 's/^.*Parent=/Parent=/g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $0"#"}' >MGEtoLTR.Chr4.list.txt
cat MGEtoLTR.Chr5.par_and_full.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"name=RT"&&$0~"name=PR"&&$0~"name=INT"&&$0~"name=RH")print $0}'|sed -e 's/^.*Parent=/Parent=/g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $0"#"}' >MGEtoLTR.Chr5.list.txt
cat MGEtoLTR.ChrX.par_and_full.tmp.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"name=RT"&&$0~"name=PR"&&$0~"name=INT"&&$0~"name=RH")print $0}'|sed -e 's/^.*Parent=/Parent=/g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $0"#"}' >MGEtoLTR.ChrX.list.txt

cat MGEtoLTR.Chr1.list.txt|awk 'BEGIN{ORS="\\\|"}{print $0}'|sed -e 's/\\|$//g' >MGEtoLTR.Chr1.grep.script.txt
cat MGEtoLTR.Chr2.list.txt|awk 'BEGIN{ORS="\\\|"}{print $0}'|sed -e 's/\\|$//g' >MGEtoLTR.Chr2.grep.script.txt
cat MGEtoLTR.Chr3.list.txt|awk 'BEGIN{ORS="\\\|"}{print $0}'|sed -e 's/\\|$//g' >MGEtoLTR.Chr3.grep.script.txt
cat MGEtoLTR.Chr4.list.txt|awk 'BEGIN{ORS="\\\|"}{print $0}'|sed -e 's/\\|$//g' >MGEtoLTR.Chr4.grep.script.txt
cat MGEtoLTR.Chr5.list.txt|awk 'BEGIN{ORS="\\\|"}{print $0}'|sed -e 's/\\|$//g' >MGEtoLTR.Chr5.grep.script.txt
cat MGEtoLTR.ChrX.list.txt|awk 'BEGIN{ORS="\\\|"}{print $0}'|sed -e 's/\\|$//g' >MGEtoLTR.ChrX.grep.script.txt

cat MGEtoLTR.Chr1.par_and_full.tmp.txt|grep -f MGEtoLTR.Chr1.grep.script.txt|tr '#' '\n' >MGEtoLTR.Chr1.full.gff.txt
cat MGEtoLTR.Chr2.par_and_full.tmp.txt|grep -f MGEtoLTR.Chr2.grep.script.txt|tr '#' '\n' >MGEtoLTR.Chr2.full.gff.txt
cat MGEtoLTR.Chr3.par_and_full.tmp.txt|grep -f MGEtoLTR.Chr3.grep.script.txt|tr '#' '\n' >MGEtoLTR.Chr3.full.gff.txt
cat MGEtoLTR.Chr4.par_and_full.tmp.txt|grep -f MGEtoLTR.Chr4.grep.script.txt|tr '#' '\n' >MGEtoLTR.Chr4.full.gff.txt
cat MGEtoLTR.Chr5.par_and_full.tmp.txt|grep -f MGEtoLTR.Chr5.grep.script.txt|tr '#' '\n' >MGEtoLTR.Chr5.full.gff.txt
cat MGEtoLTR.ChrX.par_and_full.tmp.txt|grep -f MGEtoLTR.ChrX.grep.script.txt|tr '#' '\n' >MGEtoLTR.ChrX.full.gff.txt

cat MGEtoLTR.Chr1.par_and_full.tmp.txt|grep -v -f MGEtoLTR.Chr1.grep.script.txt|tr '#' '\n' >MGEtoLTR.Chr1.par.gff.txt
cat MGEtoLTR.Chr2.par_and_full.tmp.txt|grep -v -f MGEtoLTR.Chr2.grep.script.txt|tr '#' '\n' >MGEtoLTR.Chr2.par.gff.txt
cat MGEtoLTR.Chr3.par_and_full.tmp.txt|grep -v -f MGEtoLTR.Chr3.grep.script.txt|tr '#' '\n' >MGEtoLTR.Chr3.par.gff.txt
cat MGEtoLTR.Chr4.par_and_full.tmp.txt|grep -v -f MGEtoLTR.Chr4.grep.script.txt|tr '#' '\n' >MGEtoLTR.Chr4.par.gff.txt
cat MGEtoLTR.Chr5.par_and_full.tmp.txt|grep -v -f MGEtoLTR.Chr5.grep.script.txt|tr '#' '\n' >MGEtoLTR.Chr5.par.gff.txt
cat MGEtoLTR.ChrX.par_and_full.tmp.txt|grep -v -f MGEtoLTR.ChrX.grep.script.txt|tr '#' '\n' >MGEtoLTR.ChrX.par.gff.txt

cd ..
mkdir MGE2
cd MGE2
cat ../MGE/MGEtoLTR.Chr1.full.gff.txt|sed -e 's/$/;status=full;/g' >MGEtoLTR.Chr1.full.gff.txt
cat ../MGE/MGEtoLTR.Chr2.full.gff.txt|sed -e 's/$/;status=full;/g' >MGEtoLTR.Chr2.full.gff.txt
cat ../MGE/MGEtoLTR.Chr3.full.gff.txt|sed -e 's/$/;status=full;/g' >MGEtoLTR.Chr3.full.gff.txt
cat ../MGE/MGEtoLTR.Chr4.full.gff.txt|sed -e 's/$/;status=full;/g' >MGEtoLTR.Chr4.full.gff.txt
cat ../MGE/MGEtoLTR.Chr5.full.gff.txt|sed -e 's/$/;status=full;/g' >MGEtoLTR.Chr5.full.gff.txt
cat ../MGE/MGEtoLTR.ChrX.full.gff.txt|sed -e 's/$/;status=full;/g' >MGEtoLTR.ChrX.full.gff.txt

cat ../MGE/MGEtoLTR.Chr1.par.gff.txt|sed -e 's/$/;status=par;/g' >MGEtoLTR.Chr1.par.gff.txt
cat ../MGE/MGEtoLTR.Chr2.par.gff.txt|sed -e 's/$/;status=par;/g' >MGEtoLTR.Chr2.par.gff.txt
cat ../MGE/MGEtoLTR.Chr3.par.gff.txt|sed -e 's/$/;status=par;/g' >MGEtoLTR.Chr3.par.gff.txt
cat ../MGE/MGEtoLTR.Chr4.par.gff.txt|sed -e 's/$/;status=par;/g' >MGEtoLTR.Chr4.par.gff.txt
cat ../MGE/MGEtoLTR.Chr5.par.gff.txt|sed -e 's/$/;status=par;/g' >MGEtoLTR.Chr5.par.gff.txt
cat ../MGE/MGEtoLTR.ChrX.par.gff.txt|sed -e 's/$/;status=par;/g' >MGEtoLTR.ChrX.par.gff.txt

cat ../MGE/MGEtoLTR.Chr1.solo.gff.txt|sed -e 's/$/;status=solo;/g' >MGEtoLTR.Chr1.solo.gff.txt
cat ../MGE/MGEtoLTR.Chr2.solo.gff.txt|sed -e 's/$/;status=solo;/g' >MGEtoLTR.Chr2.solo.gff.txt
cat ../MGE/MGEtoLTR.Chr3.solo.gff.txt|sed -e 's/$/;status=solo;/g' >MGEtoLTR.Chr3.solo.gff.txt
cat ../MGE/MGEtoLTR.Chr4.solo.gff.txt|sed -e 's/$/;status=solo;/g' >MGEtoLTR.Chr4.solo.gff.txt
cat ../MGE/MGEtoLTR.Chr5.solo.gff.txt|sed -e 's/$/;status=solo;/g' >MGEtoLTR.Chr5.solo.gff.txt
cat ../MGE/MGEtoLTR.ChrX.solo.gff.txt|sed -e 's/$/;status=solo;/g' >MGEtoLTR.ChrX.solo.gff.txt

cd ..
mkdir Combine
cd Combine
cat ../sep5/LTRharvest.Chr1.gff.txt ../MGE2/MGEtoLTR.Chr1.full.gff.txt ../MGE2/MGEtoLTR.Chr1.par.gff.txt ../MGE2/MGEtoLTR.Chr1.solo.gff.txt|sort -n -k4,4 >Combined.Chr1.gff.txt
cat ../sep5/LTRharvest.Chr2.gff.txt ../MGE2/MGEtoLTR.Chr2.full.gff.txt ../MGE2/MGEtoLTR.Chr2.par.gff.txt ../MGE2/MGEtoLTR.Chr2.solo.gff.txt|sort -n -k4,4 >Combined.Chr2.gff.txt
cat ../sep5/LTRharvest.Chr3.gff.txt ../MGE2/MGEtoLTR.Chr3.full.gff.txt ../MGE2/MGEtoLTR.Chr3.par.gff.txt ../MGE2/MGEtoLTR.Chr3.solo.gff.txt|sort -n -k4,4 >Combined.Chr3.gff.txt
cat ../sep5/LTRharvest.Chr4.gff.txt ../MGE2/MGEtoLTR.Chr4.full.gff.txt ../MGE2/MGEtoLTR.Chr4.par.gff.txt ../MGE2/MGEtoLTR.Chr4.solo.gff.txt|sort -n -k4,4 >Combined.Chr4.gff.txt
cat ../sep5/LTRharvest.Chr5.gff.txt ../MGE2/MGEtoLTR.Chr5.full.gff.txt ../MGE2/MGEtoLTR.Chr5.par.gff.txt ../MGE2/MGEtoLTR.Chr5.solo.gff.txt|sort -n -k4,4 >Combined.Chr5.gff.txt
cat ../sep5/LTRharvest.ChrX.gff.txt ../MGE2/MGEtoLTR.ChrX.full.gff.txt ../MGE2/MGEtoLTR.ChrX.par.gff.txt ../MGE2/MGEtoLTR.ChrX.solo.gff.txt|sort -n -k4,4 >Combined.ChrX.gff.txt

cd ..
mkdir protein_match
cd protein_match
mkdir ENV
mkdir GAG
mkdir GAGCOAT
mkdir INT
mkdir PR
mkdir Peptidase_A17
mkdir RH
mkdir RNase
mkdir RT
mkdir Retrotrans
mkdir galadriel
mkdir zf-CCHC

cd ENV
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=ENV")print}' >protein_match.Chr1.ENV.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=ENV")print}' >protein_match.Chr2.ENV.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=ENV")print}' >protein_match.Chr3.ENV.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=ENV")print}' >protein_match.Chr4.ENV.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=ENV")print}' >protein_match.Chr5.ENV.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=ENV")print}' >protein_match.ChrX.ENV.gff.txt
cd ../GAG
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAG")print}' >protein_match.Chr1.GAG.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAG")print}' >protein_match.Chr2.GAG.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAG")print}' >protein_match.Chr3.GAG.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAG")print}' >protein_match.Chr4.GAG.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAG")print}' >protein_match.Chr5.GAG.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAG")print}' >protein_match.ChrX.GAG.gff.txt
cd ../GAGCOAT
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAGCOAT")print}' >protein_match.Chr1.GAGCOAT.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAGCOAT")print}' >protein_match.Chr2.GAGCOAT.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAGCOAT")print}' >protein_match.Chr3.GAGCOAT.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAGCOAT")print}' >protein_match.Chr4.GAGCOAT.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAGCOAT")print}' >protein_match.Chr5.GAGCOAT.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=GAGCOAT")print}' >protein_match.ChrX.GAGCOAT.gff.txt
cd ../INT
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=INT")print}' >protein_match.Chr1.INT.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=INT")print}' >protein_match.Chr2.INT.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=INT")print}' >protein_match.Chr3.INT.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=INT")print}' >protein_match.Chr4.INT.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=INT")print}' >protein_match.Chr5.INT.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=INT")print}' >protein_match.ChrX.INT.gff.txt
cd ../PR
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=PR")print}' >protein_match.Chr1.PR.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=PR")print}' >protein_match.Chr2.PR.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=PR")print}' >protein_match.Chr3.PR.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=PR")print}' >protein_match.Chr4.PR.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=PR")print}' >protein_match.Chr5.PR.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=PR")print}' >protein_match.ChrX.PR.gff.txt
cd ../Peptidase_A17
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Peptidase_A17")print}' >protein_match.Chr1.Peptidase_A17.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Peptidase_A17")print}' >protein_match.Chr2.Peptidase_A17.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Peptidase_A17")print}' >protein_match.Chr3.Peptidase_A17.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Peptidase_A17")print}' >protein_match.Chr4.Peptidase_A17.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Peptidase_A17")print}' >protein_match.Chr5.Peptidase_A17.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Peptidase_A17")print}' >protein_match.ChrX.Peptidase_A17.gff.txt
cd ../RH
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RH")print}' >protein_match.Chr1.RH.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RH")print}' >protein_match.Chr2.RH.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RH")print}' >protein_match.Chr3.RH.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RH")print}' >protein_match.Chr4.RH.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RH")print}' >protein_match.Chr5.RH.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RH")print}' >protein_match.ChrX.RH.gff.txt
cd ../RNase
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RNase")print}' >protein_match.Chr1.RNase.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RNase")print}' >protein_match.Chr2.RNase.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RNase")print}' >protein_match.Chr3.RNase.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RNase")print}' >protein_match.Chr4.RNase.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RNase")print}' >protein_match.Chr5.RNase.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RNase")print}' >protein_match.ChrX.RNase.gff.txt
cd ../RT
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RT")print}' >protein_match.Chr1.RT.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RT")print}' >protein_match.Chr2.RT.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RT")print}' >protein_match.Chr3.RT.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RT")print}' >protein_match.Chr4.RT.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RT")print}' >protein_match.Chr5.RT.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=RT")print}' >protein_match.ChrX.RT.gff.txt
cd ../Retrotrans
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Retrotrans")print}' >protein_match.Chr1.Retrotrans.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Retrotrans")print}' >protein_match.Chr2.Retrotrans.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Retrotrans")print}' >protein_match.Chr3.Retrotrans.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Retrotrans")print}' >protein_match.Chr4.Retrotrans.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Retrotrans")print}' >protein_match.Chr5.Retrotrans.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=Retrotrans")print}' >protein_match.ChrX.Retrotrans.gff.txt
cd ../galadriel
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=galadriel")print}' >protein_match.Chr1.galadriel.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=galadriel")print}' >protein_match.Chr2.galadriel.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=galadriel")print}' >protein_match.Chr3.galadriel.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=galadriel")print}' >protein_match.Chr4.galadriel.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=galadriel")print}' >protein_match.Chr5.galadriel.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=galadriel")print}' >protein_match.ChrX.galadriel.gff.txt
cd ../zf-CCHC
cat ../../Combine/Combined.Chr1.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=zf-CCHC")print}' >protein_match.Chr1.zf-CCHC.gff.txt
cat ../../Combine/Combined.Chr2.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=zf-CCHC")print}' >protein_match.Chr2.zf-CCHC.gff.txt
cat ../../Combine/Combined.Chr3.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=zf-CCHC")print}' >protein_match.Chr3.zf-CCHC.gff.txt
cat ../../Combine/Combined.Chr4.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=zf-CCHC")print}' >protein_match.Chr4.zf-CCHC.gff.txt
cat ../../Combine/Combined.Chr5.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=zf-CCHC")print}' >protein_match.Chr5.zf-CCHC.gff.txt
cat ../../Combine/Combined.ChrX.gff.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$9~"name=zf-CCHC")print}' >protein_match.ChrX.zf-CCHC.gff.txt

cd ../ENV
cat protein_match.Chr1.ENV.gff.txt|grep 'MGEScan' >protein_match.Chr1.ENV.MGEScan.gff.txt
cat protein_match.Chr1.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.ENV.MGEScan.list.txt
cat protein_match.Chr1.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.ENV.MGEScan.awk.script1.txt
cat protein_match.Chr1.ENV.MGEScan.list.txt|awk -f protein_match.Chr1.ENV.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.ENV.MGEScan.awk.script2.txt
cat protein_match.Chr1.ENV.gff.txt|awk -f protein_match.Chr1.ENV.MGEScan.awk.script2.txt >protein_match.Chr1.ENV.MGEScan.min.gff.txt

cat protein_match.Chr1.ENV.gff.txt|grep 'LTRdigest' >protein_match.Chr1.ENV.LTRdigest.gff.txt
cat protein_match.Chr1.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.ENV.LTRdigest.list.txt
cat protein_match.Chr1.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.ENV.LTRdigest.awk.script1.txt
cat protein_match.Chr1.ENV.LTRdigest.list.txt|awk -f protein_match.Chr1.ENV.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.ENV.LTRdigest.awk.script2.txt
cat protein_match.Chr1.ENV.gff.txt|awk -f protein_match.Chr1.ENV.LTRdigest.awk.script2.txt >protein_match.Chr1.ENV.LTRdigest.min.gff.txt

cat protein_match.Chr1.ENV.MGEScan.min.gff.txt protein_match.Chr1.ENV.LTRdigest.min.gff.txt >../protein_match.Chr1.ENV.min.gff.txt

cat protein_match.Chr2.ENV.gff.txt|grep 'MGEScan' >protein_match.Chr2.ENV.MGEScan.gff.txt
cat protein_match.Chr2.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.ENV.MGEScan.list.txt
cat protein_match.Chr2.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.ENV.MGEScan.awk.script1.txt
cat protein_match.Chr2.ENV.MGEScan.list.txt|awk -f protein_match.Chr2.ENV.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.ENV.MGEScan.awk.script2.txt
cat protein_match.Chr2.ENV.gff.txt|awk -f protein_match.Chr2.ENV.MGEScan.awk.script2.txt >protein_match.Chr2.ENV.MGEScan.min.gff.txt

cat protein_match.Chr2.ENV.gff.txt|grep 'LTRdigest' >protein_match.Chr2.ENV.LTRdigest.gff.txt
cat protein_match.Chr2.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.ENV.LTRdigest.list.txt
cat protein_match.Chr2.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.ENV.LTRdigest.awk.script1.txt
cat protein_match.Chr2.ENV.LTRdigest.list.txt|awk -f protein_match.Chr2.ENV.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.ENV.LTRdigest.awk.script2.txt
cat protein_match.Chr2.ENV.gff.txt|awk -f protein_match.Chr2.ENV.LTRdigest.awk.script2.txt >protein_match.Chr2.ENV.LTRdigest.min.gff.txt

cat protein_match.Chr2.ENV.MGEScan.min.gff.txt protein_match.Chr2.ENV.LTRdigest.min.gff.txt >../protein_match.Chr2.ENV.min.gff.txt

cat protein_match.Chr3.ENV.gff.txt|grep 'MGEScan' >protein_match.Chr3.ENV.MGEScan.gff.txt
cat protein_match.Chr3.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.ENV.MGEScan.list.txt
cat protein_match.Chr3.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.ENV.MGEScan.awk.script1.txt
cat protein_match.Chr3.ENV.MGEScan.list.txt|awk -f protein_match.Chr3.ENV.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.ENV.MGEScan.awk.script2.txt
cat protein_match.Chr3.ENV.gff.txt|awk -f protein_match.Chr3.ENV.MGEScan.awk.script2.txt >protein_match.Chr3.ENV.MGEScan.min.gff.txt

cat protein_match.Chr3.ENV.gff.txt|grep 'LTRdigest' >protein_match.Chr3.ENV.LTRdigest.gff.txt
cat protein_match.Chr3.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.ENV.LTRdigest.list.txt
cat protein_match.Chr3.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.ENV.LTRdigest.awk.script1.txt
cat protein_match.Chr3.ENV.LTRdigest.list.txt|awk -f protein_match.Chr3.ENV.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.ENV.LTRdigest.awk.script2.txt
cat protein_match.Chr3.ENV.gff.txt|awk -f protein_match.Chr3.ENV.LTRdigest.awk.script2.txt >protein_match.Chr3.ENV.LTRdigest.min.gff.txt

cat protein_match.Chr3.ENV.MGEScan.min.gff.txt protein_match.Chr3.ENV.LTRdigest.min.gff.txt >../protein_match.Chr3.ENV.min.gff.txt

cat protein_match.Chr4.ENV.gff.txt|grep 'MGEScan' >protein_match.Chr4.ENV.MGEScan.gff.txt
cat protein_match.Chr4.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.ENV.MGEScan.list.txt
cat protein_match.Chr4.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.ENV.MGEScan.awk.script1.txt
cat protein_match.Chr4.ENV.MGEScan.list.txt|awk -f protein_match.Chr4.ENV.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.ENV.MGEScan.awk.script2.txt
cat protein_match.Chr4.ENV.gff.txt|awk -f protein_match.Chr4.ENV.MGEScan.awk.script2.txt >protein_match.Chr4.ENV.MGEScan.min.gff.txt

cat protein_match.Chr4.ENV.gff.txt|grep 'LTRdigest' >protein_match.Chr4.ENV.LTRdigest.gff.txt
cat protein_match.Chr4.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.ENV.LTRdigest.list.txt
cat protein_match.Chr4.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.ENV.LTRdigest.awk.script1.txt
cat protein_match.Chr4.ENV.LTRdigest.list.txt|awk -f protein_match.Chr4.ENV.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.ENV.LTRdigest.awk.script2.txt
cat protein_match.Chr4.ENV.gff.txt|awk -f protein_match.Chr4.ENV.LTRdigest.awk.script2.txt >protein_match.Chr4.ENV.LTRdigest.min.gff.txt

cat protein_match.Chr4.ENV.MGEScan.min.gff.txt protein_match.Chr4.ENV.LTRdigest.min.gff.txt >../protein_match.Chr4.ENV.min.gff.txt

cat protein_match.Chr5.ENV.gff.txt|grep 'MGEScan' >protein_match.Chr5.ENV.MGEScan.gff.txt
cat protein_match.Chr5.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.ENV.MGEScan.list.txt
cat protein_match.Chr5.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.ENV.MGEScan.awk.script1.txt
cat protein_match.Chr5.ENV.MGEScan.list.txt|awk -f protein_match.Chr5.ENV.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.ENV.MGEScan.awk.script2.txt
cat protein_match.Chr5.ENV.gff.txt|awk -f protein_match.Chr5.ENV.MGEScan.awk.script2.txt >protein_match.Chr5.ENV.MGEScan.min.gff.txt

cat protein_match.Chr5.ENV.gff.txt|grep 'LTRdigest' >protein_match.Chr5.ENV.LTRdigest.gff.txt
cat protein_match.Chr5.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.ENV.LTRdigest.list.txt
cat protein_match.Chr5.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.ENV.LTRdigest.awk.script1.txt
cat protein_match.Chr5.ENV.LTRdigest.list.txt|awk -f protein_match.Chr5.ENV.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.ENV.LTRdigest.awk.script2.txt
cat protein_match.Chr5.ENV.gff.txt|awk -f protein_match.Chr5.ENV.LTRdigest.awk.script2.txt >protein_match.Chr5.ENV.LTRdigest.min.gff.txt

cat protein_match.Chr5.ENV.MGEScan.min.gff.txt protein_match.Chr5.ENV.LTRdigest.min.gff.txt >../protein_match.Chr5.ENV.min.gff.txt

cat protein_match.ChrX.ENV.gff.txt|grep 'MGEScan' >protein_match.ChrX.ENV.MGEScan.gff.txt
cat protein_match.ChrX.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.ENV.MGEScan.list.txt
cat protein_match.ChrX.ENV.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.ENV.MGEScan.awk.script1.txt
cat protein_match.ChrX.ENV.MGEScan.list.txt|awk -f protein_match.ChrX.ENV.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.ENV.MGEScan.awk.script2.txt
cat protein_match.ChrX.ENV.gff.txt|awk -f protein_match.ChrX.ENV.MGEScan.awk.script2.txt >protein_match.ChrX.ENV.MGEScan.min.gff.txt

cat protein_match.ChrX.ENV.gff.txt|grep 'LTRdigest' >protein_match.ChrX.ENV.LTRdigest.gff.txt
cat protein_match.ChrX.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.ENV.LTRdigest.list.txt
cat protein_match.ChrX.ENV.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.ENV.LTRdigest.awk.script1.txt
cat protein_match.ChrX.ENV.LTRdigest.list.txt|awk -f protein_match.ChrX.ENV.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.ENV.LTRdigest.awk.script2.txt
cat protein_match.ChrX.ENV.gff.txt|awk -f protein_match.ChrX.ENV.LTRdigest.awk.script2.txt >protein_match.ChrX.ENV.LTRdigest.min.gff.txt

cat protein_match.ChrX.ENV.MGEScan.min.gff.txt protein_match.ChrX.ENV.LTRdigest.min.gff.txt >../protein_match.ChrX.ENV.min.gff.txt

cd ../GAG
cat protein_match.Chr1.GAG.gff.txt|grep 'MGEScan' >protein_match.Chr1.GAG.MGEScan.gff.txt
cat protein_match.Chr1.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.GAG.MGEScan.list.txt
cat protein_match.Chr1.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.GAG.MGEScan.awk.script1.txt
cat protein_match.Chr1.GAG.MGEScan.list.txt|awk -f protein_match.Chr1.GAG.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.GAG.MGEScan.awk.script2.txt
cat protein_match.Chr1.GAG.gff.txt|awk -f protein_match.Chr1.GAG.MGEScan.awk.script2.txt >protein_match.Chr1.GAG.MGEScan.min.gff.txt

cat protein_match.Chr1.GAG.gff.txt|grep 'LTRdigest' >protein_match.Chr1.GAG.LTRdigest.gff.txt
cat protein_match.Chr1.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.GAG.LTRdigest.list.txt
cat protein_match.Chr1.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.GAG.LTRdigest.awk.script1.txt
cat protein_match.Chr1.GAG.LTRdigest.list.txt|awk -f protein_match.Chr1.GAG.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.GAG.LTRdigest.awk.script2.txt
cat protein_match.Chr1.GAG.gff.txt|awk -f protein_match.Chr1.GAG.LTRdigest.awk.script2.txt >protein_match.Chr1.GAG.LTRdigest.min.gff.txt

cat protein_match.Chr1.GAG.MGEScan.min.gff.txt protein_match.Chr1.GAG.LTRdigest.min.gff.txt >../protein_match.Chr1.GAG.min.gff.txt

cat protein_match.Chr2.GAG.gff.txt|grep 'MGEScan' >protein_match.Chr2.GAG.MGEScan.gff.txt
cat protein_match.Chr2.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.GAG.MGEScan.list.txt
cat protein_match.Chr2.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.GAG.MGEScan.awk.script1.txt
cat protein_match.Chr2.GAG.MGEScan.list.txt|awk -f protein_match.Chr2.GAG.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.GAG.MGEScan.awk.script2.txt
cat protein_match.Chr2.GAG.gff.txt|awk -f protein_match.Chr2.GAG.MGEScan.awk.script2.txt >protein_match.Chr2.GAG.MGEScan.min.gff.txt

cat protein_match.Chr2.GAG.gff.txt|grep 'LTRdigest' >protein_match.Chr2.GAG.LTRdigest.gff.txt
cat protein_match.Chr2.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.GAG.LTRdigest.list.txt
cat protein_match.Chr2.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.GAG.LTRdigest.awk.script1.txt
cat protein_match.Chr2.GAG.LTRdigest.list.txt|awk -f protein_match.Chr2.GAG.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.GAG.LTRdigest.awk.script2.txt
cat protein_match.Chr2.GAG.gff.txt|awk -f protein_match.Chr2.GAG.LTRdigest.awk.script2.txt >protein_match.Chr2.GAG.LTRdigest.min.gff.txt

cat protein_match.Chr2.GAG.MGEScan.min.gff.txt protein_match.Chr2.GAG.LTRdigest.min.gff.txt >../protein_match.Chr2.GAG.min.gff.txt

cat protein_match.Chr3.GAG.gff.txt|grep 'MGEScan' >protein_match.Chr3.GAG.MGEScan.gff.txt
cat protein_match.Chr3.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.GAG.MGEScan.list.txt
cat protein_match.Chr3.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.GAG.MGEScan.awk.script1.txt
cat protein_match.Chr3.GAG.MGEScan.list.txt|awk -f protein_match.Chr3.GAG.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.GAG.MGEScan.awk.script2.txt
cat protein_match.Chr3.GAG.gff.txt|awk -f protein_match.Chr3.GAG.MGEScan.awk.script2.txt >protein_match.Chr3.GAG.MGEScan.min.gff.txt

cat protein_match.Chr3.GAG.gff.txt|grep 'LTRdigest' >protein_match.Chr3.GAG.LTRdigest.gff.txt
cat protein_match.Chr3.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.GAG.LTRdigest.list.txt
cat protein_match.Chr3.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.GAG.LTRdigest.awk.script1.txt
cat protein_match.Chr3.GAG.LTRdigest.list.txt|awk -f protein_match.Chr3.GAG.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.GAG.LTRdigest.awk.script2.txt
cat protein_match.Chr3.GAG.gff.txt|awk -f protein_match.Chr3.GAG.LTRdigest.awk.script2.txt >protein_match.Chr3.GAG.LTRdigest.min.gff.txt

cat protein_match.Chr3.GAG.MGEScan.min.gff.txt protein_match.Chr3.GAG.LTRdigest.min.gff.txt >../protein_match.Chr3.GAG.min.gff.txt

cat protein_match.Chr4.GAG.gff.txt|grep 'MGEScan' >protein_match.Chr4.GAG.MGEScan.gff.txt
cat protein_match.Chr4.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.GAG.MGEScan.list.txt
cat protein_match.Chr4.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.GAG.MGEScan.awk.script1.txt
cat protein_match.Chr4.GAG.MGEScan.list.txt|awk -f protein_match.Chr4.GAG.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.GAG.MGEScan.awk.script2.txt
cat protein_match.Chr4.GAG.gff.txt|awk -f protein_match.Chr4.GAG.MGEScan.awk.script2.txt >protein_match.Chr4.GAG.MGEScan.min.gff.txt

cat protein_match.Chr4.GAG.gff.txt|grep 'LTRdigest' >protein_match.Chr4.GAG.LTRdigest.gff.txt
cat protein_match.Chr4.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.GAG.LTRdigest.list.txt
cat protein_match.Chr4.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.GAG.LTRdigest.awk.script1.txt
cat protein_match.Chr4.GAG.LTRdigest.list.txt|awk -f protein_match.Chr4.GAG.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.GAG.LTRdigest.awk.script2.txt
cat protein_match.Chr4.GAG.gff.txt|awk -f protein_match.Chr4.GAG.LTRdigest.awk.script2.txt >protein_match.Chr4.GAG.LTRdigest.min.gff.txt

cat protein_match.Chr4.GAG.MGEScan.min.gff.txt protein_match.Chr4.GAG.LTRdigest.min.gff.txt >../protein_match.Chr4.GAG.min.gff.txt

cat protein_match.Chr5.GAG.gff.txt|grep 'MGEScan' >protein_match.Chr5.GAG.MGEScan.gff.txt
cat protein_match.Chr5.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.GAG.MGEScan.list.txt
cat protein_match.Chr5.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.GAG.MGEScan.awk.script1.txt
cat protein_match.Chr5.GAG.MGEScan.list.txt|awk -f protein_match.Chr5.GAG.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.GAG.MGEScan.awk.script2.txt
cat protein_match.Chr5.GAG.gff.txt|awk -f protein_match.Chr5.GAG.MGEScan.awk.script2.txt >protein_match.Chr5.GAG.MGEScan.min.gff.txt

cat protein_match.Chr5.GAG.gff.txt|grep 'LTRdigest' >protein_match.Chr5.GAG.LTRdigest.gff.txt
cat protein_match.Chr5.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.GAG.LTRdigest.list.txt
cat protein_match.Chr5.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.GAG.LTRdigest.awk.script1.txt
cat protein_match.Chr5.GAG.LTRdigest.list.txt|awk -f protein_match.Chr5.GAG.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.GAG.LTRdigest.awk.script2.txt
cat protein_match.Chr5.GAG.gff.txt|awk -f protein_match.Chr5.GAG.LTRdigest.awk.script2.txt >protein_match.Chr5.GAG.LTRdigest.min.gff.txt

cat protein_match.Chr5.GAG.MGEScan.min.gff.txt protein_match.Chr5.GAG.LTRdigest.min.gff.txt >../protein_match.Chr5.GAG.min.gff.txt

cat protein_match.ChrX.GAG.gff.txt|grep 'MGEScan' >protein_match.ChrX.GAG.MGEScan.gff.txt
cat protein_match.ChrX.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.GAG.MGEScan.list.txt
cat protein_match.ChrX.GAG.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.GAG.MGEScan.awk.script1.txt
cat protein_match.ChrX.GAG.MGEScan.list.txt|awk -f protein_match.ChrX.GAG.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.GAG.MGEScan.awk.script2.txt
cat protein_match.ChrX.GAG.gff.txt|awk -f protein_match.ChrX.GAG.MGEScan.awk.script2.txt >protein_match.ChrX.GAG.MGEScan.min.gff.txt

cat protein_match.ChrX.GAG.gff.txt|grep 'LTRdigest' >protein_match.ChrX.GAG.LTRdigest.gff.txt
cat protein_match.ChrX.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.GAG.LTRdigest.list.txt
cat protein_match.ChrX.GAG.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.GAG.LTRdigest.awk.script1.txt
cat protein_match.ChrX.GAG.LTRdigest.list.txt|awk -f protein_match.ChrX.GAG.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.GAG.LTRdigest.awk.script2.txt
cat protein_match.ChrX.GAG.gff.txt|awk -f protein_match.ChrX.GAG.LTRdigest.awk.script2.txt >protein_match.ChrX.GAG.LTRdigest.min.gff.txt

cat protein_match.ChrX.GAG.MGEScan.min.gff.txt protein_match.ChrX.GAG.LTRdigest.min.gff.txt >../protein_match.ChrX.GAG.min.gff.txt

cd ../GAGCOAT
cat protein_match.Chr1.GAGCOAT.gff.txt|grep 'MGEScan' >protein_match.Chr1.GAGCOAT.MGEScan.gff.txt
cat protein_match.Chr1.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.GAGCOAT.MGEScan.list.txt
cat protein_match.Chr1.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.GAGCOAT.MGEScan.awk.script1.txt
cat protein_match.Chr1.GAGCOAT.MGEScan.list.txt|awk -f protein_match.Chr1.GAGCOAT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.GAGCOAT.MGEScan.awk.script2.txt
cat protein_match.Chr1.GAGCOAT.gff.txt|awk -f protein_match.Chr1.GAGCOAT.MGEScan.awk.script2.txt >protein_match.Chr1.GAGCOAT.MGEScan.min.gff.txt

cat protein_match.Chr1.GAGCOAT.gff.txt|grep 'LTRdigest' >protein_match.Chr1.GAGCOAT.LTRdigest.gff.txt
cat protein_match.Chr1.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.GAGCOAT.LTRdigest.list.txt
cat protein_match.Chr1.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.GAGCOAT.LTRdigest.awk.script1.txt
cat protein_match.Chr1.GAGCOAT.LTRdigest.list.txt|awk -f protein_match.Chr1.GAGCOAT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.GAGCOAT.LTRdigest.awk.script2.txt
cat protein_match.Chr1.GAGCOAT.gff.txt|awk -f protein_match.Chr1.GAGCOAT.LTRdigest.awk.script2.txt >protein_match.Chr1.GAGCOAT.LTRdigest.min.gff.txt

cat protein_match.Chr1.GAGCOAT.MGEScan.min.gff.txt protein_match.Chr1.GAGCOAT.LTRdigest.min.gff.txt >../protein_match.Chr1.GAGCOAT.min.gff.txt

cat protein_match.Chr2.GAGCOAT.gff.txt|grep 'MGEScan' >protein_match.Chr2.GAGCOAT.MGEScan.gff.txt
cat protein_match.Chr2.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.GAGCOAT.MGEScan.list.txt
cat protein_match.Chr2.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.GAGCOAT.MGEScan.awk.script1.txt
cat protein_match.Chr2.GAGCOAT.MGEScan.list.txt|awk -f protein_match.Chr2.GAGCOAT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.GAGCOAT.MGEScan.awk.script2.txt
cat protein_match.Chr2.GAGCOAT.gff.txt|awk -f protein_match.Chr2.GAGCOAT.MGEScan.awk.script2.txt >protein_match.Chr2.GAGCOAT.MGEScan.min.gff.txt

cat protein_match.Chr2.GAGCOAT.gff.txt|grep 'LTRdigest' >protein_match.Chr2.GAGCOAT.LTRdigest.gff.txt
cat protein_match.Chr2.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.GAGCOAT.LTRdigest.list.txt
cat protein_match.Chr2.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.GAGCOAT.LTRdigest.awk.script1.txt
cat protein_match.Chr2.GAGCOAT.LTRdigest.list.txt|awk -f protein_match.Chr2.GAGCOAT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.GAGCOAT.LTRdigest.awk.script2.txt
cat protein_match.Chr2.GAGCOAT.gff.txt|awk -f protein_match.Chr2.GAGCOAT.LTRdigest.awk.script2.txt >protein_match.Chr2.GAGCOAT.LTRdigest.min.gff.txt

cat protein_match.Chr2.GAGCOAT.MGEScan.min.gff.txt protein_match.Chr2.GAGCOAT.LTRdigest.min.gff.txt >../protein_match.Chr2.GAGCOAT.min.gff.txt

cat protein_match.Chr3.GAGCOAT.gff.txt|grep 'MGEScan' >protein_match.Chr3.GAGCOAT.MGEScan.gff.txt
cat protein_match.Chr3.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.GAGCOAT.MGEScan.list.txt
cat protein_match.Chr3.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.GAGCOAT.MGEScan.awk.script1.txt
cat protein_match.Chr3.GAGCOAT.MGEScan.list.txt|awk -f protein_match.Chr3.GAGCOAT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.GAGCOAT.MGEScan.awk.script2.txt
cat protein_match.Chr3.GAGCOAT.gff.txt|awk -f protein_match.Chr3.GAGCOAT.MGEScan.awk.script2.txt >protein_match.Chr3.GAGCOAT.MGEScan.min.gff.txt

cat protein_match.Chr3.GAGCOAT.gff.txt|grep 'LTRdigest' >protein_match.Chr3.GAGCOAT.LTRdigest.gff.txt
cat protein_match.Chr3.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.GAGCOAT.LTRdigest.list.txt
cat protein_match.Chr3.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.GAGCOAT.LTRdigest.awk.script1.txt
cat protein_match.Chr3.GAGCOAT.LTRdigest.list.txt|awk -f protein_match.Chr3.GAGCOAT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.GAGCOAT.LTRdigest.awk.script2.txt
cat protein_match.Chr3.GAGCOAT.gff.txt|awk -f protein_match.Chr3.GAGCOAT.LTRdigest.awk.script2.txt >protein_match.Chr3.GAGCOAT.LTRdigest.min.gff.txt

cat protein_match.Chr3.GAGCOAT.MGEScan.min.gff.txt protein_match.Chr3.GAGCOAT.LTRdigest.min.gff.txt >../protein_match.Chr3.GAGCOAT.min.gff.txt

cat protein_match.Chr4.GAGCOAT.gff.txt|grep 'MGEScan' >protein_match.Chr4.GAGCOAT.MGEScan.gff.txt
cat protein_match.Chr4.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.GAGCOAT.MGEScan.list.txt
cat protein_match.Chr4.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.GAGCOAT.MGEScan.awk.script1.txt
cat protein_match.Chr4.GAGCOAT.MGEScan.list.txt|awk -f protein_match.Chr4.GAGCOAT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.GAGCOAT.MGEScan.awk.script2.txt
cat protein_match.Chr4.GAGCOAT.gff.txt|awk -f protein_match.Chr4.GAGCOAT.MGEScan.awk.script2.txt >protein_match.Chr4.GAGCOAT.MGEScan.min.gff.txt

cat protein_match.Chr4.GAGCOAT.gff.txt|grep 'LTRdigest' >protein_match.Chr4.GAGCOAT.LTRdigest.gff.txt
cat protein_match.Chr4.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.GAGCOAT.LTRdigest.list.txt
cat protein_match.Chr4.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.GAGCOAT.LTRdigest.awk.script1.txt
cat protein_match.Chr4.GAGCOAT.LTRdigest.list.txt|awk -f protein_match.Chr4.GAGCOAT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.GAGCOAT.LTRdigest.awk.script2.txt
cat protein_match.Chr4.GAGCOAT.gff.txt|awk -f protein_match.Chr4.GAGCOAT.LTRdigest.awk.script2.txt >protein_match.Chr4.GAGCOAT.LTRdigest.min.gff.txt

cat protein_match.Chr4.GAGCOAT.MGEScan.min.gff.txt protein_match.Chr4.GAGCOAT.LTRdigest.min.gff.txt >../protein_match.Chr4.GAGCOAT.min.gff.txt

cat protein_match.Chr5.GAGCOAT.gff.txt|grep 'MGEScan' >protein_match.Chr5.GAGCOAT.MGEScan.gff.txt
cat protein_match.Chr5.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.GAGCOAT.MGEScan.list.txt
cat protein_match.Chr5.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.GAGCOAT.MGEScan.awk.script1.txt
cat protein_match.Chr5.GAGCOAT.MGEScan.list.txt|awk -f protein_match.Chr5.GAGCOAT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.GAGCOAT.MGEScan.awk.script2.txt
cat protein_match.Chr5.GAGCOAT.gff.txt|awk -f protein_match.Chr5.GAGCOAT.MGEScan.awk.script2.txt >protein_match.Chr5.GAGCOAT.MGEScan.min.gff.txt

cat protein_match.Chr5.GAGCOAT.gff.txt|grep 'LTRdigest' >protein_match.Chr5.GAGCOAT.LTRdigest.gff.txt
cat protein_match.Chr5.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.GAGCOAT.LTRdigest.list.txt
cat protein_match.Chr5.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.GAGCOAT.LTRdigest.awk.script1.txt
cat protein_match.Chr5.GAGCOAT.LTRdigest.list.txt|awk -f protein_match.Chr5.GAGCOAT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.GAGCOAT.LTRdigest.awk.script2.txt
cat protein_match.Chr5.GAGCOAT.gff.txt|awk -f protein_match.Chr5.GAGCOAT.LTRdigest.awk.script2.txt >protein_match.Chr5.GAGCOAT.LTRdigest.min.gff.txt

cat protein_match.Chr5.GAGCOAT.MGEScan.min.gff.txt protein_match.Chr5.GAGCOAT.LTRdigest.min.gff.txt >../protein_match.Chr5.GAGCOAT.min.gff.txt

cat protein_match.ChrX.GAGCOAT.gff.txt|grep 'MGEScan' >protein_match.ChrX.GAGCOAT.MGEScan.gff.txt
cat protein_match.ChrX.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.GAGCOAT.MGEScan.list.txt
cat protein_match.ChrX.GAGCOAT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.GAGCOAT.MGEScan.awk.script1.txt
cat protein_match.ChrX.GAGCOAT.MGEScan.list.txt|awk -f protein_match.ChrX.GAGCOAT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.GAGCOAT.MGEScan.awk.script2.txt
cat protein_match.ChrX.GAGCOAT.gff.txt|awk -f protein_match.ChrX.GAGCOAT.MGEScan.awk.script2.txt >protein_match.ChrX.GAGCOAT.MGEScan.min.gff.txt

cat protein_match.ChrX.GAGCOAT.gff.txt|grep 'LTRdigest' >protein_match.ChrX.GAGCOAT.LTRdigest.gff.txt
cat protein_match.ChrX.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.GAGCOAT.LTRdigest.list.txt
cat protein_match.ChrX.GAGCOAT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.GAGCOAT.LTRdigest.awk.script1.txt
cat protein_match.ChrX.GAGCOAT.LTRdigest.list.txt|awk -f protein_match.ChrX.GAGCOAT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.GAGCOAT.LTRdigest.awk.script2.txt
cat protein_match.ChrX.GAGCOAT.gff.txt|awk -f protein_match.ChrX.GAGCOAT.LTRdigest.awk.script2.txt >protein_match.ChrX.GAGCOAT.LTRdigest.min.gff.txt

cat protein_match.ChrX.GAGCOAT.MGEScan.min.gff.txt protein_match.ChrX.GAGCOAT.LTRdigest.min.gff.txt >../protein_match.ChrX.GAGCOAT.min.gff.txt

cd ../INT
cat protein_match.Chr1.INT.gff.txt|grep 'MGEScan' >protein_match.Chr1.INT.MGEScan.gff.txt
cat protein_match.Chr1.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.INT.MGEScan.list.txt
cat protein_match.Chr1.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.INT.MGEScan.awk.script1.txt
cat protein_match.Chr1.INT.MGEScan.list.txt|awk -f protein_match.Chr1.INT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.INT.MGEScan.awk.script2.txt
cat protein_match.Chr1.INT.gff.txt|awk -f protein_match.Chr1.INT.MGEScan.awk.script2.txt >protein_match.Chr1.INT.MGEScan.min.gff.txt

cat protein_match.Chr1.INT.gff.txt|grep 'LTRdigest' >protein_match.Chr1.INT.LTRdigest.gff.txt
cat protein_match.Chr1.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.INT.LTRdigest.list.txt
cat protein_match.Chr1.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.INT.LTRdigest.awk.script1.txt
cat protein_match.Chr1.INT.LTRdigest.list.txt|awk -f protein_match.Chr1.INT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.INT.LTRdigest.awk.script2.txt
cat protein_match.Chr1.INT.gff.txt|awk -f protein_match.Chr1.INT.LTRdigest.awk.script2.txt >protein_match.Chr1.INT.LTRdigest.min.gff.txt

cat protein_match.Chr1.INT.MGEScan.min.gff.txt protein_match.Chr1.INT.LTRdigest.min.gff.txt >../protein_match.Chr1.INT.min.gff.txt

cat protein_match.Chr2.INT.gff.txt|grep 'MGEScan' >protein_match.Chr2.INT.MGEScan.gff.txt
cat protein_match.Chr2.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.INT.MGEScan.list.txt
cat protein_match.Chr2.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.INT.MGEScan.awk.script1.txt
cat protein_match.Chr2.INT.MGEScan.list.txt|awk -f protein_match.Chr2.INT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.INT.MGEScan.awk.script2.txt
cat protein_match.Chr2.INT.gff.txt|awk -f protein_match.Chr2.INT.MGEScan.awk.script2.txt >protein_match.Chr2.INT.MGEScan.min.gff.txt

cat protein_match.Chr2.INT.gff.txt|grep 'LTRdigest' >protein_match.Chr2.INT.LTRdigest.gff.txt
cat protein_match.Chr2.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.INT.LTRdigest.list.txt
cat protein_match.Chr2.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.INT.LTRdigest.awk.script1.txt
cat protein_match.Chr2.INT.LTRdigest.list.txt|awk -f protein_match.Chr2.INT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.INT.LTRdigest.awk.script2.txt
cat protein_match.Chr2.INT.gff.txt|awk -f protein_match.Chr2.INT.LTRdigest.awk.script2.txt >protein_match.Chr2.INT.LTRdigest.min.gff.txt

cat protein_match.Chr2.INT.MGEScan.min.gff.txt protein_match.Chr2.INT.LTRdigest.min.gff.txt >../protein_match.Chr2.INT.min.gff.txt

cat protein_match.Chr2.INT.gff.txt|grep 'MGEScan' >protein_match.Chr2.INT.MGEScan.gff.txt
cat protein_match.Chr2.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.INT.MGEScan.list.txt
cat protein_match.Chr2.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.INT.MGEScan.awk.script1.txt
cat protein_match.Chr2.INT.MGEScan.list.txt|awk -f protein_match.Chr2.INT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.INT.MGEScan.awk.script2.txt
cat protein_match.Chr2.INT.gff.txt|awk -f protein_match.Chr2.INT.MGEScan.awk.script2.txt >protein_match.Chr2.INT.MGEScan.min.gff.txt

cat protein_match.Chr2.INT.gff.txt|grep 'LTRdigest' >protein_match.Chr2.INT.LTRdigest.gff.txt
cat protein_match.Chr2.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.INT.LTRdigest.list.txt
cat protein_match.Chr2.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.INT.LTRdigest.awk.script1.txt
cat protein_match.Chr2.INT.LTRdigest.list.txt|awk -f protein_match.Chr2.INT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.INT.LTRdigest.awk.script2.txt
cat protein_match.Chr2.INT.gff.txt|awk -f protein_match.Chr2.INT.LTRdigest.awk.script2.txt >protein_match.Chr2.INT.LTRdigest.min.gff.txt

cat protein_match.Chr2.INT.MGEScan.min.gff.txt protein_match.Chr2.INT.LTRdigest.min.gff.txt >../protein_match.Chr2.INT.min.gff.txt

cat protein_match.Chr3.INT.gff.txt|grep 'MGEScan' >protein_match.Chr3.INT.MGEScan.gff.txt
cat protein_match.Chr3.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.INT.MGEScan.list.txt
cat protein_match.Chr3.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.INT.MGEScan.awk.script1.txt
cat protein_match.Chr3.INT.MGEScan.list.txt|awk -f protein_match.Chr3.INT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.INT.MGEScan.awk.script2.txt
cat protein_match.Chr3.INT.gff.txt|awk -f protein_match.Chr3.INT.MGEScan.awk.script2.txt >protein_match.Chr3.INT.MGEScan.min.gff.txt

cat protein_match.Chr3.INT.gff.txt|grep 'LTRdigest' >protein_match.Chr3.INT.LTRdigest.gff.txt
cat protein_match.Chr3.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.INT.LTRdigest.list.txt
cat protein_match.Chr3.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.INT.LTRdigest.awk.script1.txt
cat protein_match.Chr3.INT.LTRdigest.list.txt|awk -f protein_match.Chr3.INT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.INT.LTRdigest.awk.script2.txt
cat protein_match.Chr3.INT.gff.txt|awk -f protein_match.Chr3.INT.LTRdigest.awk.script2.txt >protein_match.Chr3.INT.LTRdigest.min.gff.txt

cat protein_match.Chr3.INT.MGEScan.min.gff.txt protein_match.Chr3.INT.LTRdigest.min.gff.txt >../protein_match.Chr3.INT.min.gff.txt

cat protein_match.Chr4.INT.gff.txt|grep 'MGEScan' >protein_match.Chr4.INT.MGEScan.gff.txt
cat protein_match.Chr4.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.INT.MGEScan.list.txt
cat protein_match.Chr4.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.INT.MGEScan.awk.script1.txt
cat protein_match.Chr4.INT.MGEScan.list.txt|awk -f protein_match.Chr4.INT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.INT.MGEScan.awk.script2.txt
cat protein_match.Chr4.INT.gff.txt|awk -f protein_match.Chr4.INT.MGEScan.awk.script2.txt >protein_match.Chr4.INT.MGEScan.min.gff.txt

cat protein_match.Chr4.INT.gff.txt|grep 'LTRdigest' >protein_match.Chr4.INT.LTRdigest.gff.txt
cat protein_match.Chr4.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.INT.LTRdigest.list.txt
cat protein_match.Chr4.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.INT.LTRdigest.awk.script1.txt
cat protein_match.Chr4.INT.LTRdigest.list.txt|awk -f protein_match.Chr4.INT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.INT.LTRdigest.awk.script2.txt
cat protein_match.Chr4.INT.gff.txt|awk -f protein_match.Chr4.INT.LTRdigest.awk.script2.txt >protein_match.Chr4.INT.LTRdigest.min.gff.txt

cat protein_match.Chr4.INT.MGEScan.min.gff.txt protein_match.Chr4.INT.LTRdigest.min.gff.txt >../protein_match.Chr4.INT.min.gff.txt

cat protein_match.Chr5.INT.gff.txt|grep 'MGEScan' >protein_match.Chr5.INT.MGEScan.gff.txt
cat protein_match.Chr5.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.INT.MGEScan.list.txt
cat protein_match.Chr5.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.INT.MGEScan.awk.script1.txt
cat protein_match.Chr5.INT.MGEScan.list.txt|awk -f protein_match.Chr5.INT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.INT.MGEScan.awk.script2.txt
cat protein_match.Chr5.INT.gff.txt|awk -f protein_match.Chr5.INT.MGEScan.awk.script2.txt >protein_match.Chr5.INT.MGEScan.min.gff.txt

cat protein_match.Chr5.INT.gff.txt|grep 'LTRdigest' >protein_match.Chr5.INT.LTRdigest.gff.txt
cat protein_match.Chr5.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.INT.LTRdigest.list.txt
cat protein_match.Chr5.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.INT.LTRdigest.awk.script1.txt
cat protein_match.Chr5.INT.LTRdigest.list.txt|awk -f protein_match.Chr5.INT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.INT.LTRdigest.awk.script2.txt
cat protein_match.Chr5.INT.gff.txt|awk -f protein_match.Chr5.INT.LTRdigest.awk.script2.txt >protein_match.Chr5.INT.LTRdigest.min.gff.txt

cat protein_match.Chr5.INT.MGEScan.min.gff.txt protein_match.Chr5.INT.LTRdigest.min.gff.txt >../protein_match.Chr5.INT.min.gff.txt

cat protein_match.ChrX.INT.gff.txt|grep 'MGEScan' >protein_match.ChrX.INT.MGEScan.gff.txt
cat protein_match.ChrX.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.INT.MGEScan.list.txt
cat protein_match.ChrX.INT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.INT.MGEScan.awk.script1.txt
cat protein_match.ChrX.INT.MGEScan.list.txt|awk -f protein_match.ChrX.INT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.INT.MGEScan.awk.script2.txt
cat protein_match.ChrX.INT.gff.txt|awk -f protein_match.ChrX.INT.MGEScan.awk.script2.txt >protein_match.ChrX.INT.MGEScan.min.gff.txt

cat protein_match.ChrX.INT.gff.txt|grep 'LTRdigest' >protein_match.ChrX.INT.LTRdigest.gff.txt
cat protein_match.ChrX.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.INT.LTRdigest.list.txt
cat protein_match.ChrX.INT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.INT.LTRdigest.awk.script1.txt
cat protein_match.ChrX.INT.LTRdigest.list.txt|awk -f protein_match.ChrX.INT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.INT.LTRdigest.awk.script2.txt
cat protein_match.ChrX.INT.gff.txt|awk -f protein_match.ChrX.INT.LTRdigest.awk.script2.txt >protein_match.ChrX.INT.LTRdigest.min.gff.txt

cat protein_match.ChrX.INT.MGEScan.min.gff.txt protein_match.ChrX.INT.LTRdigest.min.gff.txt >../protein_match.ChrX.INT.min.gff.txt

cd ../PR
cat protein_match.Chr1.PR.gff.txt|grep 'MGEScan' >protein_match.Chr1.PR.MGEScan.gff.txt
cat protein_match.Chr1.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.PR.MGEScan.list.txt
cat protein_match.Chr1.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.PR.MGEScan.awk.script1.txt
cat protein_match.Chr1.PR.MGEScan.list.txt|awk -f protein_match.Chr1.PR.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.PR.MGEScan.awk.script2.txt
cat protein_match.Chr1.PR.gff.txt|awk -f protein_match.Chr1.PR.MGEScan.awk.script2.txt >protein_match.Chr1.PR.MGEScan.min.gff.txt

cat protein_match.Chr1.PR.gff.txt|grep 'LTRdigest' >protein_match.Chr1.PR.LTRdigest.gff.txt
cat protein_match.Chr1.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.PR.LTRdigest.list.txt
cat protein_match.Chr1.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.PR.LTRdigest.awk.script1.txt
cat protein_match.Chr1.PR.LTRdigest.list.txt|awk -f protein_match.Chr1.PR.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.PR.LTRdigest.awk.script2.txt
cat protein_match.Chr1.PR.gff.txt|awk -f protein_match.Chr1.PR.LTRdigest.awk.script2.txt >protein_match.Chr1.PR.LTRdigest.min.gff.txt

cat protein_match.Chr1.PR.MGEScan.min.gff.txt protein_match.Chr1.PR.LTRdigest.min.gff.txt >../protein_match.Chr1.PR.min.gff.txt

cat protein_match.Chr2.PR.gff.txt|grep 'MGEScan' >protein_match.Chr2.PR.MGEScan.gff.txt
cat protein_match.Chr2.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.PR.MGEScan.list.txt
cat protein_match.Chr2.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.PR.MGEScan.awk.script1.txt
cat protein_match.Chr2.PR.MGEScan.list.txt|awk -f protein_match.Chr2.PR.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.PR.MGEScan.awk.script2.txt
cat protein_match.Chr2.PR.gff.txt|awk -f protein_match.Chr2.PR.MGEScan.awk.script2.txt >protein_match.Chr2.PR.MGEScan.min.gff.txt

cat protein_match.Chr2.PR.gff.txt|grep 'LTRdigest' >protein_match.Chr2.PR.LTRdigest.gff.txt
cat protein_match.Chr2.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.PR.LTRdigest.list.txt
cat protein_match.Chr2.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.PR.LTRdigest.awk.script1.txt
cat protein_match.Chr2.PR.LTRdigest.list.txt|awk -f protein_match.Chr2.PR.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.PR.LTRdigest.awk.script2.txt
cat protein_match.Chr2.PR.gff.txt|awk -f protein_match.Chr2.PR.LTRdigest.awk.script2.txt >protein_match.Chr2.PR.LTRdigest.min.gff.txt

cat protein_match.Chr2.PR.MGEScan.min.gff.txt protein_match.Chr2.PR.LTRdigest.min.gff.txt >../protein_match.Chr2.PR.min.gff.txt

cat protein_match.Chr3.PR.gff.txt|grep 'MGEScan' >protein_match.Chr3.PR.MGEScan.gff.txt
cat protein_match.Chr3.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.PR.MGEScan.list.txt
cat protein_match.Chr3.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.PR.MGEScan.awk.script1.txt
cat protein_match.Chr3.PR.MGEScan.list.txt|awk -f protein_match.Chr3.PR.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.PR.MGEScan.awk.script2.txt
cat protein_match.Chr3.PR.gff.txt|awk -f protein_match.Chr3.PR.MGEScan.awk.script2.txt >protein_match.Chr3.PR.MGEScan.min.gff.txt

cat protein_match.Chr3.PR.gff.txt|grep 'LTRdigest' >protein_match.Chr3.PR.LTRdigest.gff.txt
cat protein_match.Chr3.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.PR.LTRdigest.list.txt
cat protein_match.Chr3.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.PR.LTRdigest.awk.script1.txt
cat protein_match.Chr3.PR.LTRdigest.list.txt|awk -f protein_match.Chr3.PR.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.PR.LTRdigest.awk.script2.txt
cat protein_match.Chr3.PR.gff.txt|awk -f protein_match.Chr3.PR.LTRdigest.awk.script2.txt >protein_match.Chr3.PR.LTRdigest.min.gff.txt

cat protein_match.Chr3.PR.MGEScan.min.gff.txt protein_match.Chr3.PR.LTRdigest.min.gff.txt >../protein_match.Chr3.PR.min.gff.txt

cat protein_match.Chr4.PR.gff.txt|grep 'MGEScan' >protein_match.Chr4.PR.MGEScan.gff.txt
cat protein_match.Chr4.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.PR.MGEScan.list.txt
cat protein_match.Chr4.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.PR.MGEScan.awk.script1.txt
cat protein_match.Chr4.PR.MGEScan.list.txt|awk -f protein_match.Chr4.PR.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.PR.MGEScan.awk.script2.txt
cat protein_match.Chr4.PR.gff.txt|awk -f protein_match.Chr4.PR.MGEScan.awk.script2.txt >protein_match.Chr4.PR.MGEScan.min.gff.txt

cat protein_match.Chr4.PR.gff.txt|grep 'LTRdigest' >protein_match.Chr4.PR.LTRdigest.gff.txt
cat protein_match.Chr4.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.PR.LTRdigest.list.txt
cat protein_match.Chr4.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.PR.LTRdigest.awk.script1.txt
cat protein_match.Chr4.PR.LTRdigest.list.txt|awk -f protein_match.Chr4.PR.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.PR.LTRdigest.awk.script2.txt
cat protein_match.Chr4.PR.gff.txt|awk -f protein_match.Chr4.PR.LTRdigest.awk.script2.txt >protein_match.Chr4.PR.LTRdigest.min.gff.txt

cat protein_match.Chr4.PR.MGEScan.min.gff.txt protein_match.Chr4.PR.LTRdigest.min.gff.txt >../protein_match.Chr4.PR.min.gff.txt

cat protein_match.Chr5.PR.gff.txt|grep 'MGEScan' >protein_match.Chr5.PR.MGEScan.gff.txt
cat protein_match.Chr5.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.PR.MGEScan.list.txt
cat protein_match.Chr5.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.PR.MGEScan.awk.script1.txt
cat protein_match.Chr5.PR.MGEScan.list.txt|awk -f protein_match.Chr5.PR.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.PR.MGEScan.awk.script2.txt
cat protein_match.Chr5.PR.gff.txt|awk -f protein_match.Chr5.PR.MGEScan.awk.script2.txt >protein_match.Chr5.PR.MGEScan.min.gff.txt

cat protein_match.Chr5.PR.gff.txt|grep 'LTRdigest' >protein_match.Chr5.PR.LTRdigest.gff.txt
cat protein_match.Chr5.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.PR.LTRdigest.list.txt
cat protein_match.Chr5.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.PR.LTRdigest.awk.script1.txt
cat protein_match.Chr5.PR.LTRdigest.list.txt|awk -f protein_match.Chr5.PR.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.PR.LTRdigest.awk.script2.txt
cat protein_match.Chr5.PR.gff.txt|awk -f protein_match.Chr5.PR.LTRdigest.awk.script2.txt >protein_match.Chr5.PR.LTRdigest.min.gff.txt

cat protein_match.Chr5.PR.MGEScan.min.gff.txt protein_match.Chr5.PR.LTRdigest.min.gff.txt >../protein_match.Chr5.PR.min.gff.txt

cat protein_match.ChrX.PR.gff.txt|grep 'MGEScan' >protein_match.ChrX.PR.MGEScan.gff.txt
cat protein_match.ChrX.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.PR.MGEScan.list.txt
cat protein_match.ChrX.PR.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.PR.MGEScan.awk.script1.txt
cat protein_match.ChrX.PR.MGEScan.list.txt|awk -f protein_match.ChrX.PR.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.PR.MGEScan.awk.script2.txt
cat protein_match.ChrX.PR.gff.txt|awk -f protein_match.ChrX.PR.MGEScan.awk.script2.txt >protein_match.ChrX.PR.MGEScan.min.gff.txt

cat protein_match.ChrX.PR.gff.txt|grep 'LTRdigest' >protein_match.ChrX.PR.LTRdigest.gff.txt
cat protein_match.ChrX.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.PR.LTRdigest.list.txt
cat protein_match.ChrX.PR.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.PR.LTRdigest.awk.script1.txt
cat protein_match.ChrX.PR.LTRdigest.list.txt|awk -f protein_match.ChrX.PR.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.PR.LTRdigest.awk.script2.txt
cat protein_match.ChrX.PR.gff.txt|awk -f protein_match.ChrX.PR.LTRdigest.awk.script2.txt >protein_match.ChrX.PR.LTRdigest.min.gff.txt

cat protein_match.ChrX.PR.MGEScan.min.gff.txt protein_match.ChrX.PR.LTRdigest.min.gff.txt >../protein_match.ChrX.PR.min.gff.txt

cd ../Peptidase_A17
cat protein_match.Chr1.Peptidase_A17.gff.txt|grep 'MGEScan' >protein_match.Chr1.Peptidase_A17.MGEScan.gff.txt
cat protein_match.Chr1.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.Peptidase_A17.MGEScan.list.txt
cat protein_match.Chr1.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.Peptidase_A17.MGEScan.awk.script1.txt
cat protein_match.Chr1.Peptidase_A17.MGEScan.list.txt|awk -f protein_match.Chr1.Peptidase_A17.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.Peptidase_A17.MGEScan.awk.script2.txt
cat protein_match.Chr1.Peptidase_A17.gff.txt|awk -f protein_match.Chr1.Peptidase_A17.MGEScan.awk.script2.txt >protein_match.Chr1.Peptidase_A17.MGEScan.min.gff.txt

cat protein_match.Chr1.Peptidase_A17.gff.txt|grep 'LTRdigest' >protein_match.Chr1.Peptidase_A17.LTRdigest.gff.txt
cat protein_match.Chr1.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.Peptidase_A17.LTRdigest.list.txt
cat protein_match.Chr1.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.Peptidase_A17.LTRdigest.awk.script1.txt
cat protein_match.Chr1.Peptidase_A17.LTRdigest.list.txt|awk -f protein_match.Chr1.Peptidase_A17.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.Peptidase_A17.LTRdigest.awk.script2.txt
cat protein_match.Chr1.Peptidase_A17.gff.txt|awk -f protein_match.Chr1.Peptidase_A17.LTRdigest.awk.script2.txt >protein_match.Chr1.Peptidase_A17.LTRdigest.min.gff.txt

cat protein_match.Chr1.Peptidase_A17.MGEScan.min.gff.txt protein_match.Chr1.Peptidase_A17.LTRdigest.min.gff.txt >../protein_match.Chr1.Peptidase_A17.min.gff.txt

cat protein_match.Chr2.Peptidase_A17.gff.txt|grep 'MGEScan' >protein_match.Chr2.Peptidase_A17.MGEScan.gff.txt
cat protein_match.Chr2.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.Peptidase_A17.MGEScan.list.txt
cat protein_match.Chr2.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.Peptidase_A17.MGEScan.awk.script1.txt
cat protein_match.Chr2.Peptidase_A17.MGEScan.list.txt|awk -f protein_match.Chr2.Peptidase_A17.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.Peptidase_A17.MGEScan.awk.script2.txt
cat protein_match.Chr2.Peptidase_A17.gff.txt|awk -f protein_match.Chr2.Peptidase_A17.MGEScan.awk.script2.txt >protein_match.Chr2.Peptidase_A17.MGEScan.min.gff.txt

cat protein_match.Chr2.Peptidase_A17.gff.txt|grep 'LTRdigest' >protein_match.Chr2.Peptidase_A17.LTRdigest.gff.txt
cat protein_match.Chr2.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.Peptidase_A17.LTRdigest.list.txt
cat protein_match.Chr2.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.Peptidase_A17.LTRdigest.awk.script1.txt
cat protein_match.Chr2.Peptidase_A17.LTRdigest.list.txt|awk -f protein_match.Chr2.Peptidase_A17.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.Peptidase_A17.LTRdigest.awk.script2.txt
cat protein_match.Chr2.Peptidase_A17.gff.txt|awk -f protein_match.Chr2.Peptidase_A17.LTRdigest.awk.script2.txt >protein_match.Chr2.Peptidase_A17.LTRdigest.min.gff.txt

cat protein_match.Chr2.Peptidase_A17.MGEScan.min.gff.txt protein_match.Chr2.Peptidase_A17.LTRdigest.min.gff.txt >../protein_match.Chr2.Peptidase_A17.min.gff.txt

cat protein_match.Chr3.Peptidase_A17.gff.txt|grep 'MGEScan' >protein_match.Chr3.Peptidase_A17.MGEScan.gff.txt
cat protein_match.Chr3.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.Peptidase_A17.MGEScan.list.txt
cat protein_match.Chr3.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.Peptidase_A17.MGEScan.awk.script1.txt
cat protein_match.Chr3.Peptidase_A17.MGEScan.list.txt|awk -f protein_match.Chr3.Peptidase_A17.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.Peptidase_A17.MGEScan.awk.script2.txt
cat protein_match.Chr3.Peptidase_A17.gff.txt|awk -f protein_match.Chr3.Peptidase_A17.MGEScan.awk.script2.txt >protein_match.Chr3.Peptidase_A17.MGEScan.min.gff.txt

cat protein_match.Chr3.Peptidase_A17.gff.txt|grep 'LTRdigest' >protein_match.Chr3.Peptidase_A17.LTRdigest.gff.txt
cat protein_match.Chr3.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.Peptidase_A17.LTRdigest.list.txt
cat protein_match.Chr3.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.Peptidase_A17.LTRdigest.awk.script1.txt
cat protein_match.Chr3.Peptidase_A17.LTRdigest.list.txt|awk -f protein_match.Chr3.Peptidase_A17.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.Peptidase_A17.LTRdigest.awk.script2.txt
cat protein_match.Chr3.Peptidase_A17.gff.txt|awk -f protein_match.Chr3.Peptidase_A17.LTRdigest.awk.script2.txt >protein_match.Chr3.Peptidase_A17.LTRdigest.min.gff.txt

cat protein_match.Chr3.Peptidase_A17.MGEScan.min.gff.txt protein_match.Chr3.Peptidase_A17.LTRdigest.min.gff.txt >../protein_match.Chr3.Peptidase_A17.min.gff.txt

cat protein_match.Chr4.Peptidase_A17.gff.txt|grep 'MGEScan' >protein_match.Chr4.Peptidase_A17.MGEScan.gff.txt
cat protein_match.Chr4.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.Peptidase_A17.MGEScan.list.txt
cat protein_match.Chr4.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.Peptidase_A17.MGEScan.awk.script1.txt
cat protein_match.Chr4.Peptidase_A17.MGEScan.list.txt|awk -f protein_match.Chr4.Peptidase_A17.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.Peptidase_A17.MGEScan.awk.script2.txt
cat protein_match.Chr4.Peptidase_A17.gff.txt|awk -f protein_match.Chr4.Peptidase_A17.MGEScan.awk.script2.txt >protein_match.Chr4.Peptidase_A17.MGEScan.min.gff.txt

cat protein_match.Chr4.Peptidase_A17.gff.txt|grep 'LTRdigest' >protein_match.Chr4.Peptidase_A17.LTRdigest.gff.txt
cat protein_match.Chr4.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.Peptidase_A17.LTRdigest.list.txt
cat protein_match.Chr4.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.Peptidase_A17.LTRdigest.awk.script1.txt
cat protein_match.Chr4.Peptidase_A17.LTRdigest.list.txt|awk -f protein_match.Chr4.Peptidase_A17.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.Peptidase_A17.LTRdigest.awk.script2.txt
cat protein_match.Chr4.Peptidase_A17.gff.txt|awk -f protein_match.Chr4.Peptidase_A17.LTRdigest.awk.script2.txt >protein_match.Chr4.Peptidase_A17.LTRdigest.min.gff.txt

cat protein_match.Chr4.Peptidase_A17.MGEScan.min.gff.txt protein_match.Chr4.Peptidase_A17.LTRdigest.min.gff.txt >../protein_match.Chr4.Peptidase_A17.min.gff.txt

cat protein_match.Chr5.Peptidase_A17.gff.txt|grep 'MGEScan' >protein_match.Chr5.Peptidase_A17.MGEScan.gff.txt
cat protein_match.Chr5.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.Peptidase_A17.MGEScan.list.txt
cat protein_match.Chr5.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.Peptidase_A17.MGEScan.awk.script1.txt
cat protein_match.Chr5.Peptidase_A17.MGEScan.list.txt|awk -f protein_match.Chr5.Peptidase_A17.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.Peptidase_A17.MGEScan.awk.script2.txt
cat protein_match.Chr5.Peptidase_A17.gff.txt|awk -f protein_match.Chr5.Peptidase_A17.MGEScan.awk.script2.txt >protein_match.Chr5.Peptidase_A17.MGEScan.min.gff.txt

cat protein_match.Chr5.Peptidase_A17.gff.txt|grep 'LTRdigest' >protein_match.Chr5.Peptidase_A17.LTRdigest.gff.txt
cat protein_match.Chr5.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.Peptidase_A17.LTRdigest.list.txt
cat protein_match.Chr5.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.Peptidase_A17.LTRdigest.awk.script1.txt
cat protein_match.Chr5.Peptidase_A17.LTRdigest.list.txt|awk -f protein_match.Chr5.Peptidase_A17.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.Peptidase_A17.LTRdigest.awk.script2.txt
cat protein_match.Chr5.Peptidase_A17.gff.txt|awk -f protein_match.Chr5.Peptidase_A17.LTRdigest.awk.script2.txt >protein_match.Chr5.Peptidase_A17.LTRdigest.min.gff.txt

cat protein_match.Chr5.Peptidase_A17.MGEScan.min.gff.txt protein_match.Chr5.Peptidase_A17.LTRdigest.min.gff.txt >../protein_match.Chr5.Peptidase_A17.min.gff.txt

cat protein_match.ChrX.Peptidase_A17.gff.txt|grep 'MGEScan' >protein_match.ChrX.Peptidase_A17.MGEScan.gff.txt
cat protein_match.ChrX.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.Peptidase_A17.MGEScan.list.txt
cat protein_match.ChrX.Peptidase_A17.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.Peptidase_A17.MGEScan.awk.script1.txt
cat protein_match.ChrX.Peptidase_A17.MGEScan.list.txt|awk -f protein_match.ChrX.Peptidase_A17.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.Peptidase_A17.MGEScan.awk.script2.txt
cat protein_match.ChrX.Peptidase_A17.gff.txt|awk -f protein_match.ChrX.Peptidase_A17.MGEScan.awk.script2.txt >protein_match.ChrX.Peptidase_A17.MGEScan.min.gff.txt

cat protein_match.ChrX.Peptidase_A17.gff.txt|grep 'LTRdigest' >protein_match.ChrX.Peptidase_A17.LTRdigest.gff.txt
cat protein_match.ChrX.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.Peptidase_A17.LTRdigest.list.txt
cat protein_match.ChrX.Peptidase_A17.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.Peptidase_A17.LTRdigest.awk.script1.txt
cat protein_match.ChrX.Peptidase_A17.LTRdigest.list.txt|awk -f protein_match.ChrX.Peptidase_A17.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.Peptidase_A17.LTRdigest.awk.script2.txt
cat protein_match.ChrX.Peptidase_A17.gff.txt|awk -f protein_match.ChrX.Peptidase_A17.LTRdigest.awk.script2.txt >protein_match.ChrX.Peptidase_A17.LTRdigest.min.gff.txt

cat protein_match.ChrX.Peptidase_A17.MGEScan.min.gff.txt protein_match.ChrX.Peptidase_A17.LTRdigest.min.gff.txt >../protein_match.ChrX.Peptidase_A17.min.gff.txt

cd ../RH
cat protein_match.Chr1.RH.gff.txt|grep 'MGEScan' >protein_match.Chr1.RH.MGEScan.gff.txt
cat protein_match.Chr1.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.RH.MGEScan.list.txt
cat protein_match.Chr1.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.RH.MGEScan.awk.script1.txt
cat protein_match.Chr1.RH.MGEScan.list.txt|awk -f protein_match.Chr1.RH.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.RH.MGEScan.awk.script2.txt
cat protein_match.Chr1.RH.gff.txt|awk -f protein_match.Chr1.RH.MGEScan.awk.script2.txt >protein_match.Chr1.RH.MGEScan.min.gff.txt

cat protein_match.Chr1.RH.gff.txt|grep 'LTRdigest' >protein_match.Chr1.RH.LTRdigest.gff.txt
cat protein_match.Chr1.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.RH.LTRdigest.list.txt
cat protein_match.Chr1.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.RH.LTRdigest.awk.script1.txt
cat protein_match.Chr1.RH.LTRdigest.list.txt|awk -f protein_match.Chr1.RH.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.RH.LTRdigest.awk.script2.txt
cat protein_match.Chr1.RH.gff.txt|awk -f protein_match.Chr1.RH.LTRdigest.awk.script2.txt >protein_match.Chr1.RH.LTRdigest.min.gff.txt

cat protein_match.Chr1.RH.MGEScan.min.gff.txt protein_match.Chr1.RH.LTRdigest.min.gff.txt >../protein_match.Chr1.RH.min.gff.txt

cat protein_match.Chr2.RH.gff.txt|grep 'MGEScan' >protein_match.Chr2.RH.MGEScan.gff.txt
cat protein_match.Chr2.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.RH.MGEScan.list.txt
cat protein_match.Chr2.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.RH.MGEScan.awk.script1.txt
cat protein_match.Chr2.RH.MGEScan.list.txt|awk -f protein_match.Chr2.RH.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.RH.MGEScan.awk.script2.txt
cat protein_match.Chr2.RH.gff.txt|awk -f protein_match.Chr2.RH.MGEScan.awk.script2.txt >protein_match.Chr2.RH.MGEScan.min.gff.txt

cat protein_match.Chr2.RH.gff.txt|grep 'LTRdigest' >protein_match.Chr2.RH.LTRdigest.gff.txt
cat protein_match.Chr2.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.RH.LTRdigest.list.txt
cat protein_match.Chr2.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.RH.LTRdigest.awk.script1.txt
cat protein_match.Chr2.RH.LTRdigest.list.txt|awk -f protein_match.Chr2.RH.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.RH.LTRdigest.awk.script2.txt
cat protein_match.Chr2.RH.gff.txt|awk -f protein_match.Chr2.RH.LTRdigest.awk.script2.txt >protein_match.Chr2.RH.LTRdigest.min.gff.txt

cat protein_match.Chr2.RH.MGEScan.min.gff.txt protein_match.Chr2.RH.LTRdigest.min.gff.txt >../protein_match.Chr2.RH.min.gff.txt

cat protein_match.Chr3.RH.gff.txt|grep 'MGEScan' >protein_match.Chr3.RH.MGEScan.gff.txt
cat protein_match.Chr3.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.RH.MGEScan.list.txt
cat protein_match.Chr3.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.RH.MGEScan.awk.script1.txt
cat protein_match.Chr3.RH.MGEScan.list.txt|awk -f protein_match.Chr3.RH.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.RH.MGEScan.awk.script2.txt
cat protein_match.Chr3.RH.gff.txt|awk -f protein_match.Chr3.RH.MGEScan.awk.script2.txt >protein_match.Chr3.RH.MGEScan.min.gff.txt

cat protein_match.Chr3.RH.gff.txt|grep 'LTRdigest' >protein_match.Chr3.RH.LTRdigest.gff.txt
cat protein_match.Chr3.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.RH.LTRdigest.list.txt
cat protein_match.Chr3.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.RH.LTRdigest.awk.script1.txt
cat protein_match.Chr3.RH.LTRdigest.list.txt|awk -f protein_match.Chr3.RH.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.RH.LTRdigest.awk.script2.txt
cat protein_match.Chr3.RH.gff.txt|awk -f protein_match.Chr3.RH.LTRdigest.awk.script2.txt >protein_match.Chr3.RH.LTRdigest.min.gff.txt

cat protein_match.Chr3.RH.MGEScan.min.gff.txt protein_match.Chr3.RH.LTRdigest.min.gff.txt >../protein_match.Chr3.RH.min.gff.txt

cat protein_match.Chr4.RH.gff.txt|grep 'MGEScan' >protein_match.Chr4.RH.MGEScan.gff.txt
cat protein_match.Chr4.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.RH.MGEScan.list.txt
cat protein_match.Chr4.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.RH.MGEScan.awk.script1.txt
cat protein_match.Chr4.RH.MGEScan.list.txt|awk -f protein_match.Chr4.RH.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.RH.MGEScan.awk.script2.txt
cat protein_match.Chr4.RH.gff.txt|awk -f protein_match.Chr4.RH.MGEScan.awk.script2.txt >protein_match.Chr4.RH.MGEScan.min.gff.txt

cat protein_match.Chr4.RH.gff.txt|grep 'LTRdigest' >protein_match.Chr4.RH.LTRdigest.gff.txt
cat protein_match.Chr4.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.RH.LTRdigest.list.txt
cat protein_match.Chr4.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.RH.LTRdigest.awk.script1.txt
cat protein_match.Chr4.RH.LTRdigest.list.txt|awk -f protein_match.Chr4.RH.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.RH.LTRdigest.awk.script2.txt
cat protein_match.Chr4.RH.gff.txt|awk -f protein_match.Chr4.RH.LTRdigest.awk.script2.txt >protein_match.Chr4.RH.LTRdigest.min.gff.txt

cat protein_match.Chr4.RH.MGEScan.min.gff.txt protein_match.Chr4.RH.LTRdigest.min.gff.txt >../protein_match.Chr4.RH.min.gff.txt

cat protein_match.Chr5.RH.gff.txt|grep 'MGEScan' >protein_match.Chr5.RH.MGEScan.gff.txt
cat protein_match.Chr5.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.RH.MGEScan.list.txt
cat protein_match.Chr5.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.RH.MGEScan.awk.script1.txt
cat protein_match.Chr5.RH.MGEScan.list.txt|awk -f protein_match.Chr5.RH.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.RH.MGEScan.awk.script2.txt
cat protein_match.Chr5.RH.gff.txt|awk -f protein_match.Chr5.RH.MGEScan.awk.script2.txt >protein_match.Chr5.RH.MGEScan.min.gff.txt

cat protein_match.Chr5.RH.gff.txt|grep 'LTRdigest' >protein_match.Chr5.RH.LTRdigest.gff.txt
cat protein_match.Chr5.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.RH.LTRdigest.list.txt
cat protein_match.Chr5.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.RH.LTRdigest.awk.script1.txt
cat protein_match.Chr5.RH.LTRdigest.list.txt|awk -f protein_match.Chr5.RH.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.RH.LTRdigest.awk.script2.txt
cat protein_match.Chr5.RH.gff.txt|awk -f protein_match.Chr5.RH.LTRdigest.awk.script2.txt >protein_match.Chr5.RH.LTRdigest.min.gff.txt

cat protein_match.Chr5.RH.MGEScan.min.gff.txt protein_match.Chr5.RH.LTRdigest.min.gff.txt >../protein_match.Chr5.RH.min.gff.txt

cat protein_match.ChrX.RH.gff.txt|grep 'MGEScan' >protein_match.ChrX.RH.MGEScan.gff.txt
cat protein_match.ChrX.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.RH.MGEScan.list.txt
cat protein_match.ChrX.RH.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.RH.MGEScan.awk.script1.txt
cat protein_match.ChrX.RH.MGEScan.list.txt|awk -f protein_match.ChrX.RH.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.RH.MGEScan.awk.script2.txt
cat protein_match.ChrX.RH.gff.txt|awk -f protein_match.ChrX.RH.MGEScan.awk.script2.txt >protein_match.ChrX.RH.MGEScan.min.gff.txt

cat protein_match.ChrX.RH.gff.txt|grep 'LTRdigest' >protein_match.ChrX.RH.LTRdigest.gff.txt
cat protein_match.ChrX.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.RH.LTRdigest.list.txt
cat protein_match.ChrX.RH.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.RH.LTRdigest.awk.script1.txt
cat protein_match.ChrX.RH.LTRdigest.list.txt|awk -f protein_match.ChrX.RH.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.RH.LTRdigest.awk.script2.txt
cat protein_match.ChrX.RH.gff.txt|awk -f protein_match.ChrX.RH.LTRdigest.awk.script2.txt >protein_match.ChrX.RH.LTRdigest.min.gff.txt

cat protein_match.ChrX.RH.MGEScan.min.gff.txt protein_match.ChrX.RH.LTRdigest.min.gff.txt >../protein_match.ChrX.RH.min.gff.txt

cd ../RNase
cat protein_match.Chr1.RNase.gff.txt|grep 'MGEScan' >protein_match.Chr1.RNase.MGEScan.gff.txt
cat protein_match.Chr1.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.RNase.MGEScan.list.txt
cat protein_match.Chr1.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.RNase.MGEScan.awk.script1.txt
cat protein_match.Chr1.RNase.MGEScan.list.txt|awk -f protein_match.Chr1.RNase.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.RNase.MGEScan.awk.script2.txt
cat protein_match.Chr1.RNase.gff.txt|awk -f protein_match.Chr1.RNase.MGEScan.awk.script2.txt >protein_match.Chr1.RNase.MGEScan.min.gff.txt

cat protein_match.Chr1.RNase.gff.txt|grep 'LTRdigest' >protein_match.Chr1.RNase.LTRdigest.gff.txt
cat protein_match.Chr1.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.RNase.LTRdigest.list.txt
cat protein_match.Chr1.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.RNase.LTRdigest.awk.script1.txt
cat protein_match.Chr1.RNase.LTRdigest.list.txt|awk -f protein_match.Chr1.RNase.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.RNase.LTRdigest.awk.script2.txt
cat protein_match.Chr1.RNase.gff.txt|awk -f protein_match.Chr1.RNase.LTRdigest.awk.script2.txt >protein_match.Chr1.RNase.LTRdigest.min.gff.txt

cat protein_match.Chr1.RNase.MGEScan.min.gff.txt protein_match.Chr1.RNase.LTRdigest.min.gff.txt >../protein_match.Chr1.RNase.min.gff.txt

cat protein_match.Chr2.RNase.gff.txt|grep 'MGEScan' >protein_match.Chr2.RNase.MGEScan.gff.txt
cat protein_match.Chr2.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.RNase.MGEScan.list.txt
cat protein_match.Chr2.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.RNase.MGEScan.awk.script1.txt
cat protein_match.Chr2.RNase.MGEScan.list.txt|awk -f protein_match.Chr2.RNase.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.RNase.MGEScan.awk.script2.txt
cat protein_match.Chr2.RNase.gff.txt|awk -f protein_match.Chr2.RNase.MGEScan.awk.script2.txt >protein_match.Chr2.RNase.MGEScan.min.gff.txt

cat protein_match.Chr2.RNase.gff.txt|grep 'LTRdigest' >protein_match.Chr2.RNase.LTRdigest.gff.txt
cat protein_match.Chr2.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.RNase.LTRdigest.list.txt
cat protein_match.Chr2.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.RNase.LTRdigest.awk.script1.txt
cat protein_match.Chr2.RNase.LTRdigest.list.txt|awk -f protein_match.Chr2.RNase.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.RNase.LTRdigest.awk.script2.txt
cat protein_match.Chr2.RNase.gff.txt|awk -f protein_match.Chr2.RNase.LTRdigest.awk.script2.txt >protein_match.Chr2.RNase.LTRdigest.min.gff.txt

cat protein_match.Chr2.RNase.MGEScan.min.gff.txt protein_match.Chr2.RNase.LTRdigest.min.gff.txt >../protein_match.Chr2.RNase.min.gff.txt

cat protein_match.Chr3.RNase.gff.txt|grep 'MGEScan' >protein_match.Chr3.RNase.MGEScan.gff.txt
cat protein_match.Chr3.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.RNase.MGEScan.list.txt
cat protein_match.Chr3.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.RNase.MGEScan.awk.script1.txt
cat protein_match.Chr3.RNase.MGEScan.list.txt|awk -f protein_match.Chr3.RNase.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.RNase.MGEScan.awk.script2.txt
cat protein_match.Chr3.RNase.gff.txt|awk -f protein_match.Chr3.RNase.MGEScan.awk.script2.txt >protein_match.Chr3.RNase.MGEScan.min.gff.txt

cat protein_match.Chr3.RNase.gff.txt|grep 'LTRdigest' >protein_match.Chr3.RNase.LTRdigest.gff.txt
cat protein_match.Chr3.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.RNase.LTRdigest.list.txt
cat protein_match.Chr3.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.RNase.LTRdigest.awk.script1.txt
cat protein_match.Chr3.RNase.LTRdigest.list.txt|awk -f protein_match.Chr3.RNase.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.RNase.LTRdigest.awk.script2.txt
cat protein_match.Chr3.RNase.gff.txt|awk -f protein_match.Chr3.RNase.LTRdigest.awk.script2.txt >protein_match.Chr3.RNase.LTRdigest.min.gff.txt

cat protein_match.Chr3.RNase.MGEScan.min.gff.txt protein_match.Chr3.RNase.LTRdigest.min.gff.txt >../protein_match.Chr3.RNase.min.gff.txt

cat protein_match.Chr4.RNase.gff.txt|grep 'MGEScan' >protein_match.Chr4.RNase.MGEScan.gff.txt
cat protein_match.Chr4.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.RNase.MGEScan.list.txt
cat protein_match.Chr4.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.RNase.MGEScan.awk.script1.txt
cat protein_match.Chr4.RNase.MGEScan.list.txt|awk -f protein_match.Chr4.RNase.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.RNase.MGEScan.awk.script2.txt
cat protein_match.Chr4.RNase.gff.txt|awk -f protein_match.Chr4.RNase.MGEScan.awk.script2.txt >protein_match.Chr4.RNase.MGEScan.min.gff.txt

cat protein_match.Chr4.RNase.gff.txt|grep 'LTRdigest' >protein_match.Chr4.RNase.LTRdigest.gff.txt
cat protein_match.Chr4.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.RNase.LTRdigest.list.txt
cat protein_match.Chr4.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.RNase.LTRdigest.awk.script1.txt
cat protein_match.Chr4.RNase.LTRdigest.list.txt|awk -f protein_match.Chr4.RNase.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.RNase.LTRdigest.awk.script2.txt
cat protein_match.Chr4.RNase.gff.txt|awk -f protein_match.Chr4.RNase.LTRdigest.awk.script2.txt >protein_match.Chr4.RNase.LTRdigest.min.gff.txt

cat protein_match.Chr4.RNase.MGEScan.min.gff.txt protein_match.Chr4.RNase.LTRdigest.min.gff.txt >../protein_match.Chr4.RNase.min.gff.txt

cat protein_match.Chr5.RNase.gff.txt|grep 'MGEScan' >protein_match.Chr5.RNase.MGEScan.gff.txt
cat protein_match.Chr5.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.RNase.MGEScan.list.txt
cat protein_match.Chr5.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.RNase.MGEScan.awk.script1.txt
cat protein_match.Chr5.RNase.MGEScan.list.txt|awk -f protein_match.Chr5.RNase.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.RNase.MGEScan.awk.script2.txt
cat protein_match.Chr5.RNase.gff.txt|awk -f protein_match.Chr5.RNase.MGEScan.awk.script2.txt >protein_match.Chr5.RNase.MGEScan.min.gff.txt

cat protein_match.Chr5.RNase.gff.txt|grep 'LTRdigest' >protein_match.Chr5.RNase.LTRdigest.gff.txt
cat protein_match.Chr5.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.RNase.LTRdigest.list.txt
cat protein_match.Chr5.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.RNase.LTRdigest.awk.script1.txt
cat protein_match.Chr5.RNase.LTRdigest.list.txt|awk -f protein_match.Chr5.RNase.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.RNase.LTRdigest.awk.script2.txt
cat protein_match.Chr5.RNase.gff.txt|awk -f protein_match.Chr5.RNase.LTRdigest.awk.script2.txt >protein_match.Chr5.RNase.LTRdigest.min.gff.txt

cat protein_match.Chr5.RNase.MGEScan.min.gff.txt protein_match.Chr5.RNase.LTRdigest.min.gff.txt >../protein_match.Chr5.RNase.min.gff.txt

cat protein_match.ChrX.RNase.gff.txt|grep 'MGEScan' >protein_match.ChrX.RNase.MGEScan.gff.txt
cat protein_match.ChrX.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.RNase.MGEScan.list.txt
cat protein_match.ChrX.RNase.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.RNase.MGEScan.awk.script1.txt
cat protein_match.ChrX.RNase.MGEScan.list.txt|awk -f protein_match.ChrX.RNase.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.RNase.MGEScan.awk.script2.txt
cat protein_match.ChrX.RNase.gff.txt|awk -f protein_match.ChrX.RNase.MGEScan.awk.script2.txt >protein_match.ChrX.RNase.MGEScan.min.gff.txt

cat protein_match.ChrX.RNase.gff.txt|grep 'LTRdigest' >protein_match.ChrX.RNase.LTRdigest.gff.txt
cat protein_match.ChrX.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.RNase.LTRdigest.list.txt
cat protein_match.ChrX.RNase.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.RNase.LTRdigest.awk.script1.txt
cat protein_match.ChrX.RNase.LTRdigest.list.txt|awk -f protein_match.ChrX.RNase.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.RNase.LTRdigest.awk.script2.txt
cat protein_match.ChrX.RNase.gff.txt|awk -f protein_match.ChrX.RNase.LTRdigest.awk.script2.txt >protein_match.ChrX.RNase.LTRdigest.min.gff.txt

cat protein_match.ChrX.RNase.MGEScan.min.gff.txt protein_match.ChrX.RNase.LTRdigest.min.gff.txt >../protein_match.ChrX.RNase.min.gff.txt

cd ../RT
cat protein_match.Chr1.RT.gff.txt|grep 'MGEScan' >protein_match.Chr1.RT.MGEScan.gff.txt
cat protein_match.Chr1.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.RT.MGEScan.list.txt
cat protein_match.Chr1.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.RT.MGEScan.awk.script1.txt
cat protein_match.Chr1.RT.MGEScan.list.txt|awk -f protein_match.Chr1.RT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.RT.MGEScan.awk.script2.txt
cat protein_match.Chr1.RT.gff.txt|awk -f protein_match.Chr1.RT.MGEScan.awk.script2.txt >protein_match.Chr1.RT.MGEScan.min.gff.txt

cat protein_match.Chr1.RT.gff.txt|grep 'LTRdigest' >protein_match.Chr1.RT.LTRdigest.gff.txt
cat protein_match.Chr1.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.RT.LTRdigest.list.txt
cat protein_match.Chr1.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.RT.LTRdigest.awk.script1.txt
cat protein_match.Chr1.RT.LTRdigest.list.txt|awk -f protein_match.Chr1.RT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.RT.LTRdigest.awk.script2.txt
cat protein_match.Chr1.RT.gff.txt|awk -f protein_match.Chr1.RT.LTRdigest.awk.script2.txt >protein_match.Chr1.RT.LTRdigest.min.gff.txt

cat protein_match.Chr1.RT.MGEScan.min.gff.txt protein_match.Chr1.RT.LTRdigest.min.gff.txt >../protein_match.Chr1.RT.min.gff.txt

cat protein_match.Chr2.RT.gff.txt|grep 'MGEScan' >protein_match.Chr2.RT.MGEScan.gff.txt
cat protein_match.Chr2.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.RT.MGEScan.list.txt
cat protein_match.Chr2.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.RT.MGEScan.awk.script1.txt
cat protein_match.Chr2.RT.MGEScan.list.txt|awk -f protein_match.Chr2.RT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.RT.MGEScan.awk.script2.txt
cat protein_match.Chr2.RT.gff.txt|awk -f protein_match.Chr2.RT.MGEScan.awk.script2.txt >protein_match.Chr2.RT.MGEScan.min.gff.txt

cat protein_match.Chr2.RT.gff.txt|grep 'LTRdigest' >protein_match.Chr2.RT.LTRdigest.gff.txt
cat protein_match.Chr2.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.RT.LTRdigest.list.txt
cat protein_match.Chr2.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.RT.LTRdigest.awk.script1.txt
cat protein_match.Chr2.RT.LTRdigest.list.txt|awk -f protein_match.Chr2.RT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.RT.LTRdigest.awk.script2.txt
cat protein_match.Chr2.RT.gff.txt|awk -f protein_match.Chr2.RT.LTRdigest.awk.script2.txt >protein_match.Chr2.RT.LTRdigest.min.gff.txt

cat protein_match.Chr2.RT.MGEScan.min.gff.txt protein_match.Chr2.RT.LTRdigest.min.gff.txt >../protein_match.Chr2.RT.min.gff.txt

cat protein_match.Chr3.RT.gff.txt|grep 'MGEScan' >protein_match.Chr3.RT.MGEScan.gff.txt
cat protein_match.Chr3.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.RT.MGEScan.list.txt
cat protein_match.Chr3.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.RT.MGEScan.awk.script1.txt
cat protein_match.Chr3.RT.MGEScan.list.txt|awk -f protein_match.Chr3.RT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.RT.MGEScan.awk.script2.txt
cat protein_match.Chr3.RT.gff.txt|awk -f protein_match.Chr3.RT.MGEScan.awk.script2.txt >protein_match.Chr3.RT.MGEScan.min.gff.txt

cat protein_match.Chr3.RT.gff.txt|grep 'LTRdigest' >protein_match.Chr3.RT.LTRdigest.gff.txt
cat protein_match.Chr3.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.RT.LTRdigest.list.txt
cat protein_match.Chr3.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.RT.LTRdigest.awk.script1.txt
cat protein_match.Chr3.RT.LTRdigest.list.txt|awk -f protein_match.Chr3.RT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.RT.LTRdigest.awk.script2.txt
cat protein_match.Chr3.RT.gff.txt|awk -f protein_match.Chr3.RT.LTRdigest.awk.script2.txt >protein_match.Chr3.RT.LTRdigest.min.gff.txt

cat protein_match.Chr3.RT.MGEScan.min.gff.txt protein_match.Chr3.RT.LTRdigest.min.gff.txt >../protein_match.Chr3.RT.min.gff.txt

cat protein_match.Chr4.RT.gff.txt|grep 'MGEScan' >protein_match.Chr4.RT.MGEScan.gff.txt
cat protein_match.Chr4.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.RT.MGEScan.list.txt
cat protein_match.Chr4.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.RT.MGEScan.awk.script1.txt
cat protein_match.Chr4.RT.MGEScan.list.txt|awk -f protein_match.Chr4.RT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.RT.MGEScan.awk.script2.txt
cat protein_match.Chr4.RT.gff.txt|awk -f protein_match.Chr4.RT.MGEScan.awk.script2.txt >protein_match.Chr4.RT.MGEScan.min.gff.txt

cat protein_match.Chr4.RT.gff.txt|grep 'LTRdigest' >protein_match.Chr4.RT.LTRdigest.gff.txt
cat protein_match.Chr4.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.RT.LTRdigest.list.txt
cat protein_match.Chr4.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.RT.LTRdigest.awk.script1.txt
cat protein_match.Chr4.RT.LTRdigest.list.txt|awk -f protein_match.Chr4.RT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.RT.LTRdigest.awk.script2.txt
cat protein_match.Chr4.RT.gff.txt|awk -f protein_match.Chr4.RT.LTRdigest.awk.script2.txt >protein_match.Chr4.RT.LTRdigest.min.gff.txt

cat protein_match.Chr4.RT.MGEScan.min.gff.txt protein_match.Chr4.RT.LTRdigest.min.gff.txt >../protein_match.Chr4.RT.min.gff.txt

cat protein_match.Chr5.RT.gff.txt|grep 'MGEScan' >protein_match.Chr5.RT.MGEScan.gff.txt
cat protein_match.Chr5.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.RT.MGEScan.list.txt
cat protein_match.Chr5.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.RT.MGEScan.awk.script1.txt
cat protein_match.Chr5.RT.MGEScan.list.txt|awk -f protein_match.Chr5.RT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.RT.MGEScan.awk.script2.txt
cat protein_match.Chr5.RT.gff.txt|awk -f protein_match.Chr5.RT.MGEScan.awk.script2.txt >protein_match.Chr5.RT.MGEScan.min.gff.txt

cat protein_match.Chr5.RT.gff.txt|grep 'LTRdigest' >protein_match.Chr5.RT.LTRdigest.gff.txt
cat protein_match.Chr5.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.RT.LTRdigest.list.txt
cat protein_match.Chr5.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.RT.LTRdigest.awk.script1.txt
cat protein_match.Chr5.RT.LTRdigest.list.txt|awk -f protein_match.Chr5.RT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.RT.LTRdigest.awk.script2.txt
cat protein_match.Chr5.RT.gff.txt|awk -f protein_match.Chr5.RT.LTRdigest.awk.script2.txt >protein_match.Chr5.RT.LTRdigest.min.gff.txt

cat protein_match.Chr5.RT.MGEScan.min.gff.txt protein_match.Chr5.RT.LTRdigest.min.gff.txt >../protein_match.Chr5.RT.min.gff.txt

cat protein_match.ChrX.RT.gff.txt|grep 'MGEScan' >protein_match.ChrX.RT.MGEScan.gff.txt
cat protein_match.ChrX.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.RT.MGEScan.list.txt
cat protein_match.ChrX.RT.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.RT.MGEScan.awk.script1.txt
cat protein_match.ChrX.RT.MGEScan.list.txt|awk -f protein_match.ChrX.RT.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.RT.MGEScan.awk.script2.txt
cat protein_match.ChrX.RT.gff.txt|awk -f protein_match.ChrX.RT.MGEScan.awk.script2.txt >protein_match.ChrX.RT.MGEScan.min.gff.txt

cat protein_match.ChrX.RT.gff.txt|grep 'LTRdigest' >protein_match.ChrX.RT.LTRdigest.gff.txt
cat protein_match.ChrX.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.RT.LTRdigest.list.txt
cat protein_match.ChrX.RT.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.RT.LTRdigest.awk.script1.txt
cat protein_match.ChrX.RT.LTRdigest.list.txt|awk -f protein_match.ChrX.RT.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.RT.LTRdigest.awk.script2.txt
cat protein_match.ChrX.RT.gff.txt|awk -f protein_match.ChrX.RT.LTRdigest.awk.script2.txt >protein_match.ChrX.RT.LTRdigest.min.gff.txt

cat protein_match.ChrX.RT.MGEScan.min.gff.txt protein_match.ChrX.RT.LTRdigest.min.gff.txt >../protein_match.ChrX.RT.min.gff.txt

cd ../Retrotrans
cat protein_match.Chr1.Retrotrans.gff.txt|grep 'MGEScan' >protein_match.Chr1.Retrotrans.MGEScan.gff.txt
cat protein_match.Chr1.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.Retrotrans.MGEScan.list.txt
cat protein_match.Chr1.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.Retrotrans.MGEScan.awk.script1.txt
cat protein_match.Chr1.Retrotrans.MGEScan.list.txt|awk -f protein_match.Chr1.Retrotrans.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.Retrotrans.MGEScan.awk.script2.txt
cat protein_match.Chr1.Retrotrans.gff.txt|awk -f protein_match.Chr1.Retrotrans.MGEScan.awk.script2.txt >protein_match.Chr1.Retrotrans.MGEScan.min.gff.txt

cat protein_match.Chr1.Retrotrans.gff.txt|grep 'LTRdigest' >protein_match.Chr1.Retrotrans.LTRdigest.gff.txt
cat protein_match.Chr1.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.Retrotrans.LTRdigest.list.txt
cat protein_match.Chr1.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.Retrotrans.LTRdigest.awk.script1.txt
cat protein_match.Chr1.Retrotrans.LTRdigest.list.txt|awk -f protein_match.Chr1.Retrotrans.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.Retrotrans.LTRdigest.awk.script2.txt
cat protein_match.Chr1.Retrotrans.gff.txt|awk -f protein_match.Chr1.Retrotrans.LTRdigest.awk.script2.txt >protein_match.Chr1.Retrotrans.LTRdigest.min.gff.txt

cat protein_match.Chr1.Retrotrans.MGEScan.min.gff.txt protein_match.Chr1.Retrotrans.LTRdigest.min.gff.txt >../protein_match.Chr1.Retrotrans.min.gff.txt

cat protein_match.Chr2.Retrotrans.gff.txt|grep 'MGEScan' >protein_match.Chr2.Retrotrans.MGEScan.gff.txt
cat protein_match.Chr2.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.Retrotrans.MGEScan.list.txt
cat protein_match.Chr2.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.Retrotrans.MGEScan.awk.script1.txt
cat protein_match.Chr2.Retrotrans.MGEScan.list.txt|awk -f protein_match.Chr2.Retrotrans.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.Retrotrans.MGEScan.awk.script2.txt
cat protein_match.Chr2.Retrotrans.gff.txt|awk -f protein_match.Chr2.Retrotrans.MGEScan.awk.script2.txt >protein_match.Chr2.Retrotrans.MGEScan.min.gff.txt

cat protein_match.Chr2.Retrotrans.gff.txt|grep 'LTRdigest' >protein_match.Chr2.Retrotrans.LTRdigest.gff.txt
cat protein_match.Chr2.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.Retrotrans.LTRdigest.list.txt
cat protein_match.Chr2.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.Retrotrans.LTRdigest.awk.script1.txt
cat protein_match.Chr2.Retrotrans.LTRdigest.list.txt|awk -f protein_match.Chr2.Retrotrans.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.Retrotrans.LTRdigest.awk.script2.txt
cat protein_match.Chr2.Retrotrans.gff.txt|awk -f protein_match.Chr2.Retrotrans.LTRdigest.awk.script2.txt >protein_match.Chr2.Retrotrans.LTRdigest.min.gff.txt

cat protein_match.Chr2.Retrotrans.MGEScan.min.gff.txt protein_match.Chr2.Retrotrans.LTRdigest.min.gff.txt >../protein_match.Chr2.Retrotrans.min.gff.txt

cat protein_match.Chr3.Retrotrans.gff.txt|grep 'MGEScan' >protein_match.Chr3.Retrotrans.MGEScan.gff.txt
cat protein_match.Chr3.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.Retrotrans.MGEScan.list.txt
cat protein_match.Chr3.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.Retrotrans.MGEScan.awk.script1.txt
cat protein_match.Chr3.Retrotrans.MGEScan.list.txt|awk -f protein_match.Chr3.Retrotrans.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.Retrotrans.MGEScan.awk.script2.txt
cat protein_match.Chr3.Retrotrans.gff.txt|awk -f protein_match.Chr3.Retrotrans.MGEScan.awk.script2.txt >protein_match.Chr3.Retrotrans.MGEScan.min.gff.txt

cat protein_match.Chr3.Retrotrans.gff.txt|grep 'LTRdigest' >protein_match.Chr3.Retrotrans.LTRdigest.gff.txt
cat protein_match.Chr3.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.Retrotrans.LTRdigest.list.txt
cat protein_match.Chr3.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.Retrotrans.LTRdigest.awk.script1.txt
cat protein_match.Chr3.Retrotrans.LTRdigest.list.txt|awk -f protein_match.Chr3.Retrotrans.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.Retrotrans.LTRdigest.awk.script2.txt
cat protein_match.Chr3.Retrotrans.gff.txt|awk -f protein_match.Chr3.Retrotrans.LTRdigest.awk.script2.txt >protein_match.Chr3.Retrotrans.LTRdigest.min.gff.txt

cat protein_match.Chr3.Retrotrans.MGEScan.min.gff.txt protein_match.Chr3.Retrotrans.LTRdigest.min.gff.txt >../protein_match.Chr3.Retrotrans.min.gff.txt

cat protein_match.Chr4.Retrotrans.gff.txt|grep 'MGEScan' >protein_match.Chr4.Retrotrans.MGEScan.gff.txt
cat protein_match.Chr4.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.Retrotrans.MGEScan.list.txt
cat protein_match.Chr4.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.Retrotrans.MGEScan.awk.script1.txt
cat protein_match.Chr4.Retrotrans.MGEScan.list.txt|awk -f protein_match.Chr4.Retrotrans.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.Retrotrans.MGEScan.awk.script2.txt
cat protein_match.Chr4.Retrotrans.gff.txt|awk -f protein_match.Chr4.Retrotrans.MGEScan.awk.script2.txt >protein_match.Chr4.Retrotrans.MGEScan.min.gff.txt

cat protein_match.Chr4.Retrotrans.gff.txt|grep 'LTRdigest' >protein_match.Chr4.Retrotrans.LTRdigest.gff.txt
cat protein_match.Chr4.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.Retrotrans.LTRdigest.list.txt
cat protein_match.Chr4.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.Retrotrans.LTRdigest.awk.script1.txt
cat protein_match.Chr4.Retrotrans.LTRdigest.list.txt|awk -f protein_match.Chr4.Retrotrans.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.Retrotrans.LTRdigest.awk.script2.txt
cat protein_match.Chr4.Retrotrans.gff.txt|awk -f protein_match.Chr4.Retrotrans.LTRdigest.awk.script2.txt >protein_match.Chr4.Retrotrans.LTRdigest.min.gff.txt

cat protein_match.Chr4.Retrotrans.MGEScan.min.gff.txt protein_match.Chr4.Retrotrans.LTRdigest.min.gff.txt >../protein_match.Chr4.Retrotrans.min.gff.txt

cat protein_match.Chr5.Retrotrans.gff.txt|grep 'MGEScan' >protein_match.Chr5.Retrotrans.MGEScan.gff.txt
cat protein_match.Chr5.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.Retrotrans.MGEScan.list.txt
cat protein_match.Chr5.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.Retrotrans.MGEScan.awk.script1.txt
cat protein_match.Chr5.Retrotrans.MGEScan.list.txt|awk -f protein_match.Chr5.Retrotrans.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.Retrotrans.MGEScan.awk.script2.txt
cat protein_match.Chr5.Retrotrans.gff.txt|awk -f protein_match.Chr5.Retrotrans.MGEScan.awk.script2.txt >protein_match.Chr5.Retrotrans.MGEScan.min.gff.txt

cat protein_match.Chr5.Retrotrans.gff.txt|grep 'LTRdigest' >protein_match.Chr5.Retrotrans.LTRdigest.gff.txt
cat protein_match.Chr5.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.Retrotrans.LTRdigest.list.txt
cat protein_match.Chr5.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.Retrotrans.LTRdigest.awk.script1.txt
cat protein_match.Chr5.Retrotrans.LTRdigest.list.txt|awk -f protein_match.Chr5.Retrotrans.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.Retrotrans.LTRdigest.awk.script2.txt
cat protein_match.Chr5.Retrotrans.gff.txt|awk -f protein_match.Chr5.Retrotrans.LTRdigest.awk.script2.txt >protein_match.Chr5.Retrotrans.LTRdigest.min.gff.txt

cat protein_match.Chr5.Retrotrans.MGEScan.min.gff.txt protein_match.Chr5.Retrotrans.LTRdigest.min.gff.txt >../protein_match.Chr5.Retrotrans.min.gff.txt

cat protein_match.ChrX.Retrotrans.gff.txt|grep 'MGEScan' >protein_match.ChrX.Retrotrans.MGEScan.gff.txt
cat protein_match.ChrX.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.Retrotrans.MGEScan.list.txt
cat protein_match.ChrX.Retrotrans.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.Retrotrans.MGEScan.awk.script1.txt
cat protein_match.ChrX.Retrotrans.MGEScan.list.txt|awk -f protein_match.ChrX.Retrotrans.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.Retrotrans.MGEScan.awk.script2.txt
cat protein_match.ChrX.Retrotrans.gff.txt|awk -f protein_match.ChrX.Retrotrans.MGEScan.awk.script2.txt >protein_match.ChrX.Retrotrans.MGEScan.min.gff.txt

cat protein_match.ChrX.Retrotrans.gff.txt|grep 'LTRdigest' >protein_match.ChrX.Retrotrans.LTRdigest.gff.txt
cat protein_match.ChrX.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.Retrotrans.LTRdigest.list.txt
cat protein_match.ChrX.Retrotrans.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.Retrotrans.LTRdigest.awk.script1.txt
cat protein_match.ChrX.Retrotrans.LTRdigest.list.txt|awk -f protein_match.ChrX.Retrotrans.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.Retrotrans.LTRdigest.awk.script2.txt
cat protein_match.ChrX.Retrotrans.gff.txt|awk -f protein_match.ChrX.Retrotrans.LTRdigest.awk.script2.txt >protein_match.ChrX.Retrotrans.LTRdigest.min.gff.txt

cat protein_match.ChrX.Retrotrans.MGEScan.min.gff.txt protein_match.ChrX.Retrotrans.LTRdigest.min.gff.txt >../protein_match.ChrX.Retrotrans.min.gff.txt

cd ../galadriel
cat protein_match.Chr1.galadriel.gff.txt|grep 'MGEScan' >protein_match.Chr1.galadriel.MGEScan.gff.txt
cat protein_match.Chr1.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.galadriel.MGEScan.list.txt
cat protein_match.Chr1.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.galadriel.MGEScan.awk.script1.txt
cat protein_match.Chr1.galadriel.MGEScan.list.txt|awk -f protein_match.Chr1.galadriel.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.galadriel.MGEScan.awk.script2.txt
cat protein_match.Chr1.galadriel.gff.txt|awk -f protein_match.Chr1.galadriel.MGEScan.awk.script2.txt >protein_match.Chr1.galadriel.MGEScan.min.gff.txt

cat protein_match.Chr1.galadriel.gff.txt|grep 'LTRdigest' >protein_match.Chr1.galadriel.LTRdigest.gff.txt
cat protein_match.Chr1.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.galadriel.LTRdigest.list.txt
cat protein_match.Chr1.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.galadriel.LTRdigest.awk.script1.txt
cat protein_match.Chr1.galadriel.LTRdigest.list.txt|awk -f protein_match.Chr1.galadriel.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.galadriel.LTRdigest.awk.script2.txt
cat protein_match.Chr1.galadriel.gff.txt|awk -f protein_match.Chr1.galadriel.LTRdigest.awk.script2.txt >protein_match.Chr1.galadriel.LTRdigest.min.gff.txt

cat protein_match.Chr1.galadriel.MGEScan.min.gff.txt protein_match.Chr1.galadriel.LTRdigest.min.gff.txt >../protein_match.Chr1.galadriel.min.gff.txt

cat protein_match.Chr2.galadriel.gff.txt|grep 'MGEScan' >protein_match.Chr2.galadriel.MGEScan.gff.txt
cat protein_match.Chr2.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.galadriel.MGEScan.list.txt
cat protein_match.Chr2.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.galadriel.MGEScan.awk.script1.txt
cat protein_match.Chr2.galadriel.MGEScan.list.txt|awk -f protein_match.Chr2.galadriel.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.galadriel.MGEScan.awk.script2.txt
cat protein_match.Chr2.galadriel.gff.txt|awk -f protein_match.Chr2.galadriel.MGEScan.awk.script2.txt >protein_match.Chr2.galadriel.MGEScan.min.gff.txt

cat protein_match.Chr2.galadriel.gff.txt|grep 'LTRdigest' >protein_match.Chr2.galadriel.LTRdigest.gff.txt
cat protein_match.Chr2.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.galadriel.LTRdigest.list.txt
cat protein_match.Chr2.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.galadriel.LTRdigest.awk.script1.txt
cat protein_match.Chr2.galadriel.LTRdigest.list.txt|awk -f protein_match.Chr2.galadriel.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.galadriel.LTRdigest.awk.script2.txt
cat protein_match.Chr2.galadriel.gff.txt|awk -f protein_match.Chr2.galadriel.LTRdigest.awk.script2.txt >protein_match.Chr2.galadriel.LTRdigest.min.gff.txt

cat protein_match.Chr2.galadriel.MGEScan.min.gff.txt protein_match.Chr2.galadriel.LTRdigest.min.gff.txt >../protein_match.Chr2.galadriel.min.gff.txt

cat protein_match.Chr3.galadriel.gff.txt|grep 'MGEScan' >protein_match.Chr3.galadriel.MGEScan.gff.txt
cat protein_match.Chr3.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.galadriel.MGEScan.list.txt
cat protein_match.Chr3.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.galadriel.MGEScan.awk.script1.txt
cat protein_match.Chr3.galadriel.MGEScan.list.txt|awk -f protein_match.Chr3.galadriel.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.galadriel.MGEScan.awk.script2.txt
cat protein_match.Chr3.galadriel.gff.txt|awk -f protein_match.Chr3.galadriel.MGEScan.awk.script2.txt >protein_match.Chr3.galadriel.MGEScan.min.gff.txt

cat protein_match.Chr3.galadriel.gff.txt|grep 'LTRdigest' >protein_match.Chr3.galadriel.LTRdigest.gff.txt
cat protein_match.Chr3.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.galadriel.LTRdigest.list.txt
cat protein_match.Chr3.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.galadriel.LTRdigest.awk.script1.txt
cat protein_match.Chr3.galadriel.LTRdigest.list.txt|awk -f protein_match.Chr3.galadriel.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.galadriel.LTRdigest.awk.script2.txt
cat protein_match.Chr3.galadriel.gff.txt|awk -f protein_match.Chr3.galadriel.LTRdigest.awk.script2.txt >protein_match.Chr3.galadriel.LTRdigest.min.gff.txt

cat protein_match.Chr3.galadriel.MGEScan.min.gff.txt protein_match.Chr3.galadriel.LTRdigest.min.gff.txt >../protein_match.Chr3.galadriel.min.gff.txt

cat protein_match.Chr4.galadriel.gff.txt|grep 'MGEScan' >protein_match.Chr4.galadriel.MGEScan.gff.txt
cat protein_match.Chr4.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.galadriel.MGEScan.list.txt
cat protein_match.Chr4.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.galadriel.MGEScan.awk.script1.txt
cat protein_match.Chr4.galadriel.MGEScan.list.txt|awk -f protein_match.Chr4.galadriel.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.galadriel.MGEScan.awk.script2.txt
cat protein_match.Chr4.galadriel.gff.txt|awk -f protein_match.Chr4.galadriel.MGEScan.awk.script2.txt >protein_match.Chr4.galadriel.MGEScan.min.gff.txt

cat protein_match.Chr4.galadriel.gff.txt|grep 'LTRdigest' >protein_match.Chr4.galadriel.LTRdigest.gff.txt
cat protein_match.Chr4.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.galadriel.LTRdigest.list.txt
cat protein_match.Chr4.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.galadriel.LTRdigest.awk.script1.txt
cat protein_match.Chr4.galadriel.LTRdigest.list.txt|awk -f protein_match.Chr4.galadriel.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.galadriel.LTRdigest.awk.script2.txt
cat protein_match.Chr4.galadriel.gff.txt|awk -f protein_match.Chr4.galadriel.LTRdigest.awk.script2.txt >protein_match.Chr4.galadriel.LTRdigest.min.gff.txt

cat protein_match.Chr4.galadriel.MGEScan.min.gff.txt protein_match.Chr4.galadriel.LTRdigest.min.gff.txt >../protein_match.Chr4.galadriel.min.gff.txt

cat protein_match.Chr5.galadriel.gff.txt|grep 'MGEScan' >protein_match.Chr5.galadriel.MGEScan.gff.txt
cat protein_match.Chr5.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.galadriel.MGEScan.list.txt
cat protein_match.Chr5.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.galadriel.MGEScan.awk.script1.txt
cat protein_match.Chr5.galadriel.MGEScan.list.txt|awk -f protein_match.Chr5.galadriel.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.galadriel.MGEScan.awk.script2.txt
cat protein_match.Chr5.galadriel.gff.txt|awk -f protein_match.Chr5.galadriel.MGEScan.awk.script2.txt >protein_match.Chr5.galadriel.MGEScan.min.gff.txt

cat protein_match.Chr5.galadriel.gff.txt|grep 'LTRdigest' >protein_match.Chr5.galadriel.LTRdigest.gff.txt
cat protein_match.Chr5.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.galadriel.LTRdigest.list.txt
cat protein_match.Chr5.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.galadriel.LTRdigest.awk.script1.txt
cat protein_match.Chr5.galadriel.LTRdigest.list.txt|awk -f protein_match.Chr5.galadriel.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.galadriel.LTRdigest.awk.script2.txt
cat protein_match.Chr5.galadriel.gff.txt|awk -f protein_match.Chr5.galadriel.LTRdigest.awk.script2.txt >protein_match.Chr5.galadriel.LTRdigest.min.gff.txt

cat protein_match.Chr5.galadriel.MGEScan.min.gff.txt protein_match.Chr5.galadriel.LTRdigest.min.gff.txt >../protein_match.Chr5.galadriel.min.gff.txt

cat protein_match.ChrX.galadriel.gff.txt|grep 'MGEScan' >protein_match.ChrX.galadriel.MGEScan.gff.txt
cat protein_match.ChrX.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.galadriel.MGEScan.list.txt
cat protein_match.ChrX.galadriel.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.galadriel.MGEScan.awk.script1.txt
cat protein_match.ChrX.galadriel.MGEScan.list.txt|awk -f protein_match.ChrX.galadriel.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.galadriel.MGEScan.awk.script2.txt
cat protein_match.ChrX.galadriel.gff.txt|awk -f protein_match.ChrX.galadriel.MGEScan.awk.script2.txt >protein_match.ChrX.galadriel.MGEScan.min.gff.txt

cat protein_match.ChrX.galadriel.gff.txt|grep 'LTRdigest' >protein_match.ChrX.galadriel.LTRdigest.gff.txt
cat protein_match.ChrX.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.galadriel.LTRdigest.list.txt
cat protein_match.ChrX.galadriel.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.galadriel.LTRdigest.awk.script1.txt
cat protein_match.ChrX.galadriel.LTRdigest.list.txt|awk -f protein_match.ChrX.galadriel.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.galadriel.LTRdigest.awk.script2.txt
cat protein_match.ChrX.galadriel.gff.txt|awk -f protein_match.ChrX.galadriel.LTRdigest.awk.script2.txt >protein_match.ChrX.galadriel.LTRdigest.min.gff.txt

cat protein_match.ChrX.galadriel.MGEScan.min.gff.txt protein_match.ChrX.galadriel.LTRdigest.min.gff.txt >../protein_match.ChrX.galadriel.min.gff.txt

cd ../zf-CCHC
cat protein_match.Chr1.zf-CCHC.gff.txt|grep 'MGEScan' >protein_match.Chr1.zf-CCHC.MGEScan.gff.txt
cat protein_match.Chr1.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.zf-CCHC.MGEScan.list.txt
cat protein_match.Chr1.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.zf-CCHC.MGEScan.awk.script1.txt
cat protein_match.Chr1.zf-CCHC.MGEScan.list.txt|awk -f protein_match.Chr1.zf-CCHC.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.zf-CCHC.MGEScan.awk.script2.txt
cat protein_match.Chr1.zf-CCHC.gff.txt|awk -f protein_match.Chr1.zf-CCHC.MGEScan.awk.script2.txt >protein_match.Chr1.zf-CCHC.MGEScan.min.gff.txt

cat protein_match.Chr1.zf-CCHC.gff.txt|grep 'LTRdigest' >protein_match.Chr1.zf-CCHC.LTRdigest.gff.txt
cat protein_match.Chr1.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr1.zf-CCHC.LTRdigest.list.txt
cat protein_match.Chr1.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr1.zf-CCHC.LTRdigest.awk.script1.txt
cat protein_match.Chr1.zf-CCHC.LTRdigest.list.txt|awk -f protein_match.Chr1.zf-CCHC.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr1.zf-CCHC.LTRdigest.awk.script2.txt
cat protein_match.Chr1.zf-CCHC.gff.txt|awk -f protein_match.Chr1.zf-CCHC.LTRdigest.awk.script2.txt >protein_match.Chr1.zf-CCHC.LTRdigest.min.gff.txt

cat protein_match.Chr1.zf-CCHC.MGEScan.min.gff.txt protein_match.Chr1.zf-CCHC.LTRdigest.min.gff.txt >../protein_match.Chr1.zf-CCHC.min.gff.txt

cat protein_match.Chr2.zf-CCHC.gff.txt|grep 'MGEScan' >protein_match.Chr2.zf-CCHC.MGEScan.gff.txt
cat protein_match.Chr2.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.zf-CCHC.MGEScan.list.txt
cat protein_match.Chr2.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.zf-CCHC.MGEScan.awk.script1.txt
cat protein_match.Chr2.zf-CCHC.MGEScan.list.txt|awk -f protein_match.Chr2.zf-CCHC.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.zf-CCHC.MGEScan.awk.script2.txt
cat protein_match.Chr2.zf-CCHC.gff.txt|awk -f protein_match.Chr2.zf-CCHC.MGEScan.awk.script2.txt >protein_match.Chr2.zf-CCHC.MGEScan.min.gff.txt

cat protein_match.Chr2.zf-CCHC.gff.txt|grep 'LTRdigest' >protein_match.Chr2.zf-CCHC.LTRdigest.gff.txt
cat protein_match.Chr2.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr2.zf-CCHC.LTRdigest.list.txt
cat protein_match.Chr2.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr2.zf-CCHC.LTRdigest.awk.script1.txt
cat protein_match.Chr2.zf-CCHC.LTRdigest.list.txt|awk -f protein_match.Chr2.zf-CCHC.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr2.zf-CCHC.LTRdigest.awk.script2.txt
cat protein_match.Chr2.zf-CCHC.gff.txt|awk -f protein_match.Chr2.zf-CCHC.LTRdigest.awk.script2.txt >protein_match.Chr2.zf-CCHC.LTRdigest.min.gff.txt

cat protein_match.Chr2.zf-CCHC.MGEScan.min.gff.txt protein_match.Chr2.zf-CCHC.LTRdigest.min.gff.txt >../protein_match.Chr2.zf-CCHC.min.gff.txt

cat protein_match.Chr3.zf-CCHC.gff.txt|grep 'MGEScan' >protein_match.Chr3.zf-CCHC.MGEScan.gff.txt
cat protein_match.Chr3.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.zf-CCHC.MGEScan.list.txt
cat protein_match.Chr3.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.zf-CCHC.MGEScan.awk.script1.txt
cat protein_match.Chr3.zf-CCHC.MGEScan.list.txt|awk -f protein_match.Chr3.zf-CCHC.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.zf-CCHC.MGEScan.awk.script2.txt
cat protein_match.Chr3.zf-CCHC.gff.txt|awk -f protein_match.Chr3.zf-CCHC.MGEScan.awk.script2.txt >protein_match.Chr3.zf-CCHC.MGEScan.min.gff.txt

cat protein_match.Chr3.zf-CCHC.gff.txt|grep 'LTRdigest' >protein_match.Chr3.zf-CCHC.LTRdigest.gff.txt
cat protein_match.Chr3.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr3.zf-CCHC.LTRdigest.list.txt
cat protein_match.Chr3.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr3.zf-CCHC.LTRdigest.awk.script1.txt
cat protein_match.Chr3.zf-CCHC.LTRdigest.list.txt|awk -f protein_match.Chr3.zf-CCHC.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr3.zf-CCHC.LTRdigest.awk.script2.txt
cat protein_match.Chr3.zf-CCHC.gff.txt|awk -f protein_match.Chr3.zf-CCHC.LTRdigest.awk.script2.txt >protein_match.Chr3.zf-CCHC.LTRdigest.min.gff.txt

cat protein_match.Chr3.zf-CCHC.MGEScan.min.gff.txt protein_match.Chr3.zf-CCHC.LTRdigest.min.gff.txt >../protein_match.Chr3.zf-CCHC.min.gff.txt

cat protein_match.Chr4.zf-CCHC.gff.txt|grep 'MGEScan' >protein_match.Chr4.zf-CCHC.MGEScan.gff.txt
cat protein_match.Chr4.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.zf-CCHC.MGEScan.list.txt
cat protein_match.Chr4.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.zf-CCHC.MGEScan.awk.script1.txt
cat protein_match.Chr4.zf-CCHC.MGEScan.list.txt|awk -f protein_match.Chr4.zf-CCHC.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.zf-CCHC.MGEScan.awk.script2.txt
cat protein_match.Chr4.zf-CCHC.gff.txt|awk -f protein_match.Chr4.zf-CCHC.MGEScan.awk.script2.txt >protein_match.Chr4.zf-CCHC.MGEScan.min.gff.txt

cat protein_match.Chr4.zf-CCHC.gff.txt|grep 'LTRdigest' >protein_match.Chr4.zf-CCHC.LTRdigest.gff.txt
cat protein_match.Chr4.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr4.zf-CCHC.LTRdigest.list.txt
cat protein_match.Chr4.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr4.zf-CCHC.LTRdigest.awk.script1.txt
cat protein_match.Chr4.zf-CCHC.LTRdigest.list.txt|awk -f protein_match.Chr4.zf-CCHC.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr4.zf-CCHC.LTRdigest.awk.script2.txt
cat protein_match.Chr4.zf-CCHC.gff.txt|awk -f protein_match.Chr4.zf-CCHC.LTRdigest.awk.script2.txt >protein_match.Chr4.zf-CCHC.LTRdigest.min.gff.txt

cat protein_match.Chr4.zf-CCHC.MGEScan.min.gff.txt protein_match.Chr4.zf-CCHC.LTRdigest.min.gff.txt >../protein_match.Chr4.zf-CCHC.min.gff.txt

cat protein_match.Chr5.zf-CCHC.gff.txt|grep 'MGEScan' >protein_match.Chr5.zf-CCHC.MGEScan.gff.txt
cat protein_match.Chr5.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.zf-CCHC.MGEScan.list.txt
cat protein_match.Chr5.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.zf-CCHC.MGEScan.awk.script1.txt
cat protein_match.Chr5.zf-CCHC.MGEScan.list.txt|awk -f protein_match.Chr5.zf-CCHC.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.zf-CCHC.MGEScan.awk.script2.txt
cat protein_match.Chr5.zf-CCHC.gff.txt|awk -f protein_match.Chr5.zf-CCHC.MGEScan.awk.script2.txt >protein_match.Chr5.zf-CCHC.MGEScan.min.gff.txt

cat protein_match.Chr5.zf-CCHC.gff.txt|grep 'LTRdigest' >protein_match.Chr5.zf-CCHC.LTRdigest.gff.txt
cat protein_match.Chr5.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.Chr5.zf-CCHC.LTRdigest.list.txt
cat protein_match.Chr5.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.Chr5.zf-CCHC.LTRdigest.awk.script1.txt
cat protein_match.Chr5.zf-CCHC.LTRdigest.list.txt|awk -f protein_match.Chr5.zf-CCHC.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.Chr5.zf-CCHC.LTRdigest.awk.script2.txt
cat protein_match.Chr5.zf-CCHC.gff.txt|awk -f protein_match.Chr5.zf-CCHC.LTRdigest.awk.script2.txt >protein_match.Chr5.zf-CCHC.LTRdigest.min.gff.txt

cat protein_match.Chr5.zf-CCHC.MGEScan.min.gff.txt protein_match.Chr5.zf-CCHC.LTRdigest.min.gff.txt >../protein_match.Chr5.zf-CCHC.min.gff.txt

cat protein_match.ChrX.zf-CCHC.gff.txt|grep 'MGEScan' >protein_match.ChrX.zf-CCHC.MGEScan.gff.txt
cat protein_match.ChrX.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.zf-CCHC.MGEScan.list.txt
cat protein_match.ChrX.zf-CCHC.MGEScan.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.zf-CCHC.MGEScan.awk.script1.txt
cat protein_match.ChrX.zf-CCHC.MGEScan.list.txt|awk -f protein_match.ChrX.zf-CCHC.MGEScan.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.zf-CCHC.MGEScan.awk.script2.txt
cat protein_match.ChrX.zf-CCHC.gff.txt|awk -f protein_match.ChrX.zf-CCHC.MGEScan.awk.script2.txt >protein_match.ChrX.zf-CCHC.MGEScan.min.gff.txt

cat protein_match.ChrX.zf-CCHC.gff.txt|grep 'LTRdigest' >protein_match.ChrX.zf-CCHC.LTRdigest.gff.txt
cat protein_match.ChrX.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$6}'|sed -e 's/^Parent=LTR_retrotransposon//g' >protein_match.ChrX.zf-CCHC.LTRdigest.list.txt
cat protein_match.ChrX.zf-CCHC.LTRdigest.gff.txt|sed -e 's/;reading_frame=.*$//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print $9}'|sed -e 's/^Parent=LTR_retrotransposon//g'|sort -n|uniq|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{m"$1"=10}{if($1==\""$1"\"&&m"$1">$2)m"$1"=$2}END{print \""$1"\",m"$1"}"}' >protein_match.ChrX.zf-CCHC.LTRdigest.awk.script1.txt
cat protein_match.ChrX.zf-CCHC.LTRdigest.list.txt|awk -f protein_match.ChrX.zf-CCHC.LTRdigest.awk.script1.txt|sed -e 's/^/Parent=LTR_retrotransposon/g'|awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1";"}'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($9~\""$2"\"&&$6==\""$1"\")print}"}' >protein_match.ChrX.zf-CCHC.LTRdigest.awk.script2.txt
cat protein_match.ChrX.zf-CCHC.gff.txt|awk -f protein_match.ChrX.zf-CCHC.LTRdigest.awk.script2.txt >protein_match.ChrX.zf-CCHC.LTRdigest.min.gff.txt

cat protein_match.ChrX.zf-CCHC.MGEScan.min.gff.txt protein_match.ChrX.zf-CCHC.LTRdigest.min.gff.txt >../protein_match.ChrX.zf-CCHC.min.gff.txt

cd ..
cat ../Combine/Combined.Chr1.gff.txt|grep -v 'protein_match' >protein_match.Chr1.without.protein_match.gff.txt
cat ../Combine/Combined.Chr2.gff.txt|grep -v 'protein_match' >protein_match.Chr2.without.protein_match.gff.txt
cat ../Combine/Combined.Chr3.gff.txt|grep -v 'protein_match' >protein_match.Chr3.without.protein_match.gff.txt
cat ../Combine/Combined.Chr4.gff.txt|grep -v 'protein_match' >protein_match.Chr4.without.protein_match.gff.txt
cat ../Combine/Combined.Chr5.gff.txt|grep -v 'protein_match' >protein_match.Chr5.without.protein_match.gff.txt
cat ../Combine/Combined.ChrX.gff.txt|grep -v 'protein_match' >protein_match.ChrX.without.protein_match.gff.txt

mkdir ../complete1
cat protein_match.Chr1.*|sort -n -k4,4 >../complete1/Chr1.gff.txt
cat protein_match.Chr2.*|sort -n -k4,4 >../complete1/Chr2.gff.txt
cat protein_match.Chr3.*|sort -n -k4,4 >../complete1/Chr3.gff.txt
cat protein_match.Chr4.*|sort -n -k4,4 >../complete1/Chr4.gff.txt
cat protein_match.Chr5.*|sort -n -k4,4 >../complete1/Chr5.gff.txt
cat protein_match.ChrX.*|sort -n -k4,4 >../complete1/ChrX.gff.txt

cd ..
mkdir fasta
cd fasta
cat ../OR/Sp34.genome.v7.7.fa|awk 'BEGIN{FS="\t";OFS="\t"}{if($1~">")print $1"#";else print $0}'|tr -d '\n'|tr '#' '\n'|sed -e 's/>/\n>/g' >Sp34.genome.v7.7.Chr.fa.txt

cat Sp34.genome.v7.7.Chr.fa.txt|grep '>Sp34_Chr1' -A1|grep -v '>' >Sp34.Chr1.split.fa.txt
cat Sp34.genome.v7.7.Chr.fa.txt|grep '>Sp34_Chr2' -A1|grep -v '>' >Sp34.Chr2.split.fa.txt
cat Sp34.genome.v7.7.Chr.fa.txt|grep '>Sp34_Chr3' -A1|grep -v '>' >Sp34.Chr3.split.fa.txt
cat Sp34.genome.v7.7.Chr.fa.txt|grep '>Sp34_Chr4' -A1|grep -v '>' >Sp34.Chr4.split.fa.txt
cat Sp34.genome.v7.7.Chr.fa.txt|grep '>Sp34_Chr5' -A1|grep -v '>' >Sp34.Chr5.split.fa.txt
cat Sp34.genome.v7.7.Chr.fa.txt|grep '>Sp34_ChrX' -A1|grep -v '>' >Sp34.ChrX.split.fa.txt

cd ..
mkdir cd-hit
cd cd-hit
cat ../complete1/Chr1.gff.txt|grep 'long_terminal_repeat'|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$4,$5,$9}'|sed -e 's/;status=solo;/s/g'|sed -e 's/;status=par;/p/g'|sed -e 's/;status=full;/f/g'|sed -e 's/Scan_LTR//g'|sed -e 's/harvest//g'|sed -e 's/Parent=LTR_retrotransposon//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print "1"$1$4,$2,$3}' >Chr1.cd-hit.list.txt
cat ../complete1/Chr2.gff.txt|grep 'long_terminal_repeat'|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$4,$5,$9}'|sed -e 's/;status=solo;/s/g'|sed -e 's/;status=par;/p/g'|sed -e 's/;status=full;/f/g'|sed -e 's/Scan_LTR//g'|sed -e 's/harvest//g'|sed -e 's/Parent=LTR_retrotransposon//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print "2"$1$4,$2,$3}' >Chr2.cd-hit.list.txt
cat ../complete1/Chr3.gff.txt|grep 'long_terminal_repeat'|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$4,$5,$9}'|sed -e 's/;status=solo;/s/g'|sed -e 's/;status=par;/p/g'|sed -e 's/;status=full;/f/g'|sed -e 's/Scan_LTR//g'|sed -e 's/harvest//g'|sed -e 's/Parent=LTR_retrotransposon//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print "3"$1$4,$2,$3}' >Chr3.cd-hit.list.txt
cat ../complete1/Chr4.gff.txt|grep 'long_terminal_repeat'|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$4,$5,$9}'|sed -e 's/;status=solo;/s/g'|sed -e 's/;status=par;/p/g'|sed -e 's/;status=full;/f/g'|sed -e 's/Scan_LTR//g'|sed -e 's/harvest//g'|sed -e 's/Parent=LTR_retrotransposon//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print "4"$1$4,$2,$3}' >Chr4.cd-hit.list.txt
cat ../complete1/Chr5.gff.txt|grep 'long_terminal_repeat'|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$4,$5,$9}'|sed -e 's/;status=solo;/s/g'|sed -e 's/;status=par;/p/g'|sed -e 's/;status=full;/f/g'|sed -e 's/Scan_LTR//g'|sed -e 's/harvest//g'|sed -e 's/Parent=LTR_retrotransposon//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print "5"$1$4,$2,$3}' >Chr5.cd-hit.list.txt
cat ../complete1/ChrX.gff.txt|grep 'long_terminal_repeat'|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$4,$5,$9}'|sed -e 's/;status=solo;/s/g'|sed -e 's/;status=par;/p/g'|sed -e 's/;status=full;/f/g'|sed -e 's/Scan_LTR//g'|sed -e 's/harvest//g'|sed -e 's/Parent=LTR_retrotransposon//g'|awk 'BEGIN{FS="\t";OFS="\t"}{print "X"$1$4,$2,$3}' >ChrX.cd-hit.list.txt

cat Chr1.cd-hit.list.txt|awk 'BEGIN{FS="\t";OFS="\t"}{print "{print \">"$1"\\n\"substr($0,"$2","$3-$2")}"}' >Chr1.cd-hit.awk.script.txt
cat Chr2.cd-hit.list.txt|awk 'BEGIN{FS="\t";OFS="\t"}{print "{print \">"$1"\\n\"substr($0,"$2","$3-$2")}"}' >Chr2.cd-hit.awk.script.txt
cat Chr3.cd-hit.list.txt|awk 'BEGIN{FS="\t";OFS="\t"}{print "{print \">"$1"\\n\"substr($0,"$2","$3-$2")}"}' >Chr3.cd-hit.awk.script.txt
cat Chr4.cd-hit.list.txt|awk 'BEGIN{FS="\t";OFS="\t"}{print "{print \">"$1"\\n\"substr($0,"$2","$3-$2")}"}' >Chr4.cd-hit.awk.script.txt
cat Chr5.cd-hit.list.txt|awk 'BEGIN{FS="\t";OFS="\t"}{print "{print \">"$1"\\n\"substr($0,"$2","$3-$2")}"}' >Chr5.cd-hit.awk.script.txt
cat ChrX.cd-hit.list.txt|awk 'BEGIN{FS="\t";OFS="\t"}{print "{print \">"$1"\\n\"substr($0,"$2","$3-$2")}"}' >ChrX.cd-hit.awk.script.txt

cat ../fasta/Sp34.Chr1.split.fa.txt|awk -f Chr1.cd-hit.awk.script.txt >Sp34.Chr1.cd-hit.fa.txt
cat ../fasta/Sp34.Chr2.split.fa.txt|awk -f Chr2.cd-hit.awk.script.txt >Sp34.Chr2.cd-hit.fa.txt
cat ../fasta/Sp34.Chr3.split.fa.txt|awk -f Chr3.cd-hit.awk.script.txt >Sp34.Chr3.cd-hit.fa.txt
cat ../fasta/Sp34.Chr4.split.fa.txt|awk -f Chr4.cd-hit.awk.script.txt >Sp34.Chr4.cd-hit.fa.txt
cat ../fasta/Sp34.Chr5.split.fa.txt|awk -f Chr5.cd-hit.awk.script.txt >Sp34.Chr5.cd-hit.fa.txt
cat ../fasta/Sp34.ChrX.split.fa.txt|awk -f ChrX.cd-hit.awk.script.txt >Sp34.ChrX.cd-hit.fa.txt

cat Sp34.Chr*.cd-hit.fa.txt >Sp34.cd-hit.fa.txt

cd-hit -i Sp34.cd-hit.fa.txt -o out -c 0.8

cat ../cd-hit/out.clstr|sed -e 's/Cluster //g'|sed -e 's/^.*aa, >//g'|sed -e 's/\.\.\..*$//g'|tr '\n' ';'|tr '>' '\n' >Cluster.table1.txt
cat Cluster.table1.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($0~"f"||$0~"p")print}' >f_and_p.table.txt
cat f_and_p.table.txt |awk 'BEGIN{FS=";";OFS="\t"}{for(i=2;i<NF;i++)print $i,$1}'|sed -e 's/^/Chr/g'|sed -e 's/LTR/\tLTR\t/g'|sed -e 's/MGE/\tMGE\t/g'|sed -e 's/s\t/;\t/g'|sed -e 's/p\t/;\t/g'|sed -e 's/f\t/;\t/g'|sed -e 's/$/;/g'|awk 'BEGIN{FS="\t";OFS="\t"}{print "BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($1~\""$1"\"&&$2~\""$2"\"&&$9~\"LTR_retrotransposon"$3"\")print $0\"Cluster="$4"\"}"}'|sort|uniq >f_and_p.awk.script.txt

cat ../complete1/Chr1.gff.txt |awk -f f_and_p.awk.script.txt|sort -n -k4,4 >Chr1.gff.txt
cat ../complete1/Chr2.gff.txt |awk -f f_and_p.awk.script.txt|sort -n -k4,4 >Chr2.gff.txt
cat ../complete1/Chr3.gff.txt |awk -f f_and_p.awk.script.txt|sort -n -k4,4 >Chr3.gff.txt
cat ../complete1/Chr4.gff.txt |awk -f f_and_p.awk.script.txt|sort -n -k4,4 >Chr4.gff.txt
cat ../complete1/Chr5.gff.txt |awk -f f_and_p.awk.script.txt|sort -n -k4,4 >Chr5.gff.txt
cat ../complete1/ChrX.gff.txt |awk -f f_and_p.awk.script.txt|sort -n -k4,4 >ChrX.gff.txt
cat Chr* >Clustered.v1.gff.txt

cat Cluster.table1.txt|awk 'BEGIN{FS=";"}{if($0~"p"&&$0!~"f"&&$0!~"s")print "Cluster="$1";"}' >par_only.cluster.list.txt
cat par_only.cluster.list.txt |awk 'BEGIN{ORS="\\\|"}{print}'|sed -e 's/\\|$//g' >par_only.grep.script1.txt
cat Clustered.v1.gff.txt |grep -f par_only.grep.script1.txt|awk 'BEGIN{FS="\t";OFS="\t"}{if($3~"protein_match"&&$6<1e-20)print}'|sed -e 's/^.*Cluster=//g'|sed -e 's/^/Cluster=/g'|sort|uniq|awk 'BEGIN{ORS="\\|"}{print}'|sed -e 's/\\|$//g' >par_only.grep.script2.txt
cat par_only.cluster.list.txt|grep -v -f par_only.grep.script2.txt|awk 'BEGIN{ORS="\\|"}{print}'|sed -e 's/\\|$//g' >par_only.grep.script3.txt
cat Clustered.v1.gff.txt|grep -v -f par_only.grep.script3.txt >../Clustered.v2.gff.txt
