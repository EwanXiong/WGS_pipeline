####=================================================set.dir===================================================####
file_name=/opt/tsinghua/NuoHe/WGS_23.6.1/shell1/name.txt
input=/opt/tsinghua/NuoHe/WGS_23.6.1/rawdata
output=/opt/tsinghua/NuoHe/WGS_23.6.1/output1
#match=/opt/tsinghua/NuoHe/WES_23.6.5/shell/match.txt
align=/opt/tsinghua/cfDNApipeTest/file/hg19/bwa
#align=/opt/tsinghua/cfDNApipeTest/file/hg38/bwa
gatk_vcf=/opt/tsinghua/cfDNApipeTest/file/vcf
pon_hg19_wes=/opt/tsinghua/cfDNApipeTest/file/vcf/pon/somatic-hg19_Mutect2-exome-panel.vcf.gz
pon_hg19_wgs=/opt/tsinghua/cfDNApipeTest/file/vcf/pon/somatic-hg19_Mutect2-WGS-panel.vcf.gz
pon_hg38=/opt/tsinghua/cfDNApipeTest/file/vcf/pon/somatic-hg38_1000g_pon.hg38.vcf.gz
pon=$pon_hg19_wgs
hg=hg19
CNVkit_ref1=/opt/tsinghua/cfDNApipeTest/file/CNVkit/refFlat_hg19.txt
CNVkit_ref2=/opt/tsinghua/cfDNApipeTest/file/CNVkit/access-5kb-mappable.hg19.bed
#CNVkit_ref1=/opt/tsinghua/cfDNApipeTest/file/CNVkit/refFlat_hg38.txt
#CNVkit_ref2=/opt/tsinghua/cfDNApipeTest/file/CNVkit/access-excludes.hg38.bed
####================================================pipeline====================================================####
echo "[`date`] ===============================WGS Analysis Start====================================="
name=($(cat $file_name)) 
if [ ! -d "$output/1_qc" ]; then
		echo
		echo 
		echo "[`date`] qc"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/1_qc
		myvar=0
		for i in ${name[@]}
		do
			time /opt/tsinghua/software/fastqc/fastqc/FastQC/fastqc \
			--outdir $output/1_qc  \
			--threads 30  \
			$input/$i* &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
				myvar=0
				wait
			fi
		done
		wait
		multiqc $output/1_qc/ -o $output/1_qc/
		echo '[`date`] Run complete'
		echo '-----------------------------------------------------------------------------'
fi
if [ ! -d "$output/2_low_complexity_filter" ]; then
		echo
		echo 
		echo "[`date`] low_complexity_filter"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/2_low_complexity_filter
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/2_low_complexity_filter/"$i"
			time fastp \
			--thread 16 \
			--correction \
			--detect_adapter_for_pe \
			--json $output/2_low_complexity_filter/"$i"/"$i".json \
			--html $output/2_low_complexity_filter/"$i"/"$i".html \
			--report_title "$i" \
			-i $input/"$i"_1.fq.gz \
			-I $input/"$i"_2.fq.gz \
			-o $output/2_low_complexity_filter/"$i"/"$i".pair1.truncated.gz  \
			-O $output/2_low_complexity_filter/"$i"/"$i".pair2.truncated.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
				myvar=0
				wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '-----------------------------------------------------------------------------'
fi
if [ ! -d "$output/3_bwa/" ]; then
		echo
		echo
		echo "[`date`] bwa/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/3_bwa
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/3_bwa/"$i"
			time bwa mem \
			-t 30 -M \
			-R "@RG\tID:$i\tPL:illumina\tLB:$i\tSM:$i" \
			$align/"$hg".fa \
			$output/2_low_complexity_filter/"$i"/"$i".pair1.truncated.gz \
			$output/2_low_complexity_filter/"$i"/"$i".pair2.truncated.gz | samtools view -Sb - -o $output/3_bwa/"$i"/"$i".bam &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
				myvar=0
				wait
			fi
		done 
		wait
		
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi
if [ ! -d "$output/4_bamsort/" ]; then
		echo
		echo
		echo "[`date`] bamsort/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/4_bamsort
		myvar=0
		for i in ${name[@]}
		do
			time samtools sort \
			-@ 30 \
			-O bam \
			-o $output/4_bamsort/"$i"_sorted.bam $output/3_bwa/"$i"/"$i".bam ;\
			time samtools index \
			-b  $output/4_bamsort/"$i"_sorted.bam &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
				myvar=0
				wait
			fi
		done 
		wait
		mkdir -p $output/4_bamsort/stat 
		for i in ${name[@]}
		do
			time samtools flagstat $output/4_bamsort/"$i"_sorted.bam > $output/4_bamsort/stat/"$i".stat &
		done
		wait
		multiqc $output/4_bamsort/stat/*.stat -o $output/4_bamsort/stat/
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi
if [ ! -d "$output/5_Qualimap/" ]; then
		echo
		echo
		echo "[`date`] Qualimap/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/5_Qualimap
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/5_Qualimap/"$i"
			time qualimap bamqc \
			--java-mem-size=8G \
			-nt 30 -bam $output/4_bamsort/"$i"_sorted.bam \
			-outdir $output/5_Qualimap/"$i"/"$i" \
			-outformat PDF:HTML  &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
				myvar=0
				wait
			fi
		done 
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi
#if [ ! -d "$output/collectHSmetrics/" ]; then
#		echo
#		echo
#		echo "[`date`] collectHSmetrics/"
#		echo '-----------------------------------------------------------------------------'
#		mkdir -p $output/collectHSmetrics
#		gatk BedToIntervalList \
#		-I $bed/xgen-exome-research-panel-v2-targets-${hg}.bed \
#		-O $output/collectHSmetrics/"$hg".exon_interval.bed \
#		-SD $align/"$hg".dict
#		wait
#		myvar=0
#		for i in ${name[@]}
#		do
#			time gatk CollectHsMetrics \
#			-BI $output/collectHSmetrics/"$hg".exon_interval.bed \
#			-I $output/4_bamsort/"$i"_sorted.bam \
#			-O $output/collectHSmetrics/"$i".gatk.stat.txt \
#			-TI $output/collectHSmetrics/"$hg".exon_interval.bed &
#			myvar=$(($myvar + 1 ))
#			if [ "$myvar" = "4" ]
#			then
#				myvar=0
#				wait
#			fi
#		done
#		wait
#		multiqc $output/collectHSmetrics/*.stat.txt -o $output/collectHSmetrics/
#		echo
#		echo 'run complete'
#		echo '-----------------------------------------------------------------------------'
#fi
echo "[`date`] ======================================GATK PROCESS START====================================="
if [ ! -d "$output/6_RemoveDuplicates/" ]; then
		echo
		echo
		echo "[`date`] RemoveDuplicates/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/6_RemoveDuplicates
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/6_RemoveDuplicates/"$i"
			time gatk --java-options "-Xmx6G" MarkDuplicates \
			--REMOVE_DUPLICATES true \
			--INPUT $output/4_bamsort/"$i"_sorted.bam  \
			--METRICS_FILE $output/6_RemoveDuplicates/"$i"/"$i"-rmdup.txt \
			--OUTPUT $output/6_RemoveDuplicates/"$i"/"$i"-rmdup.bam ;\
			time samtools index \
			-b  $output/6_RemoveDuplicates/"$i"/"$i"-rmdup.bam &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
				myvar=0
				wait
			fi
		done 
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi

if [ ! -d "$output/8_BaseRecalibrator/" ]; then
		echo
		echo
		echo "[`date`] BaseRecalibrator/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/8_BaseRecalibrator
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/8_BaseRecalibrator/"$i"
			time gatk BaseRecalibrator \
			-I $output/6_RemoveDuplicates/"$i"/"$i"-rmdup.bam \
			-R $align/"$hg".fa \
			--known-sites $gatk_vcf/1000G_omni2.5."$hg".vcf \
			--known-sites $gatk_vcf/1000G_phase1.snps.high_confidence."$hg".vcf \
			--known-sites $gatk_vcf/dbsnp_146."$hg".vcf  \
			--known-sites $gatk_vcf/Mills_and_1000G_gold_standard.indels."$hg".vcf  \
			-O $output/8_BaseRecalibrator/"$i"/"$i".recal.table &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
				myvar=0
				wait
			fi
		done 
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi
if [ ! -d "$output/9_BQSR/" ]; then
		echo
		echo
		echo "[`date`]BQSR/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/9_BQSR
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/9_BQSR/"$i"
			time gatk ApplyBQSR \
			-I $output/6_RemoveDuplicates/"$i"/"$i"-rmdup.bam \
			-R $align/"$hg".fa \
			--bqsr-recal-file $output/8_BaseRecalibrator/"$i"/"$i".recal.table  \
			-O $output/9_BQSR/"$i"/"$i"-BQSR.bam &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done 
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi

if [ ! -d "$output/10_getPileup/" ]; then
		echo
		echo
		echo "[`date`] getPileup/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/10_getPileup
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/10_getPileup/"$i"
			time gatk GetPileupSummaries \
			-I $output/9_BQSR/"$i"/"$i"-BQSR.bam \
			-V $gatk_vcf/small_exac_common_3_"$hg".SNP_biallelic.vcf \
			-L $gatk_vcf/small_exac_common_3_"$hg".SNP_biallelic.vcf \
			-O $output/10_getPileup/"$i"/"$i".pileups.table &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done 
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi
if [ ! -d "$output/11_contamination/" ]; then
		echo
		echo
		echo "[`date`] contamination/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/11_contamination
		myvar=0
		for i in ${name[@]}
		do
		  tumor=$i
		  mkdir -p $output/11_contamination/"$tumor"
			time gatk CalculateContamination \
			-I $output/10_getPileup/"$tumor"/"$tumor".pileups.table \
			-O $output/11_contamination/"$tumor"/"$tumor".contamination.table &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi

if [ ! -d "$output/12_Somatic_mutect2/" ]; then
		echo
		echo
		echo "[`date`] Somatic_mutect2/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/12_Somatic_mutect2
		myvar=0
		for i in ${name[@]}
		do
		  tumor=$i
		  mkdir -p $output/12_Somatic_mutect2/"$tumor"
			time gatk Mutect2 \
			-R $align/"$hg".fa \
			-I $output/9_BQSR/"$tumor"/"$tumor"-BQSR.bam \
			--germline-resource $gatk_vcf/af-only-gnomad."$hg".vcf.gz \
			--panel-of-normals $pon \
			-O $output/12_Somatic_mutect2/"$tumor"/"$tumor".unfiltered.vcf.gz \
			--f1r2-tar-gz $output/12_Somatic_mutect2/"$tumor"/"$tumor".f1r2.tar.gz;\
			time gatk LearnReadOrientationModel \
			-I $output/12_Somatic_mutect2/"$tumor"/"$tumor".f1r2.tar.gz \
			-O $output/12_Somatic_mutect2/"$tumor"/"$tumor".read-orientation-model.tar.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi
if [ ! -d "$output/13_Somatic_filtermutectcalls/" ]; then
		echo
		echo
		echo "[`date`] Somatic_filtermutectcalls/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/13_Somatic_filtermutectcalls
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/13_Somatic_filtermutectcalls/"$i"
			time gatk FilterMutectCalls \
			-V $output/12_Somatic_mutect2/"$i"/"$i".unfiltered.vcf.gz \
			-R $align/"$hg".fa  \
			--min-median-mapping-quality 20 \
			--ob-priors $output/12_Somatic_mutect2/"$i"/"$i".read-orientation-model.tar.gz \
			--contamination-table $output/11_contamination/"$i"/"$i".contamination.table \
			-O $output/13_Somatic_filtermutectcalls/"$i"/"$i".filtered.vcf.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi

if [ ! -d "$output/13_Somatic_bcftools/" ]; then
		echo
		echo
		echo "[`date`] Somatic_bcftools/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/13_Somatic_bcftools
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/13_Somatic_bcftools/"$i"
			time bcftools view -f 'PASS' -i 'INFO/DP>10' \
			$output/13_Somatic_filtermutectcalls/"$i"/"$i".filtered.vcf.gz \
			-o $output/13_Somatic_bcftools/"$i"/"$i".somatic.vcf.gz -O z
			time gatk IndexFeatureFile --input $output/13_Somatic_bcftools/"$i"/"$i".somatic.vcf.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi

if [ ! -d "$output/Somatic_SNP_vcf/" ]; then
		echo
		echo
		echo "[`date`] Somatic_SNP_vcf/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/Somatic_SNP_vcf
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/Somatic_SNP_vcf/"$i"
			time gatk SelectVariants \
			-R $align/"$hg".fa \
			--select-type-to-include SNP \
			-V $output/13_Somatic_bcftools/"$i"/"$i".somatic.vcf.gz \
			-O $output/Somatic_SNP_vcf/"$i"/"$i".SNP.filtered.vcf.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi
if [ ! -d "$output/Somatic_INDEL_vcf/" ]; then
		echo
		echo
		echo "[`date`] Somatic_INDEL_vcf/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/Somatic_INDEL_vcf
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/Somatic_INDEL_vcf/"$i"
			time gatk SelectVariants \
			-R $align/"$hg".fa \
			--select-type-to-include INDEL \
			-V $output/13_Somatic_bcftools/"$i"/"$i".somatic.vcf.gz \
			-O $output/Somatic_INDEL_vcf/"$i"/"$i".indel.filtered.vcf.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi
if [ ! -d "$output/14_Somatic_annovar/" ]; then
		echo
		echo
		echo "[`date`] Somatic_annovar/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/14_Somatic_annovar/snp
		mkdir -p $output/14_Somatic_annovar/indel
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/14_Somatic_annovar/snp/"$i"
			mkdir -p $output/14_Somatic_annovar/indel/"$i"
			time /opt/tsinghua/cfDNApipeTest/software/annovar/annovar/table_annovar.pl \
			$output/Somatic_SNP_vcf/"$i"/"$i".SNP.filtered.vcf.gz \
			/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/humandb \
			--out $output/14_Somatic_annovar/snp/"$i"/"$i".snp.somatic \
			--buildver "$hg" \
			--protocol refGene,HGMD,omim,mimTitles,morbidmap,cytoBand,wgRna,genomicSuperDups,avsnp150,cosmic70,clinvar_20210501,gwasCatalog,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,esp6500siv2_all,exac03,dbnsfp35a,gnomad_genome,gnomad211_exome,gnomad211_genome,cancer_hotspots_INDEL_v2,cancer_hotspots_SNV_v2,icgc28,intervar_20180118,pmkb_simple \
			--operation g,r,r,f,f,r,r,r,f,f,f,r,f,f,f,f,f,f,f,f,f,r,r,r,f,f,r,r,f \
			--argument ',,,,,,,,,,,,,,,,,,,,,,,,,,,,' \
			--thread 30 \
			--remove \
			--nastring . \
			--vcfinput \
			--intronhgvs 300 ;\
			time /opt/tsinghua/cfDNApipeTest/software/annovar/annovar/table_annovar.pl \
			$output/Somatic_INDEL_vcf/"$i"/"$i".indel.filtered.vcf.gz \
			/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/humandb \
			--out $output/14_Somatic_annovar/indel/"$i"/"$i".indel.somatic \
			--buildver "$hg" \
			--protocol refGene,HGMD,omim,mimTitles,morbidmap,cytoBand,wgRna,genomicSuperDups,avsnp150,cosmic70,clinvar_20210501,gwasCatalog,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,esp6500siv2_all,exac03,dbnsfp35a,gnomad_genome,gnomad211_exome,gnomad211_genome,cancer_hotspots_INDEL_v2,cancer_hotspots_SNV_v2,icgc28,intervar_20180118,pmkb_simple \
			--operation g,r,r,f,f,r,r,r,f,f,f,r,f,f,f,f,f,f,f,f,f,r,r,r,f,f,r,r,f \
			--argument ',,,,,,,,,,,,,,,,,,,,,,,,,,,,' \
			--thread 30 \
			--remove \
			--nastring . \
			--vcfinput \
			--intronhgvs 300 &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi 
if [ ! -d "$output/15_VariantQC_Somatic/" ]; then
		echo
		echo
		echo "[`date`] VariantQC_Somatic/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/15_VariantQC_Somatic
		myvar=0
		for i in ${name[@]}
		do
			#echo "-V:"$i" $output/13_Somatic_filtermutectcalls/"$i"/"$i".filtered.vcf.gz" >> $output/15_VariantQC_Somatic/Somatic.txt &
			echo "-V:"$i" $output/13_Somatic_bcftools/"$i"/"$i".somatic.vcf.gz" >> $output/15_VariantQC_Somatic/Somatic.txt &
			echo "-V:"$i" $output/Somatic_SNP_vcf/"$i"/"$i".SNP.filtered.vcf.gz" >> $output/15_VariantQC_Somatic/Somatic_SNP.txt &
			echo "-V:"$i" $output/Somatic_INDEL_vcf/"$i"/"$i".indel.filtered.vcf.gz" >> $output/15_VariantQC_Somatic/Somatic_indel.txt &
		done
		wait
		cat $output/15_VariantQC_Somatic/Somatic.txt|tr "\n" " " > $output/15_VariantQC_Somatic/Somatic_ref.txt &
		cat $output/15_VariantQC_Somatic/Somatic_SNP.txt|tr "\n" " " > $output/15_VariantQC_Somatic/Somatic_SNP_ref.txt &
		cat $output/15_VariantQC_Somatic/Somatic_indel.txt|tr "\n" " " > $output/15_VariantQC_Somatic/Somatic_indel_ref.txt &
		wait
		somatic=($(less $output/15_VariantQC_Somatic/Somatic_ref.txt))
		snp_ref=($(less $output/15_VariantQC_Somatic/Somatic_SNP_ref.txt))
		indel_ref=($(less $output/15_VariantQC_Somatic/Somatic_indel_ref.txt))
		time java -jar /opt/tsinghua/software/DISCVRSeq-1.3.13/DISCVRSeq-1.3.13.jar VariantQC \
		-R $align/"$hg".fa \
		"${somatic[@]}" \
		--rd $output/15_VariantQC_Somatic/somatic_reports.json\
		-O $output/15_VariantQC_Somatic/somatic_report.html &
		time java -jar /opt/tsinghua/software/DISCVRSeq-1.3.13/DISCVRSeq-1.3.13.jar VariantQC \
		-R $align/"$hg".fa \
		"${snp_ref[@]}" \
		--rd $output/15_VariantQC_Somatic/snp_reports.json\
		-O $output/15_VariantQC_Somatic/snp_report.html &
		time java -jar /opt/tsinghua/software/DISCVRSeq-1.3.13/DISCVRSeq-1.3.13.jar VariantQC \
		-R $align/"$hg".fa \
		"${indel_ref[@]}" \
		--rd $output/15_VariantQC_Somatic/indel_reports.json\
		-O $output/15_VariantQC_Somatic/indel_report.html &
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi

if [ ! -d "$output/12_germline_Haplotypecaller/" ]; then
		echo
		echo
		echo "[`date`] germline_Haplotypecaller/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/12_germline_Haplotypecaller
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/12_germline_Haplotypecaller/"$i"
			time gatk --java-options "-Xmx6g" HaplotypeCaller \
			-R $align/"$hg".fa \
			-I $output/9_BQSR/"$i"/"$i"-BQSR.bam \
			-O $output/12_germline_Haplotypecaller/"$i"/"$i".unfiltered.vcf.gz \
			--bamout $output/12_germline_Haplotypecaller/"$i"/"$i".out.bam &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi
#exit
if [ ! -d "$output/13_germline_CNNScoreVariants/" ]; then
		echo
		echo
		echo "[`date`] germline_CNNScoreVariants/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/13_germline_CNNScoreVariants
		source ~/anaconda3/bin/activate gatk1
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/13_germline_CNNScoreVariants/"$i"
      #1D:
      time gatk CNNScoreVariants \
         -V $output/12_germline_Haplotypecaller/"$i"/"$i".unfiltered.vcf.gz \
         -R $align/"$hg".fa \
         -O $output/13_germline_CNNScoreVariants/"$i"/"$i".1Dannotated.vcf.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi

if [ ! -d "$output/13_germline_FilterVariantTranches/" ]; then
		echo
		echo
		echo "[`date`] germline_FilterVariantTranches/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/13_germline_FilterVariantTranches
		source ~/anaconda3/bin/activate gatk1
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/13_germline_FilterVariantTranches/"$i"
      time gatk FilterVariantTranches \
         -V $output/13_germline_CNNScoreVariants/"$i"/"$i".1Dannotated.vcf.gz \
         --resource $gatk_vcf/hapmap_3.3."$hg".vcf \
         --resource $gatk_vcf/Mills_and_1000G_gold_standard.indels."$hg".vcf \
         --info-key CNN_1D \
         --snp-tranche 99.95 \
         --indel-tranche 99.4 \
         -O $output/13_germline_FilterVariantTranches/"$i"/"$i".1D.filtered.vcf.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi


if [ ! -d "$output/germline_SNP_vcf/" ]; then
		echo
		echo
		echo "[`date`] germline_SNP_vcf/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/germline_SNP_vcf
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/germline_SNP_vcf/"$i"
			time gatk SelectVariants \
			-R $align/"$hg".fa \
			--select-type-to-include SNP \
			-V $output/13_germline_FilterVariantTranches/"$i"/"$i".1D.filtered.vcf.gz \
			-O $output/germline_SNP_vcf/"$i"/"$i".SNP.vcf.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi
if [ ! -d "$output/germline_INDEL_vcf/" ]; then
		echo
		echo
		echo "[`date`] germline_INDEL_vcf/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/germline_INDEL_vcf
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/germline_INDEL_vcf/"$i"
			time gatk SelectVariants \
			-R $align/"$hg".fa \
			--select-type-to-include INDEL \
			-V $output/13_germline_FilterVariantTranches/"$i"/"$i".1D.filtered.vcf.gz \
			-O $output/germline_INDEL_vcf/"$i"/"$i".indel.vcf.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi
if [ ! -d "$output/14_germline_annovar/" ]; then
		echo
		echo
		echo "[`date`] germline_annovar/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/14_germline_annovar/snp
		mkdir -p $output/14_germline_annovar/indel
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/14_germline_annovar/snp/"$i"
			mkdir -p $output/14_germline_annovar/indel/"$i"
			time /opt/tsinghua/cfDNApipeTest/software/annovar/annovar/table_annovar.pl \
			$output/germline_SNP_vcf/"$i"/"$i".SNP.vcf.gz \
			/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/humandb \
			--out $output/14_germline_annovar/snp/"$i"/"$i".snp.germline \
			--buildver "$hg" \
			--protocol refGene,HGMD,omim,mimTitles,morbidmap,cytoBand,wgRna,genomicSuperDups,avsnp150,cosmic70,clinvar_20210501,gwasCatalog,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,esp6500siv2_all,exac03,dbnsfp35a,gnomad_genome,gnomad211_exome,gnomad211_genome,cancer_hotspots_INDEL_v2,cancer_hotspots_SNV_v2,icgc28,intervar_20180118,pmkb_simple \
			--operation g,r,r,f,f,r,r,r,f,f,f,r,f,f,f,f,f,f,f,f,f,r,r,r,f,f,r,r,f \
			--argument ',,,,,,,,,,,,,,,,,,,,,,,,,,,,' \
			--thread 30 \
			--remove \
			--nastring . \
			--vcfinput \
			--intronhgvs 300 &
			time /opt/tsinghua/cfDNApipeTest/software/annovar/annovar/table_annovar.pl \
			$output/germline_INDEL_vcf/"$i"/"$i".indel.vcf.gz \
			/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/humandb \
			--out $output/14_germline_annovar/indel/"$i"/"$i".indel.germline \
			--buildver "$hg" \
			--protocol refGene,HGMD,omim,mimTitles,morbidmap,cytoBand,wgRna,genomicSuperDups,avsnp150,cosmic70,clinvar_20210501,gwasCatalog,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,esp6500siv2_all,exac03,dbnsfp35a,gnomad_genome,gnomad211_exome,gnomad211_genome,cancer_hotspots_INDEL_v2,cancer_hotspots_SNV_v2,icgc28,intervar_20180118,pmkb_simple \
			--operation g,r,r,f,f,r,r,r,f,f,f,r,f,f,f,f,f,f,f,f,f,r,r,r,f,f,r,r,f \
			--argument ',,,,,,,,,,,,,,,,,,,,,,,,,,,,' \
			--thread 30 \
			--remove \
			--nastring . \
			--vcfinput \
			--intronhgvs 300 &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi 
if [ ! -d "$output/15_VariantQC_germline/" ]; then
		echo
		echo
		echo "[`date`] VariantQC_germline/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/15_VariantQC_germline
		myvar=0
		for i in ${name[@]}
		do
			echo "-V:"$i" $output/13_germline_FilterVariantTranches/"$i"/"$i".1D.filtered.vcf.gz" >> $output/15_VariantQC_germline/germline.txt &
			echo "-V:"$i" $output/germline_SNP_vcf/"$i"/"$i".SNP.vcf.gz" >> $output/15_VariantQC_germline/germline_SNP.txt &
			echo "-V:"$i" $output/germline_INDEL_vcf/"$i"/"$i".indel.vcf.gz" >> $output/15_VariantQC_germline/germline_indel.txt &
		done
		wait
		cat $output/15_VariantQC_germline/germline.txt|tr "\n" " " > $output/15_VariantQC_germline/germline_ref.txt &
		cat $output/15_VariantQC_germline/germline_SNP.txt|tr "\n" " " > $output/15_VariantQC_germline/germline_SNP_ref.txt &
		cat $output/15_VariantQC_germline/germline_indel.txt|tr "\n" " " > $output/15_VariantQC_germline/germline_indel_ref.txt &
		wait
		germline=($(less $output/15_VariantQC_germline/germline_ref.txt))
		snp_ref=($(less $output/15_VariantQC_germline/germline_SNP_ref.txt))
		indel_ref=($(less $output/15_VariantQC_germline/germline_indel_ref.txt))
		time java -jar /opt/tsinghua/software/DISCVRSeq-1.3.13/DISCVRSeq-1.3.13.jar VariantQC \
		-R $align/"$hg".fa \
		"${germline[@]}" \
		--rd $output/15_VariantQC_germline/germline_reports.json \
		-O $output/15_VariantQC_germline/germline_report.html &
		time java -jar /opt/tsinghua/software/DISCVRSeq-1.3.13/DISCVRSeq-1.3.13.jar VariantQC \
		-R $align/"$hg".fa \
		"${snp_ref[@]}" \
		--rd $output/15_VariantQC_germline/snp_reports.json \
		-O $output/15_VariantQC_germline/snp_report.html &
		time java -jar /opt/tsinghua/software/DISCVRSeq-1.3.13/DISCVRSeq-1.3.13.jar VariantQC \
		-R $align/"$hg".fa \
		"${indel_ref[@]}" \
		--rd $output/15_VariantQC_germline/indel_reports.json \
		-O $output/15_VariantQC_germline/indel_report.html &
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi

echo "[`date`] ======================================CNV Start====================================="
if [ ! -d "$output/16_cnvbatch" ]; then
		echo
		echo 
		echo "[`date`] cnvbatch"
		echo '-------------------------'
		mkdir -p $output/16_cnvbatch
		bam_tumor=($(for i in ${name[@]};do ls $output/6_RemoveDuplicates/$i/*-rmdup.bam;done))
		time cnvkit.py batch \
		"${bam_tumor[@]}" \
		--normal \
		--annotate $CNVkit_ref1 \
		--access $CNVkit_ref2 \
		--fasta $align/"$hg".fa \
		-y \
		--method wgs \
		--output-reference $output/16_cnvbatch/reference.cnn \
		--output-dir $output/16_cnvbatch \
		-p 30
		wait
		echo '[`date`] Run complete'
		echo '------------------------'
fi
#-x male/female在已知的情况下可以指定，否则不加-x参数软件自动guess sex
if [ ! -d "$output/17_cnvHeatmap" ]; then
		echo
		echo
		echo "[`date`] cnvHeatmap"
		echo '-------------------------'
		mkdir -p $output/17_cnvHeatmap
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/17_cnvHeatmap/"$i"
			time cnvkit.py heatmap \
			$output/16_cnvbatch/"$i"-rmdup.cnr \
			-d \
			-y \
			-x female \
			-o $output/17_cnvHeatmap/"$i"/"$i"_heatmap.pdf  &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		#time cnvkit.py heatmap $(ls $output/16_cnvbatch/*-rmdup.cnr) -d -y -o $output/17_cnvHeatmap/heatmap.pdf &
		echo '[`date`] Run complete'
		echo '------------------------'
fi

if [ ! -d "$output/18_cnvPlot" ]; then
		echo
		echo
		echo "[`date`] cnvPlot"
		echo '-------------------------'
		mkdir -p $output/18_cnvPlot
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/18_cnvPlot/"$i"
			time cnvkit.py diagram \
			$output/16_cnvbatch/"$i"-rmdup.cnr \
			-s $output/16_cnvbatch/"$i"-rmdup.cns \
			-o $output/18_cnvPlot/"$i"/"$i"_diagram.pdf \
			--title "$i" \
			--threshold 0.5 \
			--min-probes 3 -y -x female &
			time cnvkit.py scatter \
			$output/16_cnvbatch/"$i"-rmdup.cnr \
			-s $output/16_cnvbatch/"$i"-rmdup.cns \
			-o $output/18_cnvPlot/"$i"/"$i"_scatter.pdf \
			--title "$i" &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------'
fi
if [ ! -d "$output/19_cnvTable" ]; then
		echo
		echo
		echo "[`date`] cnvTable"
		echo '-------------------------'
		mkdir -p $output/19_cnvTable
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/19_cnvTable/"$i"
			time cnvkit.py breaks $output/16_cnvbatch/"$i"-rmdup.cnr $output/16_cnvbatch/"$i"-rmdup.cns  --min-probes 1 \
			-o $output/19_cnvTable/"$i"/"$i"_breaks.txt &
			time cnvkit.py genemetrics $output/16_cnvbatch/"$i"-rmdup.cnr -s $output/16_cnvbatch/"$i"-rmdup.cns --threshold 0.1 --min-probes 3 -y -x female \
			-o $output/19_cnvTable/"$i"/"$i"_genemetrics_cnrs.txt;cnvkit.py genemetrics $output/16_cnvbatch/"$i"-rmdup.cnr --threshold 0.1 --min-probes 3 -y -x female \
			-o $output/19_cnvTable/"$i"/"$i"_genemetrics_cnr.txt;tail -n +2 \
			$output/19_cnvTable/"$i"/"$i"_genemetrics_cnr.txt| cut -f 1 |sort > $output/19_cnvTable/"$i"/"$i"_cnrs_gene.txt;tail -n +2 \
			$output/19_cnvTable/"$i"/"$i"_genemetrics_cnr.txt| cut -f 1 |sort > $output/19_cnvTable/"$i"/"$i"_cnr_gene.txt;comm -12 \
			$output/19_cnvTable/"$i"/"$i"_cnrs_gene.txt $output/19_cnvTable/"$i"/"$i"_cnr_gene.txt  > $output/19_cnvTable/"$i"/"$i"_genemetrics_gene.txt  &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------'
fi
echo "[`date`] ======================================SUMMARY COUNT START====================================="
if [ ! -d "$output/summary/filter_result" ]; then
		echo
		echo
		echo "[`date`] filter_result"
		echo '-------------------------'
		mkdir -p $output/summary
		mkdir -p $output/summary/filter_result
		myvar=0
		for i in ${name[@]}
		do
			sed -n '15p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\t"total_reads"://'|sed 's/,//' >>$output/summary/filter_result/clean_reads.txt ;\
			sed -n '19p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\t"q20_rate"://'|sed 's/,//' >>$output/summary/filter_result/Q20.txt ;\
			sed -n '20p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\t"q30_rate"://'|sed 's/,//' >>$output/summary/filter_result/Q30.txt ;\
			sed -n '23p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\t"gc_content"://' >>$output/summary/filter_result/GC_content.txt ;\
			sed -n '36p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\"rate"://'|sed 's/[ ]//'|sed 's/,//' >>$output/summary/filter_result/duplication.txt ;\
			sed -n '41p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\"peak"://'|sed 's/[ ]//'|sed 's/,//' >>$output/summary/filter_result/insert_size_peak.txt
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		paste -d '\t' $file_name $output/summary/filter_result/clean_reads.txt $output/summary/filter_result/Q20.txt $output/summary/filter_result/Q30.txt \
		$output/summary/filter_result/GC_content.txt $output/summary/filter_result/duplication.txt $output/summary/filter_result/insert_size_peak.txt|sed '1i\Sample\tclean_reads\tQ20\tQ30\tGC_content\tduplication\tinsert_size_peak' >>$output/summary/filter_result/result.txt
		echo '[`date`] Run complete'
		echo '------------------------'
fi
if [ ! -d "$output/summary/mapping_result" ]; then
		echo
		echo
		echo "[`date`] mapping_result"
		echo '-------------------------'
		mkdir -p $output/summary/mapping_result
		myvar=0
		for i in ${name[@]}
		do
			sed -n '5p' $output/4_bamsort/stat/"$i".stat |cut -d ':' -f1|cut -d '(' -f2  >>$output/summary/mapping_result/Mapping_effciency.txt ;\
			grep "number of mapped paired reads (both in pair)" $output/5_Qualimap/"$i"/"$i"/genome_results.txt|sed 's/number of mapped paired reads (both in pair) = //'|sed 's/[ ]*//' >>$output/summary/mapping_result/Paried_map.txt ;\
			grep "median insert size =" $output/5_Qualimap/"$i"/"$i"/genome_results.txt|sed 's/mean insert size = //'|sed 's/[ ]*//' >>$output/summary/mapping_result/Insert_Size.txt ;\
			grep "duplication rate =" $output/5_Qualimap/"$i"/"$i"/genome_results.txt|sed 's/duplication rate = //'|sed 's/[ ]*//' >>$output/summary/mapping_result/Duplication.txt ;\
			grep "mean coverageData =" $output/5_Qualimap/"$i"/"$i"/genome_results.txt|sed 's/mean coverageData = //'|sed 's/[ ]*//'|sed 's/X//' >>$output/summary/mapping_result/Average_Depth.txt ;\
			grep "of reference with a coverageData >= 4X" $output/5_Qualimap/"$i"/"$i"/genome_results.txt|sed 's/There is a //'|sed 's/of reference with a coverageData >= 4X//'|sed 's/[ ]*//' >>$output/summary/mapping_result/Coverage_4X.txt ;\
			grep "of reference with a coverageData >= 10X" $output/5_Qualimap/"$i"/"$i"/genome_results.txt|sed 's/There is a //'|sed 's/of reference with a coverageData >= 10X//'|sed 's/[ ]*//' >>$output/summary/mapping_result/Coverage_10X.txt ;\
			grep "of reference with a coverageData >= 20X" $output/5_Qualimap/"$i"/"$i"/genome_results.txt|sed 's/There is a //'|sed 's/of reference with a coverageData >= 20X//'|sed 's/[ ]*//' >>$output/summary/mapping_result/Coverage_20X.txt
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "1" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		paste -d '\t' $file_name $output/summary/mapping_result/Mapping_effciency.txt $output/summary/mapping_result/Paried_map.txt $output/summary/mapping_result/Insert_Size.txt \
		$output/summary/mapping_result/Duplication.txt $output/summary/mapping_result/Average_Depth.txt $output/summary/mapping_result/Coverage_4X.txt \
		$output/summary/mapping_result/Coverage_10X.txt $output/summary/mapping_result/Coverage_20X.txt|sed '1i\Sample\tMapping_effciency\tParied_map\tInsert_Size\tDuplication\tAverage_Depth\tCoverage>4X\tCoverage>10X\tCoverage>20X' >>$output/summary/mapping_result/result.txt
		echo '[`date`] Run complete'
		echo '------------------------'
fi
if [ ! -d "$output/summary/snp_result" ]; then
		echo
		echo
		echo "[`date`] snp_result"
		echo '-------------------------'
		mkdir -p $output/summary/snp_result
		for i in ${name[@]}
		do
			mkdir -p $output/summary/snp_result/"$i"
			cut -f 1,9 $output/14_germline_annovar/indel/"$i"/"$i".indel.germline."$hg"_multianno.txt |grep -v -w "\." >>$output/summary/snp_result/"$i"/"$i".germline.indel_plot.txt ;\
			cut -f 1,9 $output/14_germline_annovar/snp/"$i"/"$i".snp.germline."$hg"_multianno.txt |grep -v -w "\." >>$output/summary/snp_result/"$i"/"$i".germline.snp_plot.txt ;\
			cut -f 1,9 $output/14_Somatic_annovar/indel/"$i"/"$i".indel.somatic."$hg"_multianno.txt |grep -v -w "\." >>$output/summary/snp_result/"$i"/"$i".somatic.indel_plot.txt  ;\
			cut -f 1,9 $output/14_Somatic_annovar/snp/"$i"/"$i".snp.somatic."$hg"_multianno.txt |grep -v -w "\." >>$output/summary/snp_result/"$i"/"$i".somatic.snp_plot.txt
		done
		wait
		for i in ${name[@]}
		do
			Rscript /opt/tsinghua/pipleline_sh/WES/new_sh/SNV_plot.r --input_dir $output/summary/snp_result/"$i"
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------'
fi
#svaba配对样本有-n可以生成somatic和germline，单样本不能分
if [ ! -d "$output/17_Svaba" ]; then
    echo
    echo
    echo "[`date`] Svaba"
    echo '-------------------------'
    mkdir -p $output/17_Svaba
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/17_Svaba/"$i"
			time svaba run \
			-p 30 \
			-G $align/"$hg".fa \
			-t $output/6_RemoveDuplicates/"$i"/"$i"-rmdup.bam \
			-D $gatk_vcf/dbsnp_146."$hg".vcf \
			-a $output/17_Svaba/"$i"/"$i" \
		  --g-zip &
		done
        wait
        echo '[`date`] Run complete'
        echo '------------------------'
fi
#if [ ! -d "$output/17_Svaba_germline" ]; then
#    echo
#    echo
#    echo "[`date`] Svaba_germline"
#    echo '-------------------------'
#    mkdir -p $output/17_Svaba_germline
#		myvar=0
#		for i in ${name[@]}
#		do
#			mkdir -p $output/17_Svaba_germline/"$i"
#			time svaba run \
#			-p 30 \
#			-G $align/"$hg".fa \
#			-t $output/6_RemoveDuplicates/"$i"/"$i"-rmdup.bam \
#			-D $gatk_vcf/dbsnp_146."$hg".vcf \
#			-a $output/17_Svaba_germline/"$i"/"$i".germline_run \
#			--germline \
#		  --g-zip &
#		done
#        wait
#        echo '[`date`] Run complete'
#        echo '------------------------'
#fi
exit # A+ 分析， 去掉exit
if [ ! -d "$output/susceptibility_gene" ]; then
		echo
		echo 
		echo "[`date`] susceptibility_gene"
		echo '-------------------------'
		mkdir -p $output/susceptibility_gene
		for i in ${name[@]}
		do
			mkdir -p $output/susceptibility_gene/"$i"
			cut -f '1-10,17,29-31,121,131' $output/14_germline_annovar/indel/"$i"/"$i".indel.germline."$hg"_multianno.txt >$output/susceptibility_gene/"$i"/"$i"_indel_input.txt 
			cut -f '1-10,17,29-31,121,131' $output/14_germline_annovar/snp/"$i"/"$i".snp.germline."$hg"_multianno.txt >$output/susceptibility_gene/"$i"/"$i"_snp_input.txt 
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------'
fi
if [ ! -d "$output/Driver_gene/" ]; then
		echo
		echo
		echo "[`date`] Driver_gene/"
		echo '-----------------------------------------------------------------------------'
		mkdir -p $output/Driver_gene
		for i in ${name[@]}
		do
			mkdir -p $output/Driver_gene/"$i"
			cut -f '1-10,17,29-31,121,131' $output/14_Somatic_annovar/indel/"$i"/"$i".indel.somatic."$hg"_multianno.txt >$output/Driver_gene/"$i"/"$i"_indel_input.txt 
			cut -f '1-10,17,29-31,121,131' $output/14_Somatic_annovar/snp/"$i"/"$i".snp.somatic."$hg"_multianno.txt >$output/Driver_gene/"$i"/"$i"_snp_input.txt 
		done
		wait
		echo 
		echo 'run complete'
		echo '-----------------------------------------------------------------------------'
fi 
echo "[`date`] ======================================TASK FINISH====================================="
