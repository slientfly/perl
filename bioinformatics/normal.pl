use DBI;
use strict;
use warnings;
use feature 'say';
use Term::ANSIColor;
use Parallel::ForkManager;
use Getopt::Long qw/GetOptions/;

my ($sample,$threads,$r1,$r2,$panel,$lane,$pbs_id,$random,$email,$gender);
GetOptions(
	's=s' => \$sample,
	'r1=s' => \$r1,
	'r2=s' => \$r2,
	'c=s' => \$panel,
	'l=s' => \$lane,
	'r=s' => \$random,
	'e=s' => \$email,
	'g=s' => \$gender,
	't=s' => \$threads,
);
my $usage =<<EOF;

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

程序简介：      遗传病分析流程，在线版
使用方法：      perl target_pipeline_web.pl arg1 arg2 ...
更新日期：      2017年5月3日

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EOF
unless($sample && $r1 && $r2){
        say colored($usage, "yellow on_magenta");
        exit 0;
}

$threads = 7;
my $intervals = &intervals();

if($panel =~ /MitoChip/){
	&mtDNA;
	$panel =~ s/MitoChip[,]?//g;
	$panel =~ s/,MitoChip//g;
}
if($panel =~ /,CNV|^CNV,|^CNV$/){
	&CNV;
	$panel =~ s/,CNV//g;
	$panel =~ s/^CNV[,]?//g;
}
if($panel ne ''){
	&normal;
}


sub normal{
	my $cutadapt = '/bin/cutadapt';
	my $ref = '/local_disk/DB/hg19/Sequence/BWAIndex/genome.fa';
	my $ref_fa = '/local_disk/DB/hg19/Sequence/WholeGenomeFasta/genome.fa';
	my $picard = '/disk1/software/picard-tools-2.2.3/picard.jar';
	my $bamtools = '/disk1/software/bin/bamtools';
	my $bwa = '/disk1/software/bin/bwa';
	my $samtools = '/disk1/software/bin/samtools';
	my $bamToBed = '/disk1/software/bin/bamToBed';
	my $coverageBed = '/disk1/software/bin/coverageBed';
	my $GATK = '/disk1/software/GATK-3.7/GenomeAnalysisTK.jar';
	my $dbsnp = '/local_disk/DB/dbsnp/dbsnp_147.hg19.vcf';
	my $varscan  = '/disk1/software/VarScan.v2.3.7.jar';
	my $beddir = '/disk1/bed/';	
	my $disk = "/ssd1/mygeno/pengj/project/$sample";
	my $java = "java -Djava.io.tmpdir=$disk -Xmx15g -Xms5g -jar";
	unless(-e $disk){`mkdir -p $disk`;}
	chdir $disk;
	if ($r1=~/,/ && $r2=~/,/){
		my @reads1=split ",", $r1;
		my @reads2=split ",", $r2;
		for(my $i=0; $i<@reads1; $i++){
			&parallel(
				"cat $reads1[$i] >> ${sample}_R1.fq.gz",
				"cat $reads2[$i] >> ${sample}_R2.fq.gz"
			);
		}
	} else {
		&parallel(
			"cp $r1 ${sample}_R1.fq.gz",
			"cp $r2 ${sample}_R2.fq.gz"
		);
	}
	###&mysql_inject("fastq文件质控","正在运行",10);
	chomp(my $hardware = `zcat ${sample}_R1.fq.gz | head -1 | cut -d: -f1`);
	if($hardware =~ /NS500/){
		`$cutadapt -m 80 -q10,10 --nextseq-trim=10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
	}else{
		`$cutadapt -m 80 -q10,10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
		`$cutadapt -m 20 -q10,10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
	}
	###&mysql_inject("将fastq比对到人类参考基因组上","正在运行",20);
	`$bwa mem -M -t $threads $ref $sample.clean_R1.fq.gz $sample.clean_R2.fq.gz | $samtools fixmate -O bam - - | $samtools sort -\@ $threads -m 1G - $sample.sort`;
	`$samtools index $sample.sort.bam > $sample.sort.bam.bai`;
	###&mysql_inject("添加BAM文件头信息","正在运行",30);
	`$java $picard AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=$sample.sort.bam OUTPUT=$sample.sort.header.bam RGLB=$ref RGPL=ILLUMINA RGSM=GP1 RGPU=GRP1 SORT_ORDER=coordinate CREATE_INDEX=true`;
	###&mysql_inject("过滤BAM文件","正在运行",35);
	`$bamtools filter -isMapped true -isPaired true -isProperPair true -in $sample.sort.header.bam -out $sample.sort.flt.bam`;
	`rm -f $sample.sort.header.bam $sample.sort.header.bai`;
	if ($panel =~ /DMD|Clinican|Canmute|SLC26A4|Deaf_V2|lung8|LAMA2|PKU|LBS|MG_BIK_V1/){
		my @parallel_command=(
			"$java $picard MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.rmdup.sorted.bam M=duplication-report2.txt",
		);
		if($panel =~ /PKU/){
			push @parallel_command,"perl /disk1/software/crest/extractSClip.pl -i $sample.sort.bam --ref_genome $ref -r chr6 > $sample.extractSClip.log 2>&1";
		}else{
			for my $chr (1 .. 22,'X','Y'){
				push @parallel_command,"perl /disk1/software/crest/extractSClip.pl -i $sample.sort.bam --ref_genome $ref -r chr$chr > $sample.extractSClip.log 2>&1";
			}
		}
		###&mysql_inject("去除PCR冗余序列/提取soft-clip reads","正在运行",40);
		&parallel(@parallel_command);
		###&mysql_inject("合并$sample.sort.bam.*.cover/$sample.sort.bam.*.sclip.txt","正在运行",45);
		`cat $sample.sort.bam.*.cover > $sample.sort.bam.cover`;
		`cat $sample.sort.bam.*.sclip.txt > $sample.sort.bam.sclip.txt`;
	}else{
		###&mysql_inject("去除PCR冗余序列","正在运行",40);
		`$java $picard MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.rmdup.sorted.bam M=duplication-report2.txt`;
	}
	###&mysql_inject("标记PCR冗余序列/拷贝BAM文件到公共盘/bamtools stats/bamToBed","正在运行",45);
	&parallel(
		"$java $picard MarkDuplicates VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.sort.flt.mark.bam M=duplication-report.txt",
		"$bamtools stats -in $sample.sort.bam > $sample.sort.bam.stat",
		"$bamToBed -i $sample.rmdup.sorted.bam > $sample.rmdup.sorted.bamtobed",
	);
	`rm -f $sample.sort.flt.ba?`;
	if ($panel =~ /DMD|Clinican|Canmute|SLC26A4|Deaf_V2|lung8|LAMA2|PKU|LBS|MG_BIK_V1/){
		my @parallel_command=(
			"$java $GATK -T BaseRecalibrator -R $ref_fa -I $sample.sort.flt.mark.bam -knownSites $dbsnp -o recal.table -nct $threads $intervals"
		);
		if ($panel =~ /PKU/ || $panel =~ /LAMA2/){
			push @parallel_command,"perl /disk1/software/crest/CREST.pl -l 150 -f $sample.sort.bam.cover -d $sample.sort.bam --ref_genome $ref -t /local_disk/DB/hg19/Sequence/blat/genome.2bit -r chr6 -p $sample.sort.bam.chr6 > $sample.chr6.CREST.log 2>&1";
		}else{
			for my $chr (1 .. 22,'X','Y'){
				push @parallel_command,"perl /disk1/software/crest/CREST.pl -l 150 -f $sample.sort.bam.cover -d $sample.sort.bam --ref_genome $ref -t /local_disk/DB/hg19/Sequence/blat/genome.2bit -r chr$chr -p $sample.sort.bam.chr$chr > $sample.chr$chr.CREST.log 2>&1";
			}
		}			
		if($panel =~ /DMD/){
			push @parallel_command,"$samtools depth -aa -r chrX:31135345-33359192 -d 100000 $sample.rmdup.sorted.bam > $sample.DMD.depth";
		}elsif($panel =~ /SLC26A4/){
		    push @parallel_command,"$samtools depth -aa -r chr7:107301080-107358252 -d 100000 $sample.rmdup.sorted.bam > $sample.SLC26A4.depth";
		}elsif($panel =~ /PKU/ || $panel =~ /LAMA2/){
			push @parallel_command,"$samtools depth -aa -r chr6:129204100-129838000 -d 100000 $sample.rmdup.sorted.bam > $sample.LAMA2.depth";
		}
		###&mysql_inject("并行:校正平台误差/检测融合基因","正在运行",50);
		&parallel(@parallel_command);
		`cat *.chr*.predSV.txt > $sample.sort.bam.predSV.txt`;
		`perl /disk2/mygeno/web/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt`;
		if($panel =~ /DMD/){
			&parallel(
				"perl /disk2/mygeno/web/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt",
				"perl /disk2/mygeno/web/script/batchmendelian/siteplot_dmd.pl $sample.DMD.depth"
			);
		}elsif($panel =~ /SLC26A4/){
			&parallel(
				"perl /disk2/mygeno/web/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt",
				"perl /disk2/mygeno/web/script/batchmendelian/siteplot_SLC26A4.pl $sample.SLC26A4.depth SLC26A4 SLC26A4 chr7:107301080-107358252"
			);
		}elsif($panel =~ /PKU/ || $panel =~ /LAMA2/){
			&parallel(
                "perl /disk2/mygeno/web/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt",
                "perl /disk2/mygeno/web/script/batchmendelian/siteplot_LAMA2.pl $sample.LAMA2.depth"
            );
		}
	}
	unless (glob "${sample}*crest.xls"){
		###&mysql_inject("校正平台误差","正在运行",55);
		`$java $GATK -T BaseRecalibrator -R $ref_fa -I $sample.sort.flt.mark.bam -knownSites $dbsnp -o recal.table -nct $threads $intervals`;
	}
	###&mysql_inject("输出高质量比对序列","正在运行",60);
	`$java $GATK -T PrintReads -R $ref_fa -I $sample.sort.flt.mark.bam -BQSR recal.table -o $sample.recal.bam -nct $threads`;
	`rm -f $sample.sort.flt.mark.ba?`;
	for my $sub_chip ((split /,/,$panel)){
		next if $sub_chip eq '';
		if($sub_chip =~ /Wscancer|Canmute|Clinican|BloodV2|Immu_All|ImmuV2|MetFull|MetAdd|geneis|Can06|Lym_MyGeno|lung8/){
			#MPIleup+GATK
			###&mysql_inject("并行:mpileup低频突变/HaplotypeCaller常规突变/coverageBed","正在运行",65);
			&parallel(
				"$samtools mpileup -s -f $ref -l $beddir/$sub_chip.bed -d 10000 -L 10000 $sample.recal.bam | $java $varscan mpileup2cns --min-reads2 1 --min-coverage 1 --strand-filter 1 --output-vcf 1 --variants 1 --min-var-freq 0.0001 --p-value 1 --min-avg-qual 20 > $sample.varscan.$sub_chip.vcf",
				"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/$sub_chip.bed -ip 100 -D $dbsnp -o $sample.gatk.$sub_chip.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
				"$coverageBed -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.coverage",
				"$coverageBed -d -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.target.coverage"
			);
			###&mysql_inject("捕获效率统计","正在运行",65);
			`perl /disk2/mygeno/web/script/batchmendelian/depth_count.pl $sample.sort.$sub_chip.coverage $sample.sort.$sub_chip.target.coverage > $sample.$sub_chip.sort.bam.depth.statistics.xls`;
			`perl /disk2/mygeno/web/script/batchmendelian/read_count.pl $disk/$sample.$sub_chip.sort.bam.depth.statistics.xls $sub_chip`;
			`perl /disk2/mygeno/web/script/batchmendelian/vcf2vcf.pl $sample.gatk.$sub_chip.vcf > $sample.gatk.$sub_chip.tmp.vcf`;
			`perl /disk2/mygeno/web/script/batchmendelian/merge_vcf_for_gatk_varscan.pl $sample.gatk.$sub_chip.tmp.vcf $sample.varscan.$sub_chip.vcf > $sample.gatk.$sub_chip.split.vcf`;

		}else{
			#GATK
			if ($sub_chip eq 'Deaf_V2'){
				###&mysql_inject("并行:HaplotypeCaller(常规区域和miRNA区域)/coverageBed","正在运行",65);
				&parallel(
					"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/$sub_chip.bed -ip 100 -D $dbsnp -o $sample.gatk.$sub_chip.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
					"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/${sub_chip}_miRNA_region.bed -ip 100 -D $dbsnp -o $sample.gatk.${sub_chip}.miRNA.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
					"$coverageBed -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.coverage",
					"$coverageBed -d -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.target.coverage",
				);
				&parallel(
					"perl /disk2/mygeno/web/script/batchmendelian/depth_count.pl $sample.sort.$sub_chip.coverage $sample.sort.$sub_chip.target.coverage > $sample.$sub_chip.sort.bam.depth.statistics.xls"
				);
				`perl /disk2/mygeno/web/script/batchmendelian/read_count.pl $disk/$sample.$sub_chip.sort.bam.depth.statistics.xls $sub_chip`;
			}else{
				###&mysql_inject("并行:HaplotypeCaller/coverageBed","正在运行",65);
				&parallel(
					"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/$sub_chip.bed -ip 100 -D $dbsnp -o $sample.gatk.$sub_chip.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
					"$coverageBed -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.coverage",
					"$coverageBed -d -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.target.coverage"
				);
				`perl /disk2/mygeno/web/script/batchmendelian/depth_count.pl $sample.sort.$sub_chip.coverage $sample.sort.$sub_chip.target.coverage > $sample.$sub_chip.sort.bam.depth.statistics.xls`;
				`perl /disk2/mygeno/web/script/batchmendelian/read_count.pl $disk/$sample.$sub_chip.sort.bam.depth.statistics.xls $sub_chip`;
			}
			#不需要mpileup的通用
			###&mysql_inject("split vcf","正在运行",60);
			`perl /disk2/mygeno/web/script/batchmendelian/vcf2vcf.pl $sample.gatk.$sub_chip.vcf > $sample.gatk.$sub_chip.split.vcf`;
		}
		#所有Panel通用
		###&mysql_inject("并行:ANNOVAR功能注释/拷贝统计文件","正在运行",70);
		&parallel(
			"perl /disk1/software/annovar/table_annovar.pl $sample.gatk.$sub_chip.split.vcf /local_disk/DB/annovar_db/humandb/ -remove -buildver hg19 -protocol refGene,cytoBand,dbNsfpInterPro,Inhouse,esp6500siv2_all,1000g2015aug_all,exac03,avsnp147,cosmic74,dbnsfp30a,clinvar_20160302,spidex,mcap10,revel,dbscsnv11,dbnsfp31a_interpro,gnomad_exome,gnomad_genome,InterVar -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -argument \"-splicing_threshold 5 --hgvs,,,,,,,,,,,,,,,,,,\" -vcfinput -nastring . --thread $threads --maxgenethread $threads",
			"perl /disk2/mygeno/web/script/batchmendelian/depth_count3.pl $sample.sort.$sub_chip.target.coverage $sample $sub_chip",
			"mv $sample.${sub_chip}.count.xls Statistics.xls",
		);
		###&mysql_inject("并行:vcf2gff/拷贝捕获效率统计文件","正在运行",80);
		&parallel(
			"perl /disk2/mygeno/web/script/batchmendelian/ParseVCFHgvs.pl -v $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf -f $sample.gatk.$sub_chip.split.vcf.refGene.exonic_variant_function -g $gender > $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff",
		);
		###&mysql_inject("关联HGMD数据库","正在运行",90);
		`perl /disk2/mygeno/web/script/batchmendelian/HGMD.pl $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff $sample $sub_chip`;
		###&mysql_inject("将捕获效率统计信息追加到Statistics.xls","正在运行",99);
	
		`sed -i '1iSample\tRaw_data_bases(Mb)\tClean_data_bases(Mb)\tAligned_bases(Mb)\tAligned\tInitial bases on target\tBase covered on target\tCoverage of target region\tEffective bases on target\tFraction of effective bases on target\tAverage sequencing depth on target\tFraction of target covered with at least 4X\tFraction of target covered with at least 10X\tFraction of target covered with at least 20X\tduplication rate\tchip' Statistics.xls`;
		`sed -i '1iDescr\thgvs\tDisease\tTag\tPublication\tEffect\tSample\tGene_Symbol\tID\tRef_Transcript\tExon\tNucleotide_Changes\tAmino_Acid_Changes\tGene_Type\t1000g2015aug_all\tPathogenic_Analysis\tInhert\tDisease\tDiseasePhenotype\tDiseaseInformation\tclinvar_20160302\tref/alt(fre)\tChr\tBegin\tEnd\tRef\tAlt\tMutRatio\tRefCount\tMutCount\tDepth\tMutation_Type\tGene_Region\tCytoBand\tdbsnp147\tPathSNP\tMutInNormal\t1000Genome\tMutInDatabase\t1000g2015aug_all\tESP6500si\tInhouse\tExAC_ALL\tExAC_EAS\tSIFT\tSIFT_Predict\tPolyPhen_2\tPolyPhen_2_Predict\tMutationTaster\tMutationTaster_Predict\tGERP++\tGERP++_Predict\tSPIDEX\tAA_change\tdbNsfpInterPro\tREVEL_score\tMCAP_score\tMCAP_pred\tgnomAD_exome_ALL\tgnomAD_exome_EAS\tgnomAD_genome_ALL\tgnomAD_genome_EAS\tInterVar\tHighest-MAF' $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff.hgmd.xls`;
		`perl /disk2/mygeno/web/script/batchmendelian/MakeSelectedGeneWeb.pl  $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff.hgmd.xls ${sample}_SNP_INDEL.xlsx	$disk`;
	}
}

sub mtDNA{
	my $cutadapt = '/usr/bin/cutadapt';
	my $ref = '/local_disk/DB/hg38/chrM/chrM.fa';
	my $picard = '/disk1/software/picard-tools-2.2.3/picard.jar';
	my $bamtools = '/disk1/software/bin/bamtools';
	my $bwa = '/disk1/software/bin/bwa';
	my $samtools = '/disk1/software/bin/samtools';
	my $bamToBed = '/disk1/software/bin/bamToBed';
	my $coverageBed = '/disk1/software/bin/coverageBed';
	my $GATK = '/disk1/software/GATK-3.6/GenomeAnalysisTK.jar';
	my $dbsnp = '/local_disk/DB/dbsnp/dbsnp_144.chrM.hg38.vcf';
	my $beddir = '/disk1/bed/';
	my $disk = "/ssd1/mygeno/hanwj/project/MitoChip/$sample";
	my $java = "java -Djava.io.tmpdir=$disk -Xmx15g -Xms5g -jar";
	unless(-e $disk){`mkdir -p $disk`;}
	chdir $disk;
	if ($r1=~/,/ && $r2=~/,/){
        my @reads1=split ",", $r1;
        my @reads2=split ",", $r2;
        for(my $i=0; $i<@reads1; $i++){
            &parallel(
                "cat $reads1[$i] >> ${sample}_R1.fq.gz",
                "cat $reads2[$i] >> ${sample}_R2.fq.gz"
            );
        }
    } else {
        &parallel(
            "cp $r1 ${sample}_R1.fq.gz",
            "cp $r2 ${sample}_R2.fq.gz"
        );
    }	
	###&mysql_inject("去除低质量Reads(连续N序列、低于80bp序列、引物自连序列)","正在运行",10);
        chomp(my $hardware = `zcat ${sample}_R1.fq.gz | head -1 | cut -d: -f1`);
        if($hardware =~ /NS500/){
                `$cutadapt -m 80 -q10,10 --nextseq-trim=10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a ATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
        }else{
                `$cutadapt -m 80 -q10,10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a ATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
        }
	###&mysql_inject("BWA mem(比对到UCSC hg38)","正在运行",20);
	`$bwa mem -M -t $threads $ref $sample.clean_R1.fq.gz $sample.clean_R2.fq.gz | $samtools fixmate -O bam - - | $samtools sort -\@ $threads -m 1G - $sample.sort`;
	###&mysql_inject("samtools index","正在运行",22);
	`$samtools index $sample.sort.bam > $sample.sort.bam.bai`;
	###&mysql_inject("AddOrReplaceReadGroups","正在运行",25);
	`$java $picard AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=$sample.sort.bam OUTPUT=$sample.sort.header.bam RGLB=$ref RGPL=ILLUMINA RGSM=GP1 RGPU=GRP1 SORT_ORDER=coordinate CREATE_INDEX=true`;
	###&mysql_inject("bamtools filter","正在运行",30);
	`$bamtools filter -isMapped true -isPaired true -isProperPair true -in $sample.sort.header.bam -out $sample.sort.flt.bam`;
	###&mysql_inject("samtools index","正在运行",30);
	`$samtools index $sample.sort.flt.bam > $sample.sort.flt.bam.bai`;
	###&mysql_inject("picard MarkDuplicates(去除冗余)","正在运行",35);
	`$java $picard MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.rmdup.sorted.bam M=duplication-report2.txt`;
	###&mysql_inject("并行:picard MarkDuplicates(标记冗余)/拷贝bam文件到公共盘/bamtools stats/bamToBed","正在运行",40);
	&parallel(
		"$java $picard MarkDuplicates VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.sort.flt.mark.bam M=duplication-report.txt",
		"$bamtools stats -in $sample.sort.bam > $sample.sort.bam.stat",
		"$bamToBed -i $sample.rmdup.sorted.bam > $sample.rmdup.sorted.bamtobed",
	);
	###&mysql_inject("并行:BaseRecalibrator","正在运行",45);
	&parallel(
		"$java $GATK -T BaseRecalibrator -R $ref -I $sample.sort.flt.mark.bam -knownSites $dbsnp -o recal.table -nct $threads -L chrM:1-16569",
		"$coverageBed -a $sample.rmdup.sorted.bamtobed -b $beddir/hg38_MT_wholegenome.bed > $sample.sort.MitoChip.coverage",
		"$coverageBed -d -a $sample.rmdup.sorted.bamtobed -b $beddir/hg38_MT_wholegenome.bed > $sample.sort.MitoChip.target.coverage",
	);
	###&mysql_inject("并行:PrintReads","正在运行",50);
	&parallel(
		"$java $GATK -T PrintReads -R $ref -I $sample.sort.flt.mark.bam -BQSR recal.table -o $sample.recal.bam -nct $threads",
		"perl /disk2/mygeno/web/script/batchmendelian/depth_count.pl $sample.sort.MitoChip.coverage $sample.sort.MitoChip.target.coverage > $sample.sort.bam.depth.statistics.xls",
	);
	###&mysql_inject("并行:UnifiedGenotyper/samtools depth/统计捕获效率","正在运行",55);
	&parallel(
		"$java $GATK -T UnifiedGenotyper -R $ref -I $sample.recal.bam -L chrM:1-16569 -D $dbsnp -o $sample.gatk.raw.vcf -A Coverage -A HaplotypeScore --downsample_to_coverage 1000000 -nct $threads -glm BOTH",
		"$samtools depth -aa -r chrM:1-16569 -d 100000 $sample.recal.bam > $sample.recal.bam.depth",
		"perl /disk2/mygeno/web/script/batchmendelian/depth_count3.pl $sample.sort.MitoChip.target.coverage $sample MitoChip",
	);
	###&mysql_inject("SitePlot","正在运行",60);
	`perl /disk2/mygeno/web/script/batchmendelian/smooth_mt_depth_and_plot_r.pl $sample.recal.bam.depth`;
	###&mysql_inject("ANNOVAR功能注释","正在运行",70);
	`perl /disk1/software/annovar/table_annovar.pl $sample.gatk.raw.vcf.split.vcf /local_disk/DB/annovar_db/humandb -buildver hg38_MT -protocol refGene,snp142NonFlagged -operation g,f -nastring . -vcfinput`;
	###&mysql_inject("提取结果","正在运行",75);
	`perl /disk2/mygeno/web/script/batchmendelian/extraction_vcf_for_mtDNA.pl $sample.gatk.raw.vcf.split.vcf.hg38_MT_multianno.vcf`;
	# `sed -i '1iSample\tRaw_data_bases(Mb)\tClean_data_bases(Mb)\tAligned_bases(Mb)\tAligned\tInitial bases on target\tBase covered on target\tCoverage of target region\tEffective bases on target\tFraction of effective bases on target\tAverage sequencing depth on target\tFraction of target covered with at least 4X\tFraction of target covered with at least 10X\tFraction of target covered with at least 20X\tduplication rate\tchip' Statistics.xls`;
	# `sed -i '1iDescr\thgvs\tDisease\tTag\tPublication\tEffect\tSample\tGene_Symbol\tID\tRef_Transcript\tExon\tNucleotide_Changes\tAmino_Acid_Changes\tGene_Type\t1000g2015aug_all\tPathogenic_Analysis\tInhert\tDisease\tDiseasePhenotype\tDiseaseInformation\tclinvar_20160302\tref/alt(fre)\tChr\tBegin\tEnd\tRef\tAlt\tMutRatio\tRefCount\tMutCount\tDepth\tMutation_Type\tGene_Region\tCytoBand\tdbsnp147\tPathSNP\tMutInNormal\t1000Genome\tMutInDatabase\t1000g2015aug_all\tESP6500si\tInhouse\tExAC_ALL\tExAC_EAS\tSIFT\tSIFT_Predict\tPolyPhen_2\tPolyPhen_2_Predict\tMutationTaster\tMutationTaster_Predict\tGERP++\tGERP++_Predict\tSPIDEX\tAA_change\tdbNsfpInterPro\tREVEL_score\tMCAP_score\tMCAP_pred\tgnomAD_exome_ALL\tgnomAD_exome_EAS\tgnomAD_genome_ALL\tgnomAD_genome_EAS\tInterVar\tHighest-MAF' `;
	
}

sub CNV{
	my $cutadapt = '/usr/bin/cutadapt';
	my $samtools = '/disk1/software/bin/samtools';
	my $ref = '/local_disk/DB/hg19/Sequence/BWAIndex/genome.fa';
	my $PSCC = '/disk2/mygeno/yangh/lccnv/script/';
	my $disk = "/ssd1/mygeno/hanwj/project/CNV/$sample";
	my $java = "java -Djava.io.tmpdir=$disk -Xmx15g -Xms5g -jar";
	my $picard = '/disk1/software/picard-tools-2.2.3/picard.jar';
	unless(-e $disk){`mkdir -p $disk`;}
	chdir $disk;
	###&mysql_inject("拷贝fastq文件","正在运行",5);
	if ($r1=~/,/ && $r2=~/,/){
        my @reads1=split ",", $r1;
        my @reads2=split ",", $r2;
        for(my $i=0; $i<@reads1; $i++){
            &parallel(
                "cat $reads1[$i] >> ${sample}_R1.fq.gz",
                "cat $reads2[$i] >> ${sample}_R2.fq.gz"
            );
        }
    } else {
        &parallel(
            "cp $r1 ${sample}_R1.fq.gz",
            "cp $r2 ${sample}_R2.fq.gz"
        );
    }
	###&mysql_inject("去除低质量Reads(连续N序列、低于50bp序列、引物自连序列)","正在运行",10);
	chomp(my $hardware = `zcat ${sample}_R1.fq.gz | head -1 | cut -d: -f1`);
	if($hardware =~ /NS500/){
		`$cutadapt -m 50 -q10,10 --nextseq-trim=10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a ATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
	}else{
		`$cutadapt -m 50 -q10,10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a ATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
	}
	###&mysql_inject("并行:BWA aln(比对到UCSC hg19)","正在运行",20);
	&parallel(
		"/disk1/software/bin/bwa aln -t 4 $ref $sample.clean_R1.fq.gz > $sample.clean_R1.sai",
		"/disk1/software/bin/bwa aln -t 4 $ref $sample.clean_R2.fq.gz > $sample.clean_R2.sai"
	);
	###&mysql_inject("BWA sampe/fixmate","正在运行",25);
	`bwa sampe $ref $sample.clean_R1.sai $sample.clean_R2.sai $sample.clean_R1.fq.gz $sample.clean_R2.fq.gz > $sample.sam`;
	`$samtools fixmate -O bam $sample.sam -| $samtools sort -\@ 5 -m 1G - $sample.sort`;
	`$samtools index $sample.sort.bam > $sample.sort.bam.bai`;
	###&mysql_inject("检测插入片段长度分布","正在运行中",28);
	`$java $picard CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT I=$sample.sort.bam O=$sample.raw.insertsize.txt H=$sample.raw.insertsize.pdf M=0.5`;
	###&mysql_inject("测序数据及捕获效率统计","正在运行",29);
	#`$java $picard CollectRawWgsMetrics VALIDATION_STRINGENCY=SILENT I=$sample.sort.bam O=raw_wgs_metrics.txt R=$ref INCLUDE_BQ_HISTOGRAM=true`;
	my $depth = `perl /disk2/mygeno/web/script/batchmendelian/cnv-rawfq-Statistics.pl $sample ${sample}_R1.fq.gz ${sample}_R2.fq.gz $sample.clean_R1.fq.gz $sample.clean_R2.fq.gz`;
	###&mysql_inject("将sam文件转换为ext文件","正在运行",30);
	`perl $PSCC/Script/Soap2Ext.pl -soap $sample.sam -o $sample.ext -pe`;
	chomp(my $result = `perl /disk2/mygeno/yangh/lccnv/script/karyotype.pl $sample.sort.bam`);
	my ($gender,$karotype) = split /\t/,$result;
	if($gender eq 'M' || $gender eq 'F'){
		###&mysql_inject("生成ratio文件，性别:$gender","正在运行",40);
		`perl $PSCC/Script/Ratio.pl -window /disk2/mygeno/yangh/lccnv/script/Window/0.6X/winRatio.Tags25k.Slide1k.Hg19.gz -ext $sample.ext -gender $gender -o $sample.ratio`;
		###&mysql_inject("并行:Segmentation/SVG核型图","正在运行",50);
		&parallel(
			"perl $PSCC/Script/Segmentation.pl -w 400 -in $sample.ratio -o $sample.cnv -c $PSCC/Attachment/hg19.length -n $PSCC/Attachment/hg19.N_region -log $sample.ratio.cnv.log",
			"perl /disk2/mygeno/web/script/batchmendelian/svg_v3.pl $sample.ratio $gender",
			"perl /disk2/mygeno/yangh/lccnv/script/svg_v4.pl $sample.ratio $gender $sample.downsample.point"
		);
		###&mysql_inject("CNVFilterD","正在运行",60);
		`perl $PSCC/Script/CNVFilterD.pl -i $sample.cnv -g $gender -o $sample.cnv.filter -cyto /disk1/DB/annovar_db//humandb/hg19_cytoBand.txt -ref /disk1/DB/annovar_db//humandb/hg19_refGene.txt -sp 0.05 -pp 0.05`;
		`perl $PSCC/cnv-autofilter.pl $sample.cnv.filter $sample.$gender.$karotype.CNV.AutoFilter.xlsx $gender $karotype`;
		`python /disk2/mygeno/yangh/lccnv/script/doc.py $sample.ratio.png $sample.$gender.$karotype.draft.docx`;
	}else{
		for my $sub_gender ('M','F'){
			###&mysql_inject("生成ratio文件，性别:$gender","正在运行",40);
			`perl $PSCC/Script/Ratio.pl -window /disk2/mygeno/yangh/lccnv/script/Window/0.6X/winRatio.Tags25k.Slide1k.Hg19.gz -ext $sample.ext -gender $sub_gender -o $sample.ratio`;
			###&mysql_inject("并行:Segmentation/SVG核型图","正在运行",50);
			&parallel(
				"perl $PSCC/Script/Segmentation.pl -w 400 -in $sample.ratio -o $sample.cnv -c $PSCC/Attachment/hg19.length -n $PSCC/Attachment/hg19.N_region -log $sample.ratio.cnv.log",
				"perl /disk2/mygeno/web/script/batchmendelian/svg_v3.pl $sample.ratio $sub_gender",
				"perl /disk2/mygeno/yangh/lccnv/script/svg_v4.pl $sample.ratio $sub_gender $sample.downsample.point"
			);
			###&mysql_inject("CNVFilterD","正在运行",60);
			`perl $PSCC/Script/CNVFilterD.pl -i $sample.cnv -g $sub_gender -o $sample.cnv.filter -cyto /disk1/DB/annovar_db//humandb/hg19_cytoBand.txt -ref /disk1/DB/annovar_db//humandb/hg19_refGene.txt -sp 0.05 -pp 0.05`;
			`perl $PSCC/cnv-autofilter.pl $sample.cnv.filter $sample.$sub_gender.$karotype.CNV.AutoFilter.xlsx $sub_gender $karotype`;
			`python /disk2/mygeno/yangh/lccnv/script/doc.py $sample.ratio.png $sample.$sub_gender.$karotype.draft.docx`;
		}
	}
}

sub parallel{
	my $pm = new Parallel::ForkManager(7);
	for (@_){
		my $pid = $pm->start and next;
		&do($_);
		$pm->finish;
	}
	$pm->wait_all_children;
	sub do{
		`$_[0]`;
	}
}

sub intervals{
	my $initial_intervals = $panel;
	$initial_intervals =~ s/MitoChip,//g;
	$initial_intervals =~ s/,MitoChip//g;
	$initial_intervals =~ s/,CNV//g;
	$initial_intervals =~ s/^CNV,//g;
	my $intervals = join " ",map {$_ = "-L /disk1/bed/$_.bed"} grep {$_ ne ""} (split /,/,$initial_intervals);
	return $intervals;
}

use strict;
use warnings;
use feature 'say';
use Data::Dumper;

my %h;
while(<DATA>){
    chomp;
    next unless $_;
    next if /^#/;
    my $depth;
    my ($sample,$lane,$chip) = split /\s+/;
    my $flg = undef;
    my ($r1,$r2) = (undef,undef);
    for my $path ("/disk1/data/MG-N111/17O003184_OC017_CapNGS","/Data/FQDATA4/","/Data/FQDATA3/","/Data/FQDATA2/","/Data/FQDATA/","/Data/FASTQ/","/Data/FQDATA5/tech_211","/Data_tmp/FASTQ/"){#/Data/FASTQ
        if (-d "$path/$lane"){
            if(-e "$path/$lane/$sample/${sample}_R1.fq.gz" && -e "$path/$lane/$sample/${sample}_R2.fq.gz" && -s "$path/$lane/$sample/${sample}_R1.fq.gz" > 0 && -s "$path/$lane/$sample/${sample}_R2.fq.gz"> 0){
                $r1 = "$path/$lane/$sample/${sample}_R1.fq.gz";
                $r2 = "$path/$lane/$sample/${sample}_R2.fq.gz";
                $flg = 1;last;
            }
            if(-d "$path/$lane/Sample_$sample"){
                chomp(my $tmp1 = `ls $path/$lane/Sample_$sample/*_R1_001.fastq.gz 2>/dev/null`);
                chomp(my $tmp2 = `ls $path/$lane/Sample_$sample/*_R2_001.fastq.gz 2>/dev/null`);
                if(-e $tmp1 && -e $tmp2 && -s $tmp1 > 0 && -s $tmp2 > 0){
                    $r1 = $tmp1;$r2 = $tmp2;
                    $flg = 1;last;
                }
            }
            chomp(my $tmp1 = `ls $path/$lane/$sample*_R1_001.fastq.gz 2>/dev/null`);
            chomp(my $tmp2 = `ls $path/$lane/$sample*_R2_001.fastq.gz 2>/dev/null`);
            if(-e $tmp1 && -e $tmp2 && -s $tmp1 > 0 && -s $tmp2 > 0){
                $r1 = $tmp1;$r2 = $tmp2;
                $flg = 1;last;
            }
            chomp($tmp1 = `ls $path/$lane/$sample\_R1.fq.gz 2>/dev/null`);
            chomp($tmp2 = `ls $path/$lane/$sample\_R2.fq.gz 2>/dev/null`);
            if(-e $tmp1 && -e $tmp2 && -s $tmp1 > 0 && -s $tmp2 > 0){
                $r1 = $tmp1;$r2 = $tmp2;
                $flg = 1;last;
            }
        }
    }
    unless ($flg){
        say STDERR "找不到 $sample 样本的fastq文件";
        next;
    }
    if($r1 && $r2){
        $h{$sample}{'R1'} = $r1;
        $h{$sample}{'R2'} = $r2;
        $h{$sample}{'lane'} = $lane;
    }
    open S,">$sample.sh";
    say S "cd /ssd1/mygeno/pengj/project;perl /disk2/mygeno/pengj/pipline/normal.pl -s $sample -r1 $r1 -r2 $r2 -c $chip\n\n";
    close S;
    print my $rst = `qsub -q all.q -l mem=15gb,walltime=100:00:00,nodes=1:ppn=1 $sample.sh`;
}
#print Dumper \%h;
#lane编号  样本编号   性别（F-女性，M-男性，不知道性别用 - 代替）
__DATA__
17C012763_H001_CapNGS MG-L626 MG_BIK_V1
#17C001237_MedE003_CapNGS    MG-L645  MG_One_V2
# zy-4-352    MG-N082 PKU
# zy-4-403    MG-N082 PKU
#17C016971_KY019_CapNGS  MG-L581 PKU
#17C016970_KY019_CapNGS  MG-L581 PKU
# 17C025548_KNGS003_CapNGS    MG-L615 PKU
# 17C025549_KNGS003_CapNGS    MG-L615 PKU
# 17C025547_KNGS003_CapNGS    MG-L615 PKU
# 17C025546_KNGS003_CapNGS    MG-L615 PKU
# 17C025545_KNGS003_CapNGS    MG-L615 PKU
# 17C025544_KNGS003_CapNGS    MG-L615 PKU
# 17C025543_KNGS003_CapNGS    MG-L615 PKU
# 17C025542_KNGS003_CapNGS    MG-L615 PKU
# 17C025541_KNGS003_CapNGS    MG-L615 PKU
# 17C025540_KNGS003_CapNGS    MG-L615 PKU
# 17C025539_KNGS003_CapNGS    MG-L615 PKU
# 17C025538_KNGS003_CapNGS    MG-L615 PKU
# 17C025537_KNGS003_CapNGS    MG-L615 PKU
# 17C025536_KNGS003_CapNGS    MG-L615 PKU
# 17C025535_KNGS003_CapNGS    MG-L615 PKU
# 17C025534_KNGS003_CapNGS    MG-L615 PKU

use DBI;
use strict;
use warnings;
use feature 'say';
use Term::ANSIColor;
use Parallel::ForkManager;
use Getopt::Long qw/GetOptions/;

my ($sample,$threads,$r1,$r2,$panel,$lane,$pbs_id,$random,$email,$gender);
GetOptions(
	's=s' => \$sample,
	'r1=s' => \$r1,
	'r2=s' => \$r2,
	'c=s' => \$panel,
	'l=s' => \$lane,
	'r=s' => \$random,
	'e=s' => \$email,
	'g=s' => \$gender,
	't=s' => \$threads,
);
my $usage =<<EOF;

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

程序简介：      遗传病分析流程，在线版
使用方法：      perl target_pipeline_web.pl arg1 arg2 ...
更新日期：      2017年5月3日

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EOF
unless($sample && $r1 && $r2 && $lane && $random && $email){
        say colored($usage, "yellow on_magenta");
        exit 0;
}

#$threads = 7;
chdir "/disk2/mygeno/web/batchmendelian/web-server/$random/$lane/$sample";
my $panel_for_mail = $panel;
my $tans='/TNAS3/tech';
my $intervals = &intervals();
chomp($pbs_id = `cat ${lane}_${sample}.pbsid`);
my $start_t = time();
&initial();

&mysql_inject("创建目录","正在运行");
&check_ssd1_TNAS2;
unless(-e "/$tans/${lane}-web"){`mkdir -p "/$tans/${lane}-web"`}
unless(-e "/$tans/${lane}-web/BAM"){`mkdir -p "/$tans/${lane}-web/BAM"`}
unless(-e "/$tans/${lane}-web/Regioncount"){`mkdir -p "/$tans/${lane}-web/Regioncount"`}
unless(-e "/$tans/${lane}-web/Siteplot"){`mkdir -p "/$tans/${lane}-web/Siteplot"`}
unless(-e "/$tans/${lane}-web/SNP_INDEL"){`mkdir -p "/$tans/${lane}-web/SNP_INDEL"`}
unless(-e "/$tans/${lane}-web/Statistics"){`mkdir -p "/$tans/${lane}-web/Statistics"`}
unless(-e "/$tans/${lane}-web/SV"){`mkdir -p "/$tans/${lane}-web/SV"`}
unless(-e "/$tans/${lane}-web/VCF"){`mkdir -p "/$tans/${lane}-web/VCF"`}
unless(-e "/$tans/${lane}-web/mtDNA_HotSpot"){`mkdir -p "/$tans/${lane}-web/mtDNA_HotSpot"`}
unless(-e "/$tans/${lane}-web/CNV"){`mkdir -p "/$tans/${lane}-web/CNV"`}
unless(-e "/$tans/${lane}-web/GFF"){`mkdir -p "/$tans/${lane}-web/GFF"`}
unless(-e "/$tans/${lane}-web/RawDataSelectedGene"){`mkdir -p "/$tans/${lane}-web/RawDataSelectedGene"`}
unless(-e "/$tans/${lane}-web/这个目录是由网站自动产生的请不要删除"){`touch "/$tans/${lane}-web/这个目录是由网站自动产生的请不要删除"`}
unless(-e "/$tans/${lane}-web/网站正在分析这批数据请不要占用文件"){`touch "/$tans/${lane}-web/网站正在分析这批数据请不要占用文件"`}

if($panel =~ /MitoChip/){
	&mtDNA;
	$panel =~ s/MitoChip[,]?//g;
	$panel =~ s/,MitoChip//g;
}
if($panel =~ /,CNV|^CNV,|^CNV$/){
	&CNV;
	$panel =~ s/,CNV//g;
	$panel =~ s/^CNV[,]?//g;
}
if($panel ne ''){
	&normal;
}
&mysql_mail_manager;

sub normal{
	my $cutadapt = '/bin/cutadapt';
	my $ref = '/local_disk/DB/hg19/Sequence/BWAIndex/genome.fa';
	my $ref_fa = '/local_disk/DB/hg19/Sequence/WholeGenomeFasta/genome.fa';
	my $picard = '/disk1/software/picard-tools-2.2.3/picard.jar';
	my $bamtools = '/disk1/software/bin/bamtools';
	my $bwa = '/disk1/software/bin/bwa';
	my $samtools = '/disk1/software/bin/samtools';
	my $bamToBed = '/disk1/software/bin/bamToBed';
	my $coverageBed = '/disk1/software/bin/coverageBed';
	my $GATK = '/disk1/software/GATK-3.7/GenomeAnalysisTK.jar';
	my $dbsnp = '/local_disk/DB/dbsnp/dbsnp_147.hg19.vcf';
	my $varscan  = '/disk1/software/VarScan.v2.3.7.jar';
	my $beddir = '/disk1/bed/';
	my $local_disk = "/local_disk/tmp/web/batchmendelian/$random/$lane/$sample";	
	my $disk = "/ssd1/mygeno/web/batchmendelian/$random/$lane/$sample";
	my $java = "java -Djava.io.tmpdir=$local_disk -Xmx15g -Xms5g -jar";
	unless(-e $disk){`mkdir -p $disk`;}
	unless(-e $local_disk){`mkdir -p $local_disk`;}
	chdir $disk;
	&mysql_inject("拷贝fastq文件","正在运行",5);
	if ($r1=~/,/ && $r2=~/,/){
		my @reads1=split ",", $r1;
		my @reads2=split ",", $r2;
		for(my $i=0; $i<@reads1; $i++){
			&parallel(
				"cat $reads1[$i] >> ${sample}_R1.fq.gz",
				"cat $reads2[$i] >> ${sample}_R2.fq.gz"
			);
		}
	} else {
		&parallel(
			"cp $r1 ${sample}_R1.fq.gz",
			"cp $r2 ${sample}_R2.fq.gz"
		);
	}
	&mysql_inject("fastq文件质控","正在运行",10);
	chomp(my $hardware = `zcat ${sample}_R1.fq.gz | head -1 | cut -d: -f1`);
	if($hardware =~ /NS500/){
		`$cutadapt -m 80 -q10,10 --nextseq-trim=10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
	}else{
		`$cutadapt -m 80 -q10,10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
		#`$cutadapt -m 20 -q10,10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
	}
	&mysql_inject("将fastq比对到人类参考基因组上","正在运行",20);
	`$bwa mem -M -t $threads $ref $sample.clean_R1.fq.gz $sample.clean_R2.fq.gz | $samtools fixmate -O bam - - | $samtools sort -\@ $threads -m 1G - $sample.sort`;
	&check_ssd1_TNAS2;
	`$samtools index $sample.sort.bam > $sample.sort.bam.bai`;
	&mysql_inject("添加BAM文件头信息","正在运行",30);
	`$java $picard AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=$sample.sort.bam OUTPUT=$sample.sort.header.bam RGLB=$ref RGPL=ILLUMINA RGSM=GP1 RGPU=GRP1 SORT_ORDER=coordinate CREATE_INDEX=true`;
	&mysql_inject("过滤BAM文件","正在运行",35);
	`$bamtools filter -isMapped true -isPaired true -isProperPair true -in $sample.sort.header.bam -out $sample.sort.flt.bam`;
	`rm -f $sample.sort.header.bam $sample.sort.header.bai`;
	if ($panel =~ /DMD|Clinican|Canmute|SLC26A4|Deaf_V2|lung8|LAMA2|PKU|LBS/){
		my @parallel_command=(
			"$java $picard MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.rmdup.sorted.bam M=duplication-report2.txt",
		);
		if($panel =~ /PKU/){
			push @parallel_command,"perl /disk1/software/crest/extractSClip.pl -i $sample.sort.bam --ref_genome $ref -r chr6 > $sample.extractSClip.log 2>&1";
		}else{
			for my $chr (1 .. 22,'X','Y'){
				push @parallel_command,"perl /disk1/software/crest/extractSClip.pl -i $sample.sort.bam --ref_genome $ref -r chr$chr > $sample.extractSClip.log 2>&1";
			}
		}
		&mysql_inject("去除PCR冗余序列/提取soft-clip reads","正在运行",40);
		&parallel(@parallel_command);
		&mysql_inject("合并$sample.sort.bam.*.cover/$sample.sort.bam.*.sclip.txt","正在运行",45);
		`cat $sample.sort.bam.*.cover > $sample.sort.bam.cover`;
		`cat $sample.sort.bam.*.sclip.txt > $sample.sort.bam.sclip.txt`;
	}else{
		&mysql_inject("去除PCR冗余序列","正在运行",40);
		`$java $picard MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.rmdup.sorted.bam M=duplication-report2.txt`;
	}
	&mysql_inject("标记PCR冗余序列/拷贝BAM文件到公共盘/bamtools stats/bamToBed","正在运行",45);
	&parallel(
		"$java $picard MarkDuplicates VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.sort.flt.mark.bam M=duplication-report.txt",
		"cp $sample.rmdup.sorted.bam /$tans/${lane}-web/BAM/",
		"cp $sample.rmdup.sorted.bai /$tans/${lane}-web/BAM/$sample.rmdup.sorted.bam.bai",
		"$bamtools stats -in $sample.sort.bam > $sample.sort.bam.stat",
		"$bamToBed -i $sample.rmdup.sorted.bam > $sample.rmdup.sorted.bamtobed",
	);
	`rm -f $sample.sort.flt.ba?`;
	if ($panel =~ /DMD|Clinican|Canmute|SLC26A4|Deaf_V2|lung8|LAMA2|PKU|LBS/){
		my @parallel_command=(
			"$java $GATK -T BaseRecalibrator -R $ref_fa -I $sample.sort.flt.mark.bam -knownSites $dbsnp -o recal.table -nct $threads $intervals"
		);
		if ($panel =~ /PKU/ || $panel =~ /LAMA2/){
			push @parallel_command,"perl /disk1/software/crest/CREST.pl -l 150 -f $sample.sort.bam.cover -d $sample.sort.bam --ref_genome $ref -t /local_disk/DB/hg19/Sequence/blat/genome.2bit -r chr6 -p $sample.sort.bam.chr6 > $sample.chr6.CREST.log 2>&1";
		}else{
			for my $chr (1 .. 22,'X','Y'){
				push @parallel_command,"perl /disk1/software/crest/CREST.pl -l 150 -f $sample.sort.bam.cover -d $sample.sort.bam --ref_genome $ref -t /local_disk/DB/hg19/Sequence/blat/genome.2bit -r chr$chr -p $sample.sort.bam.chr$chr > $sample.chr$chr.CREST.log 2>&1";
			}
		}			
		if($panel =~ /DMD/){
			push @parallel_command,"$samtools depth -aa -r chrX:31135345-33359192 -d 100000 $sample.rmdup.sorted.bam > $sample.DMD.depth";
		}elsif($panel =~ /SLC26A4/){
		    push @parallel_command,"$samtools depth -aa -r chr7:107301080-107358252 -d 100000 $sample.rmdup.sorted.bam > $sample.SLC26A4.depth";
		}elsif($panel =~ /PKU/ || $panel =~ /LAMA2/){
			push @parallel_command,"$samtools depth -aa -r chr6:129204100-129838000 -d 100000 $sample.rmdup.sorted.bam > $sample.LAMA2.depth";
		}
		&mysql_inject("并行:校正平台误差/检测融合基因","正在运行",50);
		&parallel(@parallel_command);
		`cat *.chr*.predSV.txt > $sample.sort.bam.predSV.txt`;
		`perl /disk2/mygeno/web/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt`;
		if($panel =~ /DMD/){
			&parallel(
				"perl /disk2/mygeno/web/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt",
				"perl /disk2/mygeno/web/script/batchmendelian/siteplot_dmd.pl $sample.DMD.depth"
			);
		}elsif($panel =~ /SLC26A4/){
			&parallel(
				"perl /disk2/mygeno/web/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt",
				"perl /disk2/mygeno/web/script/batchmendelian/siteplot_SLC26A4.pl $sample.SLC26A4.depth SLC26A4 SLC26A4 chr7:107301080-107358252"
			);
		}elsif($panel =~ /PKU/ || $panel =~ /LAMA2/){
			&parallel(
                "perl /disk2/mygeno/web/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt",
                "perl /disk2/mygeno/web/script/batchmendelian/siteplot_LAMA2.pl $sample.LAMA2.depth"
            );
		}
		
		if($panel =~ /DMD/){
			&parallel(
				"cp *.pdf /$tans/${lane}-web/Siteplot/",
				"cp *.png /$tans/${lane}-web/Siteplot/",
				"cp $sample.smoothy_depth.txt /$tans/${lane}-web/Siteplot/$sample.DMD.depth.txt"
			);
		}elsif($panel =~ /SLC26A4/){
			&parallel(
				"cp *.pdf /$tans/${lane}-web/Siteplot/",
				"cp *.png /$tans/${lane}-web/Siteplot/",
				"cp $sample.SLC26A4.depth /$tans/${lane}-web/Siteplot/$sample.SLC26A4.depth.txt"
			);
		}elsif($panel =~ /PKU/ || $panel =~ /LAMA2/){
			&parallel(
                "cp *.pdf /$tans/${lane}-web/Siteplot/",
                "cp *.png /$tans/${lane}-web/Siteplot/",
                "cp $sample.LAMA2.depth /$tans/${lane}-web/Siteplot/$sample.LAMA2.depth.txt"
            );
		}
		if (glob "${sample}*crest.xls"){
			`cp ${sample}*crest.xls "/$tans/${lane}-web/SV"`;
		}
	}else{
		&mysql_inject("校正平台误差","正在运行",55);
		`$java $GATK -T BaseRecalibrator -R $ref_fa -I $sample.sort.flt.mark.bam -knownSites $dbsnp -o recal.table -nct $threads $intervals`;
	}
	&mysql_inject("输出高质量比对序列","正在运行",60);
	`$java $GATK -T PrintReads -R $ref_fa -I $sample.sort.flt.mark.bam -BQSR recal.table -o $sample.recal.bam -nct $threads`;
	`rm -f $sample.sort.flt.mark.ba?`;
	`rm -f $sample.rmdup.sorted.ba?`;
	for my $sub_chip ((split /,/,$panel)){
		next if $sub_chip eq '';
		if($sub_chip =~ /Wscancer|Canmute|Clinican|BloodV2|Immu_All|ImmuV2|MetFull|MetAdd|geneis|Can06|Lym_MyGeno|lung8/){
			#MPIleup+GATK
			&mysql_inject("并行:mpileup低频突变/HaplotypeCaller常规突变/coverageBed","正在运行",65);
			&parallel(
				"$samtools mpileup -s -f $ref -l $beddir/$sub_chip.bed -d 10000 -L 10000 $sample.recal.bam | $java $varscan mpileup2cns --min-reads2 1 --min-coverage 1 --strand-filter 1 --output-vcf 1 --variants 1 --min-var-freq 0.0001 --p-value 1 --min-avg-qual 20 > $sample.varscan.$sub_chip.vcf",
				"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/$sub_chip.bed -ip 100 -D $dbsnp -o $sample.gatk.$sub_chip.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
				"$coverageBed -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.coverage",
				"$coverageBed -d -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.target.coverage"
			);
			&mysql_inject("捕获效率统计","正在运行",65);
			`perl /disk2/mygeno/web/script/batchmendelian/depth_count.pl $sample.sort.$sub_chip.coverage $sample.sort.$sub_chip.target.coverage > $sample.$sub_chip.sort.bam.depth.statistics.xls`;
			`perl /disk2/mygeno/web/script/batchmendelian/read_count.pl $disk/$sample.$sub_chip.sort.bam.depth.statistics.xls $sub_chip`;
			`perl /disk2/mygeno/web/script/batchmendelian/vcf2vcf.pl $sample.gatk.$sub_chip.vcf > $sample.gatk.$sub_chip.tmp.vcf`;
			`perl /disk2/mygeno/web/script/batchmendelian/merge_vcf_for_gatk_varscan.pl $sample.gatk.$sub_chip.tmp.vcf $sample.varscan.$sub_chip.vcf > $sample.gatk.$sub_chip.split.vcf`;
			`mv $sample.varscan.$sub_chip.vcf /$tans/${lane}-web/VCF/$sample.varscan.$sub_chip.vcf`;
		}else{
			#GATK
			if ($sub_chip eq 'Deaf_V2'){
				&mysql_inject("并行:HaplotypeCaller(常规区域和miRNA区域)/coverageBed","正在运行",65);
				&parallel(
					"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/$sub_chip.bed -ip 100 -D $dbsnp -o $sample.gatk.$sub_chip.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
					"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/${sub_chip}_miRNA_region.bed -ip 100 -D $dbsnp -o $sample.gatk.${sub_chip}.miRNA.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
					"$coverageBed -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.coverage",
					"$coverageBed -d -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.target.coverage"
				);
				&parallel(
					"mv $sample.gatk.${sub_chip}.miRNA.vcf /$tans/${lane}-web/VCF/$sample.${sub_chip}.miRNA.vcf",
					"perl /disk2/mygeno/web/script/batchmendelian/depth_count.pl $sample.sort.$sub_chip.coverage $sample.sort.$sub_chip.target.coverage > $sample.$sub_chip.sort.bam.depth.statistics.xls"
				);
				`perl /disk2/mygeno/web/script/batchmendelian/read_count.pl $disk/$sample.$sub_chip.sort.bam.depth.statistics.xls $sub_chip`;
			}else{
				&mysql_inject("并行:HaplotypeCaller/coverageBed","正在运行",65);
				&parallel(
					"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/$sub_chip.bed -ip 100 -D $dbsnp -o $sample.gatk.$sub_chip.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
					"$coverageBed -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.coverage",
					"$coverageBed -d -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.target.coverage"
				);
				`perl /disk2/mygeno/web/script/batchmendelian/depth_count.pl $sample.sort.$sub_chip.coverage $sample.sort.$sub_chip.target.coverage > $sample.$sub_chip.sort.bam.depth.statistics.xls`;
				`perl /disk2/mygeno/web/script/batchmendelian/read_count.pl $disk/$sample.$sub_chip.sort.bam.depth.statistics.xls $sub_chip`;
			}
			#不需要mpileup的通用
			&mysql_inject("split vcf","正在运行",60);
			`perl /disk2/mygeno/web/script/batchmendelian/vcf2vcf.pl $sample.gatk.$sub_chip.vcf > $sample.gatk.$sub_chip.split.vcf`;
		}
		#所有Panel通用
		&mysql_inject("并行:ANNOVAR功能注释/拷贝统计文件","正在运行",70);
		&parallel(
			"perl /disk1/software/annovar/table_annovar.pl $sample.gatk.$sub_chip.split.vcf /local_disk/DB/annovar_db/humandb/ -buildver hg19 -protocol refGene,cytoBand,dbNsfpInterPro,Inhouse,esp6500siv2_all,1000g2015aug_all,exac03,avsnp147,cosmic74,dbnsfp30a,clinvar_20160302,spidex,mcap10,revel,dbscsnv11,dbnsfp31a_interpro,gnomad_exome,gnomad_genome,InterVar -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -argument \"-splicing_threshold 5 --hgvs,,,,,,,,,,,,,,,,,,\" -vcfinput -nastring . --thread $threads --maxgenethread $threads",
			"perl /disk2/mygeno/web/script/batchmendelian/depth_count3.pl $sample.sort.$sub_chip.target.coverage $sample $sub_chip",
			"mv $sample.$sub_chip.$sub_chip\_readcount.xls /$tans/${lane}-web/Regioncount/$sample.$sub_chip\_readcount.xls",
		);
		&mysql_inject("并行:vcf2gff/拷贝捕获效率统计文件","正在运行",80);
		&parallel(
			"perl /disk2/mygeno/web/script/batchmendelian/ParseVCFHgvs.pl -v $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf -f $sample.gatk.$sub_chip.split.vcf.refGene.exonic_variant_function -g $gender > $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff",
			"mv $sample.${sub_chip}.count.xls /$tans/${lane}-web/Statistics/$sample.${sub_chip}.Statistics.xls"
		);
		&mysql_inject("关联HGMD数据库","正在运行",90);
		`perl /disk2/mygeno/web/script/batchmendelian/HGMD.pl $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff $sample $sub_chip`;
		&mysql_inject("拷贝SNP/INDEL文件","正在运行",95);
		if($sub_chip eq 'MR_MUT-MR_CNV'){
			if(-e "/$tans/${lane}-web/Statistics/$sample.MR_MUT-MR_CNV.Statistics.xls"){
				`rm -f /$tans/${lane}-web/Statistics/$sample.MR_MUT-MR_CNV.Statistics.xls`;
			}
			if(-e "/$tans/${lane}-web/SNP_INDEL/${sample}.MR_MUT-MR_CNV.SNP_INDEL.xls"){
				`rm -f /$tans/${lane}-web/SNP_INDEL/${sample}.MR_MUT-MR_CNV.SNP_INDEL.xls`;
			}
			goto NEXT;
		}
		&mysql_inject("并行:拷贝经过注释和过滤的xls/拷贝原始未经过滤的gff/备份vcf和gff文件","正在运行",98);
		&parallel(
			"mv $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff.hgmd.xls /$tans/${lane}-web/SNP_INDEL/${sample}.${sub_chip}.SNP_INDEL.xls",
			"mv $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff /$tans/${lane}-web/GFF/$sample.$sub_chip.gff",
			"mv $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf /$tans/${lane}-web/VCF/$sample.$sub_chip.vcf"
		);
		&mysql_inject("将捕获效率统计信息追加到Statistics.xls","正在运行",99);
		`ls /$tans/${lane}-web/Statistics/*.xls | xargs -i sed 1d {} >> Statistics.xls`;
		`sed -i '1iSample\tRaw_data_bases(Mb)\tClean_data_bases(Mb)\tAligned_bases(Mb)\tAligned\tInitial bases on target\tBase covered on target\tCoverage of target region\tEffective bases on target\tFraction of effective bases on target\tAverage sequencing depth on target\tFraction of target covered with at least 4X\tFraction of target covered with at least 10X\tFraction of target covered with at least 20X\tduplication rate\tchip' Statistics.xls`;
		`mv Statistics.xls /$tans/${lane}-web/Statistics.xls`;
NEXT:
	}
	&mysql_inject("清理空间","正在运行",99);
	`rm -rf /local_disk/tmp/web/batchmendelian/$random/$lane/$sample`;
	`rm -rf *`;
	&mysql_inject("该任务已运行完成","已完成",100);
}

sub mtDNA{
	my $cutadapt = '/usr/bin/cutadapt';
	my $ref = '/local_disk/DB/hg38/chrM/chrM.fa';
	my $picard = '/disk1/software/picard-tools-2.2.3/picard.jar';
	my $bamtools = '/disk1/software/bin/bamtools';
	my $bwa = '/disk1/software/bin/bwa';
	my $samtools = '/disk1/software/bin/samtools';
	my $bamToBed = '/disk1/software/bin/bamToBed';
	my $coverageBed = '/disk1/software/bin/coverageBed';
	my $GATK = '/disk1/software/GATK-3.6/GenomeAnalysisTK.jar';
	my $dbsnp = '/local_disk/DB/dbsnp/dbsnp_144.chrM.hg38.vcf';
	my $beddir = '/disk1/bed/';
	my $local_disk = "/local_disk/tmp/web/batchmendelian/$random/$lane/MitoChip/$sample";
	my $disk = "/ssd1/mygeno/web/batchmendelian/$random/$lane/MitoChip/$sample";
	my $java = "java -Djava.io.tmpdir=$local_disk -Xmx15g -Xms5g -jar";
	unless(-e $disk){`mkdir -p $disk`;}
	unless(-e $local_disk){`mkdir -p $local_disk`;}
	chdir $disk;
	&mysql_inject("并行：拷贝fastq文件","正在运行",5);
	if ($r1=~/,/ && $r2=~/,/){
        my @reads1=split ",", $r1;
        my @reads2=split ",", $r2;
        for(my $i=0; $i<@reads1; $i++){
            &parallel(
                "cat $reads1[$i] >> ${sample}_R1.fq.gz",
                "cat $reads2[$i] >> ${sample}_R2.fq.gz"
            );
        }
    } else {
        &parallel(
            "cp $r1 ${sample}_R1.fq.gz",
            "cp $r2 ${sample}_R2.fq.gz"
        );
    }	
	&mysql_inject("去除低质量Reads(连续N序列、低于80bp序列、引物自连序列)","正在运行",10);
        chomp(my $hardware = `zcat ${sample}_R1.fq.gz | head -1 | cut -d: -f1`);
        if($hardware =~ /NS500/){
                `$cutadapt -m 80 -q10,10 --nextseq-trim=10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a ATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
        }else{
                `$cutadapt -m 80 -q10,10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a ATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
        }
	&mysql_inject("BWA mem(比对到UCSC hg38)","正在运行",20);
	`$bwa mem -M -t $threads $ref $sample.clean_R1.fq.gz $sample.clean_R2.fq.gz | $samtools fixmate -O bam - - | $samtools sort -\@ $threads -m 1G - $sample.sort`;
	&check_ssd1_TNAS2;
	&mysql_inject("samtools index","正在运行",22);
	`$samtools index $sample.sort.bam > $sample.sort.bam.bai`;
	&mysql_inject("AddOrReplaceReadGroups","正在运行",25);
	`$java $picard AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=$sample.sort.bam OUTPUT=$sample.sort.header.bam RGLB=$ref RGPL=ILLUMINA RGSM=GP1 RGPU=GRP1 SORT_ORDER=coordinate CREATE_INDEX=true`;
	&mysql_inject("bamtools filter","正在运行",30);
	`$bamtools filter -isMapped true -isPaired true -isProperPair true -in $sample.sort.header.bam -out $sample.sort.flt.bam`;
	&mysql_inject("samtools index","正在运行",30);
	`$samtools index $sample.sort.flt.bam > $sample.sort.flt.bam.bai`;
	&mysql_inject("picard MarkDuplicates(去除冗余)","正在运行",35);
	`$java $picard MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.rmdup.sorted.bam M=duplication-report2.txt`;
	&mysql_inject("并行:picard MarkDuplicates(标记冗余)/拷贝bam文件到公共盘/bamtools stats/bamToBed","正在运行",40);
	&parallel(
		"$java $picard MarkDuplicates VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.sort.flt.mark.bam M=duplication-report.txt",
		"cp $sample.rmdup.sorted.bam /$tans/${lane}-web/BAM/$sample.MitoChip.rmdup.sorted.bam",
		"cp $sample.rmdup.sorted.bai /$tans/${lane}-web/BAM/$sample.MitoChip.rmdup.sorted.bai",
		"$bamtools stats -in $sample.sort.bam > $sample.sort.bam.stat",
		"$bamToBed -i $sample.rmdup.sorted.bam > $sample.rmdup.sorted.bamtobed",
	);
	&mysql_inject("并行:BaseRecalibrator","正在运行",45);
	&parallel(
		"$java $GATK -T BaseRecalibrator -R $ref -I $sample.sort.flt.mark.bam -knownSites $dbsnp -o recal.table -nct $threads -L chrM:1-16569",
		"$coverageBed -a $sample.rmdup.sorted.bamtobed -b $beddir/hg38_MT_wholegenome.bed > $sample.sort.MitoChip.coverage",
		"$coverageBed -d -a $sample.rmdup.sorted.bamtobed -b $beddir/hg38_MT_wholegenome.bed > $sample.sort.MitoChip.target.coverage",
	);
	&mysql_inject("并行:PrintReads","正在运行",50);
	&parallel(
		"$java $GATK -T PrintReads -R $ref -I $sample.sort.flt.mark.bam -BQSR recal.table -o $sample.recal.bam -nct $threads",
		"perl /disk2/mygeno/web/script/batchmendelian/depth_count.pl $sample.sort.MitoChip.coverage $sample.sort.MitoChip.target.coverage > $sample.sort.bam.depth.statistics.xls",
	);
	&mysql_inject("并行:UnifiedGenotyper/samtools depth/统计捕获效率","正在运行",55);
	&parallel(
		"$java $GATK -T UnifiedGenotyper -R $ref -I $sample.recal.bam -L chrM:1-16569 -D $dbsnp -o $sample.gatk.raw.vcf -A Coverage -A HaplotypeScore --downsample_to_coverage 1000000 -nct $threads -glm BOTH",
		"$samtools depth -aa -r chrM:1-16569 -d 100000 $sample.recal.bam > $sample.recal.bam.depth",
		"perl /disk2/mygeno/web/script/batchmendelian/depth_count3.pl $sample.sort.MitoChip.target.coverage $sample MitoChip",
	);
	&mysql_inject("SitePlot","正在运行",60);
	`perl /disk2/mygeno/web/script/batchmendelian/smooth_mt_depth_and_plot_r.pl $sample.recal.bam.depth`;
	&mysql_inject("并行:拷贝一些结果文件/split VCF","正在运行",68);
	&parallel(
		"cp $sample\_whole_mtDNA_length_depth.pdf /$tans/${lane}-web/Siteplot",
		"cp $sample\_whole_mtDNA_length_depth.png /$tans/${lane}-web/Siteplot",
		"cp $sample.recal.bam.depth /$tans/${lane}-web/Siteplot/$sample.MitoChip.depth.txt",
		"perl /disk2/mygeno/web/script/batchmendelian/split_mtdna_vcf.pl $sample.gatk.raw.vcf",
		"cp $sample.gatk.raw.vcf /$tans/${lane}-web/VCF/$sample.MitoChip.vcf"
	);
	&mysql_inject("ANNOVAR功能注释","正在运行",70);
	`perl /disk1/software/annovar/table_annovar.pl $sample.gatk.raw.vcf.split.vcf /local_disk/DB/annovar_db/humandb -buildver hg38_MT -protocol refGene,snp142NonFlagged -operation g,f -nastring . -vcfinput`;
	&mysql_inject("提取结果","正在运行",75);
	`perl /disk2/mygeno/web/script/batchmendelian/extraction_vcf_for_mtDNA.pl $sample.gatk.raw.vcf.split.vcf.hg38_MT_multianno.vcf`;
	&mysql_inject("将$sample.SNP_INDEL.xls追加到mtDNA_HotSpot.xls","正在运行",80);
	if(!-e "/$tans/${lane}-web/mtDNA_HotSpot/mtDNA_HotSpot.xls"){
		`cp $sample.mt.SNP_INDEL.xls "/$tans/${lane}-web/mtDNA_HotSpot/mtDNA_HotSpot.xls"`;
	}else{
		`sed 1d $sample.mt.SNP_INDEL.xls >> "/$tans/${lane}-web/mtDNA_HotSpot/mtDNA_HotSpot.xls"`;
	}
	&mysql_inject("并行:拷贝一些结果文件","正在运行",85);
	&parallel(
		"mv $sample.MitoChip.count.xls /$tans/${lane}-web/Statistics/$sample.MitoChip.Statistics.xls",
		"mv $sample.mt.SNP_INDEL.xls /$tans/${lane}-web/SNP_INDEL/$sample.MitoChip_SNP_INDEL.xls",
	);
	`ls /$tans/${lane}-web/Statistics/*.xls | xargs -i sed 1d {} >> Statistics.xls`;
	`sed -i '1iSample\tRaw_data_bases(Mb)\tClean_data_bases(Mb)\tAligned_bases(Mb)\tAligned\tInitial bases on target\tBase covered on target\tCoverage of target region\tEffective bases on target\tFraction of effective bases on target\tAverage sequencing depth on target\tFraction of target covered with at least 4X\tFraction of target covered with at least 10X\tFraction of target covered with at least 20X\tduplication rate\tchip' Statistics.xls`;
	`mv Statistics.xls /$tans/${lane}-web/Statistics.xls`;
	&mysql_inject("清理空间","正在运行",90);
	`rm -rf /local_disk/tmp/web/batchmendelian/$random/$lane/MitoChip/$sample`;
	`rm -rf *`;
	&mysql_inject("该任务已运行完成","已完成",100);
}

sub CNV{
	my $cutadapt = '/usr/bin/cutadapt';
	my $samtools = '/disk1/software/bin/samtools';
	my $ref = '/local_disk/DB/hg19/Sequence/BWAIndex/genome.fa';
	my $PSCC = '/disk2/mygeno/yangh/lccnv/script/';
	my $local_disk = "/local_disk/tmp/web/batchmendelian/$random/$lane/CNV/$sample";
	my $disk = "/ssd1/mygeno/web/batchmendelian/$random/$lane/CNV/$sample";
	my $java = "java -Djava.io.tmpdir=$local_disk -Xmx15g -Xms5g -jar";
	my $picard = '/disk1/software/picard-tools-2.2.3/picard.jar';
	unless(-e $disk){`mkdir -p $disk`;}
	unless(-e $local_disk){`mkdir -p $local_disk`}
	chdir $disk;
	&mysql_inject("拷贝fastq文件","正在运行",5);
	if ($r1=~/,/ && $r2=~/,/){
        my @reads1=split ",", $r1;
        my @reads2=split ",", $r2;
        for(my $i=0; $i<@reads1; $i++){
            &parallel(
                "cat $reads1[$i] >> ${sample}_R1.fq.gz",
                "cat $reads2[$i] >> ${sample}_R2.fq.gz"
            );
        }
    } else {
        &parallel(
            "cp $r1 ${sample}_R1.fq.gz",
            "cp $r2 ${sample}_R2.fq.gz"
        );
    }
	&mysql_inject("去除低质量Reads(连续N序列、低于50bp序列、引物自连序列)","正在运行",10);
	chomp(my $hardware = `zcat ${sample}_R1.fq.gz | head -1 | cut -d: -f1`);
	if($hardware =~ /NS500/){
		`$cutadapt -m 50 -q10,10 --nextseq-trim=10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a ATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
	}else{
		`$cutadapt -m 50 -q10,10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a ATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
	}
	&mysql_inject("并行:BWA aln(比对到UCSC hg19)","正在运行",20);
	&parallel(
		"/disk1/software/bin/bwa aln -t 4 $ref $sample.clean_R1.fq.gz > $sample.clean_R1.sai",
		"/disk1/software/bin/bwa aln -t 4 $ref $sample.clean_R2.fq.gz > $sample.clean_R2.sai"
	);
	&mysql_inject("BWA sampe/fixmate","正在运行",25);
	`bwa sampe $ref $sample.clean_R1.sai $sample.clean_R2.sai $sample.clean_R1.fq.gz $sample.clean_R2.fq.gz > $sample.sam`;
	`$samtools fixmate -O bam $sample.sam -| $samtools sort -\@ 5 -m 1G - $sample.sort`;
	&check_ssd1_TNAS2;
	`$samtools index $sample.sort.bam > $sample.sort.bam.bai`;
	&mysql_inject("检测插入片段长度分布","正在运行中",28);
	`$java $picard CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT I=$sample.sort.bam O=$sample.raw.insertsize.txt H=$sample.raw.insertsize.pdf M=0.5`;
	&mysql_inject("测序数据及捕获效率统计","正在运行",29);
	#`$java $picard CollectRawWgsMetrics VALIDATION_STRINGENCY=SILENT I=$sample.sort.bam O=raw_wgs_metrics.txt R=$ref INCLUDE_BQ_HISTOGRAM=true`;
	my $depth = `perl /disk2/mygeno/web/script/batchmendelian/cnv-rawfq-Statistics.pl $sample ${sample}_R1.fq.gz ${sample}_R2.fq.gz $sample.clean_R1.fq.gz $sample.clean_R2.fq.gz`;
	&mysql_inject("将sam文件转换为ext文件","正在运行",30);
	`perl $PSCC/Script/Soap2Ext.pl -soap $sample.sam -o $sample.ext -pe`;
	chomp(my $result = `perl /disk2/mygeno/yangh/lccnv/script/karyotype.pl $sample.sort.bam`);
	my ($gender,$karotype) = split /\t/,$result;
	if($gender eq 'M' || $gender eq 'F'){
		unless(-e "/$tans/${lane}-web/CNV/$sample/$gender"){
			`mkdir -p /$tans/${lane}-web/CNV/$sample/$gender`;
		}
		&mysql_inject("生成ratio文件，性别:$gender","正在运行",40);
		`perl $PSCC/Script/Ratio.pl -window /disk2/mygeno/yangh/lccnv/script/Window/0.6X/winRatio.Tags25k.Slide1k.Hg19.gz -ext $sample.ext -gender $gender -o $sample.ratio`;
		&mysql_inject("并行:Segmentation/SVG核型图","正在运行",50);
		&parallel(
			"perl $PSCC/Script/Segmentation.pl -w 400 -in $sample.ratio -o $sample.cnv -c $PSCC/Attachment/hg19.length -n $PSCC/Attachment/hg19.N_region -log $sample.ratio.cnv.log",
			"perl /disk2/mygeno/web/script/batchmendelian/svg_v3.pl $sample.ratio $gender",
			"perl /disk2/mygeno/yangh/lccnv/script/svg_v4.pl $sample.ratio $gender $sample.downsample.point"
		);
		&mysql_inject("CNVFilterD","正在运行",60);
		`perl $PSCC/Script/CNVFilterD.pl -i $sample.cnv -g $gender -o $sample.cnv.filter -cyto /disk1/DB/annovar_db//humandb/hg19_cytoBand.txt -ref /disk1/DB/annovar_db//humandb/hg19_refGene.txt -sp 0.05 -pp 0.05`;
		`perl $PSCC/cnv-autofilter.pl $sample.cnv.filter $sample.$gender.$karotype.CNV.AutoFilter.xlsx $gender $karotype`;
		`python /disk2/mygeno/yangh/lccnv/script/doc.py $sample.ratio.png $sample.$gender.$karotype.draft.docx`;
		&mysql_inject("并行:拷贝xls/png","任务即将完成",70);
		&parallel(
			"cp $sample.$gender.$karotype.CNV.AutoFilter.xlsx /$tans/${lane}-web/CNV/$sample/$gender/$sample.${depth}x.$gender.$karotype.CNV.AutoFilter.xlsx ",
			"cp $sample.ratio.png /$tans/${lane}-web/CNV/$sample/$gender/$sample.${depth}x.$gender.$karotype.CNV.png",
			"cp $sample.downsample.point.png /$tans/${lane}-web/CNV/$sample/$gender/$sample.${depth}x.$gender.$karotype.downsample.point.png",
			"cp $sample.$gender.$karotype.draft.docx /$tans/${lane}-web/CNV/$sample/$gender/"
		);
	}else{
		for my $sub_gender ('M','F'){
			unless(-e "/$tans/${lane}-web/CNV/$sample/$sub_gender"){
				`mkdir -p /$tans/${lane}-web/CNV/$sample/$sub_gender`;
			}
			&mysql_inject("生成ratio文件，性别:$gender","正在运行",40);
			`perl $PSCC/Script/Ratio.pl -window /disk2/mygeno/yangh/lccnv/script/Window/0.6X/winRatio.Tags25k.Slide1k.Hg19.gz -ext $sample.ext -gender $sub_gender -o $sample.ratio`;
			&mysql_inject("并行:Segmentation/SVG核型图","正在运行",50);
			&parallel(
				"perl $PSCC/Script/Segmentation.pl -w 400 -in $sample.ratio -o $sample.cnv -c $PSCC/Attachment/hg19.length -n $PSCC/Attachment/hg19.N_region -log $sample.ratio.cnv.log",
				"perl /disk2/mygeno/web/script/batchmendelian/svg_v3.pl $sample.ratio $sub_gender",
				"perl /disk2/mygeno/yangh/lccnv/script/svg_v4.pl $sample.ratio $sub_gender $sample.downsample.point"
			);
			&mysql_inject("CNVFilterD","正在运行",60);
			`perl $PSCC/Script/CNVFilterD.pl -i $sample.cnv -g $sub_gender -o $sample.cnv.filter -cyto /disk1/DB/annovar_db//humandb/hg19_cytoBand.txt -ref /disk1/DB/annovar_db//humandb/hg19_refGene.txt -sp 0.05 -pp 0.05`;
			`perl $PSCC/cnv-autofilter.pl $sample.cnv.filter $sample.$sub_gender.$karotype.CNV.AutoFilter.xlsx $sub_gender $karotype`;
			`python /disk2/mygeno/yangh/lccnv/script/doc.py $sample.ratio.png $sample.$sub_gender.$karotype.draft.docx`;
			&mysql_inject("并行:拷贝xls/png","任务即将完成",70);
			&parallel(
				"cp $sample.$sub_gender.$karotype.CNV.AutoFilter.xlsx /$tans/${lane}-web/CNV/$sample/$sub_gender/$sample.${depth}x.$sub_gender.$karotype.CNV.AutoFilter.xlsx ",
				"cp $sample.ratio.png /$tans/${lane}-web/CNV/$sample/$sub_gender/$sample.${depth}x.$sub_gender.$karotype.CNV.png",
				"cp $sample.downsample.point.png /$tans/${lane}-web/CNV/$sample/$sub_gender/$sample.${depth}x.$sub_gender.$karotype.downsample.point.png",
				"cp $sample.$sub_gender.$karotype.draft.docx /$tans/${lane}-web/CNV/$sample/$sub_gender/"
			);
		}
	}
	
	&mysql_inject("拷贝BAM文件/拷贝统计文件","正在运行",80);
	unless(-e "/$tans/${lane}-web/CNV/$sample/BAM"){
	        `mkdir -p /$tans/${lane}-web/CNV/$sample/BAM`;
	}
	&parallel(
	    "cp $sample.Statistics.xls /$tans/${lane}-web/CNV/$sample/",
	    "cp $sample.raw.insertsize.pdf /$tans/${lane}-web/CNV/$sample/",
	    "cp $sample.sort.bam /$tans/${lane}-web/CNV/$sample/BAM/",
	    "cp $sample.sort.bam.bai /$tans/${lane}-web/CNV/$sample/BAM/"
	);
	&mysql_inject("清理空间","正在运行",90);
	`rm -rf /local_disk/tmp/web/batchmendelian/$random/$lane/CNV/$sample`;
	`rm -rf *`;
	&mysql_inject("该任务已运行完成","已完成",100);
}

sub mysql_inject{
	my $comment = shift @_;
	my $status = shift @_;
	my $progress = shift @_;
	$progress ||= 0;
	my $dbh = DBI->connect("DBI:mysql:database=db_inhouse;host=11.11.11.69","web",'web@mygeno');
	$dbh->do("SET NAMES utf8");
	chomp(my $time = `date \"+%Y-%m-%d %H:%M:%S\"`);
	my $end_t = time();
	my $elap = sprintf "%.3f",($end_t - $start_t) / 3600;
	my $rows = $dbh->do("update pipeline_batchmendelian set comment=\"$comment\",status=\"$status\",end_time=\"$time\",duration=$elap,progress=\"$progress\" where pbs_id=\"$pbs_id\"");
}

sub mysql_mail_manager{
	my $dbh = DBI->connect("DBI:mysql:database=db_inhouse;host=11.11.11.69","web",'web@mygeno');
	$dbh->do("SET NAMES utf8");
	my $sth = $dbh->prepare("select done,total from pipeline_batchmendelian_mail_manager where random_string=\"$random\"");
	$sth->execute();
	my @mysql_rst = $sth->fetchrow_array();
	my ($done_sample_count,$total_sample_count) = @mysql_rst[0,1];
	$done_sample_count += ~~ (split /,/,$panel_for_mail);
	if($total_sample_count == $done_sample_count){
		#说明这批lane的所有任务已经完成了
		`ls /$tans/${lane}-web/CNV/*/*.xls | xargs -i sed 1d {} |sed '1iSample\tRaw_data_bases(Mb)\tClean_data_bases(Mb)\tAverage sequencing depth\tQ30' >/$tans/${lane}-web/CNV/${lane}.CNV.Statistics.xls`;
		
		`perl /disk2/mygeno/web/batchmendelian/scripts/MakeClassify.pl /disk2/mygeno/web/batchmendelian/web-server/$random/output.txt /$tans/${lane}-web/`;
		`rm -rf /local_disk/tmp/web/batchmendelian/$random`;
		my $rows = $dbh->do("update pipeline_batchmendelian_mail_manager set done=\"$done_sample_count\" where random_string=\"$random\"");
		if(-e "${lane}_SNP_INDEL.xls"){
			`rm -f ${lane}_SNP_INDEL.xls`;
		}
		`ls /$tans/${lane}-web/SNP_INDEL/*.xls|grep -v Clinican|grep -v Canmute|grep -v MitoChip|grep -v Exome|xargs -i sed 1d {} >> ${lane}_SNP_INDEL.xls`;
		`sed -i '1iDescr\thgvs\tDisease\tTag\tPublication\tEffect\tSample\tGene_Symbol\tID\tRef_Transcript\tExon\tNucleotide_Changes\tAmino_Acid_Changes\tGene_Type\t1000g2015aug_all\tPathogenic_Analysis\tInhert\tDisease\tDiseasePhenotype\tDiseaseInformation\tclinvar_20160302\tref/alt(fre)\tChr\tBegin\tEnd\tRef\tAlt\tMutRatio\tRefCount\tMutCount\tDepth\tMutation_Type\tGene_Region\tCytoBand\tdbsnp147\tPathSNP\tMutInNormal\t1000Genome\tMutInDatabase\t1000g2015aug_all\tESP6500si\tInhouse\tExAC_ALL\tExAC_EAS\tSIFT\tSIFT_Predict\tPolyPhen_2\tPolyPhen_2_Predict\tMutationTaster\tMutationTaster_Predict\tGERP++\tGERP++_Predict\tSPIDEX\tAA_change\tdbNsfpInterPro\tREVEL_score\tMCAP_score\tMCAP_pred\tgnomAD_exome_ALL\tgnomAD_exome_EAS\tgnomAD_genome_ALL\tgnomAD_genome_EAS\tInterVar\tHighest-MAF' ${lane}_SNP_INDEL.xls`;
		`perl /disk2/mygeno/web/script/batchmendelian/MakeSelectedGeneWeb.pl ${lane}_SNP_INDEL.xls /$tans/${lane}-web/${lane}_SNP_INDEL.xlsx /$tans/${lane}-web/RawDataSelectedGene/`;
		`rm -f ${lane}_SNP_INDEL.xls`;
		if(-e "Statistics.xls"){
			`rm -f Statistics.xls`;
		}
		`ls /$tans/${lane}-web/Statistics/*.xls | xargs -i sed 1d {} >> Statistics.xls`;
		`sed -i '1iSample\tRaw_data_bases(Mb)\tClean_data_bases(Mb)\tAligned_bases(Mb)\tAligned\tInitial bases on target\tBase covered on target\tCoverage of target region\tEffective bases on target\tFraction of effective bases on target\tAverage sequencing depth on target\tFraction of target covered with at least 4X\tFraction of target covered with at least 10X\tFraction of target covered with at least 20X\tduplication rate\tchip' Statistics.xls`;
		`mv Statistics.xls /$tans/${lane}-web/Statistics.xls`;
		`rm -f /$tans/${lane}-web/网站正在分析这批数据请不要占用文件`;
		`rm -rf /ssd1/mygeno/web/batchmendelian/$random`;
		&sent_email;
	}else{
		#说明这批lane还有任务在运行
		my $rows = $dbh->do("update pipeline_batchmendelian_mail_manager set done=\"$done_sample_count\" where random_string=\"$random\"");
	}
}

sub sent_email{
	`echo -e "恭喜！您提交的遗传病任务已经分析完成，结果保存在 \\\\\\\\\\192.168.0.211\\\\\\tech 新公共盘，账号：tech，密码：tech147258（请不要长期占用、随意修改或删除公共盘的数据）。\n" | /usr/bin/mailx -s "遗传病任务分析完成" $email`;
}

sub parallel{
	my $pm = new Parallel::ForkManager(7);
	for (@_){
		my $pid = $pm->start and next;
		&do($_);
		$pm->finish;
	}
	$pm->wait_all_children;
	sub do{
		`$_[0]`;
	}
}

sub intervals{
	my $initial_intervals = $panel;
	$initial_intervals =~ s/MitoChip,//g;
	$initial_intervals =~ s/,MitoChip//g;
	$initial_intervals =~ s/,CNV//g;
	$initial_intervals =~ s/^CNV,//g;
	my $intervals = join " ",map {$_ = "-L /disk1/bed/$_.bed"} grep {$_ ne ""} (split /,/,$initial_intervals);
	return $intervals;
}

sub initial{
	my $dbh = DBI->connect("DBI:mysql:database=db_inhouse;host=11.11.11.69","web",'web@mygeno');
	$dbh->do("SET NAMES utf8");
	chomp(my $time = `date \"+%Y-%m-%d %H:%M:%S\"`);
	my $rows = $dbh->do("update pipeline_batchmendelian set start_time=\"$time\" where pbs_id=\"$pbs_id\"");
}

sub check_ssd1_TNAS2 {
    &mysql_inject("检测硬盘存储是否充足","正在运行",0);
    my $memory = 0;
    while ($memory < 1) {
        my ($ssd1, $TNAS2) = split "\t", `perl /disk2/mygeno/web/script/batchmendelian/check_memory.pl`;
        if ($ssd1 == 1 && $TNAS2 ==1){
            $memory = 1;
        }else {
            &mysql_inject("存储不足，任务正在等待，请联系生物信息部陈亚军！","任务正在运行中");
            sleep(30);
        }
    }
}
