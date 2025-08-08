#!/bin/bash

# Setup environment and directories
WORKDIR="../upd_trio_routine"
TOOLDIR="../tools"  # Parent directory containing existing tools
UPDIOHOME="../data803/yfan/fixedUPDio/UPDio-master/version_1.0"  # Path to existing UPDio installation

# Create necessary directories if they don't exist
mkdir -p $WORKDIR/{input,temp,output}

# Set up Perl environment for VCFtools
export PERL5LIB="../tools/vcftools-vcftools-581c231/src/perl:$PERL5LIB"

# Change to temp directory for processing
cd $WORKDIR/temp

# Read trio information from manifest file
while read probandID fatherID motherID probandvcf fathervcf mothervcf
do
    echo "Processing trio: $probandID $fatherID $motherID"
    
    # Compress VCF files if not already compressed
    if [[ ! $probandvcf =~ \.gz$ ]]; then
        bgzip -c $probandvcf > $probandID.vcf.gz
    else
        cp $probandvcf $probandID.vcf.gz
    fi

    if [[ ! $fathervcf =~ \.gz$ ]]; then
        bgzip -c $fathervcf > $fatherID.vcf.gz
    else
        cp $fathervcf $fatherID.vcf.gz
    fi

    if [[ ! $mothervcf =~ \.gz$ ]]; then
        bgzip -c $mothervcf > $motherID.vcf.gz
    else
        cp $mothervcf $motherID.vcf.gz
    fi

    # Index compressed VCF files
    tabix -p vcf $probandID.vcf.gz
    tabix -p vcf $fatherID.vcf.gz
    tabix -p vcf $motherID.vcf.gz

    # Filter chromosomes
    CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
    bcftools filter -r $CHROMS $probandID.vcf.gz -O z -o $probandID.processed.vcf.gz
    bcftools filter -r $CHROMS $fatherID.vcf.gz -O z -o $fatherID.processed.vcf.gz
    bcftools filter -r $CHROMS $motherID.vcf.gz -O z -o $motherID.processed.vcf.gz

    # Add homozygous reference sites
    perl $UPDIOHOME/scripts/add_hom_refs_to_vcf.pl \
        --polymorphic_sites $UPDIOHOME/sample_data/common_variants_within_well_covered_target_regions.txt \
        --no_homREF_vcf $probandID.processed.vcf.gz | bgzip > $probandID.homREFed.vcf.gz

    perl $UPDIOHOME/scripts/add_hom_refs_to_vcf.pl \
        --polymorphic_sites $UPDIOHOME/sample_data/common_variants_within_well_covered_target_regions.txt \
        --no_homREF_vcf $fatherID.processed.vcf.gz | bgzip > $fatherID.homREFed.vcf.gz

    perl $UPDIOHOME/scripts/add_hom_refs_to_vcf.pl \
        --polymorphic_sites $UPDIOHOME/sample_data/common_variants_within_well_covered_target_regions.txt \
        --no_homREF_vcf $motherID.processed.vcf.gz | bgzip > $motherID.homREFed.vcf.gz

    # Run UPDio analysis
    perl $UPDIOHOME/UPDio.pl \
        --child_vcf $WORKDIR/temp/$probandID.homREFed.vcf.gz \
        --dad_vcf $WORKDIR/temp/$fatherID.homREFed.vcf.gz \
        --mom_vcf $WORKDIR/temp/$motherID.homREFed.vcf.gz \
        --include_MI \
        --common_cnv_file $UPDIOHOME/sample_data/common_dels_1percent.tsv \
        --path_to_R /data/sysoft/miniconda3/bin/Rscript \
        --R_scripts_dir $UPDIOHOME/scripts \
        --output_path $WORKDIR/output

    echo "Completed processing trio: $probandID"

done < $WORKDIR/input/manifest.txt

# Clean up temporary files
rm -f $WORKDIR/temp/*.tbi $WORKDIR/temp/*.vcf.gz

echo "Analysis complete. Results are in $WORKDIR/output"