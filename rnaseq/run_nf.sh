samples="samplesheet.csv"
outdir="run00"
annotation="/dados02/jorge/israel_rnaseq/reference/gene_prediction/annot.gtf"
genome="/dados02/jorge/israel_rnaseq/reference/GCA_033110985.1_ASM3311098v1_genomic.fna"

conda activate nextflow_25.10.4 

nextflow run nf-core/rnaseq --input $samples --outdir $outdir --gtf $annotation --fasta $genome --aligner star_salmon --skip_qc false -resume -profile mamba --seq_platform ILLUMINA -c custom.config --skip_dupradar


# using annotation and genome used with diego
samples="samplesheet.csv"
outdir="run01"
annotation="/dados01/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.gtf.gz"
genome="/dados01/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.fa.gz"


nextflow run nf-core/rnaseq --input $samples --outdir $outdir --gtf $annotation --fasta $genome --aligner star_salmon --skip_qc false -resume -profile mamba --seq_platform ILLUMINA -c custom.config --skip_dupradar

extflow run nf-core/rnaseq --input samplesheet.csv --outdir run01 --gtf /dados01/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.gtf.gz --fasta /dados01/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.fa.gz --aligner star_salmon --skip_qc false -resume -profile mamba --seq_platform ILLUMINA -c custom.config --skip_dupradar
