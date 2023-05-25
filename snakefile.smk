
rule get_var:
    input:
        ref = 'ref_ind.fasta',
        reads = 'reads_aligned_sorted.bam'
    output:
        'reads.mpileup'
    shell:
        'samtools mpileup -d 0 -f {input.ref} {input.reads} > {output}'


rule index_bam:
    input:
        'reads_aligned_sorted.bam'
    shell:
        'samtools index {input}'

rule align_:
    input:
        ref = 'ref_ind.fasta',
        reads = 'reads.fastq.gz'
    output:
        'reads_aligned_sorted.bam'
    shell:
        'bwa mem {input.ref} {input.reads} | samtools view -S -b |  samtools sort -o {output}'


rule index_ref:
    input:
        ref = 'ref.fasta'
    output:
        ref_idx = 'ref_ind.fasta'
    shell:
        'mv {input.ref} {output.ref_idx} | bwa index -p {output.ref_idx} {output.ref_idx}'

rule downoload_ref:
    output:
        ref = 'ref.fasta'
    shell:
        'efetch -db nucleotide -id KF848938.1 -format fasta > {output.ref}'
















# rule inspect_data:
#     input:
#         reads = 'reads.fastq.gz',

#         dir = 'QC'
#     output:
#         'counts_read.txt'
#     shell:
#         'zcat {input.reads} | wc -l > {output} | mv {input.reads} {input.dir}'

# rule fastqc_:
#     input:
#         reads = 'reads.fastq.gz'
#     output:
#         dir = directory('QC')
#     shell:
#         'mkdir {output.dir} | fastqc {input.reads} -o {output.dir}'




# rule downoload_reads:
#     output:
#         'reads.fastq.gz'
#     shell:
#         'wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/SRR1705851.fastq.gz | mv SRR1705851.fastq.gz reads.fastq.gz'





