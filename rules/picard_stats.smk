
rule picard_InsertSizeMetrics:
   input:
       bam="reads/recalibrated/{sample}.dedup.recal.bam"
   output:
       metrics="reads/recalibrated/{sample}.dedup.recal.ismetrics.txt",
       histogram="reads/recalibrated/{sample}.dedup.recal.ismetrics.pdf"
   conda:
       "../envs/picard.yaml"
   params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
   benchmark:
       "benchmarks/picard/IsMetrics/{sample}.txt"
   log:
       "logs/picard/IsMetrics/{sample}.log"
   threads: config.get("rules").get("picard_InsertSizeMetrics").get("threads")
   shell:
       "picard {params.custom} CollectInsertSizeMetrics "
       "I={input.bam} "
       "O={output.metrics} "
       "H={output.histogram} "
       "&> {log} "

rule picard_WGSMetrics:
   input:
       bam="reads/recalibrated/{sample}.dedup.recal.bam"
   output:
       metrics="reads/recalibrated/{sample}.dedup.recal.wgsmetrics.txt"
   conda:
       "../envs/picard.yaml"
   params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        arguments=config.get("rules").get("picard_WGSMetrics").get("arguments")
   benchmark:
       "benchmarks/picard/WGSMetrics/{sample}.txt"
   log:
       "logs/picard/WGSMetrics/{sample}.log"
   threads: config.get("rules").get("picard_WGSMetrics").get("threads")
   shell:
       "picard {params.custom} CollectWgsMetrics "
       "{params.arguments} "
       "I={input.bam} "
       "O={output.metrics} "
       "R={params.genome} "
       "&> {log} "

rule gatk_DepthOfCoverage:
    input:
        cram="reads/recalibrated/{sample}.dedup.recal.cram",
        crai="reads/recalibrated/{sample}.dedup.recal.cram.crai"
    output:
        "reads/recalibrated/{sample}.sample_gene_summary"
    params:
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        cov_refseq=resolve_single_filepath(*references_abs_path(),config.get("depthofcov_refseq")),
        cov_intervals=resolve_single_filepath(*references_abs_path(),config.get("depthofcov_intervals")),
        prefix="reads/recalibrated/{sample}"
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk/DepthOfCoverage/{sample}.txt"
    log:
        "logs/gatk/DepthOfCoverage/{sample}.txt"
    threads: 4
    shell:
        "gatk DepthOfCoverage "
        "--omit-depth-output-at-each-base --omit-locus-table "
        "-R {params.genome} "
        "-O {params.prefix} "
        "-I {input.cram} "
        "-gene-list {params.cov_refseq} "
        "--summary-coverage-threshold 10 --summary-coverage-threshold 30 --summary-coverage-threshold 50 "
        "-L {params.cov_intervals} "
        ">& {log} "

