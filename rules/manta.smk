rule manta:
    input:   
        lambda wildcards: get_bams_by_set(wildcards, sets, target="input")
    output:
        vcf="variant_calling/sv/{set}/{set}.diploidSV.vcf.gz"
    conda:
        "../envs/manta.yaml"
    params:
        fasta=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        bams=lambda wildcards: get_bams_by_set(wildcards, sets, target="shell"),
        rundir="variant_calling/sv/{set}"
    benchmark:
        "benchmarks/manta/{set}.diploidSV.txt"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=config.get("rules").get("manta").get("threads"))
    shell:
        "configManta.py "
        "--referenceFasta {params.fasta} "
        "{params.bams} "
        "--runDir {params.rundir}; "
        "cd {params.rundir}; "
        "./runWorkflow.py -j {threads} --quiet -g 30; "
        "mv results/variants/diploidSV.vcf.gz ../../../{output.vcf} "

rule mantaINV:
    input:   
        rules.manta.output
    output:
        vcf="variant_calling/sv/{set}/{set}.diploidSV_INV.vcf.gz"
    conda:
       "../envs/manta.yaml"
    params:
        fasta=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    benchmark:
        "benchmarks/manta/{set}.diploidSV_INV.txt"
    threads: config.get("rules").get("mantaINV").get("threads")
    shell:
        """NEWPATH="$(dirname `which configManta.py`)/../share/manta-1.6.0-0/libexec"; """
        "export PATH=$PATH:$NEWPATH; "
        "convertInversion.py "
        """$NEWPATH/samtools """
        "{params.fasta} "       
        "{input} "
        "> {output.vcf} "
