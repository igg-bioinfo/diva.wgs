
def get_fastq(wildcards,units):
    # print(wildcards.unit)
    if units.loc[wildcards.unit,["fq2"]].isna().all():
        print("SE")
        # print(units.loc[wildcards.unit,["fq1"]].dropna()[0])
        return units.loc[wildcards.unit,["fq1"]].dropna()[0]
    else:
        print("PE")
        # print(units.loc[wildcards.unit,["fq1"]].dropna()[0],units.loc[wildcards.unit,["fq2"]].dropna()[0])
        return units.loc[wildcards.unit,["fq1"]].dropna()[0],units.loc[wildcards.unit,["fq2"]].dropna()[0]


rule pre_rename_fastq_pe:
    input:
        lambda wildcards: get_fastq(wildcards,units)
    output:
        r1=temp("reads/untrimmed/{unit}-R1.fq.gz"),
        r2=temp("reads/untrimmed/{unit}-R2.fq.gz")
    shell:
        "cp {input[0]} {output.r1} &&"
        "cp {input[1]} {output.r2} "
#        "ln -s {input[0]} {output.r1} &&"
#        "ln -s {input[1]} {output.r2} "

rule pre_rename_fastq_se:
    input:
       lambda wildcards: get_fastq(wildcards,units)
    output:
        r1=temp("reads/untrimmed/se/{unit}-R1.fq.gz")
    shell:
        "cp {input} {output.r1} "
#        "ln -s {input} {output.r1} "


rule trim_galore_pe:
    input:
        rules.pre_rename_fastq_pe.output
    output:
        read1=temp("reads/trimmed/{unit}-R1-trimmed.fq.gz"),
        read2=temp("reads/trimmed/{unit}-R2-trimmed.fq.gz"),
        report1="reads/trimmed/{unit}-R1.fq.gz_trimming_report.txt",
        report2="reads/trimmed/{unit}-R2.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_pe").get("arguments"),
        outdir=config.get("tmp_dir"),
        read1=config.get("tmp_dir")+"/{unit}-R1_val_1.fq.gz",
        read2=config.get("tmp_dir")+"/{unit}-R2_val_2.fq.gz",
        report1=config.get("tmp_dir")+"/{unit}-R1.fq.gz_trimming_report.txt",
        report2=config.get("tmp_dir")+"/{unit}-R2.fq.gz_trimming_report.txt"
    log:
        "logs/trim_galore/{unit}.log"
    benchmark:
        "benchmarks/trim_galore/{unit}.txt"
    conda:
        "../envs/trim_galore.yaml"
    threads: (conservative_cpu_count(reserve_cores=2, max_cores=config.get("rules").get("trim_galore_pe").get("threads")))/4 if (conservative_cpu_count(reserve_cores=2, max_cores=config.get("rules").get("trim_galore_pe").get("threads"))) > 4 else 1
    shell:
        "mkdir -p qc/fastqc; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input} "
        ">& {log}; "
        "mv {params.read1} {output.read1} && "
        "mv {params.read2} {output.read2} && "
        "mv {params.report1} {output.report1} && "
        "mv {params.report2} {output.report2}"
 

rule trim_galore_se:
    input:
        rules.pre_rename_fastq_se.output
    output:
        temp("reads/trimmed/se/{unit}-R1_trimmed.fq.gz"),
        "reads/trimmed/se/{unit}-R1.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_se").get("arguments"),
        outdir="reads/trimmed/se/"
    log:
        "logs/trim_galore/{unit}.log"
    benchmark:
        "benchmarks/trim_galore/{unit}.txt"
    conda:
        "../envs/trim_galore.yaml"
    threads: (conservative_cpu_count(reserve_cores=2, max_cores=config.get("rules").get("trim_galore_se").get("threads")))/2 if (conservative_cpu_count(reserve_cores=2, max_cores=config.get("rules").get("trim_galore_se").get("threads"))) > 2 else 1
    shell:
        "mkdir -p qc/fastqc; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input} "
        ">& {log}"




#rule post_rename_fastq_pe:
#    input:
#        rules.trim_galore_pe.output
#    output:
#        read1=temp("reads/trimmed/{unit}-R1-trimmed.fq.gz"),
#        read2=temp("reads/trimmed/{unit}-R2-trimmed.fq.gz"),
#        report1="reads/trimmed/{unit}-R1.fq.gz_trimming_report.txt",
#        report2="reads/trimmed/{unit}-R2.fq.gz_trimming_report.txt"
#    shell:
#        "mv {input[0]} {output.read1} && "
#        "mv {input[1]} {output.read2} && "
#        "mv {input[2]} {output.report1} && "
#        "mv {input[3]} {output.report2}"


rule post_rename_fastq_se:
    input:
        rules.trim_galore_se.output
    output:
        r1=temp("reads/se/trimmed/{unit}-R1-trimmed.fq.gz")
    shell:
        "mv {input[0]} {output.r1}"



def get_trimmed_reads(wildcards,units):
    print(wildcards.unit)
    if units.loc[wildcards.unit,["fq2"]].isna().all():
        # SE
        return rules.post_rename_fastq_se.output
    # PE
    else:
        return rules.trim_galore_pe.output.read1,rules.trim_galore_pe.output.read2
