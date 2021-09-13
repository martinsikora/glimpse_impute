import pandas as pd


## --------------------------------------------------------------------------------
## global parameters from config file

PREFIX = config["prefix"]  ## output prefix
REF = config["ref"]  ## reference genome
MQ = config["mq"]  ## minimum mapping quality
TMP_DIR = config["tmp_dir"]  ## path for temporary files

GLIMPSE_PATH = config["glimpse"]["path"]
GLIMPSE_WS = config["glimpse"]["windowsize"]
GLIMPSE_BS = config["glimpse"]["buffersize"]
GLIMPSE_THREADS = config["glimpse"]["threads"]


## --------------------------------------------------------------------------------
## helpers

CHROMS = range(1, 23)

unit_df = pd.read_table(config["units"], comment="#").set_index(["sampleId"])
SAMPLES = unit_df.index.drop_duplicates().values.tolist()


## --------------------------------------------------------------------------------
## functions

def get_bam(wildcards):
    return unit_df.loc[(wildcards.sample), "bam"]


## --------------------------------------------------------------------------------
## output file sets

vcf_gl_all = expand(
    "vcf/gl_merge/{chrom}." + PREFIX + ".merge.gl.{ext}",
    chrom=CHROMS,
    ext=["vcf.gz", "vcf.gz.tbi"],
)
vcf_impute_all = expand(
    "vcf/impute/{chrom}." + PREFIX + ".glimpse.{ext}",
    chrom=CHROMS,
    ext=["vcf.gz", "vcf.gz.tbi"],
)
vcf_phase_all = expand(
    "vcf/impute/{chrom}." + PREFIX + ".glimpse.phased.{ext}",
    chrom=CHROMS,
    ext=["vcf.gz", "vcf.gz.tbi"],
)


## --------------------------------------------------------------------------------
## targets


rule all:
    input:
        vcf_gl_all + vcf_impute_all + vcf_phase_all


## --------------------------------------------------------------------------------
## rules

rule compute_gl:
    input:
        bam=get_bam,
        alleles=config["panel"]["alleles"],
    wildcard_constraints:
        chrom="\d+",
    output:
        vcf=temp(TMP_DIR + "/{chrom}/{chrom}.{sample}." + PREFIX + ".gl.vcf.gz"),
        tbi=temp(TMP_DIR + "/{chrom}/{chrom}.{sample}." + PREFIX + ".gl.vcf.gz.tbi"),
        sample=temp(TMP_DIR + "/{chrom}/{sample}.sample.txt"),
    shell:
        """
        echo -e {wildcards.sample} > {output.sample}
        samtools merge - {input.bam} | bcftools mpileup -I -E -a FORMAT/DP,FORMAT/AD --ignore-RG -f {REF} -q {MQ} -T {input.alleles} - | bcftools reheader -s {output.sample} | bcftools call --ploidy GRCh37 -Aim -C alleles -T {input.alleles} -Oz > {output.vcf}
        bcftools index -t {output.vcf} 
        """


rule merge_gl:
    input:
        vcf=expand(
            TMP_DIR + "/{{chrom}}/{{chrom}}.{sample}." + PREFIX + ".gl.vcf.gz",
            sample=SAMPLES,
        ),
        tbi=expand(
            TMP_DIR + "/{{chrom}}/{{chrom}}.{sample}." + PREFIX + ".gl.vcf.gz.tbi",
            sample=SAMPLES,
        ),
    wildcard_constraints:
        chrom="\d+",
    output:
        vcf="vcf/gl_merge/{chrom}." + PREFIX + ".merge.gl.vcf.gz",
        tbi="vcf/gl_merge/{chrom}." + PREFIX + ".merge.gl.vcf.gz.tbi",
    shell:
        """
        echo {input.vcf} | perl -p -e 's/ /\\n/g' > vcf/gl_merge/{wildcards.chrom}.files.txt
        bcftools merge -l vcf/gl_merge/{wildcards.chrom}.files.txt -m all -Oz > {output.vcf}
        bcftools index -t {output.vcf}
        """


rule chunk_glimpse:
    input:
        vcf="vcf/gl_merge/{chrom}." + PREFIX + ".merge.gl.vcf.gz",
        vcf_ref=config["panel"]["vcf"],
    wildcard_constraints:
        chrom="\d+",
    output:
        "glimpse_chunks/{chrom}." + PREFIX + ".chunks.glimpse.txt",
    log:
        "logs/{chrom}.chunks.glimpse.log",
    shell:
        """
        {GLIMPSE_PATH}/GLIMPSE_chunk --input {input.vcf} --reference {input.vcf_ref} --region {wildcards.chrom} --window-size {GLIMPSE_WS} --buffer-size {GLIMPSE_BS} --output {output} --log {log}
        """


checkpoint get_chunks_chrom:
    input:
        "glimpse_chunks/{chrom}." + PREFIX + ".chunks.glimpse.txt",
    output:
        directory("glimpse_chunks/{chrom}/"),
    shell:
        """
        mkdir glimpse_chunks/{wildcards.chrom}
        cat {input} | awk '{{print >"glimpse_chunks/{wildcards.chrom}/"$2"_"$1".chunk"}}'
        """


rule impute_chunk_glimpse:
    input:
        vcf="vcf/gl_merge/{chrom}." + PREFIX + ".merge.gl.vcf.gz",
        vcf_ref=config["panel"]["vcf"],
        genmap=config["panel"]["genmap"],
        chunk="glimpse_chunks/{chrom}/{chunk}.chunk",
    wildcard_constraints:
        chrom="\d+",
    output:
        vcf=temp(TMP_DIR + "/{chrom}/{chunk}." + PREFIX + ".chunk.vcf.gz"),
        tbi=temp(TMP_DIR + "/{chrom}/{chunk}." + PREFIX + ".chunk.vcf.gz.tbi"),
    threads: GLIMPSE_THREADS
    log:
        "logs/{chrom}.{chunk}.impute_chunk.glimpse.log",
    benchmark:
        "benchmarks/{chrom}.{chunk}.impute_chunk.glimpse.txt"
    shell:
        """
        {GLIMPSE_PATH}/GLIMPSE_phase --thread {threads} --input {input.vcf} --reference {input.vcf_ref} --map {input.genmap} --input-region $(cat {input.chunk} | cut -f3) --output-region $(cat {input.chunk} | cut -f4) --output {output.vcf}
        bcftools index -t {output.vcf}
        """


def aggregate_chunks(wildcards):
    checkpoint_output = checkpoints.get_chunks_chrom.get(**wildcards).output[0]

    files = expand(
        TMP_DIR + "/{chrom}/{chunk}." + PREFIX + ".chunk.{ext}",
        chrom=wildcards.chrom,
        chunk=glob_wildcards(os.path.join(checkpoint_output, "{chunk}.chunk")).chunk,
        ext=["vcf.gz", "vcf.gz.tbi"],
    )
    return files


rule ligate_chunks_glimpse:
    input:
        aggregate_chunks
    output:
        vcf="vcf/impute/{chrom}." + PREFIX + ".glimpse.vcf.gz",
        tbi="vcf/impute/{chrom}." + PREFIX + ".glimpse.vcf.gz.tbi",
        lst=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".list")
    log:
        "logs/{chrom}.ligate.glimpse.log",
    shell:
        """
        echo {input} | perl -p -e 's/ /\\n/g' | grep -v .tbi > {output.lst}
        {GLIMPSE_PATH}/GLIMPSE_ligate --input {output.lst} --output {output.vcf} --log {log}
        bcftools index -t {output.vcf}
        """


rule phase_glimpse:
    input:
        vcf="vcf/impute/{chrom}." + PREFIX + ".glimpse.vcf.gz",
        tbi="vcf/impute/{chrom}." + PREFIX + ".glimpse.vcf.gz.tbi",
    output:
        vcf="vcf/impute/{chrom}." + PREFIX + ".glimpse.phased.vcf.gz",
        tbi="vcf/impute/{chrom}." + PREFIX + ".glimpse.phased.vcf.gz.tbi",
    log:
        "logs/{chrom}.phase.glimpse.log",
    shell:
        """
        {GLIMPSE_PATH}/GLIMPSE_sample --input {input.vcf} --solve --output {output.vcf} --log {log}
        bcftools index -t {output.vcf}
        """
