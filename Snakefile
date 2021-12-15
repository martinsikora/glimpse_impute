import pandas as pd


## --------------------------------------------------------------------------------
## global parameters from config file

PREFIX = config["prefix"]  ## output prefix
OUT_DIR = config["out_dir"]  ## path for output files
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

vcf_impute_all = expand(
    OUT_DIR + "/{chrom}." + PREFIX + ".glimpse.{ext}",
    chrom=CHROMS,
    ext=["vcf.gz", "vcf.gz.tbi"],
)


## --------------------------------------------------------------------------------
## targets


rule all:
    input:
        vcf_impute_all,


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
        vcf=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".gl_merge.vcf.gz"),
        tbi=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".gl_merge.vcf.gz.tbi"),
        lst=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".gl_merge.files.txt"),
    shell:
        """
        echo {input.vcf} | perl -p -e 's/ /\\n/g' > {output.lst}
        bcftools merge -l {output.lst} -m all -Oz > {output.vcf}
        bcftools index -t {output.vcf}
        """


rule chunk_glimpse:
    input:
        vcf=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".gl_merge.vcf.gz",
    wildcard_constraints:
        chrom="\d+",
    output:
        temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".chunks.txt"),
    log:
        OUT_DIR + "/logs/{chrom}.chunks.glimpse.log",
    shell:
        """
        {GLIMPSE_PATH}/GLIMPSE_chunk --input {input.vcf} --region {wildcards.chrom} --window-size {GLIMPSE_WS} --buffer-size {GLIMPSE_BS} --output {output} --log {log}
        """


checkpoint get_chunks_chrom:
    input:
        TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".chunks.txt",
    output:
        temp(directory(TMP_DIR + "/{chrom}/glimpse_chunks/")),
    shell:
        """
        mkdir {TMP_DIR}/{wildcards.chrom}/glimpse_chunks
        cat {input} | awk '{{print >"{TMP_DIR}/{wildcards.chrom}/glimpse_chunks/"$2"_"$1".chunk"}}'
        """


rule impute_chunk_glimpse:
    input:
        vcf=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".gl_merge.vcf.gz",
        vcf_ref=config["panel"]["vcf"],
        genmap=config["panel"]["genmap"],
        chunk=TMP_DIR + "/{chrom}/glimpse_chunks/{chunk}.chunk",
    wildcard_constraints:
        chrom="\d+",
    output:
        vcf=temp(TMP_DIR + "/{chrom}/{chunk}." + PREFIX + ".chunk.vcf.gz"),
        tbi=temp(TMP_DIR + "/{chrom}/{chunk}." + PREFIX + ".chunk.vcf.gz.tbi"),
    threads: GLIMPSE_THREADS
    log:
        OUT_DIR + "/logs/{chrom}.{chunk}.impute_chunk.glimpse.log",
    benchmark:
        OUT_DIR + "/benchmarks/{chrom}.{chunk}.impute_chunk.glimpse.txt"
    shell:
        """
        {GLIMPSE_PATH}/GLIMPSE_phase --thread {threads} --input {input.vcf} --reference {input.vcf_ref} --map {input.genmap} --input-region $(cat {input.chunk} | cut -f3) --output-region $(cat {input.chunk} | cut -f4) --output {output.vcf} --log {log} --seed $RANDOM
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
        aggregate_chunks,
    output:
        vcf=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".ligate.vcf.gz"),
        tbi=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".ligate.vcf.gz.tbi"),
        lst=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".list"),
    log:
        OUT_DIR + "/logs/{chrom}.ligate_chunks.glimpse.log",
    shell:
        """
        echo {input} | perl -p -e 's/ /\\n/g' | grep -v .tbi > {output.lst}
        {GLIMPSE_PATH}/GLIMPSE_ligate --input {output.lst} --output {output.vcf} --log {log}
        bcftools index -t {output.vcf}
        """


rule phase_glimpse:
    input:
        vcf=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".ligate.vcf.gz",
        tbi=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".ligate.vcf.gz.tbi",
    output:
        vcf=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".phased.vcf.gz"),
        tbi=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".phased.vcf.gz.tbi"),
    log:
        OUT_DIR + "/logs/{chrom}.phase.glimpse.log",
    shell:
        """
        {GLIMPSE_PATH}/GLIMPSE_sample --input {input.vcf} --solve --output {output.vcf} --log {log}
        bcftools index -t {output.vcf}
        """


rule annote_vcf_phase:
    input:
        vcf_ligate=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".ligate.vcf.gz",
        tbi_ligate=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".ligate.vcf.gz.tbi",
        vcf_phase=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".phased.vcf.gz",
        tbi_phase=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".phased.vcf.gz.tbi",
    output:
        vcf_ant=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".annot.vcf.gz"),
        tbi_ant=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".annot.vcf.gz.tbi"),
    shell:
        """
        bcftools annotate -a {input.vcf_phase} -c FMT/GT {input.vcf_ligate} -Oz > {output.vcf_ant}
        bcftools index -t {output.vcf_ant}
        """


rule annote_vcf_gl:
    input:
        vcf_ant=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".annot.vcf.gz",
        tbi_ant=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".annot.vcf.gz.tbi",
        vcf_gl=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".gl_merge.vcf.gz",
        tbi_gl=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".gl_merge.vcf.gz.tbi",
    output:
        vcf=OUT_DIR + "/{chrom}." + PREFIX + ".glimpse.vcf.gz",
        tbi=OUT_DIR + "/{chrom}." + PREFIX + ".glimpse.vcf.gz.tbi",
    shell:
        """
        bcftools annotate -a {input.vcf_gl} -c FMT/DP,FMT/PL -Oz {input.vcf_ant} > {output.vcf}
        bcftools index -t {output.vcf}
        """
