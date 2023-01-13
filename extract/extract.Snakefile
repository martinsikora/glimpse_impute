import pandas as pd


## --------------------------------------------------------------------------------
## global parameters from config file

PREFIX = config["prefix"]  ## output prefix
OUT_DIR = config["out_dir"]  ## path for output files
TMP_DIR = config["tmp_dir"]  ## path for temporary files


## --------------------------------------------------------------------------------
## helpers

CHROMS = range(1, 23)

unit_df = pd.read_table(config["units"], comment="#").set_index(["prefix"])
PX = unit_df.index.drop_duplicates().values.tolist()


## --------------------------------------------------------------------------------
## functions


def get_vcf(wildcards):
    f = unit_df.loc[[(wildcards.px)]].path.values[0] + "/" + wildcards.chrom + "." + wildcards.px + ".glimpse.vcf.gz"
    return f


## --------------------------------------------------------------------------------
## output file sets

vcf_all = expand(
    OUT_DIR + "/vcf/{chrom}." + PREFIX + ".glimpse.{ext}",
    chrom=CHROMS,
    ext=["vcf.gz", "vcf.gz.tbi"],
)
plink_all = expand(
    OUT_DIR + "/plink/" + PREFIX + ".glimpse.{ext}",
    ext=["bed", "bim", "fam", "target.snplist", "eigenvec", "eigenval", "pca.plot.pdf"],
)


## --------------------------------------------------------------------------------
## targets


rule all:
    input:
        vcf_all + plink_all


## --------------------------------------------------------------------------------
## rules

rule extract_samples:
    input:
        config["units"]
    output:
        sample=temp(expand(TMP_DIR + "/{px}.samples.txt", px=PX))
    shell:
        """
        cat {input} | tail -n+2 | awk '{{print $2>"{TMP_DIR}/"$1".samples.txt"}}'
        """


rule get_vcf_samples:
    input:
        vcf=get_vcf,
        sample=TMP_DIR + "/{px}.samples.txt"
    wildcard_constraints:
        chrom="\d+",
    output:
        vcf=temp(TMP_DIR + "/{chrom}/{chrom}.{px}." + PREFIX + ".extract.vcf.gz"),
        sample=temp(TMP_DIR + "/{chrom}/{px}." + PREFIX + ".samples_new.txt"),
    shell:
        """
        cat {input.sample} | awk '{{print $1".{wildcards.px}"}}' > {output.sample}
        bcftools view -S {input.sample} -Oz {input.vcf} | bcftools reheader -s {output.sample} > {output.vcf}
        """


rule merge_vcf:
    input:
        vcf=expand(TMP_DIR + "/{{chrom}}/{{chrom}}.{px}." + PREFIX + ".extract.vcf.gz", px=PX)
    wildcard_constraints:
        chrom="\d+",
    output:
        vcf=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".merge.vcf.gz"),
        tbi=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".merge.vcf.gz.tbi"),
        lst=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".merge.lst")
    threads:
        2
    priority:
        10
    shell:
        """
        echo {input.vcf} | tr ' ' '\\n' > {output.lst}
        bcftools merge --no-index --threads {threads} -Ou -l {output.lst} | bcftools +impute-info -Oz - > {output.vcf}
        bcftools index -t {output.vcf} 
        """


rule extract_outgroup:
    input:
        vcf=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".merge.vcf.gz",
        fa=config["og"]
    wildcard_constraints:
        chrom="\d+"
    output:
        bed=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".all.bed"),
        og=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".og.txt"),
    shell:
        """
        bcftools query -f '%CHROM\\t%POS0\\t%POS\\n' {input.vcf} > {output.bed}
        seqtk subseq -t {input.fa} {output.bed} | tr 'a-z' 'A-Z' | cut -f2,3 > {output.og}
        """

        
rule extract_targets:
    input:
        bed1=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".all.bed",
        bed2=config["targets"]
    wildcard_constraints:
        chrom="\d+"
    output:
        bed=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".targets.bed"),
        tg=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".targets.txt"),
        lst=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".snplist")
    shell:
        """
        tabix {input.bed2} {wildcards.chrom} > {output.bed}
        bedtools intersect -a {input.bed1} -b {input.bed2} -c | cut -f3,4 > {output.tg}
        cat {output.tg} | awk '$2==1{{print "{wildcards.chrom}:"$1}}' > {output.lst}
        """

        
rule write_snplist:
    input:
        expand(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".snplist", chrom=CHROMS)
    wildcard_constraints:
        chrom="\d+"
    output:
        OUT_DIR + "/plink/" + PREFIX + ".glimpse.target.snplist"
    shell:
        """
        cat {input} > {output}
        """

        
rule annotate_vcf:
    input:
        vcf=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".merge.vcf.gz",
        tg=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".targets.txt",
        og=TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".og.txt",
        hdr=config["hdr"]
    wildcard_constraints:
        chrom="\d+"
    output:
        vcf=OUT_DIR + "/vcf/{chrom}." + PREFIX + ".glimpse.vcf.gz",
        tbi=OUT_DIR + "/vcf/{chrom}." + PREFIX + ".glimpse.vcf.gz.tbi",
        ant=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".ant.txt.gz"),
        ant_tbi=temp(TMP_DIR + "/{chrom}/{chrom}." + PREFIX + ".ant.txt.gz.tbi"),
    threads:
        8
    priority:
        20
    shell:
        """
        join {input.og} {input.tg} | awk '{{print "{wildcards.chrom}\\t"$1"\\t"$2"\\t"$3}}' | bgzip -c > {output.ant}
        tabix -s1 -b2 -e2 {output.ant}
        bcftools annotate -a {output.ant} -c CHROM,POS,pan_troglodytes,CAP_TARGET_1240K -h {input.hdr} --threads {threads} -Oz {input.vcf} > {output.vcf} 
        bcftools index -t {output.vcf} 
        """


rule concat_vcf:
    input:
        expand(OUT_DIR + "/vcf/{chrom}." + PREFIX + ".glimpse.vcf.gz", chrom=CHROMS)
    wildcard_constraints:
        chrom="\d+"
    output:
        vcf=temp(TMP_DIR + "/" + PREFIX + ".vcf.gz"),
    shell:
        """
        bcftools concat -Ou {input} | bcftools annotate -x INFO,FORMAT -Oz > {output}
        """

        
rule get_plink:
    input:
        vcf=TMP_DIR + "/" + PREFIX + ".vcf.gz"
    output:
        bed=OUT_DIR + "/plink/" + PREFIX + ".glimpse.bed",
        bim=OUT_DIR + "/plink/" + PREFIX + ".glimpse.bim",
        fam=OUT_DIR + "/plink/" + PREFIX + ".glimpse.fam",
    params:
        px=OUT_DIR + "/plink/" + PREFIX + ".glimpse"
    threads:
        24
    shell:
        """
        plink --vcf {input.vcf} --double-id --real-ref-alleles --make-bed --threads {threads} --out {params.px}
        cat {params.px}.bim | awk '{{print $1"\\t"$1":"$4"\\t"$3"\\t"$4"\\t"$5"\\t"$6}}' > {params.px}.bim_tmp
        mv {params.px}.bim_tmp {output.bim}
        """

          
rule get_pca:
    input:
        bed=OUT_DIR + "/plink/" + PREFIX + ".glimpse.bed",
        bim=OUT_DIR + "/plink/" + PREFIX + ".glimpse.bim",
        fam=OUT_DIR + "/plink/" + PREFIX + ".glimpse.fam",
        lst=OUT_DIR + "/plink/" + PREFIX + ".glimpse.target.snplist"
    output:
        OUT_DIR + "/plink/" + PREFIX + ".glimpse.eigenvec",
        OUT_DIR + "/plink/" + PREFIX + ".glimpse.eigenval",
    params:
        px=OUT_DIR + "/plink/" + PREFIX + ".glimpse"
    threads:
        24
    shell:
        """
        plink --bfile {params.px} --pca header tabs --extract {input.lst} --threads {threads} --out {params.px}
        """


rule plot_pca:
    input:
        OUT_DIR + "/plink/" + PREFIX + ".glimpse.eigenvec",
    output:
        OUT_DIR + "/plink/" + PREFIX + ".glimpse.pca.plot.pdf",
    shell:
        """
        Rscript plotPca.R {input} {output}
        """

