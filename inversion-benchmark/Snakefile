# Copyright 2021 Ronald J. Nowling
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

## process 16 Anopheles genomes data
rule extract_ag16g_zipfile:
    input:
        ag16g="data/raw_data/doi_10.5061_dryad.f4114__v1.zip"
    output:
        vcf_tar_file="data/raw_data/VCFfile4DRYAD.tar.gz"
    threads:
        1
    shell:
        "unzip -d data/raw_data {input.ag16g} VCFfile4DRYAD.tar.gz"

rule extract_ag16g_tarfile:
    input:
        tar_file="data/raw_data/VCFfile4DRYAD.tar.gz"
    output:
        vcf_files=expand("data/raw_data/VCFfile4DRYAD/AGC_refHC_bialSNP_AC2_2DPGQ.{chrom}_V2.CHRcode2.DRYAD.vcf.gz", chrom=["2L", "2R", "3L"])
    threads:
        1
    shell:
        "tar -C data/raw_data -xzvf {input.tar_file}"

rule filter_ag16g_vcf:
    input:
        vcf="data/raw_data/VCFfile4DRYAD/AGC_refHC_bialSNP_AC2_2DPGQ.{chrom}_V2.CHRcode2.DRYAD.vcf.gz"
    output:
        vcf="data/ag16g/ag16g_{chrom}_gambiae_coluzzii.vcf.gz"
    threads:
        1
    shell:
        "vcftools --gzvcf {input.vcf} --recode --stdout --keep sample_lists/ag16g_gambiae_coluzzii_ids.txt --remove-indels --maf 0.01 | gzip -c > {output.vcf}"

rule ag16g_to_bed:
    input:
        vcf="data/ag16g/ag16g_{chrom}.vcf.gz"
    output:
        bed="data/ag16g/ag16g_{chrom}.bed"
    threads:
        1
    shell:
        "plink1.9 --vcf {input.vcf} --make-bed --allow-extra-chr --out data/ag16g/ag16g_{wildcards.chrom}"

rule ag16g_to_raw:
    input:
        bed="data/ag16g/ag16g_{chrom}.bed"
    output:
        raw="data/ag16g/ag16g_{chrom}.raw"
    threads:
        1
    shell:
        "plink1.9 --bfile data/ag16g/ag16g_{wildcards.chrom} --recode A --allow-extra-chr --out data/ag16g/ag16g_{wildcards.chrom}"

rule ag16g_to_inveRsion:
    input:
        raw="data/ag16g/ag16g_{chrom}.raw"
    output:
        inveRsion="data/ag16g/ag16g_{chrom}.inveRsion"
    threads:
        1
    shell:
        "scripts/raw_to_inveRsion --input-raw {input.raw} --output-txt {output.inveRsion}"
        
## process Drosophila Genetics Reference Panel v2 VCFs
rule filter_dgrp2_vcf:
    input:
        vcf="data/raw_data/dgrp2.vcf"
    output:
        filtered_vcf="data/dgrp2/dgrp2.biallelic.vcf"
    threads:
        1
    shell:
        "vcftools --vcf {input.vcf} --min-alleles 2 --max-alleles 2 --remove-indels --maf 0.01 --recode --stdout --remove-indv line_348 --remove-indv line_350 --remove-indv line_358 --remove-indv line_385 --remove-indv line_392 --remove-indv line_395 --remove-indv line_399 > {output.filtered_vcf}"

rule split_dgrp2_by_chrom:
    input:
        vcf="data/dgrp2/dgrp2.biallelic.vcf"
    output:
        chrom_vcf="data/dgrp2/dgrp2_{chrom}.biallelic.vcf"
    threads:
        1
    shell:
        "vcftools --vcf {input.vcf} --chr {wildcards.chrom} --recode --stdout > {output.chrom_vcf}"

rule dgrp2_to_bed:
    input:
        vcf="data/dgrp2/dgrp2_{chrom}.biallelic.vcf"
    output:
        bed="data/dgrp2/dgrp2_{chrom}.biallelic.bed"
    threads:
        1
    shell:
        "plink1.9 --vcf {input.vcf} --make-bed --allow-extra-chr --out data/dgrp2/dgrp2_{wildcards.chrom}.biallelic"

rule dgrp2_to_raw:
    input:
        bed="data/dgrp2/dgrp2_{chrom}.biallelic.bed"
    output:
        raw="data/dgrp2/dgrp2_{chrom}.biallelic.raw"
    threads:
        1
    shell:
        "plink1.9 --bfile data/dgrp2/dgrp2_{wildcards.chrom}.biallelic --recode A --allow-extra-chr --out data/dgrp2/dgrp2_{wildcards.chrom}.biallelic"

rule dgrp2_to_inveRsion:
    input:
        raw="data/dgrp2/dgrp2_{chrom}.raw"
    output:
        inveRsion="data/dgrp2/dgrp2_{chrom}.inveRsion"
    threads:
        1
    shell:
        "scripts/raw_to_inveRsion --input-raw {input.raw} --output-txt {output.inveRsion}"
        
## Process 1000 Anopheles Genomes VCFs
rule select_ag1000g_samples:
    input:
        vcf="data/raw_data/ag1000g.phase1.ar3.pass.biallelic.{chrom}.vcf.gz"
    output:
        vcf="data/ag1000g/ag1000g_{chrom}_bfaso.vcf.gz"
    threads:
        1
    shell:
        "vcftools --gzvcf {input.vcf} --recode --stdout --keep sample_lists/ag1000g_bfm_bfs_ids.txt --remove-indels --maf 0.01 | gzip -c > {output.vcf}"

rule select_ag1000g_gambiae_samples:
    input:
        vcf="data/ag1000g/ag1000g_{chrom}_bfaso.vcf.gz"
    output:
        vcf="data/ag1000g/ag1000g_{chrom}_bfaso_gambiae.vcf.gz"
    threads:
        1
    shell:
        "vcftools --gzvcf {input.vcf} --recode --stdout --keep sample_lists/ag1000g_bfs_ids.txt --maf 0.01 | gzip -c > {output.vcf}"

rule select_ag1000g_coluzzii_samples:
    input:
        vcf="data/ag1000g/ag1000g_{chrom}_bfaso.vcf.gz"
    output:
        vcf="data/ag1000g/ag1000g_{chrom}_bfaso_coluzzii.vcf.gz"
    threads:
        1
    shell:
        "vcftools --gzvcf {input.vcf} --recode --stdout --keep sample_lists/ag1000g_bfm_ids.txt --maf 0.01 | gzip -c > {output.vcf}"

rule ag1000g_to_bed:
    input:
        vcf="data/ag1000g/ag1000g_{chrom}.vcf.gz"
    output:
        bed="data/ag1000g/ag1000g_{chrom}.bed"
    threads:
        1
    shell:
        "plink1.9 --vcf {input.vcf} --set-missing-var-ids '@_#_\$1_\$2' --make-bed --allow-extra-chr --out data/ag1000g/ag1000g_{wildcards.chrom}"

rule ag1000g_to_raw:
    input:
        bed="data/ag1000g/ag1000g_{chrom}.bed"
    output:
        raw="data/ag1000g/ag1000g_{chrom}.raw"
    threads:
        1
    shell:
        "plink1.9 --bfile data/ag1000g/ag1000g_{wildcards.chrom} --recode A --allow-extra-chr --out data/ag1000g/ag1000g_{wildcards.chrom}"

rule ag1000g_to_inveRsion:
    input:
        raw="data/ag1000g/ag1000g_{chrom}.raw"
    output:
        inveRsion="data/ag1000g/ag1000g_{chrom}.inveRsion"
    threads:
        1
    shell:
        "scripts/raw_to_inveRsion --input-raw {input.raw} --output-txt {output.inveRsion}"
        
## Top-level rules
rule check_inputs:
    input:
        dgrp="data/raw_data/dgrp2.vcf",
        ag1000g=expand("data/raw_data/ag1000g.phase1.ar3.pass.biallelic.{chrom}.vcf.gz", chrom=["2L", "3L", "2R"]),
        ag16g="data/raw_data/doi_10.5061_dryad.f4114__v1.zip"

rule prepare_dgrp2:
    input:
        dgrp=expand("data/dgrp2/dgrp2_{chrom}.biallelic.{format}", chrom=["2L", "2R", "3R", "3L"], format=["bed", "raw", "vcf", "inveRsion"])

rule prepare_ag1000g:
    input:
        ag1000g_bfaso=expand("data/ag1000g/ag1000g_{chrom}_bfaso.{format}", format=["bed", "raw", "vcf.gz"], chrom=["2L", "2R", "3L"]),
        ag1000g_bfaso_gambiae=expand("data/ag1000g/ag1000g_{chrom}_bfaso_gambiae.{format}", chrom=["2L", "2R", "3L"], format=["bed", "raw", "vcf.gz", "inveRsion"]),
        ag1000g_bfaso_coluzii=expand("data/ag1000g/ag1000g_{chrom}_bfaso_coluzzii.{format}", chrom=["2L", "2R", "3L"], format=["bed", "raw", "vcf.gz", "inveRsion"])
        
rule prepare_ag16g:
    input:
        ag16g=expand("data/ag16g/ag16g_{chrom}_gambiae_coluzzii.{format}", chrom=["2L", "2R", "3L"], format=["bed", "raw", "vcf.gz", "inveRsion"])
