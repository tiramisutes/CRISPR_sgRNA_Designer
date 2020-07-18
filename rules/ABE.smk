rule ABE_BE_Designer:
    input:
        genome = config["ref"]["genome"],
        ABEexonseq = "Resources/{GeneID}_ABE_exon.fa"
    output:
        #directory("ABE_BE-Designer/{GeneID}")
        report(directory("ABE_BE-Designer/{GeneID}"), patterns=["{name}_ABE_results.txt"], caption="../report/ABE_BE-Designer.rst", category="4-1. CRISPR/ABE")
    threads: 1
    log:
        "logs/ABE_{GeneID}.log"
    params:
        workdir = config["workdir"],
        ABE = config["windows"]["ABE"],
        CasOFFinder = config["software"]["CasOFFinder"]
    shell:
        '''
        {params.CasOFFinder}
        mkdir -p {output}

        cd {params.workdir}/bin/BE-Designer
        for i in `cat {params.workdir}/{input.ABEexonseq} | xargs`
        do
        python3 be_designer.py ${{i#*^}} {input.genome} {params.ABE} -r A -m G --output {params.workdir}/{output}/${{i%^*}} --prefix {params.workdir}/{output}/${{i%^*}} >> {params.workdir}/{log} 2>&1
        cat {params.workdir}/{output}/${{i%^*}} \
        | awk -F"\\t" 'BEGIN{{OFS="\\t"}} {{if($0 ~ /CRISPR/) {{print $0,"codon_0_PH_Change";}} else if($10 !~ /R|K|H/ && $11 ~ /R|K|H/ || $10 ~ /R|K|H/ && $11 !~ /R|K|H/) {{print $0,"PH>7";}} else if($10 !~ /D|E/ && $11 ~ /D|E/ || $10 ~ /D|E/ && $11 !~ /D|E/) {{print $0,"PH<7";}} else if($10 !~ /X/ && $11 ~ /X/ || $10 ~ /X/ && $11 !~ /X/) {{print $0,"Find Stop";}} else {{print $0,"PH=7";}}}}' \
        | awk -F"\\t" -v a="${{i%^*}}" 'BEGIN{{OFS="\\t"}} {{if($1 ~ /CRISPR/) {{print "ExonID","Editor",$0;}} else {{print a,"ABE",$0;}}}}' \
        | awk -F"\\t" 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$22,$14,$15,$16,$17,$18,$19,$20,$21}}' > {params.workdir}/{output}/${{i%^*}}_ABE_results.txt
        \\rm {params.workdir}/{output}/${{i%^*}}
        done
        '''


rule ABE_Merge_Gene:
    input:
        "ABE_BE-Designer/{GeneID}"
    output:
        report("ABE_Merge_Gene/All_{GeneID}_ABE_sgRNA.csv", caption="../report/ABE_Merge_Gene.rst", category="4-2. ABE_Merge_gene")
    threads: 1
    params:
        "_ABE_results.txt"
    conda:
        "../envs/R.yaml"
    log:
        "logs/ABE_Merge_{GeneID}.log"
    script:
        "../scripts/BE_Merge_Gene.R"


rule ABE_Merge:
    input:
        expand("ABE_Merge_Gene/All_{GeneID.GeneID}_ABE_sgRNA.csv", GeneID=GeneIDs.itertuples())
    output:
        report("ABE_Merge/All_ABE_sgRNA.csv", caption="../report/ABE_Merge.rst", category="4-3. ABE_Merge")
    threads: 1
    conda:
        "../envs/R.yaml"
    log:
        "logs/ABE_Merge.log"
    script:
        "../scripts/BE_Merge.R"


rule ABE_sgRNAcas9:
    input:
        resf = "BE_sgRNAcas9_3.0.5/sgRNAcas9.report_20.s.{GeneID}.mRNA.fa"
    output:
        ABE = report("ABE_sgRNAcas9_3.0.5/{GeneID}_sgRNAcas9_ABE_report.xls", caption="../report/ABE_sgRNAcas9.rst", category="4-4. CRISPR/ABE")
    threads: 1
    log:
        "logs/ABE_sgRNAcas9_{GeneID}.log"
    params:
        file = "{GeneID}.mRNA.fa.sgRNAcas9_report.xls",
        workdir = config["workdir"]
    shell:
        '''
        bash {params.workdir}/scripts/sgRNAcas9toBE.sh ABE {params.workdir}/{input.resf}/{params.file} {params.workdir}/{output.ABE}
        '''


rule ABE_list:
    input:
        ABE_BE_Designer = "ABE_Merge/All_ABE_sgRNA.csv",
        ABE_sgRNAcas9 = expand("ABE_sgRNAcas9_3.0.5/{GeneID.GeneID}_sgRNAcas9_ABE_report.xls", GeneID=GeneIDs.itertuples())
    output:
        "ABE_list/All_ABE.txt"
    threads: 1
    shell:
        '''
        ls {input.ABE_BE_Designer} {input.ABE_sgRNAcas9} > {output}
        '''

