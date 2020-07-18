rule CBE_BE_Designer:
    input:
        genome = config["ref"]["genome"],
        CBEexonseq = "Resources/{GeneID}_CBE_exon.fa"
    output:
        #directory("CBE_BE-Designer/{GeneID}")
        report(directory("CBE_BE-Designer/{GeneID}"), patterns=["{name}_CBE_results.txt"], caption="../report/CBE_BE-Designer.rst", category="5-1. CRISPR/CBE")
    threads: 1
    log:
        "logs/CBE_{GeneID}.log"
    params:
        workdir = config["workdir"],
        CBE = config["windows"]["CBE"],
        CasOFFinder = config["software"]["CasOFFinder"]
    shell:
        '''
        {params.CasOFFinder}
        mkdir -p {output}

        cd {params.workdir}/bin/BE-Designer
        for i in `cat {params.workdir}/{input.CBEexonseq} | xargs`
        do
        python3 be_designer.py ${{i#*^}} {input.genome} {params.CBE} -r C -m T --output {params.workdir}/{output}/${{i%^*}} --prefix {params.workdir}/{output}/${{i%^*}} >> {params.workdir}/{log} 2>&1
        cat {params.workdir}/{output}/${{i%^*}} \
        | awk -F"\\t" 'BEGIN{{OFS="\\t"}} {{if($0 ~ /CRISPR/) {{print $0,"codon_0_PH_Change";}} else if($10 !~ /R|K|H/ && $11 ~ /R|K|H/ || $10 ~ /R|K|H/ && $11 !~ /R|K|H/) {{print $0,"PH>7";}} else if($10 !~ /D|E/ && $11 ~ /D|E/ || $10 ~ /D|E/ && $11 !~ /D|E/) {{print $0,"PH<7";}} else if($10 !~ /X/ && $11 ~ /X/ || $10 ~ /X/ && $11 !~ /X/) {{print $0,"Find Stop";}} else {{print $0,"PH=7";}}}}' \
        | awk -F"\\t" -v a="${{i%^*}}" 'BEGIN{{OFS="\\t"}} {{if($1 ~ /CRISPR/) {{print "ExonID","Editor",$0;}} else {{print a,"CBE",$0;}}}}' \
        | awk -F"\\t" 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$22,$14,$15,$16,$17,$18,$19,$20,$21}}' > {params.workdir}/{output}/${{i%^*}}_CBE_results.txt
        \\rm {params.workdir}/{output}/${{i%^*}}
        done
        '''


rule CBE_Merge_Gene:
    input:
        "CBE_BE-Designer/{GeneID}"
    output:
        report("CBE_Merge_Gene/All_{GeneID}_CBE_sgRNA.csv", caption="../report/CBE_Merge_Gene.rst", category="5-2. CBE_Merge_gene")
    threads: 1
    params:
        "_CBE_results.txt"
    conda:
        "../envs/R.yaml"
    log:
        "logs/CBE_Merge_{GeneID}.log"
    script:
        "../scripts/BE_Merge_Gene.R"


rule CBE_Merge:
    input:
        expand("CBE_Merge_Gene/All_{GeneID.GeneID}_CBE_sgRNA.csv", GeneID=GeneIDs.itertuples())
    output:
        report("CBE_Merge/All_CBE_sgRNA.csv", caption="../report/CBE_Merge.rst", category="5-3. CBE_Merge")
    threads: 1
    conda:
        "../envs/R.yaml"
    log:
        "logs/CBE_Merge.log"
    script:
        "../scripts/BE_Merge.R"


rule CBE_sgRNAcas9:
    input:
        resf = "BE_sgRNAcas9_3.0.5/sgRNAcas9.report_20.s.{GeneID}.mRNA.fa"
    output:
        CBE = report("CBE_sgRNAcas9_3.0.5/{GeneID}_sgRNAcas9_CBE_report.xls", caption="../report/CBE_sgRNAcas9.rst", category="5-4. CRISPR/CBE")
    threads: 1
    log:
        "logs/CBE_sgRNAcas9_{GeneID}.log"
    params:
        file = "{GeneID}.mRNA.fa.sgRNAcas9_report.xls",
        workdir = config["workdir"]
    shell:
        '''
        bash {params.workdir}/scripts/sgRNAcas9toBE.sh CBE {params.workdir}/{input.resf}/{params.file} {params.workdir}/{output.CBE}
        '''


rule CBE_list:
    input:
        CBE_BE_Designer = "CBE_Merge/All_CBE_sgRNA.csv",
        CBE_sgRNAcas9 = expand("CBE_sgRNAcas9_3.0.5/{GeneID.GeneID}_sgRNAcas9_CBE_report.xls", GeneID=GeneIDs.itertuples())
    output:
        "CBE_list/All_CBE.txt"
    threads: 1
    shell:
        '''
        ls {input.CBE_BE_Designer} {input.CBE_sgRNAcas9} > {output}
        '''

