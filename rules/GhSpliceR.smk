rule GhSpliceR:
    output:
        report("GhSpliceR/{GeneID}_NGG_sgRNA.csv", caption="../report/GhSpliceR.rst", category="6. GhSpliceR")
    threads: 1
    params:
        id = '"{GeneID}"',
        pam = "NGG",
        enzyme = '"CBE or ABE"',
        splice_site = '"splice-donors and splice-acceptors"',
        R = config["software"]["R"]
    log:
        "logs/GhSpliceR_{GeneID}.log"
    shell:
        '''
        {params.R}
        Rscript ./scripts/GhSpliceR.R {params.id} {params.pam} {params.enzyme} {params.splice_site} {output} > {log} 2>&1
        '''

rule GhSpliceR_Merge:
    input:
        expand("GhSpliceR/{GeneID.GeneID}_NGG_sgRNA.csv", GeneID=GeneIDs.itertuples())
    output:
        "GhSpliceR_Merge/All_GhSpliceR.csv"
    threads: 1
    log:
        "logs/GhSpliceR_Merge.log"
    shell:
        '''
        cat {input} > {output}
        '''