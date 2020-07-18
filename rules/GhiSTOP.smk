rule GhiSTOP:
    output:
        report("GhiSTOP/{GeneID}_iSTOP.csv", caption="../report/GhiSTOP.rst", category="7. GhiSTOP")
    threads: 1
    params:
        id = '"{GeneID}"',
        R = config["software"]["R"]
    log:
        "logs/GhiSTOP_{GeneID}.log"
    shell:
        '''
        {params.R}
        Rscript ./scripts/GhiSTOP.R {params.id} {output} > {log} 2>&1
        '''

rule GhiSTOP_Merge:
    input:
        expand("GhiSTOP/{GeneID.GeneID}_iSTOP.csv", GeneID=GeneIDs.itertuples())
    output:
        "GhiSTOP_Merge/All_GhiSTOP.csv"
    threads: 1
    log:
        "logs/GhiSTOP_Merge.log"
    shell:
        '''
        cat {input} > {output}
        '''