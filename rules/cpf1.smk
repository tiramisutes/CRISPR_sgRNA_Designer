rule cpf1_sgRNAcas9:
    input:
        genome = config["ref"]["genome"],
        geneseq = "Resources/{GeneID}.mRNA.fa"
    output:
        res = report(directory("cpf1_sgRNAcas9_3.0.5/cpf1.sgRNAcas9.report_23.b.{GeneID}.mRNA.fa"), patterns=["{name}_sgRNAcas9_report.xls"], caption="../report/cpf1_sgRNAcas9.rst", category="2-1. CRISPR/sgRNAcas9_cpf1"),
        fres = temp(directory("cpf1_sgRNAcas9_3.0.5/{GeneID}"))
    threads: 1
    log:
        "logs/cpf1_sgRNAcas9_{GeneID}.log"
    params:
        workdir = config["workdir"]
    shell:
        '''
        mkdir -p {output.fres}
        cd {params.workdir}/bin/sgRNAcas9_3.0.5
        perl cpf1_sgRNAcas9_3.0.5.pl -i {params.workdir}/{input.geneseq} -x 23 -l 40 -m 60 -g {input.genome} -o b -t s -v l -n 5 -p {params.workdir}/{output.fres} > {params.workdir}/{log} 2>&1
        cd {params.workdir}/{output.fres}
        mv ./cpf1.sgRNAcas9.report* ../
        '''


rule cpf1_CRISPOR:
    input:
        genomedb = "AddGenome",
        exonseq = "Resources/{GeneID}_exon.fa"
    output:
        targets = report("cpf1_CRISPOR/{GeneID}_targets.xls", caption="../report/cpf1_CRISPOR_targets.rst", category="2-2. CRISPR/CRISPOR_cpf1"),
        offtargets = report("cpf1_CRISPOR/{GeneID}_Offtargets.xls", caption="../report/cpf1_CRISPOR_offtargets.rst", category="2-2. CRISPR/CRISPOR_cpf1"),
        TMP = temp(directory("cpf1_CRISPOR/{GeneID}"))
    threads: 1
    log:
        "logs/cpf1_CRISPOR_{GeneID}.log"
    params:
        genome_name = config["ref"]["genome_name"],
        workdir = config["workdir"],
        R = config["software"]["R"]
    conda:
        "../envs/CRISPOR.yaml"
    shell:
        '''
        {params.R}
        mkdir -p {output.TMP}
        cd {params.workdir}/bin/crisporWebsite
        python crispor.py {params.genome_name} {params.workdir}/{input.exonseq} {params.workdir}/{output.targets} -p TTTN -o {params.workdir}/{output.offtargets} --tempDir {params.workdir}/{output.TMP} -g {params.workdir}/{input.genomedb} >> {params.workdir}/{log} 2>&1
        '''


rule cpf1_list:
    input:
        cpf1_sgRNAcas9 = expand("cpf1_sgRNAcas9_3.0.5/cpf1.sgRNAcas9.report_23.b.{GeneID.GeneID}.mRNA.fa", GeneID=GeneIDs.itertuples()),
        cpf1_CRISPOR = expand("cpf1_CRISPOR/{GeneID.GeneID}", GeneID=GeneIDs.itertuples())
    output:
        "cpf1_list/All_cpf1.txt"
    threads: 1
    shell:
        '''
        ls {input.cpf1_sgRNAcas9} {input.cpf1_CRISPOR} > {output}
        '''