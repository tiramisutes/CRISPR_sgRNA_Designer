rule Cas9_sgRNAcas9:
    input:
        genome = config["ref"]["genome"],
        geneseq = "Resources/{GeneID}.mRNA.fa"
    output:
        res = report(directory("Cas9_sgRNAcas9_3.0.5/sgRNAcas9.report_20.b.{GeneID}.mRNA.fa"), patterns=["{name}_sgRNAcas9_report.xls"], caption="../report/Cas9_sgRNAcas9.rst", category="1-1. CRISPR/Cas9"),
        fres = temp(directory("Cas9_sgRNAcas9_3.0.5/{GeneID}"))
    threads: 1
    log:
        "logs/Cas9_sgRNAcas9_{GeneID}.log"
    params:
        workdir = config["workdir"]
    shell:
        '''
        mkdir -p {output.fres}
        cd {params.workdir}/bin/sgRNAcas9_3.0.5
        perl sgRNAcas9_3.0.5.pl -i {params.workdir}/{input.geneseq} -x 20 -l 40 -m 60 -g {input.genome} -o b -t s -v l -n 5 -p {params.workdir}/{output.fres} > {params.workdir}/{log} 2>&1
        cd {params.workdir}/{output.fres}
        mv ./sgRNAcas9.report* ../
        '''


rule Cas9_CRISPOR:
    input:
        genomedb = "AddGenome",
        exonseq = "Resources/{GeneID}_exon.fa"
    output:
        targets = report("Cas9_CRISPOR/{GeneID}_targets.xls", caption="../report/Cas9_CRISPOR_targets.rst", category="1-2. CRISPR/Cas9"),
        offtargets = report("Cas9_CRISPOR/{GeneID}_Offtargets.xls", caption="../report/Cas9_CRISPOR_offtargets.rst", category="1-2. CRISPR/Cas9"),
        TMP = temp(directory("Cas9_CRISPOR/{GeneID}"))
    threads: 1
    log:
        "logs/Cas9_CRISPOR_{GeneID}.log"
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
        python crispor.py {params.genome_name} {params.workdir}/{input.exonseq} {params.workdir}/{output.targets} -p NGG -o {params.workdir}/{output.offtargets} --tempDir {params.workdir}/{output.TMP} -g {params.workdir}/{input.genomedb} > {params.workdir}/{log} 2>&1
        '''


if config["quickly"]:
    pass
else:
    rule Cas9_Cas_Designer:
        input:
            genome = config["ref"]["genome"],
            exonseq = "Resources/{GeneID}_exon.fa"
        output:
            sgRNA = report("Cas9_Cas_Designer/{GeneID}", caption="../report/Cas9_Cas_Designer.rst", category="1-3. CRISPR/Cas9"),
            config = "Cas9_Cas_Designer/{GeneID}_config.txt"
        threads: 1
        log:
            "logs/Cas9_Cas_Designer_{GeneID}.log"
        params:
            CasOFFinder = config["software"]["CasOFFinder"],
            CasDesigner = config["software"]["CasDesigner"]
        shell:
            '''
            {params.CasOFFinder}
            {params.CasDesigner}
            echo -e {input.genome}"\\n"{input.exonseq}"\\n"20"\\n"NGG"\\n"NRG"\\n"5"\\n"2"\\n"2"\\n" > {output.config}
            cas-designer {output.config} >> {log} 2>&1
            '''


if config["quickly"]:
    rule Cas9_list:
        input:
            Cas9_sgRNAcas9 = expand("Cas9_sgRNAcas9_3.0.5/sgRNAcas9.report_20.b.{GeneID.GeneID}.mRNA.fa", GeneID = GeneIDs.itertuples()),
            Cas9_CRISPOR = expand("Cas9_CRISPOR/{GeneID.GeneID}", GeneID = GeneIDs.itertuples())
        output:
            "Cas9_list/All_Cas9.txt"
        threads: 1
        shell:
            '''
            ls {input.Cas9_sgRNAcas9} {input.Cas9_CRISPOR} > {output}
            '''
else:
    rule Cas9_list:
        input:
            Cas9_sgRNAcas9 = expand("Cas9_sgRNAcas9_3.0.5/sgRNAcas9.report_20.b.{GeneID.GeneID}.mRNA.fa", GeneID = GeneIDs.itertuples()),
            Cas9_CRISPOR = expand("Cas9_CRISPOR/{GeneID.GeneID}", GeneID = GeneIDs.itertuples()),
            Cas9_Cas_Designer = expand("Cas9_Cas_Designer/{GeneID.GeneID}", GeneID = GeneIDs.itertuples())
        output:
            "Cas9_list/All_Cas9.txt"
        threads: 1
        shell:
            '''
            ls {input.Cas9_sgRNAcas9} {input.Cas9_CRISPOR} {input.Cas9_Cas_Designer} > {output}
            '''