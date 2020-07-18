rule c2c1:
    input:
        genome = config["ref"]["genome"],
        geneseq = "Resources/{GeneID}.mRNA.fa"
    output:
        res = report(directory("c2c1_sgRNAcas9_3.0.5/c2p2.sgRNAcas9.report_20.b.{GeneID}.mRNA.fa"), patterns=["{name}_sgRNAcas9_report.xls"], caption="../report/c2c1.rst", category="3. CRISPR/c2c1"),
        fres = temp(directory("c2c1_sgRNAcas9_3.0.5/{GeneID}"))
    threads: 1
    log:
        "logs/c2c1_sgRNAcas9_{GeneID}.log"
    params:
        workdir = config["workdir"]
    shell:
        '''
        mkdir -p {output.fres}
        cd {params.workdir}/bin/sgRNAcas9_3.0.5
        perl c2p2_sgRNAcas9_3.0.5.pl -i {params.workdir}/{input.geneseq} -x 20 -l 40 -m 60 -g {input.genome} -o b -t s -v l -n 5 -p {params.workdir}/{output.fres} > {params.workdir}/{log} 2>&1
        cd {params.workdir}/{output.fres}
        mv ./c2p2.sgRNAcas9.report* ../
        '''


rule c2c1_list:
    input:
        expand("c2c1_sgRNAcas9_3.0.5/c2p2.sgRNAcas9.report_20.b.{GeneID.GeneID}.mRNA.fa", GeneID=GeneIDs.itertuples())
    output:
        "c2c1_list/All_c2c1.txt"
    threads: 1
    shell:
        '''
        ls {input} > {output}
        '''