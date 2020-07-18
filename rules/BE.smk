rule BE_sgRNAcas9:
    input:
        genome = config["ref"]["genome"],
        geneseq = "Resources/{GeneID}.mRNA.fa"
    output:
        res = directory("BE_sgRNAcas9_3.0.5/sgRNAcas9.report_20.s.{GeneID}.mRNA.fa"),
        fres = temp(directory("BE_sgRNAcas9_3.0.5/{GeneID}"))
    threads: 1
    log:
        "logs/BE_sgRNAcas9_{GeneID}.log"
    params:
        workdir = config["workdir"]
    shell:
        '''
        mkdir -p {output.fres}
        cd {params.workdir}/bin/sgRNAcas9_3.0.5
        perl sgRNAcas9_3.0.5.pl -i {params.workdir}/{input.geneseq} -x 20 -l 40 -m 60 -g {input.genome} -o s -t s -v l -n 5 -p {params.workdir}/{output.fres} > {params.workdir}/{log} 2>&1
        cd {params.workdir}/{output.fres}
        mv ./sgRNAcas9.report* ../
        '''
