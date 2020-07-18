rule genome_exon:
    input:
        fasta = config["ref"]["genome"],
        gff = config["ref"]["genome_gff"]
    output:
        bed = "Resources/{}.bed".format(config["ref"]["genome_name"]),
        exon = "Resources/{}_exon.fa".format(config["ref"]["genome_name"])
    message:
        "Get genome exon sequencs"
    conda:
        "../envs/bedtools.yaml"
    cache: True
    shell:
        '''
        cat {input.gff} | convert2bed --input=gff - > {output.bed}
        cat {output.bed} | awk -F"\\t" 'BEGIN{{OFS="\\t"}} $8=="CDS" {{print $0}}' | bedtools getfasta -fi {input.fasta} -bed - -fo {output.exon} -name+
        '''


rule gene_exon_seq:
    input:
        exon = "Resources/{}_exon.fa".format(config["ref"]["genome_name"])
    output:
        exon = "Resources/{GeneID}_exon.fa",
        ABEexon = "Resources/{GeneID}_ABE_exon.fa",
        CBEexon = "Resources/{GeneID}_CBE_exon.fa"
    params:
        id = "{GeneID}",
        ABE = config["windows"]["ABE"],
        CBE = config["windows"]["CBE"]
    threads: 1
    shell:
        '''
        cat {input.exon} | \\grep -A 1 "{params.id}" | sed "s/--//g" | grep -v "^$" > {output.exon}
        ./scripts/Filter_exon.py -i {output.exon} -o {output.ABEexon} --module ABE {params.ABE}
        ./scripts/Filter_exon.py -i {output.exon} -o {output.CBEexon} --module CBE {params.CBE}
        '''

rule gene_mRNA_seq:
    input:
        mRNA = config["ref"]["genome_mRNA"]
    output:
        "Resources/{GeneID}.mRNA.fa"
    params:
        id = "{GeneID}"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    shell:
        '''
        samtools faidx {input.mRNA} {params.id} | sed "s/\\.[0-9]$//g" > {output}
        '''


rule AddGenome:
    input:
        fasta = config["ref"]["genome"],
        gff = config["ref"]["genome_gff"]
    output:
        directory("AddGenome")
    params:
        workdir = config["workdir"],
        genome_name = config["ref"]["genome_name"],
        UCSCscripts = config["software"]["UCSCscripts"]
    cache: True
    log:
        "logs/AddGenome.log"
    conda:
        "../envs/CRISPOR.yaml"
    shell:
        '''
        mkdir {output}
        {params.UCSCscripts}
        export PATH="$PATH:{params.workdir}/bin/crisporWebsite/bin"
        export PATH="$PATH:{params.workdir}/bin/crisporWebsite/bin/Linux"

        cd {params.workdir}/bin/crisporWebsite/tools
        if [[ -d "/tmp/{params.genome_name}" ]]
        then
            echo -e "###########################################################################\n"  > {params.workdir}/{log} 2>&1
            echo "***! The /tmp/{params.genome_name} directory exist and delete." >> {params.workdir}/{log} 2>&1
            echo -e "###########################################################################\n" >> {params.workdir}/{log} 2>&1
            rm -rf /tmp/{params.genome_name}
            ./crisprAddGenome --baseDir={params.workdir}/{output} fasta {input.fasta} --gff={input.gff} \
                --desc '{params.genome_name}|{params.genome_name}|{params.genome_name}_HZAU| 2019 (Gh_TM-1_HZAU)' >> {params.workdir}/{log} 2>&1
        else
            echo -e "###########################################################################" > {params.workdir}/{log} 2>&1
            echo "                (-: Now, add genome" >> {params.workdir}/{log} 2>&1
            echo -e "###########################################################################\n" >> {params.workdir}/{log} 2>&1
            ./crisprAddGenome --baseDir={params.workdir}/{output} fasta {input.fasta} --gff={input.gff} \
                --desc '{params.genome_name}|{params.genome_name}|{params.genome_name}_HZAU| 2019 (Gh_TM-1_HZAU)' >> {params.workdir}/{log} 2>&1
        fi
        '''
