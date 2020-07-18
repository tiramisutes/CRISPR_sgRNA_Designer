import os
import re
import yaml
import pathlib
from typing import List
import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


#======================================================
# cluster directory
#======================================================
logs = pathlib.Path("./logs/cluster")

if logs.is_dir():
    pass
    #print("#"*40 + "\n" + "        The logs exist." + "\n" + "#"*40 + "\n")
else:
    #print("#"*40 + "\n" + "Now, create the logs directory." + "\n")
    os.system('mkdir -p ./logs/cluster')
    #print("#"*40 + "\n")

#======================================================
# Config files
#======================================================
##### load config and GeneID sheets #####

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

GeneIDs = pd.read_table(config["GeneID"]).set_index("GeneID", drop=False)
validate(GeneIDs, schema="schemas/GeneID.schema.yaml")

##### target rules #####

rule all:
    input:
        "Cas9_list/All_Cas9.txt" if config["module"]["Cas9"] else [],
        "cpf1_list/All_cpf1.txt" if config["module"]["cpf1"] else [],
        "c2c1_list/All_c2c1.txt" if config["module"]["c2c1"] else [],
        #expand("c2c1_sgRNAcas9_3.0.5/{name}",name=GeneIDs['GeneID'].tolist()),
        "ABE_list/All_ABE.txt" if config["module"]["ABE"] else [],
        "CBE_list/All_CBE.txt" if config["module"]["CBE"] else [],
        "GhSpliceR_Merge/All_GhSpliceR.csv" if config["module"]["GhSpliceR"] else [],
        "GhiSTOP_Merge/All_GhiSTOP.csv" if config["module"]["GhiSTOP"] else []
            


##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


##### setup report #####

report: "report/workflow.rst"

## select which CRISPR to desiger
include: "rules/ref.smk"

if config["module"]["Cas9"]:
    include: "rules/Cas9.smk"
if config["module"]["cpf1"]:
    include: "rules/cpf1.smk"
if config["module"]["c2c1"]:
    include: "rules/c2c1.smk"
if config["module"]["ABE"]:
    include: "rules/BE.smk"
    include: "rules/ABE.smk"
if config["module"]["CBE"]:
    include: "rules/BE.smk"
    include: "rules/CBE.smk"
if config["module"]["GhSpliceR"]:
    include: "rules/GhSpliceR.smk"
if config["module"]["GhiSTOP"]:
    include: "rules/GhiSTOP.smk"
