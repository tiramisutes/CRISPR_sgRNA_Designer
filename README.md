# Snakemake workflow: CRISPR sgRNA Designer

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/rna-seq-star-deseq2.svg?branch=master)](https://travis-ci.org/snakemake-workflows/rna-seq-star-deseq2)
[![Snakemake-Report](https://img.shields.io/badge/snakemake-report-green.svg)](https://cdn.rawgit.com/snakemake-workflows/rna-seq-star-deseq2/master/.test/report.html)

This workflow performs sgRNA designer for CRISPR genome edit, which includes Cas9, Cas12a/cpf1, Cas12b/c2c1, CBE and ABE editors.

## Authors

* Zhongping Xu (@hopetogy), http://tiramisutes.github.io/

## Usage

### Simple

#### Step 1: Install workflow / git clone

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).

Git clone this repository using follows command 👇

    git clone --recursive https://github.com/tiramisutes/CRISPR_sgRNA_Designer.git

#### Step 2: Prepare data and software

- **genome.fasta**: the genome fasta file for study species.

- **genome.gff3**: the annotation file in gff3 format.

- **genome.cds.fa** the genome cds file.

#### Step 3: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

#### Step 4: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

or

    snakemake --use-conda --jobs 100 --cluster-config cluster.json --cluster "bsub -q {cluster.queue} -o {cluster.output} -e {cluster.error}"

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.

#### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.
An example (using some trivial test data) can be seen [here](https://cdn.rawgit.com/snakemake-workflows/rna-seq-star-deseq2/master/.test/report.html).

### Advanced

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for the concrete project/run on your machine.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).
4. Modify the config, and any necessary sheets (and probably the workflow) as needed.
5. Commit any changes and push the project-branch to your fork on github.
6. Run the analysis.
7. Optional: Merge back any valuable and generalizable changes to the [upstream repo](https://github.com/snakemake-workflows/rna-seq-star-deseq2) via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
8. Optional: Push results (plots/tables) to the remote branch on your fork.
9. Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
10. Optional: Delete the local clone/workdir to free space.
