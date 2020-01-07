This Workflow was modified from the Snakemake workflow provided by [phiweger](https://github.com/phiweger)  of [*Nanotext*](https://github.com/phiweger/nanotext), and with help from [C. Titus Brown](https://github.com/ctb) and [Luiz Irber](https://github.com/luizirber) from the [dib-lab](https://github.com/dib-lab)

# Running the Workflow

To run the workflow:

## Sequencing Reads

Create symbolic links to paired read files in `reads/`.

## Config File 

To setup HMM database:

Navigate to the directory where you wish for the database to live

Use `wget` and `hmmer v3.1b`

```bash

# download and set up HMM database
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm

```

Then, set the absolute path for the pfam database on the system in `config.json`.

## Execute

Execute the following using this folder as the working directory:

```bash

snakemake -j 100 \
--cluster-config cluster.json \
--cluster "sbatch -A {cluster.account} \
--mem={cluster.mem} -t {cluster.time} \
-c {threads}" \
--configfile config.json \
--use-conda

```

The `results/` folder will  be created with genome annotations in `data/` and a `log/` if anything goes wrong.

When running this workflow, check for concordance between the file extension of input read files and code in the Snakefile.

