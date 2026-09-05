### Interactions with polygenic background impact quantitative traits in the UK Biobank, _medRxiv_ (2025): _accompanying code_

- **Authors of the paper:** Lino A.F. Ferreira, Sile Hu, Robert W. Davies and Simon R. Myers
- **Authors of the code:** Lino A.F. Ferreira and Simon R. Myers
- **[Link to paper](https://doi.org/10.1101/2025.11.14.25340263)**

The `data/ folder` contains a skeleton folder structure indicating the external data required to run
the pipeline. `The README.md` file in that folder describes the files necessary and how they can be
obtained.

The `workflow/` folder contains a Snakemake pipeline that executes the computational analyses reported
in the paper. The rules specified in the Snakefile and in the `.smk` files within the subfolder `rules/`
define the workflow which when executed calls the scripts in the subfolder `scripts/` as described in
the rules.

The code is organised in five sections:

0. **Data setup:** preparation of phenotypes (real and simulated), sample and variant lists, and
   annotations.
1. **Standard GWAS:** running of a standard (i.e., additive) GWAS for all phenotyes.
2. **PGS construction:** building of PGSs in preparation for interaction testing.
3. **Interaction testing:** testing for SNP×PGS and SNP×SNP interactions, including fine-mapping.
4. **TF-PGS interactions:** build transcription factor-specific PGSs and test for interactions with
   them. Includes the second robustness check for TF interactions described in Supplementary
   Information §2.8.6.
5. **New UKB application:** the UKB application we originally used for this project expired before the
   project was complete and so a few final analyses were performed under a new application. The code
   in this section reruns some of the steps in the previous five sections using data from the new
   application, and then performs an additional robustness check for TF interactions (the second
   check described in Supplementary Information §2.8.6), runs the FAME and SME marginal epistasis
   methods, and builds and retests interactions using PGSs derived from FinnGen results.

The pipeline was run using Snakemake v7.22.0 on a Linux system. In addition to Snakemake, the other 
software packages required are:
- An R distribution with the packages listed in the file workflow/R-packages.txt
(https://www.r-project.org/);
- Plink 2 (https://www.cog-genomics.org/plink/2.0/);
- ANNOVAR (https://annovar.openbioinformatics.org/);
- BCFtools (https://github.com/samtools/bcftools) with the liftover plugin (https://github.com/freeseek/score);
- FAME (https://github.com/sriramlab/FAME);
- PRS-CS (https://github.com/getian107/PRScs).
