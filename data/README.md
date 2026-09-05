This folder contains a skeleton folder structure indicating the external data required to run the
pipeline. Most of these folders are empty as the data either cannot be shared directly (UK Biobank)
or is publicly available for download.

Below we describe the different folders, the files they must contain and how this data can be
obtained. Remaining data files required by the different scripts in the pipeline are generated
through the pipeline itself as long as these initial files are in place.


### `1000G/`
  - `1000G/relate/allele_ages_EUR/allele_ages_GBR.RData`
  - `1000G/relate/mut/1000GP_Phase3_mask_prene_chr{chr}.mut.gz`

  Results of running Relate (Speidel et al., Nat. Genet. 2019) on 1000 Genomes Project Phase 3 data.
  The information of interest is which allele is ancestral and which allele is derived, as well as
  allele frequencies in a sample from the GBR population, at a large number of genomic positions.
  The script download-data.sh retrieves the necessary data.


### `annotations/`
  - `annotations/annovar/gene-anno/annovar-chr{chr}.ensGene.variant_function.gz`
  - `annotations/annovar/gene-anno/annovar-chr{chr}.refGene.variant_function.gz`
  - `annotations/vep/vep-var-consq-rel113_202410.csv`
  - `annotations/vep/anno/vep-chr{chr}.most_severe.tab.gz`

  Variant annotations from the Gencode, RefSeq Gene and Ensembl VEP databases which are downloaded
  as part of the pipeline. An additional manually created file (`vep-var-consq-rel113_202410.csv`) is
  also required which contains the ranking of VEP variant consequences; this information was
  obtained from the following web page:
  https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html.


### `covars/ancestry/`
  - `covars/ancestry/v2_487409.rds`
  - `covars/ancestry/wb_ids.rds`

  Ancestry Components from [Hu _et al._, Nat. Genet. (2025)](https://doi.org/10.1038/s41588-024-02035-8) which partition the ancestry of each of the
  487,409 samples in the imputed data into 127 components corresponding to geographically meaningful
  regions. This data has been returned to the UK Biobank by Hu *et al.* for use by other researchers.
  It consists of an R data file (`v2_487409.rds`) containing a data frame with 487,409 rows and 127
  columns (corresponding to samples and components, respectively) which is then concatenated to
  sample IDs in the correct order for our applications.

  Indicator of which samples are considered 'White British'. This follows the definition in [Bycroft
  *et al.*, Nature (2018)](https://doi.org/10.1038/s41586-018-0579-z) and consists of an R data file (`wb_ids.rds`) which we obtained directly from
  the authors of that paper. This file simply contains a logical vector of length 487,409 (the
  number of samples in the imputed data) indicating which samples meet this criterion and can be
  approximately replicated by following the process described by Bycroft *et al*.
  

### `epigenome/`
  - `epigenome/list-blueprint-h3k4me1.txt`
  - `epigenome/list-blueprint-h3k4me3.txt`
  - `epigenome/bed/{dataset}.bed`

  Data on H3K4me1 and H3K4me3 methylation marks from the Blueprint ChIP-seq Consortium obtained via
  the International Human Epigenome Consortium Data Portal. The script `download-data.sh` retrieves
  the necessary data.


### `imputed-genotypes/`
  - `imputed-genotypes/ukb_imp_chr{chr}_v3.pgen`
  - `imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar`
  - `imputed-genotypes/bgen/ukb_imp_chr{chr}_v3.bgen`

  The only genetic data directly used in this project is imputed autosomal data from the UK Biobank
  (data field 22828). This data must be available in its original format (BGEN) in the subfolder
  `bgen/` (files named `ukb_imp_chr1_v3.bgen`, etc.) as well as in Plink 2 format (PGEN) in the main
  folder (`ukb_imp_chr1_v3.pgen`, etc.).

  Conversion from BGEN to PGEN can be easily performed with Plink 2 (simply load the `.bgen` files and
  corresponding `.sample` file and use the `--out` flag). Instead of the `.pvar` files directly resulting
  from the conversion, we build and use an alternative set which has the key advantage of having
  unique variant identifiers. These must be located in the subfolder `alternative_pvar_files/` and
  named `ukb_imp_POSID_INFO_chr1_v3.pvar`, etc. See the README and Bash script present in this
  subfolder for more details.


### `ld-panels/ldblk_1kg_eur/`
  - `ld-panels/ldblk_1kg_eur/ldblk_1kg_chr{chr}.hdf5`
  - `ld-panels/ldblk_1kg_eur/snpinfo_1kg_hm3`

  The HDF5 files contain LD reference panels from the EUR population in the 1000 Genomes Project
  Phase 3 dataset provided with the PRS-CS software. The file `snpinfo_1kg_hm3` contains information
  on the SNPs included in the panel. The script `download-data.sh` retrieves and extracts the data.  


### `liftOver/`
  - `liftOver/hg19ToHg38.over.chain`
  - `liftOver/hg38ToHg19.over.chain`

  Reference files to convert ('lift') genomic coordinates from hg19 to hg38 and vice versa. These
  are for use with the rtracklayer R Bioconductor package and can be downloaded using the
  `download-data.sh` script.


### `phenotypes/`
  - `phenotypes/clean/code-desc-map-all.tab`
  - `phenotypes/clean/code-desc-map-real.tab`
  - `phenotypes/clean/code-desc-map-sim.tab`
  - `phenotypes/clean/snp-pgs-interaction-traits.tab`
  - `phenotypes/ukb-raw/ukb26618.tab`
  - `phenotypes/ukb-raw/ukb45188.tab`
  - `phenotypes/ukb-raw/ukb677925.tab`

  The first three files contain lists of all traits used in our analysis (both real and simulated).

  The fourth file lists the 52 traits for which we identified at least one SNP×PGS interaction hit
  as described in the paper. This is used for running the last part of the pipeline that tests for
  interactions with partitioned polygenic scores only for SNPs found to interact with the whole PGS.

  The files `ukb26618.tab`, `ukb45188.tab` and `ukb677925.tab` are tab-delimited files containing a
  variety of phenotypes which we obtained directly from the UK Biobank. Any other UK Biobank
  phenotype file in this format containing the necessary fields (phenotypes listed in
  `clean/code-desc-map-real.tab` and covariates described in the paper) would be sufficient to run the
  pipeline.


### `relatedness/`
  - `relatedness/ukb_rel_a27960_s488224.dat`
  - `relatedness/ukb_rel_a103076_s487957.dat`

  Lists of pairs of individuals related up to the third degree for the two applications used in this
  project. These files were obtained directly from the UK Biobank and are described in Resource 531
  (https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=531).


### `sample-ids/`
  - `sample-ids/ukb22828_c1_b0_v3_s487256.sample`
  - `sample-ids/ukb-27960-imp-auto-s487256-20210614.psam`
  - `sample-ids/ukb-103076-imp-auto-s486989.psam`
  - `sample-ids/withdrawals/w27960_20210201.csv`
  - `sample-ids/withdrawals/withdraw103076_179_20240301.txt`

  The first three files are sample files to be paired with genetic data in binary format which does
  not contain sample IDs, and are specific to our applications. The first two files correspond to
  the autosomal imputed genotypes for application number 27960; the first of these matches the BGEN
  files (the original format in which this data is provided) and the second the PGEN (Plink 2
  converted) files. The third file is equivalent to the second but for application number 103076.

  The last two files list sample IDs from individuals who withdrew from the UK Biobank study and
  whose data should therefore not be used (one for each application).


### `sum-stats/finngen/`
  - `sum-stats/finngen/finngen_R12_{OMOPID}.gz`
  - `sum-stats/finngen/Kanta_labs_GWAS_results_v2_summary.txt`
  - `sum-stats/finngen/ukb-finngen-lab-values-correspondence.csv`
  - `sum-stats/finngen/ukb-finngen-lab-values-match.tab`

  The files `finngen_R12_{OMOPID}.gz` (where the `OMOPID` strings correspond to the phenotype codes in FinnGen)
  contain summary statistics for the 28 matching phenotypes that we analyse when building PGSs
  derived from FinnGen GWAS results. The file `Kanta_labs_GWAS_results_v2_summary.txt` contains
  metadata for lab results phenotypes in FinnGen.

  The file `ukb-finngen-lab-values-correspondence.csv` was manually created and contains the
  correspondence between UKB phenotypes in the blood biochemistry, blood count and urine assays
  categories and equivalent ones in FinnGen.

  The script download-data.R generates the file `ukb-finngen-lab-values-match.tab` (which contains the
  28 phenotypes to be analysed) and downloads FinnGen files.


### `tf-binding/`
  - `tf-binding/hocomoco/v13/cluster_list.tsv`
  - `tf-binding/hocomoco/v13/H13CORE_jaspar_format.txt`

  Transcription factor binding motif information from the HOCOMOCO v13 database, including list of
  clusters into which motifs can be partitioned. The script `download-data.sh` retrieves the necessary
  data.
