The .pvar files resulting from conversion with Plink 2 contain information on each variant's chromosome, position, rsID (Reference SNP Cluster ID) and reference/alternative alleles. For example, the first few lines of the .pvar for chromosome 1 are as follows:

#CHROM	POS	ID	REF	ALT
1	10177	rs367896724	A	AC
1	10235	rs540431307	T	TA
1	10352	rs201106462	T	TA
1	10505	rs548419688	A	T

There are two disadvantages to these files:
(i)  The rsIDs listed are not unique because multi-allelic SNPs with multiple alternative alleles are listed as separate variants with the same rsID; and
(ii) The files don't contain the imputation quality score ('INFO' score) provided by the UKB as a measure of the confidence in the imputation performed for each variant.

As a solution to these two issues, we created a second set of .pvar files which are contained in this folder (they are named ukb_imp_POSID_INFO_chr{chr}_v3.pvar).

To address issue (i), we replaced the rsIDs with new, custom variant IDs of the form CHR:POS:REF:ALT. Being unique, these new variant IDs allow us to create filtered lists of variant IDs that can be read directly by Plink 2 (using the --extract flag) without ambiguity (whereas with the default variant IDs it's not possible to specify which of multiple instances of the same rsID one is referring to).

To address issue (ii), we added an INFO field named IIS (for imputation information score) in VCF 4.3 format (see https://samtools.github.io/hts-specs/VCFv4.3.pdf), which allows for filtering by imputation quality score directly with Plink 2 (with, e.g., --extract-if-info "IIS >= 0.8").

These alternative .pvar files were created using the make-alternative-pvar.R script. To use these files instead of the default ones in Plink 2, use separate '--pgen', '--psam' and '--pvar' flags.
