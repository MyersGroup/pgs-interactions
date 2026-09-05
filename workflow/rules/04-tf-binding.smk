# 04 TF BINDING
rule prep_aux_data:
    input:
        mut = expand('../data/1000G/relate/mut/1000GP_Phase3_mask_prene_chr{chr}.mut.gz', chr=chromosomes)
    output:
        '../data/1000G/relate/mut/all-snps-relate.RData'
    threads: 4
    shell:
        'Rscript ./scripts/04-tf-binding/00-prepare-auxilliary-data.R'


rule find_tags:
    input:
        ld_clump = expand('../results/02-pgs/ld-clump/{phen}/p1_{p1}-r2_{r2}-kb_{kb}/wb_all-chr{chr}-clumped.tab',
                          chr=chromosomes, p1=p1_max, r2=0.1, kb=500, allow_missing = True),
        pgen     = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar     = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
                          chr=chromosomes),
        psam     = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        freqs    = expand('../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq', chr=chromosomes),
        variants = expand('../data/variant-ids/chr{chr}-var-ids.tab', chr=chromosomes),
        map_posid_rsid = expand('../data/variant-ids/map-posid-rsid/posid-rsid-chr{chr}.tab', chr=chromosomes),
        liftOver_chain = '../data/liftOver/hg19ToHg38.over.chain',
        relate = '../data/1000G/relate/mut/all-snps-relate.RData',
        orig_gwas = expand('../results/01-gwas/plink-output/{phen}/wb_all.chr{chr}.{phen}.res.cov.glm.linear.gz',
                           chr=chromosomes, allow_missing = True)
    params:
        phen = '{phen}', samples = 'wb_all'
    output:
        '../results/04-tf-binding/pgs-snps-tags/{phen}/pgs-snps-tags.RData'
    threads: 4
    shell:
        """
        module unload R/4.3.2-gfbf-2023a
        module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

         Rscript ./scripts/04-tf-binding/01-find-tags-of-pgs-snps.R \
        --phen {params.phen} \
        --samples {params.samples}
        """


rule mk_func_anno:
    input:
        pgs_snps_tags = '../results/04-tf-binding/pgs-snps-tags/{phen}/pgs-snps-tags.RData',
        h3k4me1_peaks_ls = '../data/epigenome/list-blueprint-h3k4me1.txt',
        h3k4me3_peaks_ls = '../data/epigenome/list-blueprint-h3k4me3.txt',
        vep = expand('../data/annotations/vep/anno/vep-chr{chr}.most_severe.tab.gz', chr=chromosomes)
    params:
        phen = '{phen}', samples = 'wb_all'
    output:
        '../results/04-tf-binding/coding-epigenetics/{phen}/annot/h3k4me1_3-annot.RData',
        '../results/04-tf-binding/coding-epigenetics/{phen}/annot/functional-annot.RData',
        '../results/04-tf-binding/coding-epigenetics/{phen}/coding-pgs-snps.RData',
        '../results/04-tf-binding/coding-epigenetics/{phen}/h3k4me1-pgs-snps.RData',
        '../results/04-tf-binding/coding-epigenetics/{phen}/h3k4me3-pgs-snps.RData'
    threads: 4
    shell:
        """
         Rscript ./scripts/04-tf-binding/02-make-coding-enhancer-promoter-anno.R \
        --phen {params.phen} \
        --samples {params.samples}
        """
 

rule mk_tf_anno:
    input:
        pgs_snps_tags = '../results/04-tf-binding/pgs-snps-tags/{phen}/pgs-snps-tags.RData',
        h3k4me1_3_annot = '../results/04-tf-binding/coding-epigenetics/{phen}/annot/h3k4me1_3-annot.RData',
        coding_annot = '../results/04-tf-binding/coding-epigenetics/{phen}/annot/functional-annot.RData',
        hocomoco = '../data/tf-binding/hocomoco/v13/H13CORE_jaspar_format.txt',
        relate = '../data/1000G/relate/allele_ages_EUR/allele_ages_GBR.RData',
        h3k4me1_peaks_ls = '../data/epigenome/list-blueprint-h3k4me1.txt',
        h3k4me3_peaks_ls = '../data/epigenome/list-blueprint-h3k4me3.txt',
        liftOver_chain = '../data/liftOver/hg38ToHg19.over.chain',
    params:
        phen = '{phen}', samples = 'wb_all'
    output:
        '../results/04-tf-binding/tf-binding/{phen}/neithermaxscorequantilematq150.RData',
        '../results/04-tf-binding/tf-binding/{phen}/motif-scores.RData',
        '../results/04-tf-binding/tf-binding/{phen}/tf-pgs-snps.RData'
    threads: 4
    shell:
        """
        module unload R/4.3.2-gfbf-2023a
        module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

         Rscript ./scripts/04-tf-binding/03-make-tf-binding-anno.R \
        --phen {params.phen} \
        --samples {params.samples}
        """


rule comp_pgs_anno:
    input:
        pgen     = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar     = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
                          chr=chromosomes),
        psam     = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        samples  = '../data/sample-ids/filtered/wb_all-ids.tab',
        orig_gwas = expand('../results/01-gwas/plink-output/{phen}/wb_all.chr{chr}.{phen}.res.cov.glm.linear.gz',
                           chr=chromosomes, allow_missing = True),
        coding_snps = '../results/04-tf-binding/coding-epigenetics/{phen}/coding-pgs-snps.RData',
        h3k4me1_snps = '../results/04-tf-binding/coding-epigenetics/{phen}/h3k4me1-pgs-snps.RData',
        h3k4me3_snps = '../results/04-tf-binding/coding-epigenetics/{phen}/h3k4me3-pgs-snps.RData'
    params:
        phen = '{phen}', samples = 'wb_all'
    output:
        temp(expand('../results/04-tf-binding/coding-epigenetics/{phen}/{annot}-pgs.RData',
                    annot=['coding', 'h3k4me1', 'h3k4me3'], allow_missing = True))
    threads: 2
    shell:
        """
         Rscript ./scripts/04-tf-binding/04a-compute-pgs-coding-enhancer-promoter.R \
        --phen {params.phen} \
        --samples {params.samples}
        """


rule comp_pgs_tf:
    input:
        pgen     = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar     = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
                          chr=chromosomes),
        psam     = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        samples  = '../data/sample-ids/filtered/wb_all-ids.tab',
        orig_gwas = expand('../results/01-gwas/plink-output/{phen}/wb_all.chr{chr}.{phen}.res.cov.glm.linear.gz',
                           chr=chromosomes, allow_missing = True),
        hocomoco = '../data/tf-binding/hocomoco/v13/H13CORE_jaspar_format.txt',
        tf_snps = '../results/04-tf-binding/tf-binding/{phen}/tf-pgs-snps.RData'
    params:
        phen = '{phen}', samples = 'wb_all', tf_batch = '{tf_batch}'
    output:
        temp('../results/04-tf-binding/tf-binding/{phen}/tf-pgs-{tf_batch}.RData')
    threads: 2
    shell:
        """
         Rscript ./scripts/04-tf-binding/04b-compute-pgs-tf-binding.R \
        --phen {params.phen} \
        --samples {params.samples} \
        --tf_batch {params.tf_batch}
        """


rule run_gwas:
    input:
        pgen     = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar     = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
                          chr=chromosomes),
        psam     = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        samples  = '../data/sample-ids/filtered/wb_all-ids.tab',
        phen_res = '../results/01-gwas/residuals-covars/{phen}/{phen}-wb_all-resid-covars.tab',
        indep_hits = '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-wb_all-loco-indep-hits.tab',
        hocomoco = '../data/tf-binding/hocomoco/v13/H13CORE_jaspar_format.txt',
        anno_pgs = expand('../results/04-tf-binding/coding-epigenetics/{phen}/{annot}-pgs.RData',
                          annot=['coding', 'h3k4me1', 'h3k4me3'], allow_missing = True),
        tf_pgs   = expand('../results/04-tf-binding/tf-binding/{phen}/tf-pgs-{tf_batch}.RData',
                          tf_batch=tf_ind, allow_missing = True)
    params:
        phen = '{phen}', samples = 'wb_all'
    output:
        temp('../results/04-tf-binding/gwas/{phen}/pgs-loco.RData'),
        '../results/04-tf-binding/gwas/{phen}/int-gwas-sumstats.RData'
    threads: 4
    shell:
        """
         Rscript ./scripts/04-tf-binding/05-run-gwas.R \
        --phen {params.phen} \
        --samples {params.samples}
        """


rule stepwise_tf:
    input:
        pgen     = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar     = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
                          chr=chromosomes),
        psam     = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        samples  = '../data/sample-ids/filtered/wb_all-ids.tab',
        phen_res = '../results/01-gwas/residuals-covars/{phen}/{phen}-wb_all-resid-covars.tab',
        pgs_loco = '../results/04-tf-binding/gwas/{phen}/pgs-loco.RData',
        indep_hits = '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-wb_all-loco-indep-hits.tab',
        hocomoco = '../data/tf-binding/hocomoco/v13/H13CORE_jaspar_format.txt',
        tf_snps = '../results/04-tf-binding/tf-binding/{phen}/tf-pgs-snps.RData',
        gwas = '../results/04-tf-binding/gwas/{phen}/int-gwas-sumstats.RData',
        anno_pgs = expand('../results/04-tf-binding/coding-epigenetics/{phen}/{annot}-pgs.RData',
                          annot=['coding', 'h3k4me1', 'h3k4me3'], allow_missing = True),
        tf_pgs   = expand('../results/04-tf-binding/tf-binding/{phen}/tf-pgs-{tf_batch}.RData',
                          tf_batch=tf_ind, allow_missing = True),
        coding_snps = '../results/04-tf-binding/coding-epigenetics/{phen}/coding-pgs-snps.RData',
        h3k4me1_snps = '../results/04-tf-binding/coding-epigenetics/{phen}/h3k4me1-pgs-snps.RData',
        h3k4me3_snps = '../results/04-tf-binding/coding-epigenetics/{phen}/h3k4me3-pgs-snps.RData'
    params:
        phen = '{phen}', samples = 'wb_all'
    output:
        '../results/04-tf-binding/gwas/{phen}/int-gwas-indep-hits.RData'
    threads: 4
    shell:
        """
         Rscript ./scripts/04-tf-binding/06-stepwise-reg.R \
        --phen {params.phen} \
        --samples {params.samples}
        """


rule agg_tf:
    input:
        indep_hits_loco = '../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn.tab',
        indep_hits_tf   = expand('../results/04-tf-binding/gwas/{phen}/int-gwas-indep-hits.RData',
                                 phen=traits_hits_loco_qn),
        map_posid_rsid = expand('../data/variant-ids/map-posid-rsid/posid-rsid-chr{chr}.tab',
                                chr=chromosomes),
        maf = expand('../data/imputed-wb_all-stats/minor-alleles/chr{chr}.tab',
                     chr=chromosomes),
        vep = expand('../data/annotations/vep/anno/vep-chr{chr}.most_severe.tab.gz',
                     chr=chromosomes),
        annovar_ensGene = expand('../data/annotations/annovar/gene-anno/annovar-chr{chr}.ensGene.variant_function.gz',
                                 chr=chromosomes),
        names_clean = '../data/phenotypes/clean/code-desc-map-real-qn-thesis.tab',
        cluster_list = '../data/tf-binding/hocomoco/v13/cluster_list.tsv'
    output:
        '../results/04-tf-binding/gwas/aggregate/int-gwas-indep-hits-aggregate.RData'
    threads: 2
    shell:
        "Rscript ./scripts/04-tf-binding/07-aggregate-results.R"


rule test_hits_piece:
    input:
        samples  = '../data/sample-ids/filtered/wb_all-ids.tab',
        hocomoco = '../data/tf-binding/hocomoco/v13/H13CORE_jaspar_format.txt',
        indep_hits = '../results/04-tf-binding/gwas/aggregate/int-gwas-indep-hits-aggregate.RData',
        phen_res = expand('../results/01-gwas/residuals-covars/{phen}/{phen}-wb_all-resid-covars.tab',
                          phen=traits_hits_loco_qn),
        pgs_loco = expand('../results/04-tf-binding/gwas/{phen}/pgs-loco.RData',
                          phen=traits_hits_loco_qn),
        pgs_coeff = expand('../results/02-pgs/initial/{phen}/r2_0.1-kb_500/wb_all/coeff/p1_opt/chr{chr}-coeff.tab',
                           phen=traits_hits_loco_qn, chr=chromosomes),
        p1_opt   = expand('../results/02-pgs/initial/{phen}/r2_0.1-kb_500/wb_all/pgs-all/p1_opt/{phen}-pgs0-rsq.tab',
                          phen=traits_hits_loco_qn),
        tf_pgs_snps = expand('../results/04-tf-binding/tf-binding/{phen}/tf-pgs-snps.RData',
                             phen=traits_hits_loco_qn),
        orig_gwas = expand('../results/01-gwas/plink-output/{phen}/wb_all.chr{chr}.{phen}.res.cov.glm.linear.gz',
                           phen=traits_hits_loco_qn, chr=chromosomes),
        cluster_list = '../data/tf-binding/hocomoco/v13/cluster_list.tsv'
    output:
        '../results/04-tf-binding/gwas/aggregate/int-gwas-indep-hits-aggregate-retest.RData'
    threads: 4
    shell:
        "Rscript ./scripts/04-tf-binding/08-test-hits-piece-pgs.R"
