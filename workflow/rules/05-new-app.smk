# 00 PRELIMINARY STEP
rule prepare_covars:
    input:
        '../data/phenotypes/ukb-raw/ukb-raw/ukb677925.tab',
        '../data/sample-ids/withdrawals/withdraw103076_179_20240301.txt',
        '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        '../data/covars/ancestry/v2_487409.rds'
    output:
        ('../data/covars/age-sex-batch-centre-ac-all.tab')
    threads: 2
    shell:
        'Rscript ./scripts/05-new-app/00a-prepare-covars.R'


rule sample_filtering:
    input:
        '../data/covars/age-sex-batch-centre-ac-all.tab',
        '../data/covars/ancestry/wb_ids.rds',
        '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        '../data/relatedness/ukb_rel_a103076_s487957.dat'
    output:
        (expand('../data/sample-ids/filtered/{samples}-ids.tab',
                    samples = ['all', 'female', 'male',
                               'wb_all', 'wb_female', 'wb_male',
                               'other_gb_vali_all', 'other_gb_vali_female', 'other_gb_vali_male',
                               'other_gb_test_all', 'other_gb_test_female', 'other_gb_test_male']))
    threads: 2
    shell:
        'Rscript ./scripts/05-new-app/00b-sample-filtering.R'


rule prepare_phenotypes:
    input:
        '../data/sample-ids/filtered/all-ids.tab',
        '../data/phenotypes/ukb-raw/ukb677925.tab',
        '../data/phenotypes/clean/code-desc-map-real.tab'
    output:
        (expand('../data/phenotypes/clean/{phen}.tab', phen=traits))
    threads: 3
    shell:
        'Rscript ./scripts/05-new-app/00c-prepare-phenotypes.R'





# 01 GWAS
rule regress_out_covars:
    input:
        phen_f   = '../data/phenotypes/clean/{phen}.tab',
        covars_f = '../data/covars/age-sex-batch-centre-ac-all.tab',
        samples  = '../data/sample-ids/filtered/{samples}-ids.tab'
    params:
        phen = '{phen}', samples = '{samples}'
    output:
        '../results/05-new-app/residuals-covars/{phen}/{phen}-{samples}-resid-covars.tab',
        '../results/05-new-app/residuals-covars/{phen}_qn/{phen}_qn-{samples}-resid-covars.tab'
    shell:
        """
        Rscript ./scripts/05-new-app/01-regress-out-covars.R \
        --phen {params.phen} \
        --samples {params.samples}
        """





# 02 BUILD PGS
rule compute_iter:
    input:
        pgen    = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar    = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam    = '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        samples = '../data/sample-ids/filtered/all-ids.tab',
        freqs   = expand('../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq', chr=chromosomes),
        coeff   = expand('../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/{train_sp}/coeff/final/{phen}-coeff-chr{chr}.tab',
                         chr=1, allow_missing = True),  # this is not available for all chromosomes for all samples
    params:
        phen = '{phen}', train_sp = '{train_sp}', r2 = '{r2}', kb = '{kb}'
    output:
        pgs_bychr = ('../results/05-new-app/iterative/{phen}/r2_{r2}-kb_{kb}/{train_sp}/pgs/{phen}-all-pgs-bychr.tab'),
        pgs_loco  = ('../results/05-new-app/iterative/{phen}/r2_{r2}-kb_{kb}/{train_sp}/pgs/{phen}-all-pgs-loco.tab'),
        pgs_full  = ('../results/05-new-app/iterative/{phen}/r2_{r2}-kb_{kb}/{train_sp}/pgs/{phen}-all-pgs-full.tab')
    threads: 2
    shell:
        """
        Rscript ./scripts/05-new-app/02a-compute-iter-pgs.R \
        --phen {params.phen}
        """


rule mean_correction:
    input:
        phen_f   = '../results/05-new-app/residuals-covars/{phen}/{phen}-{samples}-resid-covars.tab',
        pgs_full = '../results/05-new-app/iterative/{phen}/r2_{r2}-kb_{kb}/wb_all/pgs/{phen}-all-pgs-full.tab',
        pgs_loco = '../results/05-new-app/iterative/{phen}/r2_{r2}-kb_{kb}/wb_all/pgs/{phen}-all-pgs-loco.tab'
    params:
        phen = '{phen}', samples = '{samples}', r2 = '{r2}', kb = '{kb}'
    output:
        pgs_full_mc = ('../results/05-new-app/mean-corrected/{phen}/r2_{r2}-kb_{kb}/{phen}-{samples}-pgs-full_mc.tab'),
        pgs_loco_mc = (expand('../results/05-new-app/mean-corrected/{phen}/r2_{r2}-kb_{kb}/{phen}-{samples}-pgs-loco{chr}_mc.tab',
                              chr = chromosomes, allow_missing = True))
    shell:
        """
        Rscript ./scripts/05-new-app/02b-mean-correction.R \
        --phen {params.phen} \
        --samples {params.samples}
        """


rule reg_pgs_ac:
    input:
        phen_res_f  = '../results/05-new-app/residuals-covars/{phen}/{phen}-wb_all-resid-covars.tab',
        covars_f = '../data/covars/age-sex-batch-centre-ac-all.tab',
        pgs_full_mc = '../results/05-new-app/mean-corrected/{phen}/r2_{r2}-kb_{kb}/{phen}-wb_all-pgs-full_mc.tab',
        pgs_loco_mc = expand('../results/05-new-app/mean-corrected/{phen}/r2_{r2}-kb_{kb}/{phen}-wb_all-pgs-loco{chr}_mc.tab',
                             chr = chromosomes, allow_missing = True)
    params:
         phen = '{phen}', train_sp = 'wb_all', r2 = '{r2}', kb = '{kb}'
    output:
        ('../results/05-new-app/ancestry-interactions/{phen}/r2_{r2}-kb_{kb}/{phen}-{train_sp}-resid-covars-pgs-full_mc-ac-centre.tab'),
        (expand('../results/05-new-app/ancestry-interactions/{phen}/r2_{r2}-kb_{kb}/{phen}-{train_sp}-resid-covars-pgs-loco{chr}_mc-ac-centre.tab',
                chr = chromosomes, allow_missing = True))
    threads: 2
    shell:
        """
        Rscript ./scripts/05-new-app/02c-pgs-ancestry-interactions.R \
        --phen {params.phen}
        """





# 04 TF PGS
rule comp_pgs_anno:
    input:
        pgen     = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar     = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
                          chr=chromosomes),
        psam     = '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        samples  = '../data/sample-ids/filtered/wb_all-ids.tab',
        orig_gwas = expand('../results/01-gwas/plink-output/{phen}/wb_all.chr{chr}.{phen}.res.cov.glm.linear.gz',
                           chr=chromosomes, allow_missing = True),
        coding_snps = '../results/04-tf-binding/coding-epigenetics/{phen}/coding-pgs-snps.RData',
        h3k4me1_snps = '../results/04-tf-binding/coding-epigenetics/{phen}/h3k4me1-pgs-snps.RData',
        h3k4me3_snps = '../results/04-tf-binding/coding-epigenetics/{phen}/h3k4me3-pgs-snps.RData'
    params:
        phen = '{phen}'
    output:
        (expand('../results/05-new-app/coding-epigenetics/{phen}/{annot}-pgs.RData',
                annot=['coding', 'h3k4me1', 'h3k4me3'], allow_missing = True))
    threads: 2
    shell:
        """
         Rscript ./scripts/05-new-app/03a-compute-pgs-coding-enhancer-promoter.R \
        --phen {params.phen}
        """


rule comp_pgs_tf:
    input:
        pgen     = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar     = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
                          chr=chromosomes),
        psam     = '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        samples  = '../data/sample-ids/filtered/wb_all-ids.tab',
        orig_gwas = expand('../results/01-gwas/plink-output/{phen}/wb_all.chr{chr}.{phen}.res.cov.glm.linear.gz',
                           chr=chromosomes, allow_missing = True),
        hocomoco = '../data/tf-binding/hocomoco/v13/H13CORE_jaspar_format.txt',
        tf_snps = '../results/04-tf-binding/tf-binding/{phen}/tf-pgs-snps.RData'
    params:
        phen = '{phen}', tf_batch = '{tf_batch}'
    output:
        ('../results/05-new-app/tf-binding/{phen}/tf-pgs-{tf_batch}.RData')
    threads: 2
    shell:
        """
         Rscript ./scripts/05-new-app/03b-compute-pgs-tf-binding.R \
        --phen {params.phen} \
        --tf_batch {params.tf_batch}
        """


rule comp_pgs_init:
    input:
        pgen     = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar     = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam     = '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        samples  = '../data/sample-ids/filtered/wb_all-ids.tab'
    params:
        phen = '{phen}'
    output:
        ('../results/05-new-app/initial/{phen}/pgs-loco.RData')
    threads: 4
    shell:
        """
         Rscript ./scripts/05-new-app/03c-compute-cpt-pgs.R \
        --phen {params.phen}
        """


### set of QN traits for which we see SNP*TF-PGS hits before filtering out interactions due to pairwise hits
traits_hits_tf = ['f.30100_qn', 'f.30110_qn', 'f.30150_qn', 'f.30610_qn',
                  'f.30620_qn', 'f.30650_qn', 'f.30690_qn', 'f.30720_qn',
                  'f.30750_qn', 'f.30840_qn', 'f.30870_qn', 'f.30890_qn']

rule test_pairwise:
    input:
        pgen    = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar    = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam    = '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        samples  = '../data/sample-ids/filtered/wb_all-ids.tab',
        hocomoco = '../data/tf-binding/hocomoco/v13/H13CORE_jaspar_format.txt',
        indep_hits = '../results/04-tf-binding/gwas/aggregate/int-gwas-indep-hits-aggregate.RData',
        tf_pgs_snps = expand('../results/04-tf-binding/tf-binding/{phen}/tf-pgs-snps.RData',
                             phen=traits_hits_tf),
        orig_gwas = expand('../results/01-gwas/plink-output/{phen}/wb_all.chr{chr}.{phen}.res.cov.glm.linear.gz',
                           phen=traits_hits_tf, chr=chromosomes),
        phen_res = expand('../results/05-new-app/residuals-covars/{phen}/{phen}-{samples}-resid-covars.tab',
                          phen=traits_hits_tf, samples='wb_all'),
        pgs_cpt = expand('../results/05-new-app/initial/{phen}/pgs-loco.RData', phen=traits_hits_tf),
    output:
        '../results/05-new-app/tf-gwas/aggregate/int-gwas-indep-hits-aggregate-retest-pairwise.RData'
    threads: 8
    shell:
        'Rscript ./scripts/05-new-app/04-test-hits-remove-pairwise.R'





# 05 RUN FAME & SME
rule fame_prep:
    input:
        bed = expand('../data/genotype-calls/ukb_cal_chr{chr}_v2.bed', chr=chromosomes),
        bim = expand('../data/genotype-calls/ukb_snp_chr{chr}_v2.bim', chr=chromosomes),
        fam = '../data/sample-ids/ukb22418_c1_b0_v2_s487957.fam',
        samples = '../data/sample-ids/filtered/wb_all-ids.tab',
        sim_coeff_null = expand('../data/phenotypes/simulations-prep/coeff/sim-coeff-1k_rg-chr{chr}.tab', chr=chromosomes),
        sim_coeff_int = '../data/phenotypes/simulations-prep/coeff/sim-coeff-inter.tab',
        lhs_hits_loco_qn = '../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn-annot.tab',
        pgen  = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar  = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam  = '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        traits_sim = '../data/phenotypes/clean/code-desc-map-sim.tab'
    output:
        temp(expand('../results/05-new-app/fame/genotype-calls/ukb_array_filtered_chr{chr}.{ext}',
                    chr = chromosomes, ext = ['bed', 'bim', 'fam', 'log'])),
        temp(expand('../results/05-new-app/fame/genotype-calls/ukb_imputed_to_add_chr{chr}.{ext}',
                    chr = chromosomes, ext = ['bed', 'bim', 'fam', 'log'])),
        temp(expand('../results/05-new-app/fame/genotype-calls/ukb_array_imputed_concat_chr{chr}.{ext}',
                    chr = chromosomes, ext = ['bed', 'bim', 'fam', 'log'])),
        '../results/05-new-app/fame/genotype-calls/ukb_array_imputed_concat_all_files.tab',
        expand('../results/05-new-app/fame/genotype-calls/ukb_all_chr.{ext}',
               ext = ['bed', 'bim', 'fam', 'log'])
    threads: 4
    shell:
        'Rscript ./scripts/05-new-app/05a-fame-prepare-bed-files.R'


rule fame_run:
    input:
        bim = '../results/05-new-app/fame/genotype-calls/ukb_all_chr.bim',
        fam = '../results/05-new-app/fame/genotype-calls/ukb_all_chr.fam',
        samples = '../data/sample-ids/filtered/wb_all-ids.tab',
        phen_raw = lambda wildcards:  f'../data/phenotypes/clean/{wildcards.phen}'[:-3] + '.tab',
        ukb_phen = '../data/phenotypes/ukb-raw/ukb677925.tab'
    output:
        temp('../results/05-new-app/fame/output/{phen}/{phen}-{snp}.pheno'),
        temp('../results/05-new-app/fame/output/{phen}/{phen}-{snp}.covar'),
        '../results/05-new-app/fame/output/{phen}/annot-{snp}.tab',
        '../results/05-new-app/fame/output/{phen}/results-{snp}.tab'
    params:
        phen = '{phen}', snp = '{snp}'
    shell:
        """
        Rscript ./scripts/05-new-app/05b-fame-run.R \
        --phen {params.phen} \
        --snp {params.snp}
        """


rule sme_filter_imp:
    input:
        pgen  = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar  = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam  = '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        samples = '../data/sample-ids/filtered/wb_all-ids.tab'
    output:
        (expand('../results/05-new-app/sme/genotype-calls/ukb_imputed_filtered_chr{chr}.{ext}',
                chr=chromosomes, ext=['bed', 'bim', 'fam', 'log']))
    threads: 4
    shell:
        'Rscript ./scripts/05-new-app/05c-sme-filter-imputed-genotypes.R'


rule sme_prep_bed:
    input:
        pgen  = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar  = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam  = '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        samples = '../data/sample-ids/filtered/wb_all-ids.tab',
        sim_coeff_null = expand('../data/phenotypes/simulations-prep/coeff/sim-coeff-1k_rg-chr{chr}.tab', chr=chromosomes),
        sim_coeff_int = '../data/phenotypes/simulations-prep/coeff/sim-coeff-inter.tab',       
        coeff   = expand('../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/{train_sp}/coeff/final/{phen}-coeff-chr{chr}.tab',
                         r2=0.9, kb=500, train_sp='wb_all', chr=1, allow_missing = True),  # this is not available for all chromosomes for all samples      
        bed_imp_flt = expand('../results/05-new-app/sme/genotype-calls/ukb_imputed_filtered_chr{chr}.{ext}',
                             chr=chromosomes, ext=['bed', 'bim', 'fam'])
    output:
        temp(expand('../results/05-new-app/sme/genotype-calls/{phen}-{snp}-ukb_all_chr.{ext}',
                    ext=['bed', 'bim', 'fam', 'log'], allow_missing = True))
    params:
        phen = '{phen}', snp = '{snp}'
    threads: 4
    priority: 0
    shell:
        """
        Rscript ./scripts/05-new-app/05d-sme-prepare-bed-files.R \
        --phen {params.phen} \
        --snp {params.snp}
        """


rule sme_prep_phen:
    input:
        bed_files = expand('../results/05-new-app/sme/genotype-calls/{phen}-{snp}-ukb_all_chr.{ext}',
                           ext=['bed', 'bim', 'fam'], allow_missing = True),
        samples = '../data/sample-ids/filtered/wb_all-ids.tab',
        phen_raw = lambda wildcards:  f'../data/phenotypes/clean/{wildcards.phen}'[:-3] + '.tab',
        ukb_phen = '../data/phenotypes/ukb-raw/ukb677925.tab'
    output:
        temp('../results/05-new-app/sme/output/{phen}/{phen}-{snp}.pheno'),
        temp('../results/05-new-app/sme/output/{phen}/{phen}-{snp}.samples'),
        temp(expand('../results/05-new-app/sme/genotype-calls/{phen}-{snp}-ukb_all_flt_samples_chr.{ext}',
                    ext=['bed', 'bim', 'fam', 'log'], allow_missing = True))
    params:
        phen = '{phen}', snp = '{snp}'
    threads: 4
    priority: 50
    shell:
        """
        Rscript ./scripts/05-new-app/05e-sme-prepare-phen-covars.R \
        --phen {params.phen} \
        --snp {params.snp}
        """


rule sme_run:
    input:
        bed_files = expand('../results/05-new-app/sme/genotype-calls/{phen}-{snp}-ukb_all_flt_samples_chr.{ext}',
                           ext=['bed', 'bim', 'fam'], allow_missing = True),
        coeff   = expand('../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/{train_sp}/coeff/final/{phen}-coeff-chr{chr}.tab',
                         r2=0.9, kb=500, train_sp='wb_all', chr=1, allow_missing = True),  # this is not available for all chromosomes for all samples
        pheno = '../results/05-new-app/sme/output/{phen}/{phen}-{snp}.pheno'
    output:
        temp('../results/05-new-app/sme/output/{phen}/{phen}-{snp}.mask'),
        '../results/05-new-app/sme/output/{phen}/{phen}-{snp}-sme-result.RData'
    params:
        phen = '{phen}', snp = '{snp}'
    threads: 4
    priority: 100
    shell:
        """
        Rscript ./scripts/05-new-app/05f-sme-run.R \
        --phen {params.phen} \
        --snp {params.snp}
        """





## FINNGEN
rule fg_convert_hg:
    input:
        ukb_fg_key = '../data/sum-stats/finngen/ukb-finngen-lab-values-match.tab',
        fg_gwas = expand('../data/sum-stats/finngen/finngen_R12_{omopid}.gz', omopid=traits_fg_omopid),
        pvar = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        snpinfo_1kg = '../../../software/PRScs/ldblk_1kg_eur/snpinfo_1kg_hm3'
    output:
        '../data/sum-stats/finngen/variant-ids-hg38.tab.gz',
        '../data/sum-stats/finngen/variant-ids-hg38.vcf.gz',
        '../data/sum-stats/finngen/variant-ids-hg19.vcf.gz',
        '../data/sum-stats/finngen/variant-ids-hg38-hg19-ukb.tab.gz',
        '../results/05-new-app/finngen/bim/fg-ukb-1kg.bim'
    threads: 4
    shell:
        'Rscript ./scripts/05-new-app/06a-finngen-convert-hg38-to-hg19.R'


rule fg_run_prscs:
    input:
        lab_docs = '../data/sum-stats/finngen/Kanta_labs_GWAS_results_v2_summary.txt',
        fg_gwas = '../data/sum-stats/finngen/finngen_R12_{omopid}.gz',
        hg38_hg19_ukb = '../data/sum-stats/finngen/variant-ids-hg38-hg19-ukb.tab.gz',
        snpinfo_1kg = '../../../software/PRScs/ldblk_1kg_eur/snpinfo_1kg_hm3',
        samples = '../data/sample-ids/filtered/all-ids.tab',
        pgen    = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar    = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam    = '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        freqs   = expand('../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq', chr=chromosomes)
    output:
        temp('../results/05-new-app/finngen/sumstats/{phen}-{omopid}-sumstats.tab'),
        expand('../results/05-new-app/finngen/PRScs/{phen}/{phen}-{omopid}_pst_eff_a1_b0.5_phiauto_chr{chr}.txt',
               chr=chromosomes, allow_missing = True),
        '../results/05-new-app/finngen/PRScs/{phen}/{phen}-{omopid}-all-pgs-full.tab',
        '../results/05-new-app/finngen/PRScs/{phen}/{phen}-{omopid}-all-pgs-loco.tab'
    params:
        phen = '{phen}', omopid = '{omopid}'
    threads: 4
    shell:
        """
        module load SciPy-bundle/2023.07-gfbf-2023a
        module load h5py/3.9.0-foss-2023a

        Rscript ./scripts/05-new-app/06b-finngen-run-prscs.R \
        --phen {params.phen} \
        --omopid {params.omopid}
        """


rule fg_retest:
    input:
        samples = '../data/sample-ids/filtered/wb_all-ids.tab',
        pgen    = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar    = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam    = '../data/sample-ids/ukb-103076-imp-auto-s486989.psam',
        lhs_hits_loco_qn = '../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn.tab',
        ukb_fg_key = '../data/sum-stats/finngen/ukb-finngen-lab-values-match.tab',
        phen_train = expand('../results/05-new-app/residuals-covars/{phen}/{phen}-wb_all-resid-covars.tab',
                            phen=traits_fg_field_id),
        pgs_full_mc_train = expand('../results/05-new-app/mean-corrected/{phen}/r2_0.9-kb_500/{phen}-wb_all-pgs-full_mc.tab',
                                   phen=traits_fg_field_id),
        phen_test = expand('../results/05-new-app/residuals-covars/{phen}/{phen}-other_gb_test_all-resid-covars.tab',
                            phen=traits_fg_field_id),
        pgs_full_mc_test =  expand('../results/05-new-app/mean-corrected/{phen}/r2_0.9-kb_500/{phen}-other_gb_test_all-pgs-full_mc.tab',
                                   phen=traits_fg_field_id),
        pgs_prscs_full = expand('../results/05-new-app/finngen/PRScs/{phen}/{phen}-{omopid}-all-pgs-full.tab',
                                zip, phen=traits_fg_field_id, omopid=traits_fg_omopid),
        pgs_prscs_loco = expand('../results/05-new-app/finngen/PRScs/{phen}/{phen}-{omopid}-all-pgs-loco.tab',
                                zip, phen=traits_fg_field_id, omopid=traits_fg_omopid),
    output:
        '../results/05-new-app/finngen/finngen-retest-interactions.rds'
    threads: 2
    shell:
        'Rscript ./scripts/05-new-app/06c-finngen-retest.R'
 
