## 02 BUILD PGS
rule convert_pgen_bigsnp:
    input:
        variants    = '../data/variant-ids/chr{chr}-var-ids.tab',
        freqs       = '../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq',
        sample_file = '../data/sample-ids/ukb22828_c1_b0_v3_s487256.sample',
        bgen        = '../data/imputed-genotypes/bgen/ukb_imp_chr{chr}_v3.bgen'
    output:
        temp('../data/imputed-genotypes/bigsnp/chr{chr}.bk'),
        temp('../data/imputed-genotypes/bigsnp/chr{chr}.rds')
    threads: 4 
    shell:
        """
        Rscript ./scripts/02-pgs/00-convert-bgen-bigsnp.R \
        --chr {wildcards.chr} \
        --threads {threads}
        """


rule ld_clump:
    input:
        bigsnp_bk  = '../data/imputed-genotypes/bigsnp/chr{chr}.bk',
        bigsnp_rds = '../data/imputed-genotypes/bigsnp/chr{chr}.rds',
        gwas       = '../results/01-gwas/plink-output/{phen}/{samples}.chr{chr}.{phen}.res.cov.glm.linear.gz'
    params:
        phen = '{phen}', chr = '{chr}', samples = '{samples}',
        p1 = '{p1}', r2 = '{r2}', kb = '{kb}'
    output:
        '../results/02-pgs/ld-clump/{phen}/p1_{p1}-r2_{r2}-kb_{kb}/{samples}-chr{chr}-clumped.tab'
    threads: 4
    shell:
        """
        Rscript ./scripts/02-pgs/01-ld-clump.R \
        --phen {params.phen} \
        --chr {params.chr} \
        --train_sp {params.samples} \
        --p1 {params.p1} \
        --r2 {params.r2} \
        --kb {params.kb} \
        --threads {threads}
        """


## build initial PGS
rule make_coeff_vec_p1:
    input:
        expand('../results/02-pgs/ld-clump/{phen}/p1_{p1}-r2_{r2}-kb_{kb}/{samples}-chr{chr}-clumped.tab',
               chr=chromosomes, p1=p1_max, allow_missing = True)
    params:
        phen = '{phen}', samples='{samples}',
        p1_max = p1_max, r2 = '{r2}', kb = '{kb}'
    output:
        expand('../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{samples}/coeff/p1_grid/pgs-p1-n_snps.tab',
               allow_missing = True),
        temp(expand('../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{samples}/coeff/p1_grid/chr{chr}-coeff.tab',
                    chr=chromosomes, allow_missing = True))
    shell:
        """
        Rscript ./scripts/02-pgs/02a-make-coeff-vec-p1.R --phen {params.phen} --samples {params.samples} \
        --p1_max {params.p1_max} --r2 {params.r2} --kb {params.kb}
        """


rule compute_scores_p1_vali:
    input:
        pgen    = '../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen',
        pvar    = '../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
        psam    = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        vali_sp = '../data/sample-ids/filtered/other_gb_vali_all-ids.tab',
        freqs   = '../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq',
        coeff   = '../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{samples}/coeff/p1_grid/chr{chr}-coeff.tab'
    params:
        out_prefix = '../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{samples}/pgs-vali/p1_grid/chr{chr}'
    output:
        temp('../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{samples}/pgs-vali/p1_grid/chr{chr}.sscore'),
        temp('../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{samples}/pgs-vali/p1_grid/chr{chr}.log')
    threads: 2
    shell:
        """
        ./scripts/02-pgs/02-compute-score.sh {input.pgen} {input.pvar} {input.psam} \
        {input.vali_sp} {input.freqs} {input.coeff} {params.out_prefix} 
        """


rule choose_p1:
    input:
        expand('../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{train_sp}/coeff/p1_grid/chr{chr}-coeff.tab',
               train_sp='wb_all', chr=chromosomes, allow_missing = True),
        expand('../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{train_sp}/pgs-vali/p1_grid/chr{chr}.sscore',
               train_sp='wb_all', chr=chromosomes, allow_missing = True)
    params:
        phen = '{phen}', train_sp = 'wb_all', vali_sp = 'wb_all',
        r2 = '{r2}', kb = '{kb}'
    output:
        expand('../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{train_sp}/coeff/p1_opt/chr{chr}-coeff.tab',
               train_sp='wb_all', chr=chromosomes, allow_missing = True),
        expand('../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{train_sp}/pgs-vali/p1_opt/p1_opt-rsq_vali.tab',
               train_sp='wb_all', allow_missing = True),
        p1_plot = expand('../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{train_sp}/figs/{phen}-{vali_sp}-p1_rsq.png',
                         train_sp='wb_all', vali_sp='wb_all', allow_missing = True)
    threads: 2
    shell:
        """
        Rscript ./scripts/02-pgs/02b-choose-p1.R \
        --phen {params.phen} \
        --train_sp {params.train_sp} \
        --vali_sp {params.vali_sp} \
        --r2 {params.r2} --kb {params.kb}
        """


rule compute_scores_p1_opt_all:
    input:
        pgen    = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar    = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam    = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        samples = '../data/sample-ids/filtered/all-ids.tab',
        freqs   = expand('../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq', chr=chromosomes),
        coeff   = expand('../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{samples}/coeff/p1_opt/chr{chr}-coeff.tab',
                         chr=chromosomes, allow_missing = True),
        rsq_table = '../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{samples}/pgs-vali/p1_opt/p1_opt-rsq_vali.tab'
    params:
        phen = '{phen}', samples = '{samples}',
        r2 = '{r2}', kb = '{kb}'
    output:
        '../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{samples}/pgs-all/p1_opt/{phen}-all-pgs0.tab',
        '../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/{samples}/pgs-all/p1_opt/{phen}-pgs0-rsq.tab'
    shell:
        """
        Rscript ./scripts/02-pgs/02c-compute-p1_opt-pgs.R \
        --phen {params.phen} --samples {params.samples} \
        --r2 {params.r2} --kb {params.kb}
        """


## iterative algorithm
rule iterative_step:
    input:
        phen_f  = expand('../results/01-gwas/residuals-covars/{phen}/{phen}-{samples}-resid-covars.tab',
                         samples=['wb_all', 'other_gb_vali_all', 'other_gb_test_all'], allow_missing = True),
        pgs0    = '../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/wb_{sex}/pgs-all/p1_opt/{phen}-all-pgs0.tab',
        p1_opt0 = '../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/wb_{sex}/pgs-all/p1_opt/{phen}-pgs0-rsq.tab',
        coeff0  = expand('../results/02-pgs/initial/{phen}/r2_{r2}-kb_{kb}/wb_{sex}/coeff/p1_opt/chr{chr}-coeff.tab',
                         chr=chromosomes, allow_missing = True),
        ukb_sample = '../data/sample-ids/ukb22828_c1_b0_v3_s487256.sample',
        all_ids    = '../data/sample-ids/filtered/all-ids.tab',
        bgen       = expand('../data/imputed-genotypes/bgen/ukb_imp_chr{chr}_v3.bgen', chr=chromosomes)
    params:
        phen = '{phen}',
        train_sp = 'wb_{sex}', vali_sp = 'other_gb_vali_{sex}', test_sp = 'other_gb_test_{sex}',
        r2 = '{r2}', kb = '{kb}'
    output:
        step1_manh = '../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/wb_{sex}/figs/{phen}-step1-manh.png',
        p1_opt_comp = '../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/wb_{sex}/{phen}-pgs-performance.tab',
        coeff_cpn = expand('../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/wb_{sex}/coeff/cpn/{phen}-coeff-cpn-chr{chr}.tab',
                           chr=1, allow_missing = True),
        coeff_final = expand('../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/wb_{sex}/coeff/final/{phen}-coeff-chr{chr}.tab',
                             chr=1, allow_missing = True),
        pgs_bychr = '../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/wb_{sex}/pgs/{phen}-all-pgs-bychr.tab',
        pgs_loco = '../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/wb_{sex}/pgs/{phen}-all-pgs-loco.tab',
        pgs_full = '../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/wb_{sex}/pgs/{phen}-all-pgs-full.tab',
        step_rsq_plot = '../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/wb_{sex}/figs/{phen}-iter_step-rsq.png'
    threads: 8
    shell:
        """
        Rscript ./scripts/02-pgs/03-iterative-step.R \
        --phen {params.phen} \
        --train_sp {params.train_sp} --vali_sp {params.vali_sp} --test_sp {params.test_sp} \
        --r2 {params.r2} --kb {params.kb} \
        --threads {threads}
        """


rule mean_correction:
    input:
        phen_f   = '../results/01-gwas/residuals-covars/{phen}/{phen}-{samples}-resid-covars.tab',
        pgs_full = '../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/wb_all/pgs/{phen}-all-pgs-full.tab',
        pgs_loco = '../results/02-pgs/iterative/{phen}/r2_{r2}-kb_{kb}/wb_all/pgs/{phen}-all-pgs-loco.tab'
    params:
         phen = '{phen}',
         samples = '{samples}',
         r2 = '{r2}', kb = '{kb}'
    output:
        pgs_full_mc = '../results/02-pgs/mean-corrected/{phen}/r2_{r2}-kb_{kb}/{phen}-{samples}-pgs-full_mc.tab',
        pgs_loco_mc = expand('../results/02-pgs/mean-corrected/{phen}/r2_{r2}-kb_{kb}/{phen}-{samples}-pgs-loco{chr}_mc.tab',
                             chr = chromosomes, allow_missing = True)
    shell:
        """
        Rscript ./scripts/02-pgs/04-mean-correction.R \
        --phen {params.phen} \
        --samples {params.samples} \
        --r2 {params.r2} --kb {params.kb}
        """
 

rule reg_pgs_ac:
    input:
        phen_res_f  = '../results/01-gwas/residuals-covars/{phen}/{phen}-wb_all-resid-covars.tab',
        covars_f = '../data/covars/age-sex-batch-centre-ac-all.tab',
        pgs_full_mc = '../results/02-pgs/mean-corrected/{phen}/r2_{r2}-kb_{kb}/{phen}-wb_all-pgs-full_mc.tab',
        pgs_loco_mc = expand('../results/02-pgs/mean-corrected/{phen}/r2_{r2}-kb_{kb}/{phen}-wb_all-pgs-loco{chr}_mc.tab',
                             chr = chromosomes, allow_missing = True)
    params:
         phen = '{phen}',
         train_sp = 'wb_all',
         r2 = '{r2}', kb = '{kb}'
    output:
        '../results/02-pgs/ancestry-interactions/{phen}/r2_{r2}-kb_{kb}/{phen}-wb_all-resid-covars-pgs-full_mc-ac-centre.tab',
        '../results/02-pgs/ancestry-interactions/{phen}/r2_{r2}-kb_{kb}/{phen}-wb_all-resid_abs-covars-pgs-full_mc-ac-centre.tab',
        expand('../results/02-pgs/ancestry-interactions/{phen}/r2_{r2}-kb_{kb}/{phen}-wb_all-resid-covars-pgs-loco{chr}_mc-ac-centre.tab',
               chr = chromosomes, allow_missing = True),
        expand('../results/02-pgs/ancestry-interactions/{phen}/r2_{r2}-kb_{kb}/{phen}-wb_all-resid_abs-covars-pgs-loco{chr}_mc-ac-centre.tab',
                    chr = chromosomes, allow_missing = True)
    threads: 2
    shell:
        """
        Rscript ./scripts/02-pgs/05-pgs-ancestry-interactions.R \
        --phen {params.phen} \
        --train_sp {params.train_sp} \
        --r2 {params.r2} --kb {params.kb}
        """
