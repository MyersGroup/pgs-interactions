# 01 GWAS
rule regress_out_covars:
    input:
        phen_f   = '../data/phenotypes/clean/{phen}.tab',
        covars_f = '../data/covars/age-sex-batch-centre-ac-all.tab',
        samples  = '../data/sample-ids/filtered/{samples}-ids.tab'
    params:
        phen = '{phen}', samples = '{samples}'
    output:
        '../results/01-gwas/residuals-covars/{phen}_qn/{phen}_qn-{samples}-resid-covars.tab'
    shell:
        """
        Rscript ./scripts/01-gwas/00-regress-out-covars.R \
        --phen {params.phen} \
        --samples {params.samples}
        """


rule run_gwas_quant:
    input:
        pgen     = '../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen',
        pvar     = '../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
        psam     = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        sample_f = '../data/sample-ids/filtered/{samples}-ids.tab',
        freqs    = '../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq',
        variants = '../data/variant-ids/chr{chr}-var-ids.tab',
        trait_f  = '../results/01-gwas/residuals-covars/{trait}/{trait}-{samples}-resid-covars.tab'
    params:
        out_prefix = '../results/01-gwas/plink-output/{trait}/{samples}.chr{chr,\d+}',
        out_full   = '../results/01-gwas/plink-output/{trait}/{samples}.chr{chr,\d+}.{trait}.res.cov.glm.linear'
    output:
        '../results/01-gwas/plink-output/{trait}/{samples}.chr{chr,\d+}.{trait}.res.cov.glm.linear.gz'
    threads: 2
    shell:
        """
        ./scripts/01-gwas/01-run-gwas-quant.sh \
        {input.pgen} \
        {input.pvar} \
        {input.psam} \
        {input.sample_f} \
        {input.freqs} \
        {input.variants} \
        {input.trait_f} \
        {params.out_prefix}

        gzip {params.out_full}
        """
