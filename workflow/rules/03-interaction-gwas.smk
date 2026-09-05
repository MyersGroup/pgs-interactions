# 03 INTERACTION GWAS
rule run_inter_gwas_snp_pgs:
    input:
        pgen     = '../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen',
        pvar     = '../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
        psam     = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        samples  = '../data/sample-ids/filtered/wb_all-ids.tab',
        freqs    = '../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq',
        variants = '../data/variant-ids/chr{chr}-var-ids.tab',
        phen_res_full = '../results/02-pgs/ancestry-interactions/{phen}/r2_0.9-kb_500/{phen}-{samples}-resid-covars-pgs-full_mc-ac-centre.tab',
        pgs_full = '../results/02-pgs/mean-corrected/{phen}/r2_0.9-kb_500/{phen}-wb_all-pgs-full_mc.tab',
        pgs_loco = '../results/02-pgs/mean-corrected/{phen}/r2_0.9-kb_500/{phen}-wb_all-pgs-loco{chr}_mc.tab'
    params:
        phen = '{phen}', samples = 'wb_all',
        out_pgs_full = '../results/03-interaction-gwas/plink-output/{phen}/snp-pgs/{phen}.{samples}.full.chr{chr,\d+}.{phen}.res.cov.pgs_full_mc.ac.ctr.glm.linear',
        out_pgs_loco = '../results/03-interaction-gwas/plink-output/{phen}/snp-pgs/{phen}.{samples}.loco.chr{chr,\d+}.{phen}.res.cov.pgs_full_mc.ac.ctr.glm.linear'
    output:
        '../results/03-interaction-gwas/plink-output/{phen}/snp-pgs/{phen}.{samples}.full.chr{chr,\d+}.{phen}.res.cov.pgs_full_mc.ac.ctr.glm.linear.gz',
        '../results/03-interaction-gwas/plink-output/{phen}/snp-pgs/{phen}.{samples}.loco.chr{chr,\d+}.{phen}.res.cov.pgs_full_mc.ac.ctr.glm.linear.gz'
    threads: 2
    shell:
        """
        ./scripts/03-interaction-gwas/01-run-gwas-snp-pgs.sh {params.phen} {wildcards.chr} {input.pgen} {input.pvar} {input.psam} \
        {params.samples} {input.freqs} {input.variants}

        gzip -f {params.out_pgs_full}
        gzip -f {params.out_pgs_loco}
        """


rule stepwise_reg:
    input:
        pgen = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        phen_resid = '../results/02-pgs/ancestry-interactions/{phen}/r2_0.9-kb_500/{phen}-{samples}-resid-covars-pgs-full_mc-ac-centre.tab',
        pgs_loco = expand('../results/02-pgs/mean-corrected/{phen}/r2_0.9-kb_500/{phen}-{samples}-pgs-loco{chr}_mc.tab', chr=chromosomes, allow_missing=True),
        pgs_full = '../results/02-pgs/mean-corrected/{phen}/r2_0.9-kb_500/{phen}-{samples}-pgs-full_mc.tab',
        full_stats = expand('../results/03-interaction-gwas/plink-output/{phen}/snp-pgs/{phen}.{samples}.full.chr{chr}.{phen}.res.cov.pgs_full_mc.ac.ctr.glm.linear.gz',
                            chr=chromosomes, allow_missing=True),
        loco_stats = expand('../results/03-interaction-gwas/plink-output/{phen}/snp-pgs/{phen}.{samples}.loco.chr{chr}.{phen}.res.cov.pgs_full_mc.ac.ctr.glm.linear.gz',
                            chr=chromosomes, allow_missing=True),
        res_abs_full_stats = expand('../results/03-interaction-gwas/plink-output/{phen}/res_abs/{phen}.{samples}.full.chr{chr}.{phen}.res_abs.cov.pgs_full_mc.ac.ctr.glm.linear.gz',
                                    chr=chromosomes, allow_missing=True)
    params:
         phen = '{phen}', samples = '{samples}'
    output:
        '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-loco-indep-hits.tab',
        '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-full-indep-hits.tab'
    threads: 8
    shell:
        """
        Rscript ./scripts/03-interaction-gwas/02a-stepwise-reg.R \
        --phen {params.phen} \
        --samples {params.samples}
        """


rule agg_results:
    input:
        code_desc = '../data/phenotypes/clean/code-desc-map-real.tab',
        indep_hits_loco = expand('../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-loco-indep-hits.tab',
                                 phen=traits_qn, samples='wb_all'),
        indep_hits_full = expand('../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-full-indep-hits.tab',
                                 phen=traits_qn, samples='wb_all'),
        map_posid_rsid = expand('../data/variant-ids/map-posid-rsid/posid-rsid-chr{chr}.tab', chr=chromosomes),
        pgen = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        sample_f = expand('../data/sample-ids/filtered/{samples}-ids.tab', samples='wb_all')
    output:
        '../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn.tab',
        '../results/03-interaction-gwas/indep-hits/aggregate/hits-full-qn.tab'
    threads: 4
    shell:
        'Rscript ./scripts/03-interaction-gwas/02b-aggregate-results.R'


rule agg_null_sims:
    input:
        code_desc = '../data/phenotypes/clean/code-desc-map-sim.tab',
        gwas = expand('../results/03-interaction-gwas/plink-output/{phen}/snp-pgs/{phen}.{samples}.{pgs}.chr{chr}.{phen}.res.cov.pgs_full_mc.ac.ctr.glm.linear.gz',
                      phen=traits_sim_qn, samples='wb_all', pgs=['full', 'loco'], chr=chromosomes),
        indep_hits = expand('../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-{pgs}-indep-hits.tab',
                            phen=traits_sim_qn, samples='wb_all', pgs=['full', 'loco']),
    output:
        expand('../results/03-interaction-gwas/indep-hits/sim_null-overview/sim_null-fp-qn{a}-{h2}.tab',
               a=['', '-a5'], h2=[0.3, 0.6])
    threads: 12
    shell:
        'Rscript ./scripts/03-interaction-gwas/02c-aggregate-null-sims.R'


rule sim_int_overview:
    input:
        inter_coeff = '../data/phenotypes/simulations-prep/coeff/sim-coeff-inter.tab',
        maf = expand('../data/imputed-wb_all-stats/minor-alleles/chr{chr}.tab', chr=chromosomes),
        gwas = expand('../results/03-interaction-gwas/plink-output/{phen}/snp-pgs/{phen}.{samples}.{pgs}.chr{chr}.{phen}.res.cov.pgs_full_mc.ac.ctr.glm.linear.gz',
                          phen=traits_sim_int + traits_sim_int_qn, samples='wb_all', pgs=['full', 'loco'], chr=chromosomes),
        ind_hits = expand('../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-{pgs}-indep-hits.tab',
                          phen=traits_sim_int + traits_sim_int_qn, samples='wb_all', pgs=['full', 'loco'], chr=chromosomes),
    output:
        expand('../results/03-interaction-gwas/indep-hits/sim_int-overview/sim_int{a}-lhs-sumstats.tab',
               a=['', '_a5'])
    threads: 8
    shell:
        'Rscript ./scripts/03-interaction-gwas/02d-sim_int-overview.R'
 

# GWAS by interacting SNP genotype: LOCO PGS
checkpoint check_stepwise_loco:
        output:
            mydir = directory("../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/loco/")
        shell:
            """
            touch ../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/loco/.check.done
            """

rule gwas_int_loco:
    input:
        '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-loco-indep-hits.tab',
        pgen     = '../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen',
        pvar     = '../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
        psam     = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        samples  = '../data/sample-ids/filtered/{samples}-ids.tab',
        freqs    = '../data/imputed-{samples}-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq',
        variants = '../data/variant-ids/chr{chr}-var-ids.tab',
        trait_f  = '../results/01-gwas/residuals-covars/{phen}/{phen}-{samples}-resid-covars.tab'
    params:
        phen = '{phen}', samples = 'wb_all', chr = '{chr}', snp = '{snp}'
    output:
        '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/gwas/{phen}.{samples}.{snp}.int.loco.chr{chr}.{phen}.res.cov.glm.linear.gz'
    threads: 2
    shell:
        """
        ./scripts/03-interaction-gwas/03a-gwas-by-int-geno-loco.sh \
        {params.phen} \
        {params.chr} \
        {input.pgen} \
        {input.pvar} \
        {input.psam} \
        {params.samples} \
        {input.freqs} \
        {input.variants} \
        {params.snp}
        """


rule stepwise_pair_loco:
    input:
        maf   = expand('../data/imputed-wb_all-stats/minor-alleles/chr{chr}.tab', chr=chromosomes),
        gwas  = expand('../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/gwas/{phen}.{train_sp}.{snp}.int.loco.chr{chr}.{phen}.res.cov.glm.linear.gz',
                       chr=chromosomes, allow_missing=True),
        clump = expand('../results/02-pgs/ld-clump/{phen}/p1_0.05-r2_0.1-kb_500/{train_sp}-chr{chr}-clumped.tab',
                       chr=chromosomes, allow_missing=True),
        resid = expand('../results/01-gwas/residuals-covars/{phen}/{phen}-{samples}-resid-covars.tab',
                       samples=['wb_all', 'other_gb_test_all'], allow_missing=True), 
        pgen  = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar  = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam  = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
    params:
        phen = '{phen}', train_sp = 'wb_all', test_sp = 'other_gb_test_all', snp = '{snp}'
    output:
        '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{train_sp}-loco-{snp}-pairwise-hits.tab'
    threads: 12
    run:
        if params.phen in traits_sim_sc + traits_sim_sc_qn:
            shell("touch ../results/03-interaction-gwas/indep-hits/{params.phen}/snp-pgs/hits/{params.phen}-{params.train_sp}-loco-{params.snp}-pairwise-hits.tab")
        else:
            shell("Rscript ./scripts/03-interaction-gwas/03b-stepwise-reg-loco.R \
            --phen {params.phen} \
            --train_sp {params.train_sp} \
            --test_sp {params.test_sp} \
            --snp {params.snp}")


def gwas_to_run_loco(wildcards):
        checkpoint_output = checkpoints.check_stepwise_loco.get(**wildcards).output.mydir
        
        ivals = glob_wildcards(os.path.join(checkpoint_output,
                                            "{snp}-gen1.tab")).snp
        print("ivals={}".format(ivals))
        return (expand('../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/gwas/{phen}.{samples}.{snp}.int.loco.chr{chr}.{phen}.res.cov.glm.linear.gz',
                       snp=ivals, chr=chromosomes, allow_missing=True) +
                expand('../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-loco-{snp}-pairwise-hits.tab',
                       snp=ivals, allow_missing=True))


rule aggregate_loco:
    output: "../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/.{phen}.{samples}.loco.done.txt",
    input:
        '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-loco-indep-hits.tab',
        gwas_to_run_loco
    shell:
        """
        touch {output}
        """


# GWAS by interacting SNP genotype: full PGS
checkpoint check_stepwise_full:
        output:
            mydir = directory("../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/full/")
        shell:
            """
            touch ../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/full/.check.done
            """

rule gwas_int_full:
    input:
        '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-full-indep-hits.tab',
        pgen     = '../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen',
        pvar     = '../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
        psam     = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        samples  = '../data/sample-ids/filtered/{samples}-ids.tab',
        freqs    = '../data/imputed-{samples}-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq',
        variants = '../data/variant-ids/chr{chr}-var-ids.tab',
        trait_f  = '../results/01-gwas/residuals-covars/{phen}/{phen}-{samples}-resid-covars.tab'
    params:
        phen = '{phen}', samples = 'wb_all', chr = '{chr}', snp = '{snp}'
    output:
        '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/gwas/{phen}.{samples}.{snp}.int.full.chr{chr}.{phen}.res.cov.glm.linear.gz'
    threads: 2
    shell:
        """
        ./scripts/03-interaction-gwas/03a-gwas-by-int-geno-full.sh \
        {params.phen} \
        {params.chr} \
        {input.pgen} \
        {input.pvar} \
        {input.psam} \
        {params.samples} \
        {input.freqs} \
        {input.variants} \
        {params.snp}
        """


rule stepwise_pair_full:
    input:
        maf   = expand('../data/imputed-wb_all-stats/minor-alleles/chr{chr}.tab', chr=chromosomes),
        gwas  = expand('../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/gwas/{phen}.{train_sp}.{snp}.int.full.chr{chr}.{phen}.res.cov.glm.linear.gz',
                       chr=chromosomes, allow_missing=True),
        clump = expand('../results/02-pgs/ld-clump/{phen}/p1_0.05-r2_0.1-kb_500/{train_sp}-chr{chr}-clumped.tab',
                       chr=chromosomes, allow_missing=True),
        resid = expand('../results/01-gwas/residuals-covars/{phen}/{phen}-{samples}-resid-covars.tab',
                       samples=['wb_all', 'other_gb_test_all'], allow_missing=True), 
        pgen  = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar  = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam  = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
    params:
        phen = '{phen}', train_sp = 'wb_all', test_sp = 'other_gb_test_all', snp = '{snp}'
    output:
        '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{train_sp}-full-{snp}-pairwise-hits.tab'
    threads: 36
    run:
        if params.phen in traits_sim_sc + traits_sim_sc_qn:
            shell("touch ../results/03-interaction-gwas/indep-hits/{params.phen}/snp-pgs/hits/{params.phen}-{params.train_sp}-full-{params.snp}-pairwise-hits.tab")
        else:
            shell("Rscript ./scripts/03-interaction-gwas/03b-stepwise-reg-full.R \
            --phen {params.phen} \
            --train_sp {params.train_sp} \
            --test_sp {params.test_sp} \
            --snp {params.snp}")


def gwas_to_run_full(wildcards):
        checkpoint_output = checkpoints.check_stepwise_full.get(**wildcards).output.mydir
        
        ivals = glob_wildcards(os.path.join(checkpoint_output,
                                            "{snp}-gen1.tab")).snp
        print("ivals={}".format(ivals))
        return (expand("../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/gwas/{phen}.{samples}.{snp}.int.full.chr{chr}.{phen}.res.cov.glm.linear.gz",
                       snp=ivals, chr=chromosomes, allow_missing=True) +
                expand('../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-full-{snp}-pairwise-hits.tab',
                       snp=ivals, allow_missing=True))


rule aggregate_full:
    output: "../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/.{phen}.{samples}.full.done.txt",
    input:
        '../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-full-indep-hits.tab',
        gwas_to_run_full
    shell:
        """
        touch {output}
        """


rule aggregate_pairwise:
    input:
        psam = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        code_desc = '../data/phenotypes/clean/code-desc-map-real.tab',
        indep_hits = expand('../results/03-interaction-gwas/indep-hits/{phen}/snp-pgs/hits/{phen}-{samples}-{pgs}-indep-hits.tab',
                                 phen=traits_qn, samples='wb_all', pgs=['loco', 'full']),
        gwas = expand('../results/01-gwas/plink-output/{phen}/{samples}.chr{chr}.{phen}.res.cov.glm.linear.gz',
                      phen=traits_qn, samples='wb_all', chr=chromosomes),
        resid_f = expand('../results/01-gwas/residuals-covars/{phen}/{phen}-{samples}-resid-covars.tab',
                         phen=traits_qn, samples='wb_all'),
        map_posid_rsid = expand('../data/variant-ids/map-posid-rsid/posid-rsid-chr{chr}.tab', chr=chromosomes),
        vep = expand('../data/annotations/vep/anno/vep-chr{chr}.most_severe.tab.gz', chr=chromosomes),
        annovar_refGene = expand('../data/annotations/annovar/gene-anno/annovar-chr{chr}.refGene.variant_function.gz', chr=chromosomes),
        annovar_ensGene = expand('../data/annotations/annovar/gene-anno/annovar-chr{chr}.ensGene.variant_function.gz', chr=chromosomes),
    output:
        '../results/03-interaction-gwas/indep-hits/aggregate/pairwise-hits-all.tab',
        '../results/03-interaction-gwas/indep-hits/aggregate/pairwise-hits-high_conf.tab',
        '../results/03-interaction-gwas/indep-hits/aggregate/pairwise-hits-others.tab'
    threads: 8
    shell:
        'Rscript ./scripts/03-interaction-gwas/03c-aggregate-results.R'


rule find_indep_hits_tags:
    input:
        lhs_hits_loco_qn = '../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn.tab',
        pairwise_hits = '../results/03-interaction-gwas/indep-hits/aggregate/pairwise-hits-all.RData',
        map_posid_rsid = expand('../data/variant-ids/map-posid-rsid/posid-rsid-chr{chr}.tab', chr=chromosomes),
        maf = expand('../data/imputed-wb_all-stats/minor-alleles/chr{chr}.tab', chr=chromosomes),
        annovar_refGene = expand('../data/annotations/annovar/gene-anno/annovar-chr{chr}.refGene.variant_function.gz', chr=chromosomes),
        annovar_ensGene = expand('../data/annotations/annovar/gene-anno/annovar-chr{chr}.ensGene.variant_function.gz', chr=chromosomes),
        vep = expand('../data/annotations/vep/anno/vep-chr{chr}.most_severe.tab.gz', chr=chromosomes),
        vep_consq_rank = '../data/annotations/vep/vep-var-consq-rel113_202410.csv'
    output:
        '../results/03-interaction-gwas/indep-hits/aggregate/tags/hits-loco-qn-tags-0.8.tab',
        '../results/03-interaction-gwas/indep-hits/aggregate/tags/hits-loco-qn-tags-0.8-most-severe.tab'
    threads: 4
    shell:
         'Rscript ./scripts/03-interaction-gwas/04-find-lhs-rhs-tags.R'
