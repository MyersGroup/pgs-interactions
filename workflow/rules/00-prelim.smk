# 00 PRELIMINARY STEP
rule prepare_covars:
    input:
        '../data/phenotypes/ukb-raw/ukb45188.tab',
        '../data/sample-ids/withdrawals/w27960_20210201.csv',
        '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        '../data/covars/ancestry/v2_487409.rds'
    output:
        '../data/covars/age-sex-batch-centre-ac-all.tab'
    threads: 2
    shell:
        'Rscript ./scripts/00-prelim/01-prepare-covars.R'


rule sample_filtering:
    input:
        '../data/covars/age-sex-batch-centre-ac-all.tab',
        '../data/covars/ancestry/wb_ids.rds',
        '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        '../data/relatedness/ukb_rel_a27960_s488224.dat'
    output:
        expand('../data/sample-ids/filtered/{samples}-ids.tab',
               samples = ['all', 'wb_all', 'other_gb_vali_all', 'other_gb_test_all'])
    threads: 2
    shell:
        'Rscript ./scripts/00-prelim/02-sample-filtering.R'


rule prepare_phenotypes:
    input:
        '../data/sample-ids/filtered/all-ids.tab',
        '../data/phenotypes/ukb-raw/ukb45188.tab',
        '../data/phenotypes/ukb-raw/ukb26618.tab',
        '../data/phenotypes/clean/code-desc-map-real.tab',
        '../data/covars/age-sex-batch-centre-ac-all.tab'
    output:
        expand('../data/phenotypes/clean/{phen}.tab', phen=traits)
    threads: 3
    shell:
        'Rscript ./scripts/00-prelim/03-prepare-phenotypes.R'


rule variant_stats:
    input:
        pgen = '../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen',
        pvar = '../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
        psam = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        samples = '../data/sample-ids/filtered/wb_all-ids.tab'
    params:
        out_path   = '../data/imputed-wb_all-stats',
        out_prefix = 'ukb_imp_chr{chr}_v3'
    output:
        '../data/imputed-wb_all-stats/allele-counts/ukb_imp_chr{chr,\d+}_v3.acount',
        '../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr{chr,\d+}_v3.afreq',
        '../data/imputed-wb_all-stats/missing-rates/ukb_imp_chr{chr,\d+}_v3.vmiss',
        '../data/imputed-wb_all-stats/hardy-eq/ukb_imp_chr{chr,\d+}_v3.hardy'
    threads: 4
    shell:
        './scripts/00-prelim/04-variant-stats.sh {input.pgen} {input.pvar} {input.psam} \
        {input.samples} {params.out_path} {params.out_prefix}'


rule variant_filtering:
    input:
        expand('../data/imputed-wb_all-stats/allele-counts/ukb_imp_chr{chr}_v3.acount', chr=chromosomes),
        expand('../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq', chr=chromosomes),
        expand('../data/imputed-wb_all-stats/missing-rates/ukb_imp_chr{chr}_v3.vmiss', chr=chromosomes),
        expand('../data/imputed-wb_all-stats/hardy-eq/ukb_imp_chr{chr}_v3.hardy', chr=chromosomes),
        expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes)
    output:
        expand('../data/variant-ids/chr{chr}-var-ids.tab', chr=chromosomes),
        expand('../data/variant-ids/map-posid-rsid/posid-rsid-chr{chr}.tab', chr=chromosomes),
        expand('../data/imputed-wb_all-stats/minor-alleles/chr{chr}.tab', chr=chromosomes)
    threads: 4
    shell:
        'Rscript ./scripts/00-prelim/05-variant-filtering.R'


rule get_anno_annovar:
    input:
        var_list = expand('../data/variant-ids/map-posid-rsid/posid-rsid-chr{chr}.tab', chr=chromosomes)
    output:
        annovar_input    = expand('../data/annotations/annovar/inputs/chr{chr}.avinput', chr=chromosomes),
        var_fn_ensGene   = expand('../data/annotations/annovar/gene-anno/annovar-chr{chr}.ensGene.variant_function.gz', chr=chromosomes),
        var_fn_refGene   = expand('../data/annotations/annovar/gene-anno/annovar-chr{chr}.refGene.variant_function.gz', chr=chromosomes),
        var_fn_multianno = expand('../data/annotations/annovar/gene-anno/annovar-chr{chr}.hg19_multianno.txt.gz', chr=chromosomes)
    threads: 2
    shell:
        'Rscript ./scripts/00-prelim/06a-get-annotations-annovar.R'


rule get_anno_vep:
    input:
        var_list = '../data/variant-ids/map-posid-rsid/posid-rsid-chr{chr}.tab'
    params:
        chr = '{chr}'
    output:
        vep_input       = '../data/annotations/vep/inputs/chr{chr}.tab.gz',
        anno_most_sev   = '../data/annotations/vep/anno/vep-chr{chr}.most_severe.tab.gz',
        anno_everything = '../data/annotations/vep/anno/vep-chr{chr}.everything.tab.gz'
    threads: 2
    shell:
        """
        module load Anaconda3/2024.02-1
        eval "$(conda shell.bash hook)"
        conda activate pgs-interactions 

        Rscript ./scripts/00-prelim/06b-get-annotations-vep.R \
        --chr {params.chr}
        """
 

rule simul_make_coeff:
    input:
        expand('../data/variant-ids/chr{chr}-var-ids.tab', chr=chromosomes),
        expand('../data/imputed-wb_all-stats/minor-alleles/chr{chr}.tab', chr=chromosomes)
    output:
        expand('../data/phenotypes/simulations-prep/coeff/sim-coeff-{ncs}_{arch}-chr{chr}.tab',
               ncs=['100', '1k', '10k'], arch=['rg', 'cl'], chr=chromosomes)
    threads: 2
    shell:
        'Rscript ./scripts/00-prelim/07a-simul-make-coeff-vec.R'


rule simul_comp_scores:
    input:
        pgen    = '../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen',
        pvar    = '../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar',
        psam    = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        samples = '../data/sample-ids/filtered/all-ids.tab',
        freqs   = '../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr{chr}_v3.afreq',
        coeff   = '../data/phenotypes/simulations-prep/coeff/sim-coeff-{ncs}_{arch}-chr{chr}.tab'
    params:
        out_prefix = '../data/phenotypes/simulations-prep/pgs/sim-pgs-{ncs}_{arch}-chr{chr}'
    output:
        temp('../data/phenotypes/simulations-prep/pgs/sim-pgs-{ncs}_{arch}-chr{chr}.sscore'),
        temp('../data/phenotypes/simulations-prep/pgs/sim-pgs-{ncs}_{arch}-chr{chr}.log')
    threads: 2
    shell:
        './scripts/00-prelim/07b-simul-compute-score.sh {input.pgen} {input.pvar} {input.psam} \
        {input.samples} {input.freqs} {input.coeff} {params.out_prefix}' 


rule simul_make_phen:
    input:
        expand('../data/phenotypes/simulations-prep/pgs/sim-pgs-{ncs}_{arch}-chr{chr}.sscore',
               ncs=['100', '1k', '10k'], arch=['rg', 'cl'], chr=chromosomes)
    output:
        expand('../data/phenotypes/clean/f.sim_{ncs}_{arch}_{tr}_{h2}.tab',
               ncs=['100', '1k', '10k'], arch=['rg', 'cl'], tr=['nt', 'sgbn', 'sgan', 'sgscbn', 'sgscan'], h2=[0.3, 0.6]),
        expand('../data/phenotypes/clean/f.sim_{ncs}_{arch}_{tr}_a5_{h2}.tab',
               ncs=['100', '1k', '10k'], arch=['rg', 'cl'], tr=['nt', 'sgscbn', 'sgscan'], h2=[0.3, 0.6])
    threads: 2
    shell:
        'Rscript ./scripts/00-prelim/07c-simul-make-phen.R'


rule simul_make_int:
    input:
        null_sim_coeff = expand('../data/phenotypes/simulations-prep/coeff/sim-coeff-{ncs}_{arch}-chr{chr}.tab',
                                ncs='1k', arch='rg', chr=chromosomes),
        all_vars_list = expand('../data/variant-ids/chr{chr}-var-ids.tab', chr=chromosomes),
        minor_alleles = expand('../data/imputed-wb_all-stats/minor-alleles/chr{chr}.tab', chr=chromosomes),
        pgen          = expand('../data/imputed-genotypes/ukb_imp_chr{chr}_v3.pgen', chr=chromosomes),
        pvar          = expand('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr{chr}_v3.pvar', chr=chromosomes),
        psam          = '../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam',
        all_ids       = '../data/sample-ids/filtered/all-ids.tab',
    output:
        coeff    = '../data/phenotypes/simulations-prep/coeff/sim-coeff-inter.tab',
        sim_add  = '../data/phenotypes/simulations-prep/f.sim_int_1k_add.tab',
        sim_phen = expand('../data/phenotypes/clean/f.sim_int_1k_e{effect_sz}_{inter}_{alpha}0.6.tab',
                          effect_sz=['001', '010', '050', '100'], inter=['all', 'no100'], alpha=['', 'a5_'])
    threads: 8
    shell:
         'Rscript ./scripts/00-prelim/07d-positive-sims.R'
