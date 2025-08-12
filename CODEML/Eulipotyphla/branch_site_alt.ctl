seqfile = codon_alignment.phy
treefile = tree_Condylura_cristata.nwk
outfile = Condylura_cristata_branch_site_alt.out

noisy = 9
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
       model  = 2       * Branch-site model
      NSsites = 2       * Model A
        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0       * Estimate omega > 1
        omega = 1
    cleandata = 1
