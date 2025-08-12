seqfile = codon_alignment.phy
treefile = tree_Condylura_cristata.nwk
outfile = Condylura_cristata_branch_site_null.out

noisy = 9
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
       model  = 2
      NSsites = 2
        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 1       * Fix omega = 1 (null model)
        omega = 1
    cleandata = 1
