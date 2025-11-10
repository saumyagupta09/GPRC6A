seqfile = codon_alignment.phy        * same alignment file
     treefile = Afrotheria_Xenarthra.nwk     * same tree file
      outfile = codeml_one_ratio.out  * different output file

        noisy = 9
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2

        clock = 0
       aaDist = 0

       model = 0                * one-ratio model
      NSsites = 0

        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 1

    cleandata = 1
