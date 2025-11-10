seqfile = codon_alignment.phy
     treefile = Rodentia_Lagomorpha.nwk
      outfile = codeml_null_model.out

        noisy = 9
      verbose = 1
      runmode = 0

      seqtype = 1      * 1 = codons
    CodonFreq = 2      * F3x4
        clock = 0
       aaDist = 0
       model  = 0      * One-ratio model
      NSsites = 0
        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.5    * Initial guess
    cleandata = 1

