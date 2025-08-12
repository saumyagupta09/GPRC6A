seqfile = codon_alignment.phy         * codon sequence alignment file
     treefile =  Carnivora_Pholidota.nwk       * Newick tree file
      outfile = codeml_freeratio_results.out      * main result file

        noisy = 9               * verbosity (0–9)
      verbose = 1               * detailed output
      runmode = 0               * user-defined tree

      seqtype = 1               * codons
    CodonFreq = 2               * F3x4 codon frequency model

        clock = 0               * no clock
       aaDist = 0

       model = 1                * branch model (model=1 = free-ratio)
      NSsites = 0               * no site models

        icode = 0               * universal genetic code
    fix_kappa = 0
        kappa = 2               * initial kappa (transition/transversion ratio)
    fix_omega = 0
        omega = 1               * initial omega (dN/dS)
    
    cleandata = 1               * remove sites with gaps/ambiguity

