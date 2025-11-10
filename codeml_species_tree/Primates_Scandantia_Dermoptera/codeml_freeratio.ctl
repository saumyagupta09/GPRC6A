      seqfile = codon_alignment.phy              * Input alignment file (PHYLIP format)
     treefile = Primates_Scandentia_Dermoptera.nwk                   * Input tree file (Newick format)
      outfile = mlc                     * Output file name

        noisy = 9                       * 0–9: Verbosity of screen output
      verbose = 1                       * 1 = detailed output
      runmode = 0                       * 0 = use provided tree

      seqtype = 1                       * 1 = codon sequences
    CodonFreq = 2                       * 2 = F3x4 codon frequency model

        clock = 0                       * 0 = no molecular clock
        model = 1                       * 1 = free-ratio model (ω varies by branch)

      NSsites = 0                       * 0 = no site model
        icode = 0                       * 0 = universal genetic code

    fix_kappa = 0                       * 0 = estimate kappa
        kappa = 2                       * initial kappa guess

    fix_omega = 0                       * 0 = estimate omega
        omega = 1                       * initial omega guess

    cleandata = 1                       * 1 = remove sites with ambiguous data

       getSE = 0                        * 0 = don't calculate standard errors
 RateAncestor = 0                       * 0 = no ancestral state reconstruction
