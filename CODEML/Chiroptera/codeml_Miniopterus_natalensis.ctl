 seqfile = codon_alignment.phy     * input PHYLIP alignment file
     treefile = tree_Miniopterus_natalensis.nwk         * tree file with one branch labeled as #1
      outfile = branchmodel_Miniopterus_natalensis_results.txt * output file

        noisy = 9                       * verbosity level
      verbose = 1                       * detailed output
      runmode = 0                       * user tree

      seqtype = 1                       * 1:codons
    CodonFreq = 2                       * F3x4 codon frequencies

        clock = 0                       * 0: no clock
       aaDist = 0                       * 0: equal aa distances
       model = 2                        * branch model
     NSsites = 0                       * 0: one site class

     icode = 0                          * genetic code (0 = universal)
 fix_kappa = 0                          * estimate kappa
     kappa = 2                          * initial value of kappa

 fix_omega = 0                          * estimate omega
     omega = 0.5                        * initial omega

    cleandata = 1                       * remove sites with ambiguity (1: yes)
