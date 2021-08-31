          seed = -1
       seqfile = fus_proteins_dating10_mcmctree_short.phy
      treefile = fus_proteins_62T_iqtree_genepart.treefile_dating
      mcmcfile = mcmc.txt
       outfile = mcmctree_step2_out.txt
         ndata = 1
       seqtype = 2    	 * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2    	 * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = 2    	 * 1: global clock; 2: independent rates; 3: correlated rates
     cleandata = 0    	 * remove sites with ambiguity data (1:yes, 0:no)?
       BDparas = 1 1 0   * birth, death, sampling
   rgene_gamma = 1 5   	 * gamma prior for overall rates for genes
  sigma2_gamma = 1 10    * gamma prior for sigma^2     (for clock=2 or 3)
         print = 1
        burnin = 2000
      sampfreq = 10
       nsample = 20000
