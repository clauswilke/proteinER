RequireVersion("2.3.3");

/* Execute a site-wise dN/dS inference using FEL in HyPhy version >= 2.3. 
   
   Provided options, in order, represent:
        + 0: "Universal" : Select the universal genetic code for analysis. Note that this will argument should be changed if non-universal codes (i.e. mitochondrial DNA).
        + 1: "/Users/sjspielman/codon.fasta": Full path to FASTA-formatted alignment
        + 2: "/Users/sjspielman/codon.tre": Full path to newick-formatted phylogeny corresponding to alignment
        + 3: "All": Test All branches for selection (as opposed to Internal, Leaves)
        + 4: "No": Do **not** use synonymous rate variation when inferring site rates (i.e. dS=1 at all sites)
        + 5: "0.1": P-value threshold for calling sites as positively selected in markdown screen output produced during analysis
*/


LoadFunctionLibrary("SelectionAnalyses/FEL.bf", {
    "0": "Universal",
    "1": "/Users/dariyasydykova/Desktop/projects/proteinER/measuring_dNdS/HRH1_aligned_codon.fasta",
    "2": "/Users/dariyasydykova/Desktop/projects/proteinER/measuring_dNdS/RAxML_bestTree.HRH1_tree",
    "3": "All",
    "4": "No", 
    "5": "0.1"
});