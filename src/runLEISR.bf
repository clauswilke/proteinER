RequireVersion("2.3.7");

/* Execute a LEISR inference in HyPhy version >= 2.3.7 
   
   Provided options, in order, represent:
        + 0: "/Users/username/aa.fasta": Full path to FASTA-formatted protein alignment.
        + 1: "/Users/username/aa.tre": Full path to newick-formatted phylogeny corresponding to alignment.
        + 2: "Protein": Specifies that a Protein (as opposed to "Nucleotide") analysis should be executed
        + 3: "WAG": Evolutionary model to use when fitting rates. Options include JC69, WAG, LG, JTT,  and specialist models cpREV, mtMAM, HIVb/w, and AB.
        + 4: "No": Use a model of rate heterogeneity when optimizing branch lengths. Options include "No", "GDD" (four-category general discrete distribution), and "Gamma"
*/

LoadFunctionLibrary("LEISR.bf", {
                    "0": "/path/to/HRH1_aligned.fasta",
                    "1": "/path/to/RAxML_bestTree.HRH1_tree",
                    "2": "Protein",
                    "3": "WAG",    
                    "4": "No"});
