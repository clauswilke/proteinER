RequireVersion("2.3.4");

/* Execute a relative protein rates inference (analogous to Rate4Site approach) in HyPhy version >= 2.3.4 
   
   Provided options, in order, represent:
        + 0: "/Users/sjspielman/aa.fasta": Full path to FASTA-formatted protein alignment
        + 1: "/Users/sjspielman/aa.tre": Full path to newick-formatted phylogeny corresponding to alignment
        + 2: "WAG": Evolutionary model to use when fitting rates. Options include JC, WAG, LG, JTT.
        + 3: "Yes": Use the +F version of amino acid model. Options Yes/No.
*/

LoadFunctionLibrary("ProteinAnalyses/relative_prot_rates.bf", {
                    "0": "/path/to/HRH1_aligned.fasta",
                    "1": "/path/to/RAxML_bestTree.HRH1_tree",
                    "2": "WAG",    
                    "3": "Yes"});

            