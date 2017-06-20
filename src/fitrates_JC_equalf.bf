LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;
OPTIMIZATION_PRECSION = 0.0000000001;
ACCEPT_BRANCH_LENGTHS=0;

global r_global;
global t_global;
#include "setup_helper.ibf"; // code to read setup

// Setup is read from default from file "setup.txt".
// If this file doesn't exist, specify setup file from
// command line like this:
//     HYPHYMP FEL.bf <<< "mysetup.txt"
setup = readSetup();
infile = setup["INFILE"];
outfile = setup["OUTFILE"];
site_dupl = 1;
#include "JC_aa.mdl"; 


DataSet raw_data  = ReadDataFile(infile);
DataSetFilter filtered_data = CreateFilter(raw_data, 1);

fprintf(stdout, "Step 1: Global optimization of branch lengths.\n");
Model JCfull = (JC69_t, JC_freqs, 1);
UseModel(USE_NO_MODEL);
UseModel(JCfull);
Tree full_tree = DATAFILE_TREE;

LikelihoodFunction full_ln_likfn = (filtered_data, full_tree);
ReplicateConstraint ("this1.?.t:=t_global", full_tree);

Optimize(full_res, full_ln_likfn);

fprintf(stdout, "Step 2: Site-specific rate scalar computation.\n");
header = "site,rate,t,lnL\n";
fprintf(outfile, CLEAR_FILE, header);
fprintf(stdout, header);

nsites = filtered_data.sites;
for (global site_count = 0; site_count < nsites; site_count = site_count+site_dupl)
{
       
    // single amino acid
    filter_string = "";
    filter_string = filter_string + (site_count) + "-" + (site_count+(site_dupl-1));
    DataSetFilter site_filter = CreateFilter(filtered_data, 1, filter_string, "", "");
	
    Model JCsite = (JC69_rt, JC_freqs, 1);
    UseModel(USE_NO_MODEL);
    UseModel(JCsite);
    Tree site_tree = DATAFILE_TREE;
    LikelihoodFunction site_ln_likfn = (site_filter, site_tree);
    
    ReplicateConstraint("this1.?.r:=r_global", site_tree);
    ReplicateConstraint("this1.?.t:=this2.?.t__", site_tree, full_tree);
     
    Optimize(site_res, site_ln_likfn);
    site_lnlk = site_res[1][0];
    line = Format(site_count+1, 0, 0) + "," + Format(r_global, 5, 10) + "," + Format(t_global, 5, 10) + "," + Format(site_lnlk, 5, 10) + "\n";
    fprintf(stdout, line);
    fprintf(outfile, line);
}
