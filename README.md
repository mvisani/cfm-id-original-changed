# README #

## What is CFM-ID

CFM-ID provides a method for accurately and efficiently identifying metabolites in spectra generated by electrospray tandem mass spectrometry (ESI-MS/MS). The program uses Competitive Fragmentation Modeling to produce a probabilistic generative model for the MS/MS fragmentation process and machine learning techniques to adapt the model parameters from data.  A web service can be found at (<http://cfmid4.wishartlab.com/>)

## What CFM-ID can do?

* **Spectra Prediction**, given chemical structure. This task predicts low/10V, medium/20V, and high/40V energy ESI-MS/MS spectra for an input structure provided in SMILES or InChI format.
* **Peak Assignment**: Annotating the peaks in set of spectra given a known chemical structure. This task takes a set of three input spectra (for ESI spectra, low/10V, medium/20V, and high/40V energy levels) in peak list format and a chemical structure in SMILES or InChI format, then assigns a putative fragment annotation to the peaks in each spectrum.
* **Compound Identification**: Predicted ranking of possible candidate structures for a target spectrum. This task takes a set of three input spectra (for ESI spectra, low/10V, medium/20V, and high/40V energy levels) in peak list format, and ranks a list of candidate structures according to how well they match the input spectra. This candidate list need to be provided by the user. Chemical classes are predicted for each candidate molecule. The original similarity score used in the ranking was computed (Dice or DotProduct) by comparing the predicted spectra of a candidate compound with the input spectra.

## What is this repository for?

* This is a repository for CFM-ID 4 MSML (mass spectrum machine learning)
* CFM-ID MSRB (mass spectrum rule based) can be found at (<https://bitbucket.org/wishartlab/msrb-fragmenter/src/master/>)
* Original source code from CFM-ID 2.0 Can be found at (<https://sourceforge.net/p/cfm-id/wiki/Home/>)

## License

This is project is using  GNU LESSER GENERAL PUBLIC LICENSE 2.1, for detail refers to LICENSE.md file.

## What include in this tool set for? ##

* cfm-predict
* cfm-id
* cfm-id-precomputed
* cfm-annotate
* cfm-train
* fraggraph-gen
  
## How do I get set up? ##

### Install from source code ###

* Please check cfm/INSTALL.md
* Note, Only Insatllation on linux and Mac has been verified, while install on Windows from source code is possible

### Use Pre Build Docker Image ###

Docker Image can be found at <https://hub.docker.com/r/wishartlab/cfmid>.

    docker pull wishartlab/cfmid:latest

Docker Image for NPS-MS

    docker pull wishartlab/cfmid:nps-ms_1.0.0

## cfm-predict ##

**cfm-predict** predicts spectra for an input molecule given a pre-trained CFM model. It can also work in batch mode, predicting spectra for a list of molecules in an input file.

### cfm-predict usage ###

    cfm-predict <smiles_or_inchi_or_file> <prob_thresh> <param_file> <config_file> <annotate_fragments> <output_file_or_dir> <apply_postproc> <suppress_exceptions>

* **smiles_or_inchi_or_file**: The smiles or inchi string of the structure whose spectra you want to predict. Or alternatively a .txt file containing a list of space-separated (id, smiles_or_inchi) pairs one per line. e.g.

    Molecule1 CCCNNNC(O)O  
    Molecule2 InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3  

* **prob_thresh**: (optional) The probability below which to prune unlikely fragmentations during fragmentation graph generation (default 0.001).

* **param_file**: (optional) The filename where the parameters of a trained cfm model can be found (if not given, assumes param_output.log in current directory). This file is the output of cfm-train. Pre-trained models as used in the above publication can be found in the supplementary data for that paper stored within the source tree of this project.

* **config_file**: (optional) The filename where the configuration parameters of the cfm model can be found (if not given, assumes param_config.txt in current directory). This needs to match the file passed to cfm-train during training. See cfm-train documentation below for further details.

* **annotate_fragments**:(optional) Whether to include fragment information in the output spectra (0 = NO (DEFAULT), 1 = YES). Note: ignored for msp/mfg output.

* **output_file_or_dir**: (optional) The filename of the output spectra file to write to (if not given, prints to stdout). In case of batch mode using file input above, this is used to specify the name of a directory where the output files (<id>.log) will be written (if not given, uses current directory), OR an msp or mgf file.

* **apply_postproc**: (optional) Whether or not to post-process predicted spectra to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks (whichever comes first) (0 = OFF, 1 = ON (default)). If turned off, will output a peak for every possible fragment of the input molecule, as long as the prob_thresh argument above is set to 0.0.

* **suppress_exceptions**: (optional) Suppress most exceptions so that the program returns normally even when it fails to produce a result (0 = OFF (default), 1 = ON).

### Expected cfm-predict output ###

    cfm-predict CCCC 0.001 metab_ce_cfm/param_output0.log metab_ce_cfm/param_config.txt

    #ESI-MS/MS [M+H]+ Spectra
    #PREDICTED BY CFM-ID 4.0.4
    energy0
    41.00219107 0.6506146774 3 (0.36197)
    43.01784114 100 2 (55.635)
    44.99710569 0.1661752377 1 (0.092452)
    59.01275576 1.510217985 4 (0.84021)
    61.02840582 77.41540439 0 (43.07)
    energy1
    41.00219107 0.5370744473 3 (0.41541)
    43.01784114 100 2 (77.347)
    44.99710569 2.439734481 1 (1.8871)
    59.01275576 0.6513231026 4 (0.50378)
    61.02840582 25.65990247 0 (19.847)
    energy2
    41.00219107 2.298055961 3 (1.631)
    43.01784114 100 2 (70.973)
    44.99710569 16.85736078 1 (11.964)
    59.01275576 1.810322795 4 (1.2848)
    61.02840582 19.93368222 0 (14.147)

    0 61.02840582 CC(O)=[OH+]
    1 44.99710569 O=C=[OH+]
    2 43.01784114 C#C[OH2+]
    3 41.00219107 [C+]#CO
    4 59.01275576 [CH+]=C(O)O

### Run cfm-predict in a docker container ###

Assuming your current directory is in the working directory.

#### For Linux Bash ####

To predict ESI-MS/MS [M+H]+ spectra  

    sudo docker run --rm=true -v $(pwd):/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cfm-predict 'CC(C)NCC(O)COC1=CC=C(CCOCC2CC2)C=C1' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 /cfmid/public/[M+H]+/myout"

To predict ESI-MS/MS [M-H]- spectra  

    sudo docker run --rm=true -v $(pwd):/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cfm-predict 'CC(C)NCC(O)COC1=CC=C(CCOCC2CC2)C=C1' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 /cfmid/public/[M-H]-/myout"

#### For Windows PoweShell ####

To predict ESI-MS/MS [M+H]+ spectra  

    docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cfm-predict 'CC(C)NCC(O)COC1=CC=C(CCOCC2CC2)C=C1' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 /cfmid/public/[M+H]+/myout"

To predict ESI-MS/MS [M-H]- spectra  

    docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cfm-predict 'CC(C)NCC(O)COC1=CC=C(CCOCC2CC2)C=C1' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 /cfmid/public/[M-H]-/myout"

#### For NPS-MS 

To predict ESI-MS/MS [M+H]+ spectra on Linux  

    docker run --rm=true -v $(pwd):/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cfm-predict 'CC(C)NCC(O)COC1=CC=C(CCOCC2CC2)C=C1' 0.001 /trained_models_npsms/[M+H]+/param_output.log /trained_models_npsms/[M+H]+/param_config.txt 1 /cfmid/public/[M+H]+/myout"

To predict ESI-MS/MS [M+H]+ spectra on Windows   

    docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cfm-predict 'CC(C)NCC(O)COC1=CC=C(CCOCC2CC2)C=C1' 0.001 /trained_models_npsms/[M+H]+/param_output.log /trained_models_npsms/[M+H]+/param_config.txt 1 /cfmid/public/[M+H]+/myout"
## cfm-id ##

Given an input spectrum and a list of candidate smiles (or inchi) strings, **cfm-id** computes the predicted spectrum for each candidate and compares it to the input spectrum. It returns a ranking of the candidates according to how closely they match. The spectrum prediction is done using a pre-trained CFM model.

### cfm-id usage ###

    cfm-id <spectrum_file> <id> <candidate_file> <num_highest> <ppm_mass_tol> <abs_mass_tol> <prob_thresh> <param_file> <config_file> <score_type> <apply_postprocessing> <output_file> <output_msp_or_mgf>

* **spectrum_file**: The filename where the input spectra can be found. This can be a .msp file in which the desired spectrum is listed under a corresponding id (next arg). Or it could be a single file with a list of peaks 'mass intensity' delimited by lines, with either 'low','med' and 'high' lines beginning spectra of different energy levels, or 'energy0', 'energy1', etc. e.g.

    energy0
    65.02 40.0
    86.11 60.0
    energy1
    65.02 100.0 ... etc

* **id**: An identifier for the target molecule (Used to retrieve input spectrum from msp (if used). Otherwise not used but printed to output, in case of multiple concatenated results)

* **candidate_file**: The filename where the input list of candidate structures can be found as line separated 'id smiles_or_inchi' pairs.

* **num_highest** (optional): The number of (ranked) candidates to return or -1 for all (if not given, returns all in ranked order).

* **ppm_mass_tol**: (optional) The mass tolerance in ppm to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs (if not given defaults to 10ppm).

* **abs_mass_tol**: (optional) The mass tolerance in abs Da to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs ( if not given defaults to 0.01Da).

* **prob_thresh**: (optional) The probability below which to prune unlikely fragmentations (default 0.001)

* **param_file**: (optional) The filename where the parameters of a trained cfm model can be found (if not given, assumes param_output.log in current directory). This file is the output of cfm-train. Pre-trained models as used in the above publication can be found in the supplementary data for that paper stored within the source tree of this project.

* **config_file**: (optional) The filename where the configuration parameters of the cfm model can be found (if not given, assumes param_config.txt in current directory). This needs to match the file passed to cfm-train during training. See cfm-train documentation below for further details.

* **score_type**: (optional) The type of scoring function to use when comparing spectra. Options: Dice(default), DotProduct.

* **apply_postprocessing**: (optional) Whether or not to post-process predicted spectra to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks (whichever comes first) (0 = OFF (default for EI-MS), 1 = ON (default for ESI-MS/MS)).

* **output_file**: (optional) The filename of the output spectra file to write to (if not given, prints to stdout)

* **output_msp_or_mgf**: (optional) The filename for an output msp or mgf file to record predicted candidate spectra (if not given, doesn't save predicted spectra)

### Expected cfm-id output ###

    cfm-id example_spec.txt AN_ID example_candidates.txt 5 10.0 0.01 0.001 metab_ce_cfm/param_output0.log metab_ce_cfm/param_config.txt DotProduct

    TARGET ID: NO_ID
    1 0.38798085 18232127 NC(Cc1ccc(O)cc1)C(=O)NC(CO)C(=O)NC(CC(=O)O)C(=O)O    //Rank, Score, Id, Smiles
    2 0.37921759 18224136 NC(CO)C(=O)NC(Cc1ccc(O)cc1)C(=O)NC(CC(=O)O)C(=O)O
    3 0.20393876 18231916 NC(Cc1ccc(O)cc1)C(=O)NC(CC(=O)O)C(=O)NC(CO)C(=O)O
    ...

### Run cfm-id in a docker container ###

Assuming your current directory is in the working directory.
To identify ESI-MS/MS [M+H]+ spectra  

    sudo docker run --rm=true -v $(pwd):/root -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-id ./<spectrum_file> <id> root/<candidate_file>  10 0.001 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt Dice 1 output"

To identify ESI-MS/MS [M-H]- spectra  

    sudo docker run --rm=true -v $(pwd):/root -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-id ./<spectrum_file> <id> root/<candidate_file>  10 0.001 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt Dice 1 output"

## cfm-id-precomputed ##

Given an input spectrum and a set of candidate spectra, **cfm-id-precomputed** compares candidate spectra to the input spectrum. It returns a ranking of the candidates according to how closely they match. The spectrum prediction is done using a pre-trained CFM model.

### cfm-id-precomputed usage ###

    cfm-id-precomputed <spectrum_file> <id> <candidate_file> <num_highest> <ppm_mass_tol> <abs_mass_tol> <score_type> <output_file>

* **spectrum_file**: The filename where the input spectra can be found. This can be a .msp file in which the desired spectrum is listed under a corresponding id (next arg). Or it could be a single file with a list of peaks 'mass intensity' delimited by lines, with either 'low','med' and 'high' lines beginning spectra of different energy levels, or 'energy0', 'energy1', etc. e.g.

    65.02 40.0
    86.11 60.0
    energy1
    65.02 100.0 ... etc

* **id**: An identifier for the target molecule (Used to retrieve input spectrum from msp (if used). Otherwise not used but printed to output, in case of multiple concatenated results)

* **candidate_file**: The filename where the input list of candidate structures can be found as line separated 'id smiles_or_inchi spectrum_file' triples. i.e. each entry also specifies a file where the precomputed spectrum should be read from.

* **num_highest** (optional): The number of (ranked) candidates to return or -1 for all (if not given, returns all in ranked order).

* **ppm_mass_tol**: (optional) The mass tolerance in ppm to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs (if not given defaults to 10ppm).

* **score_type**: (optional) The type of scoring function to use when comparing spectra. Options: Dice(default), DotProduct.

* **output_file**: (optional) The filename of the output spectra file to write to (if not given, prints to stdout)

### Run cfm-id-precomputed in a docker container ###

Assuming your current directory is in the working directory.

    sudo docker run --rm=true -v $(pwd):/root -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-id-precomputed <spectrum_file> <id> <candidate_file>  10 0.001 ./output"

## cfm-annotate ##

 **cfm-annotate** annotates the peaks in a provided set of spectra given a known molecule. It computes the complete fragmentation graph for the provided molecule, and then performs inference within a CFM model to determine the reduced graph that likely occurred. Each peak in the spectrum is then assigned the ids of any fragments in that graph with corresponding mass, and these are listed in order from most likely to least likely.

The output contains the original spectra in the input format, but with fragment id values appended to any annotated peaks. Following an empty line, the reduced fragment graph is then printed.

### cfm-annotate Usage ###

    cfm-annotate <smiles_or_inchi> <spectrum_file> **<id>** <ppm_mass_tol> <abs_mass_tol> <param_file> <config_file> <output_file>

* **smiles_or_inchi**: The smiles or Inchi string for the input molecule

* **spectrum_file**: The filename where the input spectra can be found. This can be a .msp file in which the desired spectrum is listed under a corresponding id (next arg). Or it could be a single file with a list of peaks 'mass intensity' delimited by lines, with either 'low','med' and 'high' lines beginning spectra of different energy levels, or 'energy0', 'energy1', etc. e.g.

    65.02 40.0
    86.11 60.0
    energy1
    65.02 100.0 ... etc

* **id**: An identifier for the target molecule (Used to retrieve input spectrum from msp (if used). Otherwise not used but printed to output, in case of multiple concatenated results)

* **ppm_mass_tol**: (optional) The mass tolerance in ppm to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs (if not given defaults to 10ppm).

* **config_file**: (optional) The filename where the configuration parameters of the cfm model can be found (if not given, assumes param_config.txt in current directory). This needs to match the file passed to cfm-train during training. See cfm-train documentation below for further details.

* **output_file**: (optional) The filename of the output spectra file to write to (if not given, prints to stdout)

### Expected cfm-annotate output ###

    cfm-annotate Oc1ccc(CC(NC(=O)C(N)CO)C(=O)NC(CC(O)=O)C(O)=O)cc1 example_spec.txt AN_ID 10.0 0.01 none metab_ce_cfm/param_config.txt

    TARGET_ID: AN_ID
    energy0
    87.054687 4.071272337 81 82 (2.0394 2.0319)
    105.069174 0.9636028163 9 (0.9636) //Peak at 105.07 mass, explained by Fragment of id 9
    136.07616 7.037977857 86 (7.038)
    160.076289 1.197298221 80 (1.1973)
    178.084616 2.86173976
    ...
    energy2
    42.033909 1.244230912 89 (1.2442)
    60.043746 10.82864669 20 (10.829)
    70.027268 1.291256596 85 (1.2913)
    ...

    94
    0 384.1401411 NC(CO)C(=O)[NH2+]C(Cc1ccc(O)cc1)C(=O)NC(CC(=O)O)C(=O)O
    1 366.1295764 N=C(C=O)C(O)=[NH+]C(=CC1=CC=CCC1)C(=O)N=C(CC(O)O)C(O)O
    2 290.0982763 C=C([NH+]=C(O)C(=N)C=O)C(=O)NC(CC(O)O)C(O)O
    ...
    8 109.0647913 C=C1C=CC(=[OH+])CC1
    9 105.0658539 NC(CO)C(=[NH2+])O     //i.e. This fragment
    ...
    93 223.1077188 C#CC(=C=C(CO)[NH+]=C(O)C(=N)CO)CC

    0 1 O                                        //These transitions explain how each fragment fits
    0 2 O=C1C=CC=CC1                             //in to the overall graph.
       ...
    78 20 C#CC#CC(=NCO)C(O)NC(CC(O)O)C(O)O

### Run cfm-annotate in a docker container ###

Assuming your current directory is in the working directory.
To annotate ESI-MS/MS [M+H]+ spectra  

    sudo docker run --rm=true -v $(pwd):/root -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-annotate "Oc1ccc(CC(NC(=O)C(N)CO)C(=O)NC(CC(O)=O)C(O)=O)cc1" <spectrum_file> <id> 10 0.001 /trained_models_cfmid4.0/[M+H]+/param_config.txt output"

To annotate ESI-MS/MS [M-H]- spectra  

    sudo docker run --rm=true -v $(pwd):/root -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-annotate "Oc1ccc(CC(NC(=O)C(N)CO)C(=O)NC(CC(O)=O)C(O)=O)cc1" <spectrum_file> <id> 10 0.001 /trained_models_cfmid4.0/[M-H]-/param_config.txt output"

## cfm-train ##

 **cfm-train** trains the parameters for a CFM model using a list of input molecules and their corresponding spectra

 Note that if there is a lot of input data, it can take a long time to run. For this reason, it has been implemented so that it can exploit parallel processors using the MPI framework. To run on multiple processors, a version of MPI must be installed (e.g. mpich2) and the cfm-train executable should be called via mpirun or equivalent. It can also be run on a single processor, without MPI, directly from the command line, but MPI is required for compilation of the source code.

### cfm-train usage ###

Detailed Instruction for 4.0 :

    cfm-train -i <input_filename> -f <feature_filename> -c <config_filename> -p <spec_dir> -g <group> -l <status_filename> 

**input_filename** Text file with number of mols on first line, then id smiles_or_inchi cross_validation_group on each line after that.

    3
    YMHOBZXQZVXHBM-UHFFFAOYSA-N BrC1=CC(=C(C=C1OC)CCN)OC
    PHFAGYYTDLITTB-UHFFFAOYSA-N O=C1C(NC)(C2=CC=CC=C2F)CCCC1
    ORWQBKPSGDRPPA-UHFFFAOYSA-N OC1=CC=CC2=C1C(CCN(CC)C)=CN2

**feature_filename** Text file with list of feature names to include, line separated. List of options is contained in Features.cpp. e.g.

    BreakAtomPair
    BrokenOrigBondType
    GasteigerCharges
    HydrogenMovement
    NLRootMatrixFPN10
    IonRootMatrixFPN10

**config_filename** Text file listing configuration parameters. Line separated 'name value'. Options are listed in Config.cpp. At a minimum, this file must list the weight and depth of the spectrum configuration (note that the weights are not used, but need to be there anyway so set them to 1)

    lambda 0.0
    em_converge_thresh 0.01
    ga_converge_thresh 0.001
    model_depth 2
    spectrum_depth 2
    spectrum_weight 1
    spectrum_depth 2
    spectrum_weight 1
    spectrum_depth 2
    spectrum_weight 1
    allow_frag_detours 1
    num_em_restarts 3
    em_no_progress_count 2
    em_max_iterations 100
    #####config for adam#######
    ga_method 2
    starting_step_size 0.001
    ending_step_size 0.00025
    ga_adam_beta_1 0.9
    ga_adam_beta_2 0.999
    ###########################
    ga_max_iterations 10
    ionization_mode 1
    ga_sampling_method 2
    ga_minibatch_nth_size 3
    collected_all_used_idx 1
    ########NNConfig###########
    theta_function 2
    theta_nn_hlayer_num_nodes 128
    theta_nn_hlayer_num_nodes 128
    theta_nn_layer_act_func_ids 2
    theta_nn_layer_act_func_ids 2
    theta_nn_layer_act_func_ids 0
    nn_layer_dropout_probs 0.1
    nn_layer_dropout_probs 0.1
    nn_layer_dropout_probs 0

**peakfile_dir_or_msp** Input MSP file, with ID fields corresponding to id fields in input_file OR Directory containing files with spectra. Each file should be called <id>.txt, where <id> is the id specified in the input file, and contains a list of peaks 'mass intensity' on each line, with either 'low','med' and 'high' lines beginning spectra of different energy levels, or 'energy0', 'energy1', etc. e.g.</id></id>
  
    energy0
    65.02 40.0
    86.11 60.0
    energy1
    65.02 100.0 ... etc

**group** Cross validation group to run. Otherwise will assume 10 groups and run all of them.

## fraggraph-gen ##

 **fraggraph-gen** produces a complete fragmentation graph or list of feasible fragments for an input molecule. It systematically breaks bonds within the molecule and checks for valid resulting fragments.

### fraggraph-gen Usage ###

    fraggraph-gen <smiles or inchi> <max depth> <ionization mode> <fullgraph or fragonly> <output file>

* **smiles_or_inchi**: The smiles or Inchi string for the input molecule

* **ionization mode**: Whether to generate fragments using positive ESI or EI, or negative ESI ionization. **+** for positive mode **ESI [M+H]+**, **-** for negative mode **ESI [M-H]-**, **\*** for positive mode **EI [M+]**.

* **fullgraph or fragonly**: (optional) This specifies the type of output. fragonly will return a list of unique feasible fragments with their masses. fullgraph (default) will also return a list of the connections between fragments and their corresponding neutral losses.

* **output_file**: (optional) The filename of the output spectra file to write to (if not given, prints to stdout)

### Expected fraggraph-gen output ###

    fraggraph-gen CC 2 + fullgraph

    4                            //The number of fragments
    0 31.05422664 C[CH4+]        //id mass smiles  - the fragments
    1 15.02292652 [CH3+]         //id mass smiles
    ...

    0 1 C                        //from to neutral_loss - the transitions
    0 2 [HH]                     //from to neutral_loss
    2 3 [HH]                     //from to neutral_loss
    

    fraggraph-gen.exe CC 2 * fragonly

    0 30.04640161 C[CH3+]        //id mass smiles  - the fragments
    1 15.02292652 [CH3+]         //id mass smiles
    2 28.03075155 [CH2][CH2+]    //...etc
   ...

### Run fraggraph-gen in a docker container ###

Assuming your current directory is in the working directory.

    sudo docker run --rm=true -v $(pwd):/root -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; fraggraph-gen "Oc1ccc(CC(NC(=O)C(N)CO)C(=O)NC(CC(O)=O)C(O)=O)cc1"  2 + fullgraph output"

### Others ###

### Running cfm-predict in a Singularity container ###

Build or obtain your CFM-ID Docker image, and convert it to a Singularity image using instructions as given [[Compute_Canada_High-Performance_Computing#Converting_Docker_images_to_Singularity_images|here]].

Run the Singularity container from the SIF file:

  singularity exec --bind output:/out cfmid_2.0.0.1.sif cfm-predict \''CC(C)NCC(O)COC1=CC=C(CCOCC2CC2)C=C1'\' 0.001 /out/param_output0.log /out/param_config.txt 1 /out/positive/myout.txt

### CFM-ID 2 Wiki ###

For now ,please refer to <https://sourceforge.net/p/cfm-id/wiki/Home/>

# Build Docker Image #
## Dev Build ##
Go to top directory
```  docker build -t cfmid:dev -f .\docker\DevBuild_Dockerfile . ```

``` docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cfm-predict 'CC(C)NCC(O)COC1=CC=C(CCOCC2CC2)C=C1' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 ./public/[M+H]+/myout" ``` 