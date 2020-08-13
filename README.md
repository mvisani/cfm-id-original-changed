# README #

## What is CFM-ID ##

CFM-ID provides a method for accurately and efficiently identifying metabolites in spectra generated by electrospray tandem mass spectrometry (ESI-MS/MS). The program uses Competitive Fragmentation Modeling to produce a probabilistic generative model for the MS/MS fragmentation process and machine learning techniques to adapt the model parameters from data.  A web service can be found at (<http://cfmid4.wishartlab.com/>)

## What CFM-ID can do ##

* **Spectra Predion**, given chemical structure. This task predicts low/10V, medium/20V, and high/40V energy ESI-MS/MS spectra for an input structure provided in SMILES or InChI format.
* **Peak Assignment**: Annotating the peaks in set of spectra given a known chemical structure. This task takes a set of three input spectra (for ESI spectra, low/10V, medium/20V, and high/40V energy levels) in peak list format and a chemical structure in SMILES or InChI format, then assigns a putative fragment annotation to the peaks in each spectrum.
* **Compound Identification**: Predicted ranking of possible candidate structures for a target spectrum. This task takes a set of three input spectra (for ESI spectra, low/10V, medium/20V, and high/40V energy levels) in peak list format, and ranks a list of candidate structures according to how well they match the input spectra. This candidate list need to be provided by the user. Chemical classes are predicted for each candidate molecule. The original similarity score used in the ranking was computed (Dice or DotProduct) by comparing the predicted spectra of a candidate compound with the input spectra.

## What is this repository for? ##

* This is a repository for CFM-ID 4 MSML (mass spectrum machine learning)
* CFM-ID MSRB (mass spectrum rule based) can be found at (<https://bitbucket.org/wishartlab/msrb-fragmenter/src/master/>)
* Original source code from CFM-ID 2.0 Can be found at (<https://sourceforge.net/p/cfm-id/wiki/Home/>)

## What include in this tool set for? ##

* cfm-predict
* cfm-id
* cfm-id-precomputed
* cfm-annotate
* cfm-train
* fraggraph-gen
  
## How do I get set up? ##

### Install from source code ###

* Please check INSTALL FILE
* Note Only Insatll on linux and Mac has been verified, while install on Windows from source code is possible 

### Use Pre Build Docker ###

    docker push wishartlab/cfmid:latest


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

* **output_file_or_dir**: (optional) The filename of the output spectra file to write to (if not given, prints to stdout). In case of batch mode using file input above, this is used to specify the name of a directory where the output files (<id>.log) will be written (if not given, uses current directory), OR an msp or mgf file.</id>

* **apply_postproc**: (optional) Whether or not to post-process predicted spectra to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks (whichever comes first) (0 = OFF, 1 = ON (default)). If turned off, will output a peak for every possible fragment of the input molecule, as long as the prob_thresh argument above is set to 0.0.

* **suppress_exceptions**: (optional) Suppress most exceptions so that the program returns normally even when it fails to produce a result (0 = OFF (default), 1 = ON).

### Expected cfm-predict output ###

    cfm-predict CCCC 0.001 metab_ce_cfm/param_output0.log metab_ce_cfm/param_config.txt

    energy0
    15.02292652 0.03877135094
    27.02292652 0.0004516638069
    29.03857658 0.1823948415
    31.05422664 0.1285812238
    43.05422664 0.54978044
    59.08552677 99.10002048
    energy1
    15.02292652 0.2014284022
    27.02292652 0.006349994177
    29.03857658 0.9254281728
    31.05422664 0.3201026529
    43.05422664 1.755347781
    59.08552677 96.791343
    energy2
    15.02292652 0.6774027078
    27.02292652 0.2170199999
    29.03857658 2.333980325
    31.05422664 0.9058884643
    43.05422664 27.56483288
    59.08552677 68.30087562

### Run cfm-predict in a docker container ###

Assuming your home directory is `/home/ubuntu/`：
To predict ESI-MS/MS [M+H]+ spectra  

    sudo docker run --rm=true -v /home/ubuntu/cfm_id/cfmid/output:/root -i cfmid:latest sh -c "cd /root/; cfm-predict 'CC(C)NCC(O)COC1=CC=C(CCOCC2CC2)C=C1' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /cfmid_trained_models/[M+H]+/param_config.txt 1 /root/[M+H]+/myout"

To predict ESI-MS/MS [M-H]- spectra  

    sudo docker run --rm=true -v /home/ubuntu/cfm_id/cfmid/output:/root -i cfmid:latest sh -c "cd /root/; cfm-predict 'CC(C)NCC(O)COC1=CC=C(CCOCC2CC2)C=C1' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /cfmid_trained_models/[M-H]-/param_config.txt 1 /root/[M-H]-/myout"

## cfm-id ##

Given an input spectrum and a list of candidate smiles (or inchi) strings, **cfm-id** computes the predicted spectrum for each candidate and compares it to the input spectrum. It returns a ranking of the candidates according to how closely they match. The spectrum prediction is done using a pre-trained CFM model.

### cfm-id usage ###

    cfm-id <spectrum_file> <id> <candidate_file> <num_highest> <ppm_mass_tol> <abs_mass_tol> <prob_thresh> <param_file> <config_file> <score_type> <apply_postprocessing> <output_file> <output_msp_or_mgf>

* **spectrum_file**: The filename where the input spectra can be found. This can be a .msp file in which the desired spectrum is listed under a corresponding id (next arg). Or it could be a single file with a list of peaks 'mass intensity' delimited by lines, with either 'low','med' and 'high' lines beginning spectra of different energy levels, or 'energy0', 'energy1', etc. e.g.

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
    4 0.16378009 59444507 Cc1cc(CN(CC(=O)O)CC(=O)O)nc(CN(CC(=O)O)CC(=O)O)c1
    5 0.13664102 18219720 NC(CC(=O)O)C(=O)NC(Cc1ccc(O)cc1)C(=O)NC(CO)C(=O)O

### Run cfm-id in a docker container ###

Assuming your home directory is `/home/ubuntu/`：
To identify ESI-MS/MS [M+H]+ spectra  

    sudo docker run --rm=true -v /home/ubuntu/cfm_id/cfmid/output:/root -i cfmid:latest sh -c "cd /root/; cfm-id /root/<spectrum_file> <id> root/<candidate_file>  10 0.001 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /cfmid_trained_models/[M+H]+/param_config.txt Dice 1 /root/myout"

To identify ESI-MS/MS [M-H]- spectra  

    sudo docker run --rm=true -v /home/ubuntu/cfm_id/cfmid/output:/root -i cfmid:latest sh -c "cd /root/; cfm-id /root/<spectrum_file> <id> root/<candidate_file>  10 0.001 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /cfmid_trained_models/[M-H]-/param_config.txt Dice 1 /root/myout"


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

Assuming your home directory is `/home/ubuntu/`：

    sudo docker run --rm=true -v /home/ubuntu/cfm_id/cfmid/output:/root -i cfmid:latest sh -c "cd /root/; cfm-id-precomputed /root/<spectrum_file> <id> root/<candidate_file>  10 0.001 /root/myout"

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
    178.084616 2.861739768
    223.106608 53.80100032 32 93 92 (36.58 9.2055 8.0159)
    251.10173 21.90932756 38 87 88 91 (9.7098 5.3256 5.1959 1.678)
    297.107567 2.122976713 16 90 (1.2417 0.88124)  
    384.140384 6.034804405 0 (6.0348)
    ...
    energy2
    42.033909 1.244230912 89 (1.2442)
    60.043746 10.82864669 20 (10.829)
    70.027268 1.291256596 85 (1.2913)
    87.056272 7.489320919 81 82 (3.7496 3.7398)
    91.054494 9.60202642 79 (9.602)
    119.04828 6.415043123 83 (6.415)
    121.063402 2.97004533 84 (2.97)
    133.06551 2.057893243
    135.066238 1.480563861
    136.074907 40.8315392 86 (40.832)
    160.074409 10.80320864 80 (10.803)
    178.085454 4.986225072

    94
    0 384.1401411 NC(CO)C(=O)[NH2+]C(Cc1ccc(O)cc1)C(=O)NC(CC(=O)O)C(=O)O
    1 366.1295764 N=C(C=O)C(O)=[NH+]C(=CC1=CC=CCC1)C(=O)N=C(CC(O)O)C(O)O
    2 290.0982763 C=C([NH+]=C(O)C(=N)C=O)C(=O)NC(CC(O)O)C(O)O
    3 288.0826262 C=C([NH+]=C(O)C(=N)C=O)C(=O)N=C(CC(O)O)C(O)O
    4 278.0982763 N=C(C=O)C(=O)[NH+]=CC(O)NC(CC(O)O)C(O)O
    5 276.0826262 N=C(C=O)C(=O)[NH+]=C=C(O)NC(CC(O)O)C(O)O
    6 274.0669761 N=C(C=O)C(=O)[NH+]=C=C(O)N=C(CC(O)O)C(O)O
    7 272.0513261 N=C(C=O)C(=O)[NH+]=C=C(O)N=C(C=C(O)O)C(O)O
    8 109.0647913 C=C1C=CC(=[OH+])CC1
    9 105.0658539 NC(CO)C(=[NH2+])O     //i.e. This fragment
    ...
    89 42.03382555 CC#[NH+]
    90 297.1081127 C#CC(=C=C([NH+]=C=O)C(=O)NC(CC(O)O)C(O)O)CC
    91 251.1026334 NC(O)C(=C=C1C=CC(=O)CC1)[NH+]=C(O)C=CO
    92 223.1077188 NCC(O)=[NH+]C(=C=C1C=CC(=O)CC1)CO
    93 223.1077188 C#CC(=C=C(CO)[NH+]=C(O)C(=N)CO)CC

    0 1 O                                        //These transitions explain how each fragment fits
    0 2 O=C1C=CC=CC1                             //in to the overall graph.
    0 3 O=C1C=CCCC1
    0 4 C=C1C=CC(=O)C=C1
    0 5 C=C1C=CC(=O)CC1
    0 6 CC1C=CC(=O)CC1
    0 7 CC1CCC(=O)CC1
    0 8 N=C(C=O)C(=O)N=C=C(O)NC(CC(O)O)C(O)O
    0 9 O=C1C=CC(=C=C=C(O)N=C(C=C(O)O)C(O)O)CC1   
    0 2 C=C1C=CC(=O)CC1                        
    0 3 CC1C=CC(=O)CC1
    ...
    78 2 C=CC
    78 3 CCC
    78 4 C#CCC
    78 5 C=CCC
    78 6 CCCC
    78 9 C#CC#CC=C(O)NC(CC(O)O)C(O)O
    78 13 C#CC#CC(N)C(O)NC(CC(O)O)C(O)O
    78 20 C#CC#CC(=NCO)C(O)NC(CC(O)O)C(O)O

### Run cfm-annotate in a docker container ###

Assuming your home directory is `/home/ubuntu/`：
To annotate ESI-MS/MS [M+H]+ spectra  

    sudo docker run --rm=true -v /home/ubuntu/cfm_id/cfmid/output:/root -i cfmid:latest sh -c "cd /root/; cfm-annotate "Oc1ccc(CC(NC(=O)C(N)CO)C(=O)NC(CC(O)=O)C(O)=O)cc1" /root/<spectrum_file> <id> 10 0.001 /cfmid_trained_models/[M+H]+/param_config.txt /root/myout"

To annotate ESI-MS/MS [M-H]- spectra  

    sudo docker run --rm=true -v /home/ubuntu/cfm_id/cfmid/output:/root -i cfmid:latest sh -c "cd /root/; cfm-annotate "Oc1ccc(CC(NC(=O)C(N)CO)C(=O)NC(CC(O)=O)C(O)=O)cc1" /root/<spectrum_file> <id> 10 0.001 /cfmid_trained_models/[M-H]-/param_config.txt /root/myout"

## cfm-train ##

 **cfm-train** trains the parameters for a CFM model using a list of input molecules and their corresponding spectra

### cfm-train usage ###

Detailed Instruction for 4.0 : TBD.

For mow, please check: 

    cfm-train --help

And <https://sourceforge.net/p/cfm-id/code/HEAD/tree/supplementary_material/cfm-train_example/>. 

### Run cfm-train in a docker container ###

For now, training via docker is not recommended

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
    2 29.03912516 C=[CH3+]       //id mass smiles
    3 27.02292652 C#[CH2+]       //id mass smiles

    0 1 C                        //from to neutral_loss - the transitions
    0 2 [HH]                     //from to neutral_loss
    2 3 [HH]                     //from to neutral_loss
    

    fraggraph-gen.exe CC 2 * fragonly

    0 30.04640161 C[CH3+]        //id mass smiles  - the fragments
    1 15.02292652 [CH3+]         //id mass smiles
    2 28.03075155 [CH2][CH2+]    //...etc
    3 26.01510148 [CH]=[CH+]
    4 27.02292652 C#[CH2+]
    5 29.03857658 C=[CH3+]

### Run fraggraph-gen in a docker container ###

Assuming your home directory is `/home/ubuntu/`：

    sudo docker run --rm=true -v /home/ubuntu/cfm_id/cfmid/output:/root -i cfmid:latest sh -c "cd /root/; fraggraph-gen "Oc1ccc(CC(NC(=O)C(N)CO)C(=O)NC(CC(O)=O)C(O)=O)cc1"  2 + fullgraph /root/myout"

### Others ###

### Running cfm-predict in a Singularity container ###

Build or obtain your CFM-ID Docker image, and convert it to a Singularity image using instructions as given [[Compute_Canada_High-Performance_Computing#Converting_Docker_images_to_Singularity_images|here]].

Run the Singularity container from the SIF file:

  singularity exec --bind output:/out cfmid_2.0.0.1.sif cfm-predict \''CC(C)NCC(O)COC1=CC=C(CCOCC2CC2)C=C1'\' 0.001 /out/param_output0.log /out/param_config.txt 1 /out/positive/myout.txt 

### CFM-ID 2 Wiki ###

For now ,please refer to https://sourceforge.net/p/cfm-id/wiki/Home/