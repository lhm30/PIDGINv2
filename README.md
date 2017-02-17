PIDGIN Version 2: Prediction IncluDinG INactivity Version 2
===========
[![DOI](https://zenodo.org/badge/10824/lhm30/PIDGIN.svg)](http://dx.doi.org/10.5281/zenodo.15984)


Author : Lewis Mervin, lhm30@cam.ac.uk

Supervisor : Dr. A. Bender

Platt Scaled Random Forest Protein Target Prediction Tool trained on SARs from PubChem (Mined 21/06/16) and ChEMBL21

![](https://pubchem.ncbi.nlm.nih.gov/images/pubchemlogob.gif) ![](http://upload.wikimedia.org/wikipedia/commons/a/a1/Chembl_logo.png)

Molecular Descriptors : 2048bit Morgan Binary Fingerprints (Rdkit) - ECFP4

* Algorithm: Random Forest with 100 trees, class weight = 'balanced', sample weight = ratio Inactive:Active
* Scaling: Platt-scaling ref: http://scikit-learn.org/stable/auto_examples/calibration/plot_calibration_curve.html
* Targets: 3394 Targets (with activities over 10)
* Bioactivities: 13,918,879 (PubChem & ChEMBL21)
* Actives:	2,089,404
* Inactives:	11,829,475
* Pathway information from NCBI BioSystems (http://www.ncbi.nlm.nih.gov/Structure/biosystems/docs/biosystems_about.html)
* Disease information from DisGeNET (http://www.disgenet.org/web/DisGeNET/menu/help)

All rights reserved 2016

![](http://ibi.imim.es/wp-content/uploads/2012/10/DisGeNET_logo_roundedEdges.png)
![](http://www.ncbi.nlm.nih.gov/Structure/IMG/banner_graphics/biosystems_entrez3.png) ![](http://www.genome.jp/Fig/kegg128.gif) ![](http://biocyc.org/BioCyc.gif) ![](http://blog.openhelix.eu/wp-content/uploads/2011/01/Reactome_logo.jpg) ![](http://i.picresize.com/images/2015/04/29/oAE7h.png) ![](https://s-media-cache-ak0.pinimg.com/216x146/e3/71/2d/e3712dd81b80c17e24d4fb529f6bafab.jpg) ![](http://www.wikipathways.org/skins/common/images/earth-or-pathway_text3_beta.png)
![](http://www.rdkit.org/Images/logo.png)
![](http://scikit-learn.org/stable/_static/scikit-learn-logo-small.png) ![](http://upload.wikimedia.org/wikipedia/ru/c/cc/Numpylogo.png)
![](https://dnasu.org/DNASU/image/Uniprot300.jpg)
![](http://www.chemaxon.com/images/powered_100px.gif)

INSTALLATION
==========================================================================================

Follow these steps on Linux/OSX:
 
1. ```Download and install Anaconda2 for Python 2.7 from https://www.continuum.io/downloads```
2. Open terminal in Mac/Linux and run ```conda install -c https://conda.anaconda.org/rdkit rdkit``` 
3. Now run: ```conda install scikit-learn=0.17``` (PIDGINv2 uses Scikit-learn v17)
4. Navigate the directory you wish to install PIDGINv2 and in Mac/Linux terminal run ```git clone https://github.com/lhm30/PIDGINv2/``` (recommended) or download/extract the zip from GitHub webpage
5. Download and unzip models.zip into the PIDGINv2 directory from here: ```http://tinyurl.com/zhv3m5n``` (leave .pkl.zip files compressed)
6. If you would like decision tree capabilities, open terminal in Mac/Linux and run ```conda install pydot graphviz``` 

* N.B Step 5 may take up to 10 minutes


IMPORTANT
==========================================================================================

*	You MUST download the models before running!
*	The program recognises line-separated SMILES in .csv format
*	Molecules Should be standardized before running models
*	ChemAxon Standardizer should be used for structure canonicalization and is free for academic use at (http://www.chemaxon.com)
*	Protocol used to standardise these molecules is provided: StandMoleProt.xml
*	Do not modify the 'models', 'bg_predictions.txt' etc. names or directories 
*	cytotox_library.csv and nontoxic_background.csv are included for use as example dataset for testing

INSTRUCTIONS
==========================================================================================

1. ```predict_raw.py filename.csv N_cores```
    This script outputs the Platt-scaled (sigmoid) probabilities of the Random Forest classifier for the compounds in a matrix.
    
    Example of how to run the code:

    ```
    python predict_raw.py input.csv 4
    ```
	
	Output is a matrix of Platt-scaled probabilities for an input list of compounds calculated on a machine using 4 cores


2. ```predict_binary.py filename.csv N_cores tpr_threshold```
    This script generates binary predictions for the models after application of a user-specified predicted true-positive rate (TPR) threshold.
    
    The choice of required TPR applies a given confidence when binarizing predictions (i.e. an acceptable True Positive (TP) rate)
   
    Example of how to run the code:
    
    ```
    python predict_binary.py input.csv 30 0.5
    ```
    
    where 30 cores are used to produce predictions, and 0.5 would apply a 50% TPR confidence threshold
  
    
3. ```predict_enriched.py filename.csv N_cores tpr_threshold DisGeNET_threshold organism```
    This script enriched targets, NCBI Biosystems pathways and DisGeNET diseases for a library of compounds, when compared to a precomputed target predictions from a background set of 2,000,000 compounds from PubChem (bg_predictions.txt). 

    The protocol corrects for promiscuous models / biases in training data and to which targets are statistically associated with compounds in filename.csv.
    
    Target predictions for filename.csv are compared against PubChem predictions using the Prediction Ratio (ref. http://tinyurl.com/predratio), Odd's Ratio and Fishers Test p-values.
    
    For tables with large numbers, the (inexact) chi-square test implemented in the function chi2 test should be used. Pathways and DisGeNET predictions are compared against PubChem predictions using the Prediction Ratio, Odd's Ratio and Chi-square test of independence p-values.

    bg_predictions.txt contains rows of target models with corresponding columns for the number of background compounds from PubChem at a given TPR threshold (to 2DP).
    
    DisGeNET_diseases.txt contains disease data used to annotate target predictions. DisGeNET gene-disease score takes into account the number and type of sources (level of curation, organisms), and the number of publications supporting the association. The score ranges from 0 to 1 to give confidence for annotations. A DisGeNET_threshold can be supplied at runtime when annotating predictions with diseases (0.06 threshold applied by default, which includes associations from curated sources/animal models supporting them or reported in 20-200 papers). More info on the score here: http://disgenet.org/web/DisGeNET/menu/dbinfo#score 

	If using "Organism" it must be as specified in the classes_in_model.txt and enclosed by quotes ("")

    Example of how to run the code:

    ```
    python predict_enriched.py input.csv 4 0.5 0.06 "Homo sapiens (Human)"
    ```
    
    The output is a ranked list of targets that are more statistically associated with the input compounds. A low Prediction Ratio, Odd's Ratio and p-value metric indicates a higher enrichment for a target/pathway/disease when compared to the background rate
    
    
7. ```predict_enriched_two_libraries.py input_active_library.csv input_inactive_library.csv threshold organism```
    This script calculates enriched targets, NCBI BioSystems pathways and DisGeNET for two compound libraries (e.g could be phenotypically active compounds and to phenotypically inactive compounds).

    The protocol corrects for promiscuous models / biases in training data and to which targets are statistically associated with compounds in input_active_library.csv.
    
    Target predictions for input_active_library.csv are compared against input_inactive_library.csv predictions using the Prediction Ratio (ref. http://tinyurl.com/predratio), Odd's Ratio and Fishers Test p-values.
    
    For tables with large numbers, the (inexact) chi-square test implemented in the function chi2 test should be used. Pathways and DisGeNET predictions are compared against PubChem predictions using the Prediction Ratio, Odd's Ratio and Chi-square test of independence p-values.

	Organism must be as specified in the classes_in_model.txt and enclosed by quotes ("")
	
    Example of how to run the code:
	
    ```
    python predict_enriched_two_libraries.py filename_1.csv filename_2.csv 10 0.9 0.3 "Homo sapiens (Human)"
    ```
    
    The output is a ranked list of targets that are more statistically associated with the input compounds. A low Prediction Ratio, Odd's Ratio and p-value metric indicates a higher enrichment for a target/pathway/disease when compared to the inactive compound set.
    
    
8. ```predict_per_comp.py filename_1.csv N_cores threshold DisGeNET_threshold organism```
    This script calculates target, pathway and disease hits per compound and represents them in a matrix. The DisGeNET threshold and organism are optional. Organism must be as specified in the classes_in_model.txt and enclosed by quotes ("")
    
    Example of how to run the code:
	
    ```
    python predict_fingerprints.py input.csv 30 0.5 0.3 "Homo sapiens (Human)"
    ```
    
    
9. ```predict_target_fingerprints.py filename_1.csv N_cores organism```
    This script calculates target probabilities per compound in a transposed (columns are targets), simplified a matrix. These can be used as a fingerprint/descriptor for biological space. Organism filter is optional. If filtering predictions by organism, this must be as specified in the classes_in_model.txt and enclosed by quotes ("")
    
    Example of how to run the code:
	
    ```
    python predict_target_fingerprints.py input.csv 30 0.5 0.3 "Homo sapiens (Human)"
    ```
    
    
10. ```predict_enriched_two_libraries_decision_tree.py filename_1.csv filename_2.csv N_cores threshold DisGeNET_threshold organism minimum_sample_split minimum_leaf_split max_depth```
    This script calculates target, pathway and disease hits enrichment and visualises the target predictions in a decision tree (jpg file). The DisGeNET threshold and organism are optional. As always, organism must be enclosed by quotes ("")
    
    Example of how to run the code:
	
    ```
    python predict_enriched_two_libraries_decision_tree.py cytotox_library.csv nontoxic_background.csv 10 0.5 0.5 "Homo sapiens (Human) 2 2 5"
    ```
    
   
11. ```predict_enriched_decision_tree.py filename_1.csv N_cores threshold DisGeNET_threshold organism minimum_sample_split minimum_leaf_split max_depth no_kmeans_clusters```
    This script calculates target, pathway and disease hits enrichment and visualises the target predictions in a decision tree (jpg file). This code uses kmeans clustering to cluster predictions within the input dataset, as a method to split input data into hypothetical modes-of-action. The number of clusters is therefore subjective and unsupervised. The DisGeNET threshold and organism are optional. As always, organism must be enclosed by quotes ("")
    
    Example of how to run the code:
	
    ```
    python predict_enriched_two_libraries_decision_tree.py cytotox_library.csv nontoxic_background.csv 10 0.5 0.5 "Homo sapiens (Human) 2 2 5 5"
    ```

==========================================================================================
