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
* Targets: 3394 Targets (activities over 10)
* Bioactivities: 13,918,879
* Actives:	2,089,404
* Inactives:	11,829,475

Pathway information from NCBI BioSystems

![](http://ibi.imim.es/wp-content/uploads/2012/10/DisGeNET_logo_roundedEdges.png)
![](http://www.ncbi.nlm.nih.gov/Structure/IMG/banner_graphics/biosystems_entrez3.png) ![](http://www.genome.jp/Fig/kegg128.gif) ![](http://biocyc.org/BioCyc.gif) ![](http://blog.openhelix.eu/wp-content/uploads/2011/01/Reactome_logo.jpg) ![](http://i.picresize.com/images/2015/04/29/oAE7h.png) ![](https://s-media-cache-ak0.pinimg.com/216x146/e3/71/2d/e3712dd81b80c17e24d4fb529f6bafab.jpg) ![](http://www.wikipathways.org/skins/common/images/earth-or-pathway_text3_beta.png)

Dependencies : rdkit, sklearn, numpy

![](http://www.rdkit.org/Images/logo.png) ![](http://scikit-learn.org/stable/_static/scikit-learn-logo-small.png) ![](http://upload.wikimedia.org/wikipedia/ru/c/cc/Numpylogo.png)

ChemAxon Standardizer was used for structure canonicalization and transformation, JChem 6.0.2.

![](http://www.chemaxon.com/images/powered_100px.gif)  http://www.chemaxon.com

![](https://dnasu.org/DNASU/image/Uniprot300.jpg)


All rights reserved 2014
==========================================================================================


Dependencies: 

Requires Python 2.7, Scikit-learn [1], Numpy [2] and Rdkit [3] to be installed on system.

Follow these steps on Linux:
 
1. ```Download and install Anaconda2 for Python 2.7 from https://www.continuum.io/downloads```
2. Open terminal in Mac/Linux and run ```conda install -c https://conda.anaconda.org/rdkit rdkit``` 
3. ```git clone https://github.com/lhm30/PIDGINv2/;cd PIDGINv2;curl -L -o temp.zip https://www.dropbox.com/s/1jjatrzt2gvqzo0/model.zip?dl=0;unzip temp.zip;rm temp.zip```

N.B Models are 60GB
N.B Step 3 may take up to 10 minutes

==========================================================================================


Instructions:

IMPORTANT:
*	You MUST run ```curl -L -o temp.zip https://www.dropbox.com/s/1jjatrzt2gvqzo0/model.zip?dl=0;unzip temp.zip;rm temp.zip``` in terminal/cmd prompt to download the models before first run!
*	The program recognises line-separated SMILES in .csv format
*	Molecules Should be standardized before running models
*	ChemAxon Standardizer should be used for structure canonicalization and is free for academic use at (http://www.chemaxon.com)
*	Protocol used to standardise these molecules is provided: StandMoleProt.xml
*	Do not modify the 'models', 'bg_predictions.txt' etc. names or directories 
*	cytotox_library.csv and nontoxic_background.csv are included for use as example dataset for testing

==========================================================================================


Scripts:

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
    
3. ```predict_enriched.py filename.csv N_cores tpr_threshold DisGeNET_threshold```
    This script enriched targets, NCBI Biosystems pathways and DisGeNET diseases for a library of compounds, when compared to a precomputed target predictions from a background set of 2,000,000 compounds from PubChem (bg_predictions.txt). 
    The protocol corrects for promiscuous models / biases in training data and to which targets are statistically associated with compounds in filename.csv.
    Target predictions for filename.csv are compared against PubChem predictions using the Prediction Ratio (ref. https://www.repository.cam.ac.uk/bitstream/handle/1810/246122/Liggi%20et%20al%202014%20Future%20Medicinal%20Chemistry.pdf?sequence=3), Odd's Ratio and Fishers Test p-values.
    For tables with large numbers, the (inexact) chi-square test implemented in the function chi2 test should be used. Pathways and DisGeNET predictions are compared against PubChem predictions using the Prediction Ratio (ref. https://www.repository.cam.ac.uk/bitstream/handle/1810/246122/Liggi%20et%20al%202014%20Future%20Medicinal%20Chemistry.pdf?sequence=3), Odd's Ratio and Chi-square test of independence p-values.

    'bg_predictions.txt' contains rows of target models with corresponding columns for the number of background compounds from PubChem at a given TPR threshold (to 2DP).
    'DisGeNET_diseases.txt' contains disease data used to annotate target predictions. DisGeNET assigns a gene-disease association score to give confidence for annotations, and a DisGeNET_threshold can be supplied at runtime when annotating predictions with diseases (no threshold applied by default).
        
    Example of how to run the code:

    ```
    python predict_enriched.py input.csv 4 0.5 0.25
    ```
    
    The output is a ranked list of targets that are more statistically associated with the input compounds. A low Prediction Ratio, Odd's Ratio and p-value metric indicates a higher enrichment for a target/pathway/disease when compared to the background rate
    
7. ```predict_enriched_two_libraries.py input_active_library.csv input_inactive_library.csv threshold```
    This script calculates enriched targets, NCBI BioSystems pathways and DisGeNET for two compound libraries (e.g could be phenotypically active compounds and to phenotypically inactive compounds).
    The protocol corrects for promiscuous models / biases in training data and to which targets are statistically associated with compounds in input_active_library.csv.
    Target predictions for input_active_library.csv are compared against input_inactive_library.csv predictions using the Prediction Ratio (ref. https://www.repository.cam.ac.uk/bitstream/handle/1810/246122/Liggi%20et%20al%202014%20Future%20Medicinal%20Chemistry.pdf?sequence=3), Odd's Ratio and Fishers Test p-values.
    For tables with large numbers, the (inexact) chi-square test implemented in the function chi2 test should be used. Pathways and DisGeNET predictions are compared against PubChem predictions using the Prediction Ratio (ref. https://www.repository.cam.ac.uk/bitstream/handle/1810/246122/Liggi%20et%20al%202014%20Future%20Medicinal%20Chemistry.pdf?sequence=3), Odd's Ratio and Chi-square test of independence p-values.

    Example of how to run the code:

    ```
    python predict_enriched_two_libraries.py filename_1.csv filename_2.csv N_cores threshold
    ```
    
    The output is a ranked list of targets that are more statistically associated with the input compounds. A low Prediction Ratio, Odd's Ratio and p-value metric indicates a higher enrichment for a target/pathway/disease when compared to the inactive compound set.
    
8. ```predict_fingerprints.py threshold filename.csv```
    This script calculates target and pathway hits and represents them as binary fingerprints in a matrix.
    
    Example of how to run the code:

    ```
    python predict_fingerprints.py input.csv N_cores threshold
    ```

==========================================================================================

 [1] http://scikit-learn.org/stable/
 [2] http://www.numpy.org
 [3] http://www.rdkit.org
