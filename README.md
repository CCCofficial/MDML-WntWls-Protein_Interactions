# <center>Machine Learning Characterization of Wnt-WLS Protein Binding Dynamics</center>
***  

## Objective
***
The Wnt protein family is essential for cell development, with each Wnt protein interacting uniquely with the WLS membrane protein through varying binding residues. This study employs molecular dynamics (MD) simulations and supervised machine learning (ML) to analyze the binding differences among four Wnt proteins. Using both crystal structures of Wnt3a and Wnt8a and homology models of Wnt1 and Wnt5a, the MD simulations were run for 1.5 microseconds producing 75,000 frames per protein. Features were extracted as Wnt-WLS residue pairs within a 12Å contact distance (n=985) and autocorrelation functions were used to define uncorrelated training (n=160,000) and testing (n=70,000) sets. Subclustering within regions of known significance reduced the feature set to 210 residue pairs. Using these data we trained a Random Forest multiclass classification model, optimized through hyperparameter tuning and 10-fold cross-validation, which achieved a test accuracy of 95.9%. Key residue pairs distinguishing each WNT system were then identified using permutation feature importance. This methodology not only highlights crucial contact pairs but also corroborates previously studied residues influencing enzyme activity, demonstrating the potential of ML in understanding complex protein interactions.

<br>

<br>

This repository provides tools for interacting with and processing molecular structures and data from molecular dynamic simulations. It will also include code and workflows for building and training a variety of robust machine learning models.

<br>

***
## Code
***
This section provides information on the code used to conduct the study described above. For access to data, please send an email to the authors (see the **Contact** section below).

### Environment Set-Up  
The code has been developed for `Python 3.12`. We recommend using a virtual Python envrionment, please see [this](https://docs.python.org/3/library/venv.html) documentation for details on how do this using built-in Python modules. Once the environment has been activated, please install the required Python modules by running the following from the terminal:
```bash
pip install -r requirements.txt
```
<br> 

### Scripts
The `scripts/` directory is organizes as follows:

```md
scripts/
  |-- main.ipynb
  |-- 00_map/
  |-- 00_sequence_align/
  |-- 00_visualize_contacts/
  |-- 01_get_contact_matrix/
  |-- 02_apply_threshold/
  |-- 03_finalize_dataset/
  |-- 04_map/
  |-- 05_combine/
  |-- 06_acf/
  |-- 07_preprocess/
  |-- 08_model_build/
  |-- 09_model_eval/
```

The code within the directories can be broken into two parts:  
1. **MD Data Preprocessing:** Includes all subdirectories starting with `scripts/00_*`.  
2. **ML Pipeline:** Includes code to run all subdirectories `scripts/01_*` - `scripts/05_*` as well as includes instructions for how to run the code in subdirectories `scripts/06_*` - `scripts/09_*`.

<br> 

#### MD Data Preprocessing
##### 00_map  
The code in this directory generates residue index mappings from Wnt1, wnt3a, and Wnt5a to Wnt8. To generate index mapping, start by running ```vmd -e rmsd_copy.tcl``` (vmd script that reads in the pdb files from ```input``` and generates the RMSD alignment for each possible indices - found at ```output``` as ```.out``` files).

Once the ```.out``` files are generated, run the Jupyter notebook ```mapping.ipynb``` which will generate ```*8a.csv``` files. These are three column files with the following format:
- Column 1: Wnt index to be aligned
- Column 2: Wnt8a that best fits to column 1
- Column 3: The RMSD of the fit

```
| WNT1.       | WNT8a.      | RMSD_min. |
| ----------- | ----------- |-----------| 
|     61      |     22      |   1.76    |
|     62      |     23      |   1.72    |
```

#### 00_visualize_contacts
Three scripts required to generate the following figures:  
- **Feature Importance Bonds:** ```draw_feature_importance_bonds.tcl``` from -```../07_preprocess/output/feature_importances_region_*.txt```  
  - Figure: `plots/draw_final`  
- **Initial Bonds**: ```draw_initial_bonds.tcl``` from ```output/WNT1_threshhold12_labels.txt``` (this is the same as ```03_finalize_dataset/output``` with `AA` removed)  
  - Figures: `plots/initial_wnt1`, `plots/initial_wnt3a`, `plots/initial_wnt5a`, `plots/initial_wnt8a`
- **ML Input Bonds**:```draw_MLInput_bonds.tcl``` -> Draw the 210 features from ```../07_preprocess/output/region_*.txt```  
  -Figure: `plots/MLInput_regions`  

<br> 

#### ML Pipeline
The ```main.ipynb``` includes code to run all subdirectories `scripts/01_*` - `scripts/05_*` as well as includes instructions for how to run the code in subdirectories `scripts/06_*` - `scripts/09_*`. See the notebook for specific details on the code included within each subdirectory. The first five steps of the pipeline include:  
1. Obtaining the residue pair contact distance matrix (`01_get_contact_matrix`)  
2. Apply an Ansgtrom-based threshold to each Wnt such that only contact pairs within this threshold within at least 1 frame are included  (`02_apply_threshold`)   
3. Finalize the contact matrices for each Wnt by adding the contact pair labels  (`03_finalize_dataset`)  
4. Normalize each Wnt's contact pair labels by mapping them to Wnt8a (reference Wnt) (`04_map`)  
5. Merge each Wnt's contact matrices to create a single input for ML (`05_combine`).  
  - This directory also includes a Jupyter Notebook (`merge.ipynb`) that generates a Venn Diagram of contact pairs between all Wnts.  

Each folder has a similar structure, for example:
- ```plots```: Plots from ```output/```  
- ```output```: Output of ```run.py```

Additional information for `scripts/06_*` - `scripts/09_*` is provided below.  

<br>

##### 06_acf  
This directory contains code that applies autocorrelation functions (ACF) to identify the optimal location (i.e., frame) within each Wnt's trajectory to use for deriving training and testing splits. There are two scripts:  
- ```run_acf.py```: runs the autocorrelation function for each of the 980 features individually. Store the result in ``output```
- ```plot_acf.py```: takes the ```output/*.npy```, averages over each of the features, and plots the average ACF with the frames cut-points. 

##### 07_preprocess
This directory contains code for performing feature selection, which includes two stages of clustering. The first stage clusters contact pairs by known Wnt8a regions. The second stage, performs subclustering within the stage 1 clusters, which are then randomly sampled to identify a reduced set of contact pairs.  
- ```preprocess.ipynb```: runs the preprocessing steps explained above.
- ```model/```: Contains the train + test sets
- ```output/```: 
	- ```clusterbyindex_df.pkl```: Data structure containing each of the sub clusters from the 22 initial clusters & the elements within
	- ```subcluster_*.txt```: Contact pairs within the 22 clusters
	- ```feature_importances_region_*.txt```: Important features within each region
- ```plots/```: All of the clustering plots
  - Wnt8a Region Colors
      - Red = Hairpin 1
      - Green = PAM
      - Blue = Hairpin 2
      - Purple = Hairpin 3
      - Gray = N-Term

##### 08_model_build 
This directory contains code to build the Random Forest classifier, including tuning hyperparameter via grid search with 10-fold cross-validation. It also contains code to determine the final model performance.
```model_train.py```: Grid search over parameter grid with K-Fold
```model/```: The GridSearch result from the three iterations
	- ```gridsearch_final_stratified_0.pkl```: 1st iteration of grid search
	- ```gridsearch_final_stratified_1.pkl```: 2nd iteration of grid search
	- ```gridsearch_final_stratified_2.pkl```: 3rd iteration of grid search

##### 09_model_eval 
This directory contains code to generate a learning curve and performs permutation feature importance. 
- ```compute_and_plot_learning_curve.py```: Get the learning curve with optimized parameters
- ```compute_feature_importances.ipynb```: Compute permutation feature importance  
- ```model/final_model.pkl```: Final model from learning curve(Depth = 3; # of Trees = 10)
- ```output/```: 
  - Learning Curve Output: `train_scores_finalmodel.txt` and `test_scores_finalmodel.txt` 
  - Permutation Feature Importance Output: `feature_importances.txt`
- ```plots/```: learning curve and feature importance plots

<br><br>

***
## 📨 Contact 
***  
This project is a collaboration between researchers at IBM Research, Arizona State University, and University of Illinois Urbana-Champaign. We'd love to hear from you!
To get in touch with us, please [create an issue](https://github.ibm.com/IBM-Research-AI/MDML-WntWLS/issues) or <a href="mailto:sara.capponi@ibm.com, sara.capponi@ibm.com">send us an email</a>.

[//]: # (| <img src="https://media.github.ibm.com/user/430879/files/9e327004-2b50-4a84-b638-ac5c584fb33d" width=150> | <img src="https://media.github.ibm.com/user/430879/files/b05f557f-e590-47be-b52f-12eae38a8d23" width=150> | <img src="https://media.github.ibm.com/user/430879/files/5f1e5a7d-107d-4d05-b104-1013dcf18bc5" width=150> |)

[//]: # (|:---------------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------------------------------------:|)

[//]: # (|        [**Kevin Cheng**]&#40;https://www.linkedin.com/in/kevin-jose-cheng1&#41;<br>*UIUC*<br>*IBM Intern*         |             [**Tiffany Callahan**]&#40;https://research.ibm.com/people/tiffany-callahan&#41;<br>*IBM*             |                 [**Sara Capponi**]&#40;https://research.ibm.com/people/sara-capponi&#41;<br>*IBM*                 |)

<br>

***  
## ✨ Acknowledgements and Disclaimers
*** 
This material is based upon work supported by the National Science Foundation (NSF) under Grant `No.DBI-1548297` awarded to the [Center for Cellular Construction](https://centerforcellularconstruction.org/about/) at IBM Research.  
Any opinions, findings and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation or IBM Research.

<img src="https://media.github.ibm.com/user/430879/files/fb31b0d4-849a-4d06-8743-a3571bae867f" align="left" width=195>
<img src="https://media.github.ibm.com/user/430879/files/a1e5f43d-9497-494c-bb2a-32d9c65bc589" align="left" width=75>
