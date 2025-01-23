# <center>Machine Learning Characterization of Wnt-WLS Protein Binding Dynamics</center>
***  

## 🎯 Objective
***
<img src="https://media.github.ibm.com/user/430879/files/a882b7f7-73fc-4990-97ce-e9c98719b3e8" align="right" width=300>
The Wnt protein family is essential for cell development, with each Wnt protein interacting uniquely with the WLS membrane protein through varying binding residues. This study employs molecular dynamics (MD) simulations and supervised machine learning (ML) to analyze the binding differences among four Wnt proteins. Using both crystal structures of Wnt3a and Wnt8a and homology models of Wnt1 and Wnt5a, the MD simulations were run for 1.5 microseconds producing 75,000 frames per protein. Features were extracted as Wnt-WLS residue pairs (n=985) and autocorrelation functions were used to define uncorrelated training (n=160,000) and testing (n=70,000) sets. Subclustering within regions of known significance reduced the feature set to 210 residue pairs. Using these data, we trained a Random Forest multiclass classification model, optimized through hyperparameter tuning and 10-fold cross-validation, which achieved a test accuracy of 95.9%. Key residue pairs distinguishing each WNT system were then identified using permutation feature importance. This methodology not only highlights crucial contact pairs but also corroborates previously studied residues influencing enzyme activity, demonstrating the potential of ML in understanding complex protein interactions.
<br><br>
This repository provides tools for interacting with and processing molecular structures and data from molecular dynamic simulations. It will also include code and workflows for building and training a variety of robust machine learning models.

<br><br><br>

***
## 💻 Code
***
This section provides information on the code used to conduct the study described above. For access to data, please email the authors (see the **Contact** section below).

### 🌎 Environment Set-Up  
***
The code has been developed for `Python 3.12`.
We recommend using a virtual Python environment,
please see [this](https://docs.python.org/3/library/venv.html) documentation for details
on how to do this using built-in Python modules.
Once the environment has been activated,
install the required Python modules by running the following from the terminal:
```bash
pip install -r requirements.txt
```

<br>

### 📄 Scripts  
***
The `scripts/` directory is organized as follows:

```md
scripts/
  |-- 00_sequence_align/	(not required for machine learning pipeline)
  |-- 00_unwrap_trajectory/	(not required for machine learning pipeline)
  |-- main.ipynb
  |-- 00_map/
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

<br><br>

#### ⚛️ Molecular Dynamics Data Preprocessing  
***

<br>

#### 📁 00_sequence_align
Run the Jupyter Notebook `scripts/00_sequence_align/fasta.ipynb`, which uses the `.pdb` files in `scripts/00_map/input`
to generate a `wnt.fasta` file of the aligned sequences for each Wnt.
Not required to run the Machine Learning Pipeline.

#### 📁 00_unwrap_trajectory
Run the Jupyter Notebook `scripts/00_unwrap_trajectory/unwrap.ipynb`,
which uses the raw `.dcd` files (along with other files written to `scripts/00_unwrap_trajectory/output` --
see next sentence) to unwrap and concat the trajectories for visualization.
Note that prior to running this Notebook,
the `scripts/00_unwrap_trajectory/read_raw_dcd.tcl` and `scripts/00_unwrap_trajectory/reload.tcl` files must be run
as they generate the required input files (written to `scripts/00_unwrap_trajectory/output`).
Not required to run the Machine Learning Pipeline.

<br><br>

#### 📁 00_map
The code in this directory generates residue index mappings from Wnt1, wnt3a, and Wnt5a to Wnt8 (reference Wnt). Only those residues on *Wnt* need alignment.

**STEP 1**.
Generate mappings by running ```vmd -e 01_rmsd.tcl```, a vmd script that reads pdb files from ```scripts/00_map/input```
and generates an RMSD alignment for each possible index on the Wnt residues.
Results are written as ```.out```  to the ```scripts/00_map/output``` directory.

**STEP 2**.
Once the ```.out``` files are generated,
run the Jupyter notebook ```scripts/00_map/02_generate_mappings.ipynb```,
which takes in the ```.out``` files and computes the best RMSD fits between WNT1/3a/5a and WnT8a using the following steps:
- To identify structural correspondences between WNT proteins (WNT1, WNT3a, WNT5a, and WNT8a), 41-amino-acid segments (`[i-20, i+20]`) were analyzed around each residue `i`. This range captures entire secondary structures (e.g., alpha-helices and beta-sheets) while minimizing overlaps. Corresponding segments in WNT8a were identified based on the lowest root-mean-square deviation (RMSD).
- A two-step RMSD calculation was used:
  1. Global alignment including the conserved WLS protein.
  2. Local alignment of the 41-amino-acid segments.

A cutoff RMSD of 10 Å was applied to exclude poorly fitting matches, ensuring high-quality residue correspondences. For rmsd values greater than the cutoff, the Wnt8a mapping idx is set to -1 (indicating no fit).  

Results are written to ```*8a.csv``` files in the ```scripts/00_map/output``` directory. Each ```.csv```  file has the following format:
- Column 1: Wnt index to be aligned
- Column 2: Wnt8a that best fits to column 1
- Column 3: The RMSD of the fit

```
| WNT1.       | WNT8a.      | RMSD_min. |
| ----------- | ----------- |-----------| 
|     61      |     22      |   1.76    |
|     62      |     23      |   1.72    |
```

**STEP 3**. To enable easier access to the mapping, run ```scripts/00_map/03_generate_map_dict.ipynb```, which uses the ```.pdb``` files in ```scripts/00_map/input``` and the ```*8a.csv``` files in the ```scripts/00_map/output``` to generate a dictionary which provides the residue-level mapping between the reference Wnt (Wnt8a) and the other wnts (Wnt1/3a/5). The dictionary is saved as a JSON file (```scripts/00_map/output/wnt8a_residue_alignment_map.json```). See below for a high-level overview of the dictionary structure.  

```python
{
    'wnt1': {
        'original_residues': {
            'residue_ids': ['32', '71', '180', '242', ..., '369'],
            'residue_labels': ['GLY32', 'LEU71','PHE180', 'VAL242', ..., 'CYS369']},
        'wnt8a_alignment': {'32': '-1', '71': '32', '180': '142', '242': '204', ...}
    },
    'wnt3a': {
        'original_residues':  {
            'residue_ids': ['19', '57', '127', '229', ..., '352'],
            'residue_labels': ['SER19', 'ARG57', 'SER127', 'ASP229', ..., 'LYS352']},
        'wnt8a_alignment': {'19': '-1', '57': '34', '127': '104', '229': '208', ...} 
        },
    'wnt5a': {
        'original_residues': {
            'residue_ids': ['44', '61', '169', '211', ..., '380'],
            'residue_labels': ['ASN44', 'LEU61', 'ARG169', 'TYR211', ..., 'LYS380']},
        'wnt8a_alignment': {'44': '-1', '61': '22', '169': '119', '211': '153', ...} 
    },
    'wnt8a': {
        'original_residues': {
            'residue_ids': ['22', '72', '107', '269', ..., '337'],
            'residue_labels': ['ALA22', 'LEU72', 'GLY107', 'GLY269', ..., 'ALA337']},
        'wnt8a_alignment': {'22': '22', '72': '72', '107': '107', '269': '269', ...} 
    }
}
```

<br>

#### 📁 00_visualize_contacts
The code in this directory is designed
to generate different figures using [VMD (Visual Molecular Dynamics)](https://www.ks.uiuc.edu/Research/vmd/).
Note that it requires running code from the **Machine Learning Pipeline** section.
The following three scripts (`scripts/00_visualize_contacts/`)
are required to generate the visualizations described below (plots written to `scripts/00_visualize_contacts/plots/`):  
- **Feature Importance Bonds:** `draw_feature_importance_bonds.tcl`, which requires data from `scripts/07_preprocess/output/feature_importances_region_*.txt`  
  - Generated Figure: `draw_final`  
- **Initial Bonds**: `draw_initial_bonds.tcl`, which requires data from `scripts/00_visualize_contacts/output/WNT1_threshhold12_labels.txt` (this is the same as `scripts/03_finalize_dataset/output` with `AA` removed).  
  - Generated Figures:  
    - `initial_wnt1`  
    - `initial_wnt3a`  
    - `initial_wnt5a`  
    - `initial_wnt8a`
- **ML Input Bonds**:`draw_MLInput_bonds.tcl`, which draws the final feature set generated from `scripts/07_preprocess/output/region_*.txt`  
  - Generated Figure: `MLInput_regions`  

<br><br> 

#### 🧠 Machine Learning Pipeline  
***
The `scripts/main.ipynb` Jupyter Notebook includes the code to run all subdirectories `scripts/01_*`-`scripts/05_*`
as well as includes instructions for how to run the code in subdirectories `scripts/06_*`-`scripts/09_*`.

Important Notes:  
- The code for this section assumes access to the raw MD simulation files (`dcd` and `psf` files)
and assumes that the code in the **Molecular Dynamics Data Preprocessing** section has been run.  
- The scripts in each directory listed within this section leverage helper scripts located in `scripts/util`.  
- Each directory in this section has a similar structure, for example:
  - `plots`: Plots from `output/`  
  - `output`: Output of `run.py`

<br>

##### 📁 01_get_contact_matrix 
The first step is to obtain the residue pair contact distance matrix,
which is generated for *every frame* of the trajectory.

The Notebook runs the `scripts/01_get_contact_matrix/run.py` script in this directory,
which for each trajectory, computes the contact distance matrix between all C-α atoms.
Reshape the matrix from (`nFrames x nCalpha x nCAlpha`) to have dimensions (`nFrames x nCAlpha*nCAlpha`). 

**Inputs:** requires access to the raw `dcd` and `psf` files for each system.

**Outputs:** Results are output to the `scripts/01_get_contact_matrix/output` directory.
For each trajectory, a `npy` file of the format `WNT{wnt_protein_name}_distances_iter{chunk_num}.npy` is generated.

<br>

##### 📁 02_apply_threshold 
The second step is to apply an Angstrom-based threshold to each Wnt such that
only contact pairs within this threshold (`12Å`) within at least one frame are included.

The Notebook runs the `scripts/02_apply_threshold/run.py` script in this directory,
which iterates over the second half of each trajectory and records the indices
corresponding to atom pairs that are within the cutoff at least <u>ONCE</u>.

**Inputs:** requires the `npy` files generated from the prior step which have been written to `01_get_contact_matrix/output`.

**Outputs:** Results are output to the `scripts/02_apply_threshold/output` directory.
For each trajectory, a `txt` file in the form `output/WNT{wnt_protein_name}_idx_thresh{threshhold}.txt` is generated. This file contains the column indices of contact pairs within the `12Å` threshold.

<br>

##### 📁 03_finalize_dataset 
The third step finalizes the contact matrices for each Wnt by adding the contact pair labels.

The Notebook runs the `scripts/03_finalize_dataset/run.py` script in this directory,
which given the indices of the C-α atom pairs that are within the distance threshold,
parses the original dataset generated from step 1 (`scripts/01_get_contact_matrix`) - keeping only the atom pairs within the distance threshold.

**Inputs:** requires the `npy` files from `scripts/01_get_contact_matrix/output` and the `txt` files from `scripts/02_apply_threshhold/output`.

**Outputs:** Results are output to the `scripts/03_finalize_dataset/output` directory.
1. `txt` file of the form `scripts/03_finalize_dataset/output/WNT{wnt_protein_name}_threshhold{threshhold}_labels.txt`
   (a list of column names from the passed matrix).
2. `npy` file of the form `scripts/03_finalize_dataset/output/WNT{wnt_protein_name}_threshhold{threshhold}_matrix.npy` (the matrix containing the contact pairs within the `12Å` threshold).  

<br>

##### 📁 04_map 
The fourth step normalizes each Wnt's contact pair labels by mapping them to Wnt8a (reference Wnt).

The Notebook runs the `scripts/04_map/run.py` script in this directory.

**Inputs:** The maps at `scripts/00_map/output/{wnt_from}_to_{wnt_to}.csv` (where indices are mapped from wnt_from to wnt_to). This structure is a 3-column csv. First column contains the indices of `wnt_from` and the second column contains the indices of `wnt_to` that is the best map.

**Outputs:** New labels at `scripts/04_map/output/WNT{SYSTEM}_threshhold{distance_threshhold}_labels.txt`. These essentially replace the labels located in `scripts/03_finalize_dataset/output`.

<br>

##### 📁 05_combine 
The fifth step merges each Wnt's contact matrices to create a single input for machine learning
(only contact pairs identified in all Wnts are kept; n=`985`).
This directory also includes a Jupyter Notebook (`scripts/05_combine/merge.ipynb`)
that generates a Venn Diagram of contact pairs between all Wnts.  

The Notebook runs the `scripts/05_combine/run.py` script in this directory,

**Inputs:**  
- the `txt` labels generated in `scripts/04_map/output`
- the `npy` matrix containing the contact pairs within the `12Å` threshold generated in `scripts/03_finalize_dataset/output`.

**Outputs:** `scripts/05_combine/output/ML_input.csv`, the final input containing contact distances for each frame # (rows) and each contact pair within the threshold (column).

<br>

##### 📁 06_acf  
This directory contains code that applies autocorrelation functions (ACF) from the statsmodels library ([`statsmodels.tsa.stattools.acovf`](https://www.statsmodels.org/dev/generated/statsmodels.tsa.stattools.acovf.html)) to identify the optimal location (i.e., frame) within each Wnt's trajectory to use for deriving training and testing splits. There are two scripts:  
- `scripts/06_acf/run_acf.py`: runs the autocorrelation function for each of the `980` features individually from `scripts/05_combine/output/ML_input.csv`. Stores the result in `scripts/06_acf/output`.
- `scripts/06_acf/plot_acf.py`: takes the `scripts/06_acf/output/*.npy` generated by the prior step, averages over each of the features, and plots the average ACF with the frame cut-points. 

<br>

##### 📁 07_preprocess
This directory contains code for performing feature selection, which includes two stages of clustering:  
1. The first stage clusters contact pairs by known Wnt8a regions. 
2. The second stage, performs sub-clustering within the stage 1 clusters,
which are then randomly sampled to identify a reduced set of contact pairs.

The `scripts/07_preprocess/` directory, contains the following:
- `preprocess.ipynb`: runs the following steps:  
  - Read in the input data from `scripts/05_combine/output/ML_input.csv`
  - **Step 1:** Splits the data into training and test sets based on the tau value calculated from `scripts/06_acf/run_acf.py` (Wnt1: `16760`, Wnt3a: `19451`, Wnt5: `16689`, Wnt8a: `16677`).  
    - `test_set` contains frames from `0 - tau` (n=`69,577`).  
    - `training_set` contains frames from `2*tau - nFramesPerSystem` (n=`160,846`). 
  - **Step 2:** Derive the spearman correlation matrix from the training set 
    - The spearman correlation matrix is defined using `scipy.stats.spearmanr()`  
    - Ensure matrix is symmetric with diagonal as 1 (integer). 
  - **Step 3:** Performs two-rounds of hierarchical clustering using [scikit-learn `AgglomerativeClustering(metric='euclidean', linkage='ward')`](https://scikit-learn.org/1.5/modules/generated/sklearn.cluster.AgglomerativeClustering.html).
    1. Performs initial clustering using only the Wnt8a contact pairs indices within the `12Å` cut-off (`04_map/output/WNT8a_threshhold12_labels.txt`) 
       - Identify 22 clusters using only the Wnt8a contact pairs index (without contact distance data), determined using the highest Silouette Score (`ss=0.76`)
       - Then the contact pair cluster indexes are assigned to the 22 clusters
    2. Performs hierarchical subclustering (using [scikit-learn `AgglomerativeClustering(metric='euclidean', linkage='ward')`](https://scikit-learn.org/1.5/modules/generated/sklearn.cluster.AgglomerativeClustering.html)) on the distance matrix (1-spearman correlation matrix) containing only training data.
       - For each of the 22 clusters, do the following (note that only 18 clusters contain contact pairs in all Wnts):
         - If: the cluster contains 4 or more contact pairs, then perform subclustering and randomly select 1 contact pair from each subcluster
         - Elif: the cluster contains 2-3 contact pairs, randomly select a contact pair from it
         - Else: the cluster contains 1 contact pair select it
       - Results in 210 features (ss `average=0.316`, `range=0.22-0.42`)
- `model/`: Contains the train + test sets
- `output/`: 
	- `clusterbyindex_df.pkl`: Data structure containing each of the sub-clusters from the 22 initial clusters & the elements within
	- `subcluster_*.txt`: Contact pairs within the 22 clusters
	- `feature_importances_region_*.txt`: Important features within each region
- `plots/`: All of the clustering plots
  - Wnt8a Region Colors
      - 🔴 Red = Hairpin 1
      - 🟢 Green = PAM
      - 🔵 Blue = Hairpin 2
      - 🟣 Purple = Hairpin 3
      - 🌑 Gray = N-Term

<br>

##### 📁 08_model_build 
This directory contains code to build the one-versus-one Random Forest classifier (i.e., pairwise comparisons between all class pairs),
including tuning hyperparameter via 3 iterations of grid search with 10-fold cross-validation. It also contains code to determine the final model performance.

The following hyperparameter values were examined:
- Number of trees in forest (`n_estimators`): [10, 25, 50, 75, 100, 200, 1000]
- Maximum depth of the tree (`max_depth`): [1, 2, 3, 4, 5, 10, None]
- Minimum number of samples required to be at a leaf node (`min_samples_leaf`): [1, 1608, 4201]

The `scripts/08_model_build/` directory, contains the following:
`model_train.py`: `3` iterations of grid search with `10`-Fold cross-validation
`model/`: The GridSearch result from the three iterations.
  - `gridsearch_final_stratified_0.pkl`: 1st iteration of grid search
  - `gridsearch_final_stratified_1.pkl`: 2nd iteration of grid search
  - `gridsearch_final_stratified_2.pkl`: 3rd iteration of grid search

<br>

##### 📁 09_model_eval 
This `scripts/09_model_eval/` directory contains code to generate a learning curve and performs permutation feature importance. Permutation feature importance is computed 100 times over the final tuned model using the test set.

This directory contains code to generate a learning curve and performs permutation feature importance. 
- `compute_and_plot_learning_curve.py`: Get the learning curve with optimized parameters
- `compute_feature_importances.ipynb`: Compute permutation feature importance using [scikit-learn `permutation_importance`](https://scikit-learn.org/1.5/modules/generated/sklearn.inspection.permutation_importance.html)  
- `model/final_model.pkl`: Final model from learning curve (`max_depth`=3; `n_estimators`=10, `min_samples_leaf`=4021)
- `output/`: 
  - Learning Curve Output: `train_scores_finalmodel.txt` and `test_scores_finalmodel.txt` 
  - Permutation Feature Importance Output: `feature_importances.txt`
- `plots/`: learning curve and feature importance plots

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

<center><img src="https://media.github.ibm.com/user/430879/files/fb31b0d4-849a-4d06-8743-a3571bae867f" align="center" width=200>
<img src="https://media.github.ibm.com/user/430879/files/a1e5f43d-9497-494c-bb2a-32d9c65bc589" align="center" width=85></center>
