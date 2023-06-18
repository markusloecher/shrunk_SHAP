# Debiasing SHAP scores in random forests


This repository provides code to reproduce results in the paper "Debiasing SHAP scores in random forests" by Markus Loecher (2023)

The folder structure is as follows:

* The actual code to run simulations and the functions to create figures are all in the `src` folder.
* The resulting pdfs are in the `figures` folder.
* Some (most) simulations take a long time to run, so we have saved the data outputs from those runs in the `data` folder.
* Figures from the paper:
    - We do not include source code for Figures 1 and Figure 6, since these are legacy figures from the work leading to this [paper](https://kdd.isti.cnr.it/xkdd2022/papers/XKDD_2022_paper_1418.pdf) and were created prior to this research.
    - There are dedicated Rmd files to create Figures 2 and 3. The figures in the Appendix are created by `makeFigs_Appendix.R`.
    - Figures 4 and 5 were created by python code and requite installing the python module `TreeModelsFromScratch`, directions for which can be found on the [original github page](https://github.com/Heity94/AugmentedHierarchicalShrinkage). Deviating from the general pattern, even the plots are created in the `data/titanic/` folder.
* Table 1: We imported the AUC scores for *SHAP*, SHAP<sub>*o**o**b*</sub>, *MDA*, *MDI* from a [previous paper](https://link.springer.com/chapter/10.1007/978-3-031-14463-9_8) and only recreated the entry for $\widehat{\text{SHAP}}^{shrunk}_{in}$. The relevant files for this are `AUC_run_simulations.R` and `src/AUC_simulations_functions.R`. (Again, these simulations take a long time to run, so we have saved the data for your convenience)


I should note the following: this code base is more complex than strictly necessary as a result of having evolved over the years by contributions from various students both at the Bachelor and Master level. Roughly speaking, there are four separate parts, all of which led to the final insights in the paper:

- The initial research ideas were tested in python (sklearn) where we successfully separated out-of-bag from inbag data for trees/forests.
- Inspired by the need for more control on the features and attributes of trees, we then developed [our own random forest library](https://github.com/Heity94/AugmentedHierarchicalShrinkage/tree/main/TreeModelsFromScratch) (still in python) which we used to create Figures 4 and 5.
- In order to replicate our results in different settings, we eventually switched from the original [shap module](https://github.com/slundberg/shap) to [treeshap](https://github.com/ModelOriented/treeshap) in R. Most of the initial functions are found in the files `StroblData_ASTA2022.R` and `helperFuns.R`. However, those functions work only for the train/test methodology.
- The actual separation of out-of-bag from inbag data for trees/forests was achieved in a separate effort and can be found in the files `treewise_shap_simulation.R` and `sim_utils.R`.
