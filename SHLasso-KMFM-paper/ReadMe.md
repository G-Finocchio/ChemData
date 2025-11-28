# SHLasso-KMFM-paper

This repository contains all materials related to the **Modelling and estimation of chemical reaction yields from high-throughput experiments** paper, including the Python package implementation, analysis scripts, datasets, and supplementary results.  
The project is organized into two main components:

- **KMFM_main_paper/** — Core analyses and package implementation  
- **Supplementary_results/** — Auxiliary results and tutorials  

---

## 📂 Folder Structure Overview

###  KMFM_main_paper/
Contains the main computational framework, data, and analyses for the paper.



#### • Analysis/
Contains analysis notebooks for each HTE dataset:
- **Analysis_Buchwald_Hartwig_amination/**
  - `Analysis_shim_main_init.ipynb` is a Jupyter notebook for model fitting and coefficient extraction. It also saves an Excel table with coefficients and their names.
- **Analysis_Deoxyfluorination/**
  - `Analysis_shim_main_init.ipynb` is a Jupyter notebook for model fitting and coefficient extraction. It also saves an Excel table with coefficients and their names.
- **Analysis_Suzuki_Miyaura/**
  - `Analysis_shim_main_init.ipynb` is a Jupyter notebook for model fitting and coefficient extraction. It uses our fast SHL for 5 factors and saves an Excel file with the coefficients and their names.

#### • Data/
Datasets used in the analyses. Each subfolder corresponds to one reaction dataset:
- **Data_Buchwald_Hartwig_amination/**
- **Data_Deoxyfluorination/**
- **Data_Suzuki_Miyaura/**



### • Python_package/

#### • Strong_hierarchical_lasso/
The **Python package implementation** of the Strong Hierarchical Lasso (**SHL**) algorithm.  
Throughout the code, the model is also referred to as Strong Hierarchical Interaction Model (**SHIM**) .  
This package includes the full source code for hierarchical penalization models and related utility functions.

Submodules:
- `libs_2way_cb/`, `libs_3way_cb/`, `libs_4way_cb/`, `libs_5way_cb/` — Libraries useful for implementing SHL for 2-, 3-, 4-, and 5-factor models for continuous Bernoulli response.
- `libs_cb_lasso/` — Implementation of the continuous Bernoulli lasso and a grid search method.
- `libs_normal/` — Core implementations of SHL for 2-, 3-, 4-, and 5-factor models for normal response.
- `Helpers_SHIM.py` — Auxiliary functions for SHL model including coefficients recovery using the identifiability conditions.
- `SHIM.py` - Contains the main class that selects and fits the appropriate SHL model given the number of levels and the penalty and a grid search method for SHL.
- `__init__.py` — Enables import as a Python package.
- `Model_card.md` — Brief description of the model 

#### • Configuration and setup files
- `pyproject.toml` — Defines the package dependencies and metadata.
- `RUN_TO_INSTALL_SHL_package.ipynb` — Notebook that installs the SHL package locally using the command `!pip install .`.
- The **Strong Hierarchical Lasso (SHL)** package is currently distributed as part of this repository for research and review purposes.  
Following the completion of the paper’s peer-review process, the package will be officially published on [PyPI](https://pypi.org) to allow installation via: pip install Strong_hierarchical_lasso.

#### • (Optional) Supporting metadata and temporary build folder:


###  Supplementary_results/
Contains supporting material (test prediction performance on each dataset) and basic tutorials for using our SHL package.

#### • SHL_library_basic_tutorial/
Tutorial notebooks showing how to use the Strong Hierarchical Lasso library.
- `Example_CB_SHL.ipynb` — Tutorial for using SHL with continuous Bernoulli response.
- `Example_normal_SHL.ipynb` — Tutorial for using SHL with normal response.
- `Create_synthetic_datasets.py` — Script for generating example synthetic datasets.


#### • Test_predicition/
Validation notebooks for comparing our model with other benchmark models in terms of test prediction performance for each dataset.
- **Test_Buchwald_Hartwig_amination/**
- **Test_Deoxyfluorination/**
- **Test_Suzuki_Miyaura/**

---



## ⚙️ Installation

To install the SHL package locally, open the notebook `RUN_TO_INSTALL_SHL_package.ipynb` located in KMFM_main_paper folder and run its single cell. It just runs `!pip install .`.

## 🧩 CB SHL — Summary

- The **Strong Hierarchical Lasso (SHL)** model performs penalized hierarchical estimation for factorial problems with up to **five factors**.  
  The implementation supports both **normal** and **continuous Bernoulli** response types.

- For **2-, 3-, and 4-factor continuous Bernoulli models**, optimization is performed using a **coordinate descent** approach, where all model components are minimized sequentially and individually.

- For the **5-factor model**, a faster implementation is used: all components of the same order (e.g., all main effects, all two-way interactions, etc.) are minimized simultaneously.

- The **normal-response models** employ a similar optimization strategy to the 5-factor continuous Bernoulli model, where main effects and higher-order interactions are updated in grouped steps.
