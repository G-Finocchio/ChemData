---
# For reference on model card metadata, see the spec: https://github.com/huggingface/hub-docs/blob/main/modelcard.md?plain=1
# Doc / guide: https://huggingface.co/docs/hub/model-cards
{{ card_data }}
---

# Model Card for Strong_hierarchical_lasso

The model is a penalized GLM with hierarchical constraints.

## Model Details

### Model Description

The `Strong_hierarchical_lasso` Python package provides tools to fit **hierarchical models with two to five factors** and an arbitrary number of factor levels, supporting both **normal** and **continuous Bernoulli** response types.
Throughout the code, the model is also referred to as Strong Hierarchical Interaction Model (**SHIM**).  It has the following submodules:

Submodules:
- `libs_2way_cb/`, `libs_3way_cb/`, `libs_4way_cb/`, `libs_5way_cb/` — Libraries useful for implementing SHL for 2-, 3-, 4-, and 5-factor models for continuous Bernoulli response.
- `libs_cb_lasso/` — Implementation of the continuous Bernoulli Lasso and a grid search method.
- `libs_normal/` — Core implementations of SHL for 2-, 3-, 4-, and 5-factor models for normal response.
- `Helpers_SHIM.py` — Auxiliary functions for SHL model including coefficients recovery using the identifiability conditions.
- `SHIM.py` - Contains the main class that selects and fits the appropriate SHL model given the number of levels and the penalty and a grid search method for SHL.
- `__init__.py` — Enables import as a Python package.

Normal response model:  
  The main effects are updated simultaneously, followed by the two-way interactions under the hierarchical constraints.  
  Higher-order interactions (three-way and above) are then updated in the same grouped manner.

Continuous Bernoulli response model: 
  A coordinate descent approach is employed.  
  Each model component is updated individually for designs with up to four factors.  
  To ensure numerical stability, each coefficient update step compares the new loss to the previous one. If there is no improvement, the previous estimate is kept. This guarantees that the algorithm is non-increasing and converges.  
  For the five-factor model, a faster grouped update scheme (analogous to the one used for the normal response) is implemented to improve computational efficiency.





- **Developed by:**  Tatyana Krivobokova, Razvan-Andrei Morariu, Gianluca Finocchio
- **Model type:** Supervised learning: Hierarchical generalized linear model
- **Language:** Python
- **License:** MIT

### Model Sources [optional]

<!-- Provide the basic links for the model. -->

- **Repository:** https://github.com/G-Finocchio/ChemData/tree/main/SHLasso-KMFM-paper/KMFM_main_paper/Strong_hierarchical_lasso
- **Paper:** https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/6818bb8ae561f77ed4227799/original/modelling-and-estimation-of-chemical-reaction-yields-from-high-throughput-experiments.pdf
- **Demo:** https://github.com/G-Finocchio/ChemData/tree/main/SHLasso-KMFM-paper/Supplementary_results/SHL_library_basic_tutorial

## Uses

The model is designed to be used on factorial data up to 5 factors and normal or continuous Bernoulli response. For details consult the associated paper.

## How to Get Started with the Model

The **Demo** is a good way to firstly use the model. More complex cases where the model is used on real data can be found at [Demo](https://github.com/G-Finocchio/ChemData/tree/main/SHLasso-KMFM-paper).

## Training Details


### Training Procedure

Normal response and CB response with 5 factors: Lasso is applied on main effects and then on interactions under heredity constraints.  
Continuous Bernoulli response up to 4 factors: A coordinate descent approach is employed. For further details, consult the main paper and its *Supplementary Material*.  


#### Speeds

 The running times for our model were measured on an Intel® Core™ Ultra 7 155U 1.70 GHz
 processor. All computations were performed in python (version 3.12.7), without parallelization.
 Training the hierarchical model from scratch on the main dataset (Buchwald-Hartwig amination) takes between 0.5 and 5 hours, depending
 on the amount of regularization and the convergence tolerance, which varies between 5 × 10−2 and 1 × 10−3.
 Higher regularization levels lead to faster convergence due to the model’s hierarchical constraints.
 Training the model on the synthetically generated data from the demo is completed in seconds.

## Evaluation



###  Data & Metrics

#### Data

The real model was used on three real datasets that can be found at [Data](https://github.com/G-Finocchio/ChemData/tree/main/SHLasso-KMFM-paper/KMFM_main_paper/Data).


#### Metrics

The model is evaluated using the **Root Mean Squared Error (RMSE)** and the **R² score**.  
However, it should be noted that the primary goal of the model is **accurate estimation and interpretability**, rather than solely maximizing predictive performance.

### Results

See https://github.com/G-Finocchio/ChemData/tree/main/SHLasso-KMFM-paper.

## Environmental Impact

<!-- Total emissions (in grams of CO2eq) and additional considerations, such as electricity usage, go here. Edit the suggested text below accordingly -->

Carbon emissions can be estimated using the [Machine Learning Impact calculator](https://mlco2.github.io/impact#compute) presented in [Lacoste et al. (2019)](https://arxiv.org/abs/1910.09700).

- **Hardware Type:** Personal laptop (Intel® Core™ Ultra 7 155U, 1.70 GHz)
- **Hours used:** Approximately 0.5–5 hours per full training run (depending on regularization)
- **Cloud Provider:**  Not applicable (local computation)
- **Compute Region:**  Not applicable (local computation)
- **Carbon Emitted:** Not formally estimated; expected to be minimal


## Model Card Contact

For questions or feedback regarding this model, please contact:  
- **Tatyana Krivobokova** — [tatyana.krivobokova@univie.ac.at](mailto:tatyana.krivobokova@univie.ac.at)  
- **Razvan-Andrei Morariu** — [razvan-andrei.morariu@univie.ac.at](mailto:razvan-andrei.morariu@univie.ac.at)

