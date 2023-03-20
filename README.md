## Skagit Chinook IPM: Integrated Population Models with application to Skagit River Chinook Salmon Recovery Evaluation 

The code here will allow anyone to reproduce results from the final project report for EPA Near-Term Action (NTA) 2018-0697. The goal of this project was to develop, apply, and evaluate a Bayesian state-space life cycle IPM for Skagit River Chinook Salmon (_Oncorhynchus tshawytscha_) to enable inference about population dynamics from multiple sources of life stage specific data including age composition, smolt abundance, escapement, and harvest.

We used non-least squares (NLS) models to test effects of density dependence on freshwater and marine productivity and determined that freshwater productivity was best described by a non-density dependent function, while marine productivity was best described by a Ricker model ([code](/analysis/SkagitChinook_TwoStage_NLS_DensDepModelSel.R)). To determine the ability of the IPM to detect influences of [environmental covariates](/data/
skagit_chinook_covars_all.csv) on productivity, we compiled data on temperature, river flow, weather, and ocean condition indices and used NLS models to test best combinations of correlated variables ([code](/analysis/SkagitChinook_TwoStage_NLS_CovModelSel.R)).

We developed a state-space life cycle IPM that incorporated freshwater and marine life stages with environmental covariates, and included functions for random process and observation error ((code)[SkagitChinook_TwoStage_IPM.R]). We evaluated if our model successfully described and linked freshwater and marine life stages of Chinook Salmon using model convergence and fit metrics, comparison to parallel (non-least squares models)[SkagitChinook_TwoStage_NLS.R], and examined the ability of the model to correctly identify parameters from simulated datasets.

Future work on this model will involve improving parameter identifiability, further exploring the effects of covariates on life cycle stages, and adding additional data.
