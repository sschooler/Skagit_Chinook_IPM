## Integrated population models and habitat evaluation for Chinook Salmon recovery in the Skagit River System
**Authors: Eric Beamer, Sarah L. Schooler, Michael LeMoine** <br/>
**Collaborators: Casey Ruff, Catherine Austin** <br/>
**Organization: [Skagit River System Cooperative](http://skagitcoop.org/)** <br/>
Grant program: [EPA National Estuary Program](https://www.epa.gov/puget-sound/funding-and-grants-puget-sound)/[Puget Sound Partnership Habitat Strategic Initiative](https://pugetsoundestuary.wa.gov/habitat-strategic-initiative/) <br/>
Grant number: Near-Term Action (NTA) 2018-0697 <br/>
Contact: sschooler@skagitcoop.org

### Skagit Chinook IPM: Integrated Population Models with application to Skagit River Chinook Salmon
This code allows reproduction of results from the final project report for EPA Near-Term Action (NTA) 2018-0697. The goal of this project was to develop, apply, and evaluate a Bayesian state-space life cycle IPM for Skagit River Chinook Salmon (_Oncorhynchus tshawytscha_) to enable inference about population dynamics from multiple sources of life stage specific data including age composition, smolt abundance, escapement, and harvest.

We used [non-least squares (NLS) models](/analysis/SkagitChinook_TwoStage_NLS_DensDepModelSel.R) to test effects of density dependence on freshwater and marine productivity. To determine the ability of the IPM to detect influences of [environmental covariates](/data/skagit_chinook_covars_all.csv) on productivity, we compiled data on temperature, river flow, weather, and ocean condition indices and used [NLS models](/analysis/SkagitChinook_TwoStage_NLS_CovModelSel.R) to test best combinations of correlated variables.

We developed a [state-space life cycle IPM](SkagitChinook_TwoStage_IPM.R) based on NLS results that incorporated freshwater and marine life stages with selected environmental covariates, and included functions for random process and observation error. We evaluated if our model successfully described and linked freshwater and marine life stages of Chinook Salmon using model convergence and fit metrics, comparison to parallel [NLS models](SkagitChinook_TwoStage_NLS.R), and examined the ability of the model to correctly identify parameters from simulated datasets.

Future work on this model will involve improving parameter identifiability, further exploring the effects of covariates on life cycle stages, and adding additional data.

### Habitat status and trends results for evaluation of recovery actions for Skagit River Chinook Salmon
This project included identification of 108 restoration projects from satellite and aerial imagery and their impacts on potential habitat for spawning, freshwater rearing, and estuarine rearing for Skagit River Chinook Salmon. [Habitat trend data](/data/habitat) for landslides, freshwater habitat, and estuarine habitat indices are provided.
