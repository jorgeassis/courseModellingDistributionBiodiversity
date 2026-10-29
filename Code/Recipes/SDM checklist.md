## Workflow (all steps)
- [ ] Confirm R 4.4 or later and the course packages (terra, biooracler, bDSSDMTools) load before starting.
- [ ] Keep each model in one project folder with `Data/`, `Code/`, `Outputs/` and `Report/`; never overwrite raw data.
- [ ] Number the scripts in the order of the loop and comment each step.
- [ ] Use relative paths only; check that a peer can run the code without editing it.
- [ ] Record data sources with download date or version.
- [ ] Deposit code and data in a public repository (GitHub with a Zenodo release, or Figshare).

## 1 Conceptualization
- [ ] State the purpose of the model: explanation of drivers, or prediction (e.g., extrapolation to new regions or times).
- [ ] Name the species, the global change threat and the realm (coastal, shelf, deep sea, pelagic, intertidal).
- [ ] Choose the scenarios and periods to compare (at least one low and one high emission SSP; present and end of century).
- [ ] State the assumption of niche conservatism if you project to other times or regions.
- [ ] Write the research question, hypothesis and conceptual model.

## 2 Data preparation
- [ ] Download records from all relevant sources (OBIS, GBIF, literature, curated datasets) and keep the source and citation of each.
- [ ] Clean in order: missing coordinates, duplicated coordinates, records on land, outside the known range, outside the known depth range; record how many each step removes.
- [ ] Save the cleaned records with a note of every decision, for the Methods.
- [ ] Choose predictors that follow from the conceptual model, including variables beyond temperature where a mechanism exists.
- [ ] Choose the depth layer (surface or benthic) that matches the species' habitat.
- [ ] Choose the statistic (mean, long-term maximum or minimum) that matches the mechanism.
- [ ] Download present and future layers with identical variables, statistics and layers, and save them with clear names.
- [ ] Record the data source, version, scenarios and decades for the Methods.
- [ ] Measure collinearity (pairwise r and VIF) in the study area; flag pairs with |r| > 0.7.
- [ ] Drop correlated variables by ecological reasoning, not at random; write the reason for each.
- [ ] Define the calibration area as the accessible area: dispersal distance (and connectivity where known) and depth range.
- [ ] State the projection area separately from the calibration area.
- [ ] Choose pseudo-absences or background by algorithm, and their number following published guidelines.
- [ ] Measure the spatial autocorrelation range of the predictors and record it for spatial cross-validation.

## 3 Model fitting
- [ ] Assemble the model dataset (response and predictor values) and check the number of presences and absences in each set.
- [ ] Split the data into training, validation and test sets with spatial blocks of the autocorrelation range. Use spatial (or environmental) block cross-validation with blocks larger than the autocorrelation range; report the block size and why. 
- [ ] Choose and report the algorithm and every setting used. For MaxEnt, select the features and the regularisation multiplier; for tree ensembles, the depth, learning rate and number of trees or rounds.
- [ ] Never accept default settings without testing alternatives: tune the main hyperparameters with a grid search.
- [ ] Where theory gives a clear direction, apply monotonic constraints and record them with their sources.
- [ ] Check that the most important predictors are those the conceptual model expected.
- [ ] Inspect the response curve of every predictor with the records alongside, and flag shapes that contradict the species' biology.
- [ ] Compare the tipping points of the response curves with known physiological tolerances. Search the literature for experimental tolerances of the species (thermal, salinity, oxygen) and record their sources.
- [ ] Check whether the projection reaches conditions beyond those of the training data, where a correlative model cannot know the response.

## 4 Assessment
- [ ] Report a suite of metrics matched to the purpose (sensitivity, AUC-PR or AUC-PRG, Boyce index), not ROC-AUC alone.
- [ ] Choose the threshold rule for the purpose (e.g., maximum TSS, 10% training omission), compute it and report it; never a fixed 0.5. 
- [ ] Map novelty (MESS and the most dissimilar variable, or the area of applicability) for every projection, and state the baseline used.
- [ ] If you ensemble, tune each member first, state how models are combined (unweighted or weighted, and by which metric), and report the members.
- [ ] Map the disagreement among models next to the ensemble.

## 5 Predictions
- [ ] Project the same fitted model (or ensemble) to scenarios and state the scenarios, the period and the climate data source.
- [ ] Binarise present and future maps with the same threshold.
- [ ] Map and report loss, gain and stable areas; report gains separately and state the dispersal assumption.
- [ ] Define refugia explicitly (e.g., stable under both scenarios) and report their area.
- [ ] Report the share of the present range and of refugia inside protected areas, by protection level (all MPAs, no-take). (optional)
- [ ] Combine climate projections with local pressures before ranking areas. (optional)
