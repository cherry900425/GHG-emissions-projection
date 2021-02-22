GHG emissions projections

This code is designed for generating probabilistic GHG emissions projections in paper “Moving in the right direction: consumption-based projection shows an earlier and lower peak in global anthropogenic GHG emissions”

Code Usage

This code is mostly meant to be run sequentially.

Users must download and install JAGS software before running this code.

For use of this code, one should load all required packages in the first 4 lines, and import all data from an excel file “data_original” from line 6 to 25, including emission intensities, consumption structure, final demand structure, final-demand-to-GDP ratio, GDP per capita, and population. 

GDP projections are conducted from line 89 to 211. Population projections are conducted from line 213 to 335. Final-demand-to-GDP ratio projections are conducted from line 337 to 458. Projections for reference levels of emission intensities and consumption structure are conducted from line 460 to 930. Projections of emission intensities and consumption structure are conducted from line 933 to 3337. Projections of emissions for global emissions are conducted from line 3339 to 6075. Codes between line 6077 to 6349 are for creating median value and prediction intervals of results in the paper, which include:

	Yearly GHG Emissions Predictions (Global, Regions, 3 Final Demand Components, and 7 Consumption Items).

	LMDI Decomposition Projections

	Projections of changes in emission intensities (Global and Regions)

The results for creating plots in the paper can be found in an excel file "data_projection".
