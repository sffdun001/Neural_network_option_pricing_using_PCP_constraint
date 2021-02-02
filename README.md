This Git repository contains all the code used in the production of the dissertation in partial fulfilment of requirement for the Masters of Science in Advanced Analytics. The code presented here is catorgorised into the following files and contains the following:

Data:
  1. Data Extraction: Code used to extract data from Bloomberg terminal and explores the extracted data generating relevent figures referred to in my dissertation
  2. Data Simulation: Code used to simulate data under BSM, Black and Heston models used in the training and testing of the networks
  
Neural Networks
  1. HLP: Test and train networks based on the data simulated under BSM
  2. Black: Test and train networks based on the data simulated under Black
  3. Heston: Test and train networks based on the data simulated under Heston
  4. Market: Test and train networks based on the market data extracted from bloomberg 


This github repository should be used in conjunction with the Zivahub reposition at: https://doi.org/10.25375/uct.11770740. The zivbub repository provides all the data sets used in the production of this research. Specifically,
1) The historical price information regaring options on FTSE/JSE Top 40 futures contracts from Dec 2017 - Sept 2019.
2) Synthetic datasets simulated using BSM, Black and Heston Models
3) Networks outputs trained on data sets in 1) and 2) used to generate of the dissertation. 

Although 2) and 3) can be generated from the code in this repositoy, we include it has the some data simulation algorithms and the training of networks are very time consuming.

The first section of each code files calls required packages and imports relevant data. Please ensure you adjust the code accordingly to set the correct working directorty and install necessary packages. 


