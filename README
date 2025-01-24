# Multiple Myeloma Staging and Progression Prediction
This repository contains the code and resources associated with our study on enhancing multiple myeloma (MM) staging accuracy and predicting the progression of asymptomatic cases (MGUS) to symptomatic MM using machine learning approaches. Below is an overview of the repository structure, datasets, and analysis methods.
Publication link: https://www.mdpi.com/2072-6694/17/2/332

# Abstract
Background: The accurate staging of multiple myeloma (MM) is essential for optimizing treatment strategies, while predicting the progression of asymptomatic patients, also referred to as monoclonal gammopathy of undetermined significance (MGUS), to symptomatic MM remains a significant challenge due to limited data. This study aimed to develop machine learning models to enhance MM staging accuracy and stratify asymptomatic patients by their risk of progression. Methods: We utilized gene expression microarray datasets to develop machine learning models, combined with various data transformations. For multiple myeloma staging, models were trained on a single dataset and validated across five independent datasets, with performance evaluated using multiclass area under the curve (AUC) metrics. To predict progression in asymptomatic patients, we employed two approaches: (1) training models on a dataset comprising asymptomatic patients who either progressed or remained stable without progressing to multiple myeloma, and (2) training models on multiple datasets combining asymptomatic and multiple myeloma samples and then testing their ability to distinguish between asymptomatic and asymptomatic that progressed. We performed feature selection and enrichment analyses to identify key signaling pathways underlying disease stages and progression. Results: Multiple myeloma staging models demonstrated high efficacy, with ElasticNet achieving consistent multiclass AUC values of 0.9 across datasets and transformations, demonstrating robust generalizability. For asymptomatic progression, both modeling approaches yielded similar results, with AUC values exceeding 0.8 across datasets and algorithms (ElasticNet, Boosting, and Support Vector Machines), underscoring their potential in identifying progression risk. Enrichment analyses revealed key pathways, including PI3K-Akt, MAPK, Wnt, and mTOR, as central to MM pathogenesis. Conclusions: To the best of our knowledge, this is the first study to utilize gene expression datasets for classifying patients across different stages of multiple myeloma and to integrate multiple myeloma with asymptomatic cases to predict disease progression, offering a novel methodology with potential clinical applications in patient monitoring and early intervention.

# Repository Structure
`src/`: Python and R scripts for preprocessing, modeling, and enrichment analysis.\
`data/`: a data folder needs to be created for the data to be downloaded
`results/`: a results folder needs to be created to save the results

## Scripts for data donwload, preprocessing, normalization and transformations
`src/R/download_data/`: R scripts to download the data from GEO and ArrayExpress. \
`src/R/preprocessing/`: R scripts to preprocess the data and perform the quantile normalization in a train - test fashion. \
`src/python/`: Python scripts for the following data transformations: binary, ranking, ratios.

## Scripts for MM staging
`src/R/caret/disease_stages/`: R scripts to train, predict and evaluate the models.\
`src/R/caret/disease_stages/Interpretation/`: R scripts to interpret the models and perform enrichment analysis.

## Scripts for MGUS progression
`src/R/caret/progression/`: R scripts to train, predict and evaluate the models \
`src/R/caret/progression/Interpretation/`: R scripts to interpret the models and perform enrichment analysis.

# Requirement
- python 3.12.8
- R version 4.4.2
- Required R libraries:
    - data.table
    - GEOquery
    - affy
    - oligo
    - ggplot2
    - cowplot
    - caret
    - randomForest
    - glmnet
    - gbm
    - e1071
    - kernlab
    - pROC
    - doMC
    - clusterProfiler
    - ReactomePA
    - org.Hs.eg.db
    - annotate
    - hgu133plus2.db
    - pheatmap
    - dplyr

