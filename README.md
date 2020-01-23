# pFDA-BioMarker
PrecisionFDA BioMarker Challenge

Use ./Split.sh PATH-INPUT PATH-SPLIT 5  to split the input files (Outcome, Phenotype, and FeatureMatrix files in PATH-INPUT) to training/testing data sets (80/20 split), and further partition the training set into 5 folds for cross valiation. The resulting files are saved in PATH-SPLIT

For sc1 model training: 
./Sentieon_BioM_SP.sh ./BioM PATH-SPLIT/sc1_Phase1_GE_train PATH-SPLIT/sc1_Phase1_GE_test OUTPUT/sc1 5 OUTPUT/LOGFILE PATH-INPUT/sc1_Phase1_GE 1.0 2 0 50 0 0.8 10 4 2 4 10.0 1 0.1 1 10 5 30
(sc1_Phase1_GE_train and sc1_Phase1_GE_test are the prefixes for the cross validation and training/testing split files)

For sc2 model training: 
./Sentieon_BioM_SP.sh ./BioM PATH-SPLIT/sc2_Phase1_CN_train PATH-SPLIT/sc2_Phase1_CN_test OUTPUT/sc2 5 OUTPUT/LOGFILE PATH-INPUT/sc2_Phase1_CN 0.5 1.5 0 50 0 0.6 5 3 0.5 3 5.5 1 0.1 1 10 0 20
(sc2_Phase1_CN_train and sc2_Phase1_CN_test are the prefixes for the cross validation and training/testing split files)

For sc3 model training: 
./Sentieon_BioM_SP.sh ./BioM PATH-SPLIT/sc3_Phase1_CN_GE_train PATH-SPLIT/sc3_Phase1_CN_GE_test OUTPUT/sc3 5 OUTPUT/LOGFILE PATH-INPUT/sc3_Phase1_CN_GE 0.5 1.5 0 50 0 0.6 5 3 2 4 7.5 1 0.1 1 10 5 30
(sc3_Phase1_CN_GE_train and sc3_Phase1_CN_GE_test are the prefixes for the cross validation and training/testing split files)
