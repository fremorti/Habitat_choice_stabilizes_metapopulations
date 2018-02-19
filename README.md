# Habitat_choice_stabilizes_metapopulations

All scripts used for the study submited as "Habitat choice stabilizes metapopulation dynamics through increased ecological specialisation." in Proceedings of the Royal Society B

- __adapt.py__ contains all classes used in the individual-based model: defines the individual, the metapopulation, the landscape and the sequence of events each timestep. 
- __analysis.py__ calls the classes in adapt.py to create a metapopulation object and run through the desired number of generations. Afterwards, the desired metrics of each run are saved for further analysis
