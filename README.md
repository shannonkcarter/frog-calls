# frog-calls
project measuring shifts in phenological distributions of frog calling-- published in Ecology Letters

QUESTION: 1) do phenological shifts cause long-term trends in interaction potential among competitors? 2) Are changes in interaction potential caused by uniform or non-uniform phenological shifts?

DATA: Frog calling data collected on audio recorders from 8 ponds in East Texas (Stephen F. Austin Experimental Forest and Davy Crockett National Forest). Recordings happened at 6 1-minute intervals daily (9pm, 10pm, 11pm, 12am, 1am, 2am) from May 2000 to present (used data  through Dec. 2015 for this project). Recordings were processed manually and calls of 12 species (Hyla versicolor, Hyla cineria, Bufo valliceps, Bufo woodhouseii, Rana catesbeiana, Rana clamitans, Rana sphenocephala, Gastrophryne carolinensis, Pseudacris crucifer, Pseudacris triseriata, Acris crepitans, and Rana palustris) were distinguished.

METHODS:
1) quantify phenological distribution for each species-year-pond by smoothing time series with loess function
2) for each distribution, calculate first, median, and last calling date
3) for all pairwise species combinations at each site and year, calculate a) overlap in the phenological distribution b) days-difference in first calling date c) days difference in median calling date
4) exclude species that don't overlap
5) Q1: fit linear mixed models for each species pair: overlap ~ year + (1 | pond)
6) Q2: fit linear mixed models for each species pair: overlap ~ (first/median calling date sp A - first/median calling date sp B) + (1 | pond)

RESULTS: 1) yes 2) non-uniform
