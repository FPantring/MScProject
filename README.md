# MSc Project
## Inference of somatic mutations in cancer samples from ambiguously mapped reads

The scripts listed in this directory were used in the analysis for my MSc Thesis. The project investigated whether variant calls in MuTect2 filtered by FilterMutectCalls due to low mapping quality can be recovered. To quantify the number of true variants among the filtered variant calls, a WGS simulation was analysed by applying the scripts available here. 

The script [Pileup_Processing_py](/Pileup_Processing.py) can be applied to real data, for the format of input files read the instructions in the script.

All other scripts were run dependent on files specific to the simulation. Because the simulation is as yet unpublished the results cannot be replicated.
