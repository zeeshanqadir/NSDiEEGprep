Follow order of code in main.m

1. preprocNSDMef
	- subjectConfig: which tasks and runs (copies) to load? Which out of these are bad to skip?
	- metadata/sub_sessions.txt -- loaded and added as column in events to use as possible predictor
	- load mef data
	- create events structure that only has stim and test trials. Rest are "pre_duration"
		- Keep only trials that are good DURING trial and in rest before trial
		- Also remove trials with high interictals during or in rest before trial
		- random effects predictors: task number run, session
	- ieeg_highpass data
	- convert data to match events, on interval [-2, 2)s. Remove any trials missing this full range
	- re-referencing by CARVariance, with 'var' option and bottom 25%ile. badChs are SOZ and "bad" chs
	- concatenate all trials across all runs
	- badChs unique to each run stored in taskruns variable
	- save data as .mat file
	- save eventsST as text file
	- currently, EKG channel is not saved

2. analyzePreprocScalarBB
	- Remove all trials corresponding to bad tasks/runs
	- bipolar-reference electrodes, using manually written lead segment data in metadata/.
	- for each run, find which bipolar electrodes are bad. Create masks (chs x trs) that equal 1 when the channel is bad in the run in which the trial is located
	- for electrode localization, I use getXyzsInf to inflate electrodes. Checks 4mm of white, and also 2mm of pial, which is assigned if none found in proximity to white.
	- aCAR electrode analysis
		- PSDs: hann(srate/4), pwelch on 0-0.5s, normalized by -0.5-0s.
		- scalar broadband power = geometric mean of PSD from 70 to 170 Hz in 1 Hz bins, ignoring 115-125 range
		- rsqs are calculated using power. Keep only channels with at least 5 good runs. Saved as elecsBB1.tsv -> do this in log power as well
		- Noise ceiling SNR and permutation p-value calculated with estimateNCSNR on log-transformed scalar broadband power. Saved as elecsBB2.tsv
	- The same thing as above is done for the bipolar-referenced electrodes