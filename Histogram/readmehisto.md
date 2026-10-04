These are the current histogram analysis programs: 
SM_MakeHistogram_new.fig (Matlab figure for SM_MakeHistogram_new.m)
SM_MakeHistogram_new.m (RUN THIS; The bulk of the histogram code; currently very verbose and needs a lot of fine tuning)
bobaFRET.m (individual bootstrap code that contains the analysis; called by SM_MakeHistogram_new.m)
bobaFRET_run.m (separate driver for the bootstrap code)
moleculeSpread.m (contains the molecule spread function; called by SM_MakeHistogram_new.m)
