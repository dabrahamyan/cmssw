########################################################################################################
# Runs davidNtuplePlot.C, my own version of NtuplePlot.C from Louise with some extra stuff and changes
#
# Grabs Ntuples that I store in my eos and saves all the plots in something called output_(filename)
# Also saves info on integrated efficiency into /NtuplePlotOutput
#
# By David Abrahamyan Oct 13, 2023 
########################################################################################################

# Where to get the NtupleRoot files from 
dirName="../GitCheckoutChecks/"

# File name to process (without .root)
fileName="newkfDeug_1nmatch_SingleMuon_thomasRecs_10-16-2023_1"

# Open root, run davidNtuplePlot, quit root
root -l -b -q "davidNtuplePlot_1nmatch.C(\"${fileName}\", \"${dirName}\", \"\", 0)" | tail -n 22 > NtuplePlotOutput/${fileName}_1nmatch.out 
