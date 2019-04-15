#!/bin/sh

# ------------------------------------------------------------------------ #
# Used to run tests for evolving SEMs to classify DNA Sequences
#
# Modify the values for the 3 hyperparameters as per tuning requirememts
#   - population size
#   - max number of mutations
#   - number of states
#
# The resulting executables will run in the background
# Remove & from last command to run one after another
# ------------------------------------------------------------------------ #

# popsize
for i in 5 #500 2000 #10 100 1000
do
  # MNM
  for j in 3 7 9
  do
    # states
    for k in 3 10 30 #6 18 24
    do
      g++ -lm -O3 SemDistEvo.cpp sem.cpp stat.cpp -Dpopsize=$i -DMNM=$j -DSTATES=$k -o run_${i}_${j}_${k}.out
      ./run_${i}_${j}_${k}.out $1 &
    done
  done
done
