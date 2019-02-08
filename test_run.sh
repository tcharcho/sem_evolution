#!/bin/sh

# popsize
for i in 1000 100 10
do
  # MNM
  for j in 3 7 9
  do
    # states
    for k in 6 18 24
    do
      g++ -lm -O3 SemDistEvo.cpp sem.cpp stat.cpp -Dpopsize=$i -DMNM=$j -DSTATES=$k
      # g++ test.cpp -Dpopsize=$i -DMNM=$j -Dstates=$k
      ./a.out
    done
  done
done