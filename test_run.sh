#!/bin/sh

# popsize
for i in 5 2000 5000
do
  # MNM
  for j in 3 7 9
  do
    # states
    for k in 3 50 70
    do
      g++ -lm -O3 SemDistEvo.cpp sem.cpp stat.cpp -Dpopsize=$i -DMNM=$j -DSTATES=$k
      ./a.out
    done
  done
done
