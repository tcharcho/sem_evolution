POPSIZE=10
MNM=3
STATES=6

target:
	g++ -lm -O3 SemDistEvo.cpp sem.cpp stat.cpp -Dpopsize=${POPSIZE} -DMNM=${MNM} -DSTATES=${STATES}
	./a.out

run:
	bash test_run.sh $(dir)
	bash parse.sh $(dir)
