
PGM ?= ./midas.monitor.Abs

midas.monitor.Abs: monitor.c
	gcc monitor.c -o midas.monitor.Abs

test: midas.monitor.Abs
	./test.sh ${PGM}

install:
	cp midas.monitor.Abs ${PGM}
