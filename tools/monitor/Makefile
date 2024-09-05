
PGM ?= ./midas.monitor.Abs

midas.monitor.Abs: monitor.c
	gcc monitor.c -o midas.monitor.Abs

test: midas.monitor.Abs
	./test.sh ${PGM}

install: midas.monitor.Abs
	cp midas.monitor.Abs ${PGM}
