GCC=gcc
CODE=FRESEAN
CFLAGS=
SRCDIR=src
SRCS=$(wildcard $(SRCDIR)/*.c)
OBJS=$(patsubst $(SRCDIR)/%.c, $(SRCDIR)/%.o, $(SRCS))

all: ${OBJS} avgBinMatrixList convertGMXTOP gen_list gen-modes_omp eigen extract traj_coarse projDispTrjOnSelFreqModes projDispTrjOnPLUMEDModes

${SRCDIR}/%.o: ${SRCDIR}/%.c
	${GCC} ${CFLAGS} -c -o $@ $<

projDispTrjOnPLUMEDModes: build/projDispTrjOnPLUMEDModes.o
	gcc ${SRCDIR}/fatal.o ${SRCDIR}/GMXtrrio.o ${SRCDIR}/io.o ${SRCDIR}/geo.o ${SRCDIR}/top.o ${SRCDIR}/alloc.o ${SRCDIR}/mol.o ${SRCDIR}/job.o ${SRCDIR}/select.o ${SRCDIR}/qsort.o ${SRCDIR}/align.o ${SRCDIR}/grps.o build/projDispTrjOnPLUMEDModes.o -o bin/projDispTrjOnPLUMEDModes -lm

build/projDispTrjOnPLUMEDModes.o: ${CODE}/projDispTrjOnPLUMEDModes.c
	gcc -c ${CODE}/projDispTrjOnPLUMEDModes.c -o build/projDispTrjOnPLUMEDModes.o

projDispTrjOnSelFreqModes: build/projDispTrjOnSelFreqModes.o
	gcc ${SRCDIR}/fatal.o ${SRCDIR}/GMXtrrio.o ${SRCDIR}/io.o ${SRCDIR}/geo.o ${SRCDIR}/top.o ${SRCDIR}/alloc.o ${SRCDIR}/mol.o ${SRCDIR}/job.o ${SRCDIR}/select.o ${SRCDIR}/qsort.o ${SRCDIR}/align.o ${SRCDIR}/grps.o build/projDispTrjOnSelFreqModes.o -o bin/projDispTrjOnSelFreqModes -lm

build/projDispTrjOnSelFreqModes.o: ${CODE}/projDispTrjOnSelFreqModes.c
	gcc -c ${CODE}/projDispTrjOnSelFreqModes.c -o build/projDispTrjOnSelFreqModes.o

convertGMXTOP: ${SRCDIR}/gmxtop.o 
	gcc ${SRCDIR}/fatal.o ${SRCDIR}/GMXtrrio.o ${SRCDIR}/io.o ${SRCDIR}/geo.o ${SRCDIR}/alloc.o ${SRCDIR}/gmxtop.o -o bin/convertGMXTOP -lm -lfftw3

avgBinMatrixList: build/avgBinMatrixList.o
	gcc ${SRCDIR}/fatal.o ${SRCDIR}/GMXtrrio.o ${SRCDIR}/io.o ${SRCDIR}/geo.o ${SRCDIR}/alloc.o ${SRCDIR}/progress.o build/avgBinMatrixList.o -o bin/avgBinMatrixList -lm

build/avgBinMatrixList.o:${CODE}/avgBinMatrixList.c
	gcc -c ${CODE}/avgBinMatrixList.c -o build/avgBinMatrixList.o

traj_coarse: build/traj_coarse.o
	gcc ${SRCDIR}/fatal.o ${SRCDIR}/GMXtrrio.o ${SRCDIR}/io.o ${SRCDIR}/geo.o ${SRCDIR}/top.o ${SRCDIR}/alloc.o ${SRCDIR}/mol.o ${SRCDIR}/job.o ${SRCDIR}/select.o ${SRCDIR}/qsort.o ${SRCDIR}/align.o ${SRCDIR}/grps.o build/traj_coarse.o -o bin/traj_coarse -lfftw3 -lm

build/traj_coarse.o: ${CODE}/traj_coarse.c
	gcc -c ${CODE}/traj_coarse.c -o build/traj_coarse.o

gen_list: build/gen_list.o
	gcc build/gen_list.o -o bin/gen_list -lm
build/gen_list.o : ${CODE}/gen_list.c
	gcc -c ${CODE}/gen_list.c -o build/gen_list.o 

gen-modes_omp: build/gen-modes_omp.o
	${GCC} ${SRCDIR}/fatal.o ${SRCDIR}/GMXtrrio.o ${SRCDIR}/io.o ${SRCDIR}/geo.o ${SRCDIR}/top.o ${SRCDIR}/alloc.o ${SRCDIR}/mol.o ${SRCDIR}/job.o ${SRCDIR}/select.o ${SRCDIR}/qsort.o ${SRCDIR}/align.o ${SRCDIR}/grps.o ${SRCDIR}/progress.o build/gen-modes_omp.o -o bin/gen-modes_omp -fopenmp -lfftw3 -lm

build/gen-modes_omp.o: ${CODE}/gen-modes_omp.c
	${GCC} -c ${CODE}/gen-modes_omp.c -o build/gen-modes_omp.o -fopenmp -lfftw3

eigen: build/eigen.o
	${GCC} ${SRCDIR}/fatal.o ${SRCDIR}/GMXtrrio.o ${SRCDIR}/io.o ${SRCDIR}/geo.o ${SRCDIR}/top.o ${SRCDIR}/alloc.o ${SRCDIR}/mol.o ${SRCDIR}/job.o ${SRCDIR}/select.o ${SRCDIR}/qsort.o ${SRCDIR}/align.o ${SRCDIR}/grps.o build/eigen.o -lgsl -lgslcblas -fopenmp -lm -o bin/eigen

build/eigen.o: ${CODE}/eigen.c
	${GCC} -c ${CODE}/eigen.c -o build/eigen.o -fopenmp -lfft3

extract: build/extractEigVec.o
	${GCC} ${SRCDIR}/fatal.o ${SRCDIR}/GMXtrrio.o ${SRCDIR}/io.o ${SRCDIR}/geo.o ${SRCDIR}/alloc.o build/extractEigVec.o -lm -o bin/extract

build/extractEigVec.o: ${CODE}/extractEigVec.c
	${GCC} -c ${CODE}/extractEigVec.c -o build/extractEigVec.o -lfftw3

clean:
	rm ${SRCDIR}/*.o build/*

install:
	chmod +x bin/*
	@echo 'export PATH="$(shell pwd)/bin:$$PATH"' >> ~/.bashrc
	@echo Added command line tools
