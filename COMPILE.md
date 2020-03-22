First, compile the [LANL ANSI version of SPRAL](https://github.com/lanl-ansi/spral) using the [compilation suggestions described therein](https://github.com/lanl-ansi/spral/blob/master/COMPILE.md).

To rebuild configuration files, execute
```bash
git clone git@github.com:coin-or-tools/BuildTools.git
export AUTOTOOLS_DIR="${HOME}/local2"
./BuildTools/run_autotools
```
To compile Ipopt with SPRAL, specify `${SPRALDIR}` as the directory containing `libspral.a`, then execute
```bash
mkdir build && cd build
../configure --prefix=${PWD} --with-spral="-L${SPRALDIR} -L${METISDIR} \
    -lspral -lgfortran -lhwloc -lm -lmetis -lopenblas -lstdc++ -fopenmp" \
    --with-lapack-lflags="-llapack -lopenblas"
make && make install
```
Ensure the following environment variables are set when using the SPRAL library.
```bash
export OMP_CANCELLATION=TRUE
export OMP_NESTED=TRUE
export OMP_PROC_BIND=TRUE
```
