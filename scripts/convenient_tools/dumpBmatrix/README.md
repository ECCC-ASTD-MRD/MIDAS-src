# Bmatrix.bin dump tool

## Purpose

To obtain covariance values of the 1D $`\mathbf{B}`$ matrix.
This is especially useful when the test `extractBmatrixFor1Dvar/globalBnmcLand` fails, since `Bmatrix.bin` is a binary file:
```sh
diff -y <(./dumpBmatrix live) <(./dumpBmatrix rfrn )
```

## Compilation

1. `source ./config.dot.sh`
2. `make`
