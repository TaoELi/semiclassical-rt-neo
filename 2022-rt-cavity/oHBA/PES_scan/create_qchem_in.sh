#!/bin/bash

ratio=$1
outputname=$2

Ox=2.682115
Oy=-0.219090
Hx=1.844187
Hy=1.135044
OHx=0.837928
OHy=-1.354134

Hx_new=$(echo "print($Hx + $OHx * $ratio)" | python3)
Hy_new=$(echo "print($Hy + $OHy * $ratio)" | python3)

echo "Hx = $Hx_new, Hy = $Hy_new"

cat > $outputname << EOF
\$molecule
0 1
C -1.310008   1.258755   0.000000
C  0.019289   0.780580   0.000000
C  0.322586  -0.621636   0.000000
C -0.761283  -1.465586   0.000000
C -2.125342  -0.985468   0.000000
C -2.398971   0.362640   0.000000
O  1.008708   1.658363   0.000000
C  1.725826  -1.062678   0.000000
O  2.682115  -0.219090   0.000000
H -3.414260   0.733314   0.000000
H -0.596742  -2.537818   0.000000
H -2.926679  -1.713458   0.000000
H -1.459388   2.331114   0.000000
H  1.924136  -2.138192   0.000000
H  $Hx_new   $Hy_new   0.000000
\$end

\$rem
sym_ignore = 1
input_bohr = false
method = b3lyp
SCF_ALGORITHM = diis
MAX_SCF_CYCLES=200
basis = mixed
PURECART 1
SCF_CONVERGENCE = 8
mem_total = 7000
basis = cc-pvdz
rpa = true
set_roots = 30
\$end
EOF
