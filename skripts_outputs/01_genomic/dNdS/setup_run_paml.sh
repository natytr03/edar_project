#!/bin/bash

PREFIX=${1:-edar}

INPUTDIR="$HOME/$PREFIX/msa/codonmsa"
WORKDIR="$INPUTDIR/paml"

ALIGNMENT="species_list_${PREFIX}_codon_5.final_aligned.phy"
TREE="${PREFIX}_paml_tree_5.final.nwk"

cd "$WORKDIR" || exit 1

#create directories for each model
mkdir -p m0 m1a m2a m7 m8 free

#copy alignment and tree into each directory
for dir in m0 m1a m2a m7 m8 free; do
    cp "$INPUTDIR/$ALIGNMENT" "$dir/"
    cp "$INPUTDIR/$TREE" "$dir/"
done

if [[ ! -f "$INPUTDIR/$ALIGNMENT" ]]; then
    echo "ERROR: $INPUTDIR/$ALIGNMENT not found!"
    exit 1
fi

if [[ ! -f "$INPUTDIR/$TREE" ]]; then
    echo "ERROR: $INPUTDIR/$TREE not found!"
    exit 1
fi

#setup one-ratio model M0
cat > m0/m0.ctl <<EOF
      seqfile = $ALIGNMENT
     treefile = $TREE
      outfile = m0_output.txt

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        model = 0
      NSsites = 0

        icode = 0
        clock = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 1
EOF

#setup nearly neutral model M1a
cat > m1a/m1a.ctl <<EOF
      seqfile = $ALIGNMENT
     treefile = $TREE
      outfile = m1a_output.txt

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        model = 0
      NSsites = 1

        icode = 0
        clock = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.5
EOF

#setup counterpart positive selection model M2a
cat > m2a/m2a.ctl <<EOF
      seqfile = $ALIGNMENT
     treefile = $TREE
      outfile = m2a_output.txt

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        model = 0
      NSsites = 2

        icode = 0
        clock = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 1
EOF

#setup beta-distributed model M7
cat > m7/m7.ctl <<EOF
      seqfile = $ALIGNMENT
     treefile = $TREE
      outfile = m7_output.txt

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        model = 0
      NSsites = 7

        icode = 0
        clock = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.5
EOF

#setup counterpart positive selection model M8
cat > m8/m8.ctl <<EOF
      seqfile = $ALIGNMENT
     treefile = $TREE
      outfile = m8_output.txt

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        model = 0
      NSsites = 8

        icode = 0
        clock = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 1
EOF

#setup Free-ratio model
cat > free/free.ctl <<EOF
      seqfile = $ALIGNMENT
     treefile = $TREE
      outfile = free_output.txt

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        model = 1
      NSsites = 0

        icode = 0
        clock = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 1
EOF

echo "All control files created for dataset: $PREFIX"
echo "Running PAML models"

#run each model
for model in m0 m1a m2a m7 m8 free; do
    echo "Running $model"

    cd "$WORKDIR/$model" || exit 1
    codeml "${model}.ctl"

    cd "$WORKDIR" || exit 1
done

echo "All PAML analyses completed!"

