## 000

000-oxygen-nitrigen-ref-pseudos-convergence-verification

-> Timo
- depend on null

This is for SI "The oxygen pseudopotential for further verification" section.
Scripts to run EOS for oxygen and nitrigen (do not needed anymore) and compare with AE reference.
The result will be the nu values for all oxygen pseudos.

## 001
001-import-libraries/

-> Michail
- depend on null

Run this at very beginning to import all pseudos as UPFdata in to dedicate groups.
This cannot replaced by the aiida-pseudo, because here we don't provide the cutoffs.

## 002
002-wfc-convergence-test/

-> Michail
- depend on 001

Run the convergence test for a dedicate property (e.g. cohesive energy) and a dedicate structure configuration (e.g. bcc).

## 003
003-transferability-test/

-> Michail
- depend on 001

Run the EOS (at wavefunction cutoff 200 Ry) (also bands, but bands are rerun (at the recommended cutoffs) by Michail already so that part is redundant).

## 004
004-extract-convergence-results/

-> Michail
- depend on 002

After 002 the convergence test, this is to extract the results into a h5 file.

## 005
005-extract-transferability-results/

-> Michail
- depend on 003

After 002 the eos transferability test, this is to extract the results into a h5 file.

## 006
006-plot-element-summary/

-> Timo
- depend on 004/005

After results data of convergence test and transferability test are in the h5 files, this script is to plot the summary plots for pseudos of a element.

## 007
007-rho-convergence-test/

-> Michail
- depend on 001 + manually pseudos selection (with Nicola) into a JSON wavefunction cutoffs file.

After the summary plots from 007, use that plots to discuss with nicola to see which pseudo to pick and which cutoff to pick.
Then I manually prepare a JSON file with elements as the keys and the agreed recommended cutoffs as value.

The script run the convergence test on charge density cutoff.

## 008
008-rho-convergence-extract/

-> Michail
- depend on 007

After 007, this extract the result into h5 files, and to plot (in the same folder).

## 010
010-extract-eff-lib/

-> If any pseudos changed in the future
- depend on null

Simply copy the files into a folder.

## 011
011-extract-prec-lib/

-> If any pseudos changed in the future
- depend on null

Simply copy the files into a folder.

## 012
012-band-extract/

-> Edan
- depend on 002 (only all bands structure convergence test for all bcc/fcc/sc/dc)

This will extract from the 200 Ry bands distance convergence test.
It will generate a JSON with all the mappings and a JSON with all the bands distances results that is ready for plot.

