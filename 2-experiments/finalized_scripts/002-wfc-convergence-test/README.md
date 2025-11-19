```bash
while true; do for p in bands phonon_frequencies eos pressure cohesive_energy; do ./run_converge_test.py --protocol standard --library nc-dojo-v0.4.1-str  --concurrent 20 --property ${p} --ncpus 32 --configuration BCC; done; sleep 600; done
```
