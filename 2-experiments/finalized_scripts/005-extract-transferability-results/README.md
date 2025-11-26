run `extract-to-h5.py --group eos`

which will extract data into h5 file.

This command can run incrementally add data to the h5 file by run 
`extract-to-h5.py --group eos --override`

Then run `extract-eos.py` to get `eos.json` which is required by the following plot processes
