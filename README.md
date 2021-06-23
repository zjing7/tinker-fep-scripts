# tinker-fep-scripts
A collection of scripts for protein-ligand binding simulations using Tinker

## examples
### convert xyz to pdb
1. copy protein library file
2. generate library file for each ligand and append to the library file
3. convert xyz to pdb
```
cd example/pdb2xyz/
cp ../../scripts/dat/lib.lipid.txt lib.user.txt
python ../../scripts/pdbxyz.py -m xyz2lib -x r101.xyz -k box.key -t lib.user.txt --app
python ../../scripts/pdbxyz.py -m xyz2pdb -x box.xyz -p box.pdb -t lib.user.txt
```

