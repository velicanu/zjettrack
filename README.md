# zjettrack

## Plotting
The plot macros should work out of the box on hidsk0002 , but will run faster on a local computer that has the data files:

```bash
root makegdphi.C
root makezdphi.C
```

## Making skims
Before you can make skims you need some input files for tracking efficiency correction and others, for now just copy them from my directory:

```bash
cp -r /net/hisrv0001/home/dav2105/git/zjettrack/TrkCorr_Mar15_Iterative_* .
cp /net/hisrv0001/home/dav2105/git/zjettrack/L2L3VsPtEtaBinned_ak3PF.root .
```
Next compile the code and run:
```bash
```
