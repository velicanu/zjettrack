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
g++ zjetSkim.C $(root-config --cflags --libs)  -Werror -Wall -O2 -o zjetSkim.exe
g++ gammajetSkim.C $(root-config --cflags --libs)  -Werror -Wall -O2 -o gammajetSkim.exe

./zjetSkim.exe /mnt/hadoop/cms/store/user/rbi/azsigmon-HIRun2015E-PromptReco-AOD-DielectronSkim-ElePt8-v3_forest_csjet_v1_3/0.root /export/d00/scratch/dav2105/ztrees/test.azsigmon-HIRun2015E-PromptReco-AOD-DielectronSkim-ElePt8-v3_forest_csjet_v1_3.root akPu4CaloJetAnalyzer 0
./gammajetSkim.exe /mnt/hadoop/cms/store/user/rbi/merged/Pythia8_Photon30_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1_forest_v1/0.root /export/d00/scratch/dav2105/ztrees/test.Pythia8_Photon30_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1_forest_v1.root akPu3PFJetAnalyzer 0
```
These will make output files that can the plotting macros take as input. 
