# Package to validate basic timing performance for showers
Average time, resolution and efficiency, considering the recHits within a distance rho from the shower axis

```shell
cmsrel CMSSW_11_2_0_pre8
cd CMSSW_10_2_15/src
cmsenv
git clone https://github.com/amartelli/HGCTimingAnalysis.git
cd HGCTimingAnalysis/HGCTiming
git checkout TimingAnalysis_from11_2_0_pre8_on
```

In the test folder a template to run the analysis and a prototype for splitting jobs

In the python folder, the configuration of low level parameters to perform detailed studies (should be cross-checked in case)