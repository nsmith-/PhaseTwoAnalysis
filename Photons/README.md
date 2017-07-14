Recipes for photon selection in analysis of simulated events with the CMS phase 2 detector.
=========================

One should run:
```bash
cmsRun ConfFile_cfg.py outFilename=FilteredEvents.root inputFormat=RECO/PAT
```

whether the input file format is RECO or miniAOD.

The main producers are:
   * `plugins/PatPhotonFilter.cc` -- to run over PAT events 
   * `plugins/RecoPhotonFilter.cc` -- to run over RECO events 

Details on the object definitions are given in the `implementation` section.

For each ID quality, a vector of muons and a vector of double corresponding to the muon relative isolation are added: 
   * `doubles_photonfilter_LoosePhotonRelIso_PhotonFilter`
   * `doubles_photonfilter_MediumPhotonRelIso_PhotonFilter`
   * `doubles_photonfilter_TightPhotonRelIso_PhotonFilter`
   * `[recoGsf|pat]Photons_Photonfilter_LoosePhotons_PhotonFilter`
   * `[recoGsf|pat]Photons_Photonfilter_MediumPhotons_PhotonFilter`
   * `[recoGsf|pat]Photons_Photonfilter_TightPhotons_PhotonFilter`

The initial vector of Photons is dropped to avoid any confusion.
