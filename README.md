# RecoE2E
End-to-end ML-based CMS particle reconstruction and tagging.  

The End-to-End Framework (*E2EFW*) consists of three main packages **_DataFormats_**, **_FrameProducers_** and **_Taggers_**.    
  
  
The **_DataFormats_** package consists of all the objects and classes needed for running the **E2EFW** modules and for storing any output back into EDM-format files. These consist of convenience classes for handling inputs, association maps with other pertinent collections are defined here.  
  
  
The **_FrameProducers_** is primarily responsible for extracting detector data—either as whole event data or as object-level data and auxilliary functions aiding in  this regard. A number of modules are further provided that interface with the output of detector data producers for creating localized windows or crops around the coordinates of a desired reconstructed physics object.  
  
  
The **_Taggers_** package is included to facilitate more complex analysis to be performed on the reconstructed objects like running a deep learning model inference for particle tagging, reconstruction and classification. 

## Instructions to exucute **_E2EFW_**:  
  
Setup the CMSSW environment:   
> $ cmsrel CMSSW_10_6_X  
> $ cmsenv

Clone the repository:  
> $ git clone https://github.com/cms-e2e/RecoE2E.git  

Clean build the code:  
> $ scram build vclean  
> $ scram b
  
For performing the deep learning inference, include the model graph protobuf file (preferably created with tensorflow version less than or equal to tf 1.13 ) inside the tfModels directory. Please ensure that appropriate input-output node names are provided to the Tagger modules. Once the object-level images are generated FrameProducer, a user can set the appropriate flags (doECAL, doHBHE, doBPIX, etc.) to include the required channels and has the flexibility to reorder channels (setChannelOrder). Following is a sample execution of the deep learning inference on top jets using a demo ResNet model for 8 input channels (pT, ECAL, HCAL, d0, dz, 3 BPIX layers).  
> $ cmsRun RecoE2E/TopTagger/python/TopInference_cfg.py inputFiles=file:*[input EDM-format file location]* doTracksAtECALadjPt=False maxEvents=30 setChannelOrder="0,1,2,3,4,5,6,7"
