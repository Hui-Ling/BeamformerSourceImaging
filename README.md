# Beamformer-based Source Imaging
Brainstorm plug-ins for MEG and EEG source imaging, including (1) maximum contrast beamformer (MCB) for the use in localization of brain sources [1], (2) spatiotemporal imaging of linearly-related source components (SILSC) for the use in identification of source component linearly correlated to the activity of a seed [2], and (3) beamformer-based imaging of phase-amplitude coupling (BIPAC) for the use in identification of sources coupled to the specified seed [3], which has been awarded as the open finalist of student paper competition in IEEE EMBC 2015. MCB is a beamformer-based source localization method, whereas SILSC and BIPAC are seed-based methods for correlation imaging and PAC imaging. The accuracy of BIPAC and the application to investigate the brain circuit of emotional prosody processing in women with menstrual pain were described in [4]. 

## Installation
The programs of MCB, SILSC, and BIPAC are plug-ins of the Matlab Software Brainstorm. It is necessary to install Brainstorm before using the plug-ins. Please download Brainstorm in http://neuroimage.usc.edu/brainstorm/ and our plug-ins in Github: https://github.com/Hui-Ling/BeamformerSourceImaging. After finishing the installation of Brainstorm, the plug-ins (.m files) should be placed in the the following folder: $HOME\\.brainstorm\process

Please see the the tutorials of MCB, SILSC, and BIPAC in https://github.com/Hui-Ling/BeamformerSourceImaging/wiki.

## Reference
[1] Chen, Y.-S., C.-Y. Cheng, J.-C. Hsieh, and L.-F. Chen (2006), 'Maximum contrast beamformer for electromagnetic mapping of brain activity', IEEE Trans. Biomed., 53 (9), 1765-74. DOI: https://doi.org/10.1109/TBME.2006.878115

[2] Chan, H.-L., L.-F. Chen, I.-T. Chen, and Y.-S. Chen (2015), 'Beamformer-based spatiotemporal imaging of linearly-related source components using electromagnetic neural signals', NeuroImage, 114, 1-17. DOI: https://doi.org/10.1016/j.neuroimage.2015.03.038

[3] Chan, H.-L., Y.-S. Chen, L.-F. Chen, and S. Baillet (2015), 'Beamformer-based imaging of phase-amplitude coupling using electromagnetic brain activity', 37th Conf. Proc. IEEE Eng. Med. Biol. Soc. (Milan, Italy), 7558-61. DOI: https://doi.org/10.1109/EMBC.2015.7320141

[4]  Chan, H.-L., Low, I., Chen, L.-F, Chen Y.-S., Chu, I.-T., Hsieh, J.-C. (2021), A novel beamformer-based imaging of phase-amplitude coupling (BIPAC) unveiling the inter-regional connectivity of emotional prosody processing in women with primary dysmenorrhea', Journal of Neural Engineering. 
DOI: https://doi.org/10.1088/1741-2552/abed83


