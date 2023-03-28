# Harmonic/percussive source separation with MaD TwinNet and phase recovery GitHub Repository

Here, you will find the code related harmonic/percussive source separation (HPSS) with MaD TwinNet and phase recovery.

If you use any of the things existing in this repository, please cite the [corresponding paper](https://hal.archives-ouvertes.fr/hal-01812225). 


## How to use

### Dataset set-up

To reproduce the experiments conducted in our paper, you will need to download the [Dexmixing Secret Database (DSD100)](http://www.sisec17.audiolabs-erlangen.de) and to place its content in the `dataset`.

If you use this dataset, you will end up with the proper directory structure and file names, as used in the `functions/getdata_hpss` function.

If you want to use a different dataset, then you have two options: 
- either you format your file names and directory structure to match the one from DSD100;
- or you modify the file reading function `getdata_hpss` to suit your needs.


### Magnitude spectrograms

Be sure that you first have estimates of the magnitude spectra. Here, we use the MaD TwinNet architecture, which you can obtain on the corresponding [GitHub repository](https://github.com/dr-costas/mad-twinnet).

Place the singing voice magnitude estimates in the `magnitude_spectrograms/settingX` directory, where X is the index of the setting: in our paper, we used two different STFT settings so X=1 or 2.
The naming convention is `magn_spec_ind.mat` where `ind` is the index number of the song (e.g., it ranges from `0` to `49` for the DSD100 test databset). You can change this naming convention by modifying the `functions/getdata_hpss` function.


### Phase recovery

The experiments conducted in the paper compare 3 techniques. The first technique is the kernel additive model approach, whose corresponding function is `functions/kam_hpss/kam_hpss_mono`. The two other techniques use the MaD TwinNet magnitude estimates. One approach is a baseline phase recovery algorithm (= using the mixture's phase) and the other, PU-HPSS uses an iterative phase recovery procedure. This function is `functions/pu_hpss`.

The script to reproduce the experiments are placed in the `scripts` folder. They will notably record audio files in the `audio_files` folder, and some metrics (SDR, SIR and SAR) in the `metrics` folder.


## Acknowledgements

- Part of this research was funded by from the European Research Council under the European Union’s H2020 Framework Programme through ERC Grant Agreement 637422 EVERYSOUND.
- P. Magron is supported by the Academy of Finland, project no. 290190.
- S.-I. Mimilakis is supported by the European Union’s H2020  Framework  Programme (H2020-MSCA-ITN-2014) under grant agreement no 642685 MacSeNet.
- Part of the computations leading to these results was performed  on  a  TITAN-X GPU  donated  by  NVIDIA  to  K. Drossos.
- P. Magron, K.  Drossos  and  T.  Virtanen  wish  to  acknowledge  CSC-IT  Center  for  Science, Finland,  for  computational  resources.
