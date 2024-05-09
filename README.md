# Python scripts to combine Glide HTVS and AutoDockGPU docking with Deep-Docking

The scripts found here are useful python scripts that combine Glide and AutoDockGPU based docking with Deep-Docking developed by Gentile et al. in a fully automated manner. (Gentile F, Agrawal V, Hsing M, Ton A-T, Ban F, Norinder U, et al. *Deep Docking: A Deep Learning Platform for Augmentation of Structure Based Drug Discovery.* ACS Cent Sci 2020:acscentsci.0c00229.) Deep docking (DD) is a deep learning-based tool developed to accelerate docking-based virtual screening. Using a docking program of choice, one can screen extensive chemical libraries like ZINC15 (containing > 1.3 billion molecules) 50 times faster than typical docking. The scripts required to perform the Deep-Docking wokflow are found in the developers github repository at https://github.com/jamesgleave/DD_protocol. Please copy these files into the the folder of your choice and set the path in the Glide_DeepDocking_script.py or AutoDockGPU_DeepDocking_script.py python scripts to that specific folder.

## Requirements
* rdkit
* tensorflow >= 1.14.0 (1.15 GPU version recommended. If you are using cuda11, please use [nvidia-tensorflow](https://developer.nvidia.com/blog/accelerating-tensorflow-on-a100-gpus/))
* pandas
* numpy
* keras
* matplotlib
* scikit-learn
* Glide or AutoDockGPU
* openbabel
* DimorphiteDL
* multiprocessing
