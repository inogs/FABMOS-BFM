# FABMOS-BFM
Python scripts to perform global ocean offline simulations with BFM.

## INSTALLING
To be sure to use the last version of the BFM chechkout the master branch of the BFM repository 

```
cd fabmos/extern/fabm/extern/ogs
git checkout master
```

Back to fabmos

```
cd ../../../../
```

Create a conda environment with the necessary packages

```
conda env create -f environment.yml
conda activate fabmos
```

Install FABMOS

```
pip install .
```

## RUN SIMULATIONS

To run a simulation a transport matrix is needed for the circulation. Several matrices are available at [Samar Khatiwala](http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/). Download one in the testcases directory and unzip it.

For example the MITgcm 2.8° global configuration

```
cd testcases
wget http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/MITgcm_2.8deg.tar
tar xopf MITgcm_2.8deg.tar
```

Use the directory of the transport matrix as working directory to launch simulations

```
cd MITgcm_2.8deg
```

Link a configuration file from the BFM repository

```
ln -s ../../fabmos/extern/fabm/extern/ogs/fabm_monospectral_2xDetritus.yaml fabm.yaml
```

Link the python script to run the simulation

```
ln -s ../BFM.py
```

Launch the simulation

```
python BFM.py
```

