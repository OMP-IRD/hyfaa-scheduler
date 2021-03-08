# Operational Niger GSMAP configuration

This version is made to be run with operational parameters and pre-processing steps which acquire data from external sources. It notably requires internet access.

## Get Niger configuration files

Download and extract in this folder : https://drive.google.com/file/d/1drCTR7Wv_pcoukbUMRMhpWaZZAdsZwSA/view?usp=sharing

## Start at an already initialized state

This step is necessary for now, because hyfaa_preprocessing_assimilation routine is not yet functionnal for Niger.
Since this task responsible for the creation and filling of the assimilation database is skipped, an assimilation database must exist at init.

### Configuration 1 : simple init with only an empty assimilation database
Download and extract in this folder : https://drive.google.com/file/d/1p1QqBh5PBZPWH-_s0qr6nPWMgF6FSLUr/view?usp=sharing

### Configuration 2 : GSMAP forcing data up to 2021/02/16 + empty assimilation DB
Download and extract in this folder : https://drive.google.com/file/d/1WAZWczoiCg-hM70c4P7nTWucnGjoiXmH/view?usp=sharing

### Configuration 3 : GSMAP forcing data up to 2021/02/16 + empty assimilation DB + hydrological states and post-processing up to 2021/02/16
Download and extract in this folder : https://drive.google.com/file/d/1POfZbej01JNmye9yvb0IAAVijhtqjZRj/view?usp=sharing


## Run

Standard execution :
- `./run.sh`

