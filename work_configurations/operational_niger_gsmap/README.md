# Operational Niger GSMAP configuration

This version is made to be run with operational parameters and pre-processing steps which acquire data from external sources. It notably requires internet access.

## Get Niger configuration files

Download and extract in this folder (new from 2021/03/18) : https://drive.google.com/file/d/1l_uVyPicZAibbKxZA4N0nhFECew4MCs9/view?usp=sharing

Or using commandline:
```
FILENAME=input_data_niger.tar.gz
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1l_uVyPicZAibbKxZA4N0nhFECew4MCs9' -O $FILENAME
tar zxvf $FILENAME
```

## Start with an existing forcing database (GSMAP forcing data up to 2021/02/16, to avoid long download at init)

Download and extract in this folder : https://drive.google.com/file/d/1WAZWczoiCg-hM70c4P7nTWucnGjoiXmH/view?usp=sharing

Or using commandline:
```
FILENAME=databases_nigergsmap_noassim_noprecalc.tar.gz
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1WAZWczoiCg-hM70c4P7nTWucnGjoiXmH' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1WAZWczoiCg-hM70c4P7nTWucnGjoiXmH" -O $FILENAME && rm -rf /tmp/cookies.txt
tar zxvf $FILENAME
```

## Run

Standard execution for the portal with run for all  :
- `./run.sh`

