# Operational Guyane GSMAP configuration

This version is made to be run with operational parameters and pre-processing 
steps which acquire data from external sources. It notably requires internet access.

## Get Guyane configuration files

Download and extract in this folder : 
https://drive.google.com/file/d/1hyHFPbFz06PL4rhkAFlY3Q77OVUEvYfQ/view?usp=sharing

Or using commandline:
```
FILENAME=input_data_guyane.tar.gz
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1hyHFPbFz06PL4rhkAFlY3Q77OVUEvYfQ' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1hyHFPbFz06PL4rhkAFlY3Q77OVUEvYfQ" -O $FILENAME && rm -rf /tmp/cookies.txt
tar zxvf $FILENAME
```

## Start with an existing forcing database (GSMAP forcing data up to 2021/02/16, 
to avoid long download at init)

Download and extract in this folder : 
https://drive.google.com/file/d/18TL91zVMtnzLWzCX2teqYNe8dSYoPgsL/view?usp=sharing

Or using commandline:
```
FILENAME=databases_guyanegsmap_noassim_noprecalc.tar.gz
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=18TL91zVMtnzLWzCX2teqYNe8dSYoPgsL' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=18TL91zVMtnzLWzCX2teqYNe8dSYoPgsL" -O $FILENAME && rm -rf /tmp/cookies.txt
tar zxvf $FILENAME
```

## Run

Standard execution for the portal with run for all  :
- `./run.sh`
