## Simple example of a pipeline that generates viral consensus sequences from reads

(teaching example -- not meant for production)

### Setup

```
git clone https://github.com/neherlab/CompBio2023.git

cd CompBio2023

conda env create -f conda.yml

conda activate viral_assembly
```


### Run the pipeline on example data

```
snakemake --cores 4 -p
```




