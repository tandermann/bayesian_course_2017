# SSR extraction

## Create conda environment

conda create -n ssr_tool spyder
conda activate ssr_tool
conda install -c bioconda lastz
conda install numpy
conda install pandas
conda install fuzzywuzzy
conda install scikit-learn
conda install biopython
conda install -c bioconda blast
conda install -c bioconda mafft
