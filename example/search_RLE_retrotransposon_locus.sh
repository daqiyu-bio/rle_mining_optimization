set -xe

script=../src/RLE.mining/mining.search_RLE_locus
DB=./data
outpath=./result
mkdir -p $outpath/locus
mkdir -p $outpath/log

python $script/search_RLE_retrotransposon_locus.py  \
-genome $DB/ncbi_genome.example.csv  \
-hmm $DB/domain_hmm_model/reference-RVT_1-RLE.hmm  \
-log $outpath/log  \
-output $outpath/locus
