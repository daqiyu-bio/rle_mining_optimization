set -xe

script=../src/RLE.mining/mining.build_consensus
outpath=./result
mkdir -p $outpath/MSA
mkdir -p $outpath/consensus

# Multiple sequence alignments (MSAs) were constructed using MAFFT.
# R4-2_SLi, the focal non-LTR retrotransposon characterized in this study, was 
# originally identified as GCA_022664705.1_nonLTR_g249 during initial screening.
locus_file=$outpath/locus/GCA_022664705.1.output.csv
grep -w GCA_022664705.1_nonLTR_g249 ${locus_file} | cut -f7,21 |sed -r 's/^/>/g;s/\s/\n/g' > ${outpath}/MSA/GCA_022664705.1_nonLTR_g249.fas
grep -w GCA_022664705.1_nonLTR_g249 ${locus_file} | cut -f7,9 |sed -r 's/^/>/g;s/\s/\n/g' > ${outpath}/MSA/GCA_022664705.1_nonLTR_g249.pep.fas
mafft-linsi ${outpath}/MSA/GCA_022664705.1_nonLTR_g249.fas > ${outpath}/MSA/GCA_022664705.1_nonLTR_g249.mafft.fas

# Build consensus sequence.
echo "GCA_022664705.1_nonLTR_g249 ${outpath}/MSA/GCA_022664705.1_nonLTR_g249.mafft.fas" > $outpath/MSA.list
python $script/build_consensus_for_multicopy_RLE.py  \
-msa-list $outpath/MSA.list  \
-output $outpath/consensus
