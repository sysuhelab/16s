#!/bin/bash

### user data
raw_read1=$1
raw_read2=$2
mapping_table=$3
grpbase=$4

### 64 bit usearch
usearch='/share/tools/usearch11.0.667_i86linux64'
### seqtk_utils from biostack
demux='/share/tools/seqtk_utils/seqtk_demultiplex'
### csvtk from shenwei
csvtk='/share/tools/csvtk'
### mafft for fast tree building
mafft='/usr/local/bin/mafft'
### trimal for trim sequence alignment
trimal='/usr/local/bin/trimal'
### fasttree for building tree
fasttree='/usr/local/bin/FastTreeMP'
### treebender
# treebender='/usr/local/bin/treebender'

### set shell env
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

ulimit -n 100000
threads=24

####### Read preparation ########

echo -e "${GREEN}Starting ${RED}DEMULTIPLEX ${GREEN}sequence......${NC}"
${demux} -b <(${csvtk} cut -t -f SampleID,BarcodeF,BarcodeR ${mapping_table} | sed 1d) -1 ${raw_read1} -2  ${raw_read2} -l 6 -d ${grpbase}_demux
echo -e "${GREEN}Finished ${RED}DEMULTIPLEX ${GREEN}sequence......${NC}\n\n\n"

echo -e "${GREEN}Starting ${RED}MERGING ${GREEN}sequence......${NC}"
${usearch} -fastq_mergepairs ${grpbase}_demux/*_R1.fastq -relabel @ -fastqout ${grpbase}_merged.fq -threads ${threads}
echo -e "${GREEN}Finished ${RED}MERGING ${GREEN}sequence......${NC}\n\n\n"

echo -e "${GREEN}Starting ${RED}STRIPING ${GREEN}sequence......${NC}"
read -r primerF primerR <<< $(${csvtk} cut -t -f PrimerF,PrimerR ${mapping_table} | sed 1d | ${csvtk} -H uniq)
${usearch} -search_pcr2 ${grpbase}_merged.fq -fwdprimer ${primerF} -revprimer ${primerR} -minamp 300 -maxamp 600 -strand both -fastqout ${grpbase}_amplicons.fq -threads ${threads}
echo -e "${GREEN}Finished ${RED}STRIPING ${GREEN}sequence......${NC}\n\n\n"

echo -e "${GREEN}Starting ${RED}FILTERING ${GREEN}sequence......${NC}"
${usearch} -fastq_filter ${grpbase}_amplicons.fq -fastq_maxee 1.0 -relabel Filt -fastaout ${grpbase}_filtered.fa -threads ${threads}
echo -e "${GREEN}Finished ${RED}FILTERING ${GREEN}sequence......${NC}\n\n\n"

echo -e "${GREEN}Starting ${RED}DEREPLICATING ${GREEN}sequence......${NC}"
${usearch} -fastx_uniques ${grpbase}_filtered.fa -fastaout ${grpbase}_uniques.fa -minuniquesize 1 -sizeout -relabel Uniq -threads ${threads}
echo -e "${GREEN}Finished ${RED}DEREPLICATING ${GREEN}sequence......${NC}\n\n\n"


####### OTU generation ########

echo -e "${GREEN}Starting ${RED}CLUSTERING ${GREEN}sequence......${NC}"
${usearch} -cluster_otus ${grpbase}_uniques.fa -otus ${grpbase}_otus.fa -relabel Otu
echo -e "${GREEN}Finished ${RED}CLUSTERING ${GREEN}sequence......${NC}\n\n\n"

echo -e "${GREEN}Starting ${RED}GENERATING ${GREEN}otu table......${NC}"
${usearch} -otutab ${grpbase}_merged.fq -otus ${grpbase}_otus.fa -otutabout ${grpbase}_otutab.txt -mapout ${grpbase}_map.txt -notmatchedfq ${grpbase}_unmapped.fq -dbmatched ${grpbase}_otus_with_sizes.fa -sizeout -threads ${threads}
echo -e "${GREEN}Finished ${RED}GENERATING ${GREEN}sequence......${NC}\n\n\n"

# Quality control for OTU sequences
# ....

####### Taxonomy annotation ########

echo -e "${GREEN}Starting ${RED}UTAXING ${GREEN}sequence......${NC}"
if [ -f ./rdp16s.fa ]; then
   echo "RDP ref File exists."
else
  ##download rdp_v16.fa.tar.gz
  wget https://drive5.com/sintax/rdp_16s_v16.fa.gz -O rdp_16s_v16.fa.gz
  gunzip -c rdp_16s_v16.fa.gz > rdp_16s.fa
fi
${usearch} -makeudb_usearch rdp_16s.fa -output rdp_16s.udb
${usearch} -sintax ${grpbase}_otus.fa -db rdp_16s.udb -tabbedout ${grpbase}_otus.sintax -strand both -sintax_cutoff 0.8 -threads ${threads}
echo -e "${GREEN}Finished ${RED}UTAXING ${GREEN}sequence......${NC}\n\n\n"

####### Phylogeny construction ########

echo -e "${GREEN}Starting ${RED}PHYLOGENY ${GREEN}sequence......${NC}"
${mafft} --thread ${threads} --adjustdirectionaccurately ${grpbase}_otus.fa > ${grpbase}_otus.mafft.fa
${trimal} -in ${grpbase}_otus.mafft.fa -out ${grpbase}_otus.trimal.fa  -automated1
${fasttree} -nt ${grpbase}_otus.trimal.fa > ${grpbase}_otus.nwk
# ${treebender} --mid_point_root <${grpbase}_otus.nwk >${grpbase}_rooted.nwk
echo -e "${GREEN}Finished ${RED}PHYLOGENY ${GREEN}sequence......${NC}\n\n\n"
