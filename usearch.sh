#!/bin/bash
# my usearch8 is 64 bit
# but usearch9 is 32 bit only
usearch8="path to your usearch8"
usearch9="path to your usearch9"

### raw fastq files are in folder ${demulfq}
demulfq="path to your demultiplexed fastq file"
### sample group name for output file prefix, ybmouse for example
grpbase=$1

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${GREEN}Starting ${RED}MERGING ${GREEN}sequence......${NC}"
${usearch9} -fastq_mergepairs ${demulfq}/*_R1.fastq -relabel @ -fastqout ${grpbase}_merged.fq
echo -e "${GREEN}Finished ${RED}MERGING ${GREEN}sequence......${NC}\n\n\n"

echo -e "${GREEN}Starting ${RED}FILTERING ${GREEN}sequence......${NC}"
${usearch8} -fastq_filter ${grpbase}_merged.fq -fastq_maxee 1.0 -relabel Filt -fastaout ${grpbase}_filtered.fa
echo -e "${GREEN}Finished ${RED}FILTERING ${GREEN}sequence......${NC}\n\n\n"

echo -e "${GREEN}Starting ${RED}UNIQUIFYING ${GREEN}sequence......${NC}"
${usearch8} -derep_fulllength ${grpbase}_filtered.fa -strand both --minuniquesize 2 -relabel Uniq -sizeout -fastaout ${grpbase}_uniques.fa
echo -e "${GREEN}Finished ${RED}UNIQUIFYING ${GREEN}sequence......${NC}\n\n\n"

echo -e "${GREEN}Starting ${RED}REMOVING CHIMERA ${GREEN}sequence......${NC}"
${usearch8} -uchime_denovo ${grpbase}_uniques.fa -nonchimeras ${grpbase}_ncuniq.fa
echo -e "${GREEN}Finished ${RED}REMOVING CHIMERA ${GREEN}sequence......${NC}\n\n\n"

echo -e "${GREEN}Starting ${RED}CLUSTERING ${GREEN}sequence......${NC}"
${usearch8} -cluster_otus ${grpbase}_ncuniq.fa -otus ${grpbase}_otus.fa -relabel Otu
echo -e "${GREEN}Finished ${RED}CLUSTERING ${GREEN}sequence......${NC}\n\n\n"

echo -e "${GREEN}Starting ${RED}SEARCHING OTUS ${GREEN}sequence......${NC}"
${usearch8} -makeudb_usearch ${grpbase}_otus.fa -output ${grpbase}_otus.udb 
${usearch8} -usearch_global ${grpbase}_merged.fq -db ${grpbase}_otus.udb -strand both -id 0.97 -uc ${grpbase}_otutab.uc
echo -e "${GREEN}Finished ${RED}SEARCHING OTUS ${GREEN}sequence......${NC}\n\n\n"

# taxonomy
echo -e "${GREEN}Starting ${RED}UTAXING ${GREEN}sequence......${NC}"
##download rdp_v16.fa.tar.gz
${usearch9} -makeudb_utax rdp_v16.fa -output rdp_v16.udb -taxconfsin taxconfs/500.tc
${usearch9} -utax ${grpbase}_otus.fa -db rdp_v16.udb -utaxout ${grpbase}_otus.utax -strand both
echo -e "${GREEN}Finished ${RED}UTAXING ${GREEN}sequence......${NC}\n\n\n"

# phylogeny
echo -e "${GREEN}Starting ${RED}PHYLOGENY ${GREEN}sequence......${NC}"
mafft --thread 44 --adjustdirectionaccurately ${grpbase}_otus.fa > ${grpbase}_otus.mafft.fa
echo -e "${GREEN}Finished ${RED}PHYLOGENY ${GREEN}sequence......${NC}\n\n\n"
