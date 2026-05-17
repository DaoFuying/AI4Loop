#!/bin/bash
# The input data should have been filtered against blacklisted regions

usage()
{
  echo "$(basename "$0") [-h] INTERS GENE TFPEAKS NAME DATADIR"
  echo "-- Progam to preprocess the interactions and generate negative samples."

  echo "where:"
  echo "-h           show this help text"
  echo "INTERS       Interaction file in BEDPE format"
  echo "GENE        Human genes in BED format"
  echo "DIST         The distance to merge the anchors"
  echo "NAME         The prefix for the sample."
  echo "DATADIR      Location of the output directory."
}


if [ "$1" != "" ]; then
    case $1 in
        -h | --help )           usage
                                exit
                                ;;
    esac
fi

if [ $# -lt 5 ]; then
  usage
  exit
fi

infile=${1}
gene=${2}
dist=${3}
name=${4}
datadir=${5}

dir=$(dirname "$0")

merged_anchors=${datadir}/${name}_merged_anchors.bed
merged_anchors_both_gene=${datadir}/${name}_merged_anchors.both_gene.bed
cluster_inters=${datadir}/${name}.clustered_interactions.bedpe
cluster_inters_both_gene=${datadir}/${name}.clustered_interactions.both_gene.bedpe

cat ${infile} \
    | awk 'BEGIN{FS=OFS="\t"}{printf("%s\t%s\t%s\n%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6)}' \
    | sort -k1,1 -k2,2n -k3,3n \
    | mergeBed -d ${dist}  > ${merged_anchors}

pairToBed -a ${infile} -b ${merged_anchors} -type both \
    | python ${dir}/groupBed.py -g 1 2 3 4 5 6 7 -c 8 9 10 -o collapse  \
    | cut -f 7- \
    | sed 's/,/ /g' \
    | awk 'BEGIN{OFS="\t"}!($4==$5 && $6==$7){if($4 < $5){print $2,$4,$6,$3,$5,$7,$1}else{print $3,$5,$7,$2,$4,$6,$1}}' \
    | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n \
    | awk '$1==$4' \
    | python ${dir}/groupBed.py -g 1 2 3 4 5 6 -c 7 \
    | sed 's/\..*$//g' > ${cluster_inters}

pairToBed -a ${cluster_inters} -b ${gene}  -type both \
    | cut -f 1-7 \
    | uniq > ${cluster_inters_both_gene}

cat ${cluster_inters_both_gene} \
    | awk 'BEGIN{FS=OFS="\t"}{printf("%s\t%s\t%s\n%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6)}' \
    | sort -u > ${merged_anchors_both_gene}