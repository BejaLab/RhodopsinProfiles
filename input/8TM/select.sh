ls -l /data/uniprot/2021_02/uniref50.fasta
awk '{print"UniRef50_"$1}' selected.txt | xargs seqkit faidx -f /data/uniprot/2021_02/uniref50.fasta | seqkit replace -p UniRef50_ -o references.fasta
