docker run --rm -v $PWD:/mnt -w /mnt ncbi/blast makeblastdb -dbtype prot -in plasmidProt.faa -input_type fasta -title plasmidProt
docker run --rm -v $PWD:/mnt -w /mnt ncbi/blast makeblastdb -dbtype nucl -in plasmid_originOfReplication.nr.fasta -input_type fasta -title plasmidORIG
