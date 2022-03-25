time docker run --rm -it -v $PWD:/mnt -w /mnt mhoelzer/hmmscan:0.1 hmmscan --notextw -o assembly.hmmer.out --domtblout assembly.hmmer.domtblout --tblout assembly.hmmer.tblout AntiFam.hmm assembly.faa
#pythonw hmmer.py --fasta assembly.faa --hmm AntiFam.hmm --output temp
