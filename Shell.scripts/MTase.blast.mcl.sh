mt=$pj.MTase.faa
cat $mt $bac.MTase.faa > MTase.faa
makeblastdb -in MTase.faa -out db/MTase -dbtype prot
blastp -query MTase.faa -db db/MTase \
    -out 00.MTase.blast \
    -evalue 1e-5 -num_threads 48 \
    -outfmt "6 qseqid sseqid pident length slen qstart qend sstart send evalue bitscore"
awk 'BEGIN {FS="\t"} $1!=$2&&($4/$5)>0.75 {print $1 "\t" $2 "\t" ($3*($4/$5))/100}' 00.MTase.blast > 01.filtered.blast

mcl  01.filtered.blast --abc -o 02.mcl.out
python get_MT.py -i 02.mcl.out -o MT.csv
grep ">" MTase.faa | sed 's/>//g' |awk 'BEGIN {FS="\t"} {print $1}'|awk 'BEGIN {FS=" "} {print $1}' > all.contig.list
awk 'BEGIN {FS=","} {print $2}' MT.csv > MT_contig.list
grep -vw -Ff MT_contig.list all.contig.list > single.list
sed -i 's/\t/_/g' single.list
cat 02.mcl.out single.list > mcl.sin.output
rm MT.csv MT_contig.list single.list all.contig.list
python get_MT.py -i mcl.sin.output -o MT.csv
rm mcl.sin.output

grep  "MGYG" MT.csv|cut -f 1 -d ","|sort|uniq > mt.b
grep -v "MGYG" MT.csv|cut -f 1 -d ","|sort|uniq > mt.v
comm mt.b mt.v -1 -2 > mt.s
wc *.s
comm mt.b mt.v -2 -3|wc
comm mt.b mt.v -1 -3|wc
grep -wFf mt.s MT.csv|grep -v "MGYG"|wc
cut -f 1 MT.hq.csv -d ","|sort|uniq|wc
