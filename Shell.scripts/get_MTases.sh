mkdir -p 03_rpsblast
  domain/REase.pn 
  domain/methyltransferase.pn
  grep "DNA.*methyltransferase" domain.title|cut -f 1 -d : > domain/DNA_MTase.pn

  makeprofiledb -title DNA_MTase -in domain/DNA_MTase.pn -out ${work_dir}/db/domain/DNA_MTase -threshold 9.82 -scale 100.0 -dbtype rps -index true

  db1=/mnt/raid8/sunchuqing/2021_Virome_base_modification/MTase.REase.ana/db/domain/DNA_MTase

  pj=UHGP
  faa=/mnt/raid8/sunchuqing/Data/UHGG/Protein/uhgp-95/uhgp-95_hq.faa
  pj=HGV
  faa=HGV/HGV.faa

  /mnt/raid8/sunchuqing/Softwares/blast2.12.0/ncbi-blast-2.12.0+/bin/rpsblast  \
      -query $faa \
      -outfmt 6 -num_threads 16\
      -evalue 1E-5 -out 03_rpsblast/$pj.DNA.MTase.blast -db $db1

  faa=HGV/HGV.line.faa
  cut -f 1 03_rpsblast/$pj.DNA.MTase.blast|sort|uniq > $pj/$pj.DNA.MTase.list
  seqkit grep -f $pj/$pj.DNA.MTase.list  $faa> $pj/$pj.MTase.faa

  UHGP=/mnt/raid8/sunchuqing/Data/UHGG/Protein/uhgp-95/uhgp-95_hq.line.faa
  cut -f 1 03_rpsblast/UHGP.DNA.MTase.blast|sort|uniq > UHGP/UHGP.DNA.MTase.list
  seqkit grep -f $pj/$pj.DNA.MTase.list  $UHGP> $pj/$pj.MTase.faa
