path=${path_to_smrtlink}/smrtlink/userdata/jobs_root/cromwell-executions/pb_basemods
cd $path
ls */call-gather_basemods_gff/execution/basemods.gff |cut -f 1 -d "/" > basemods.list
while read -r line
do
    cp $line/call-gather_basemods_gff/execution/basemods.gff ${work_dir}/01_Modified_gff/$line.basemods.gff
done < "basemods.list"
cat 01_Modified_gff/*.basemods.gff |grep -v "##"|grep "m4C" |cut -f 1 |sort|uniq >> m4C.list
cat 01_Modified_gff/*.basemods.gff |grep -v "##"|grep "m6A" |cut -f 1 |sort|uniq >> m6A.list
cat 01_Modified_gff/*.basemods.gff |grep -v "##"|grep "modified_base" |cut -f 1 |sort|uniq >> modified.list
type=m4C
type=m6A

cat 01_Modified_gff/*.basemods.gff |grep $type | awk 'BEGIN {FS="[\t=]"} {print ">" $1 "_" $4 "\n" $11 }' |sed 's/;IPDRatio//g' > $type.fa
seqkit rmdup --by-seq ./$type.fa  -j 16|sed 's/_/_n_/g'> ./$type.rmdup.fa
rm -r $type.rmdup.fa.split
seqkit split -p 6 ./$type.rmdup.fa
for i in `ls $type.rmdup.fa.split`
do
  mkdir -p ./cat_${type}/$i
  ${Software}/meme-5.4.1/bin/meme-chip $type.rmdup.fa.split/$i \
  -meme-p 40 -dna -oc ./cat_${type}/$i  \
  -meme-mod zoops --meme-nmotifs 10000 -minw 3 -maxw 20 
  echo -e "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence"> ./cat_${type}/$i/${type}.fimo.tsv
  cat ./cat_${type}/$i/fimo_out*/fimo.tsv  |grep -v "#" |grep -v "motif_id"|awk 'BEGIN {FS="\t"} ($4<21 && $5 >21)||($4>21 && $5<21)'>> ./cat_${type}/$i/${type}.fimo.tsv
done

mkdir -p 03_cat_fimo
echo -e "motif_id\tmotif_alt_id\tsequence_name\tmotif_site\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence"> 03_cat_fimo/m4C.fimo.tsv
cat 02_motif/*/m4C.fimo.tsv |\
    sed 's/_n_/_/g' |\
    sed 's/^\w*-//g'|\
    grep -v "matched_sequence" |\
    sed -e 's/\(.*\)_\(.*\)/\1\t\2/'  >> 03_cat_fimo/m4C.fimo.tsv

echo -e "motif_id\tmotif_alt_id\tsequence_name\tmotif_site\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence"> 03_cat_fimo/m6A.fimo.tsv
cat 02_motif/*/m6A.fimo.tsv |\
    sed 's/_n_/_/g' |\
    sed 's/^\w*-//g'|\
    grep -v "matched_sequence" |\
    sed -e 's/\(.*\)_\(.*\)/\1\t\2/' >> 03_cat_fimo/m6A.fimo.tsv

cat 01_Modified_gff/*.basemods.gff |\
    grep -v "##"|\
    grep -v "modified_base"|\
    awk 'BEGIN {FS="[\t=]"} {print  $1 "\t" $3"\t" $4  }' |sed 's/;IPDRatio//g'|\
    sort|uniq >05_modified.bases.tsv

grep -wFf HQ.Genome.tsv 05_modified.bases.tsv > 06.selected.modified.bases.tsv
