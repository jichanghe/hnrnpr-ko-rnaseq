
python DATTSS_main.py -b /mnt/c/Users/User/Desktop/allbamfiles.txt -anno /mnt/c/Users/User/Desktop/hg38_first_exon_annotation_forDATTSS.txt -p 10 -o /mnt/c/Users/User/Desktop/DATTSS_output1.txt


python DATTSS_compare.py -b /mnt/c/Users/User/Desktop/HNRNPR/allbamfiles.txt -anno /mnt/c/Users/User/Desktop/hg38_first_exon_annotation_forDATTSS.txt -p 10 -r /mnt/c/Users/User/Desktop/hg38_GEOCODE_ExonRegion_annotation.txt -d /mnt/c/Users/User/Desktop/HNRNPR/exonCount/ -o /mnt/c/Users/User/Desktop/HNRNPR/DATTSS_output.txt


Rscript Infer_DU_tandemTSS.R \
  -b /mnt/c/Users/User/Desktop/HNRNPR/allbamfiles.txt \
  -I /mnt/c/Users/User/Desktop/HNRNPR/DATTSS_output.txt \
  -d /mnt/c/Users/User/Desktop/HNRNPR/exonCount \
  -o /mnt/c/Users/User/Desktop/HNRNPR/DATTSS_tandem_DU.txt

