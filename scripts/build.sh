
result_file="./result/as_fisher_m_hv/output.rms.csv" #change here
rm $result_file
make
echo "Image(W x H) \t Domain Density \tRMS Threshold \t Time Taken(s) \t Entropy Threshold \t Variance Threshold \t Image Entropy \t Image Variance \t Compression Ratio \t BPP \t SSIM \t PSNR " > $result_file

for file in ./pgm/*
do
  for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
  do
  
  	# echo $file
  	outfile=$(echo $file | cut -d'/' -f 3 )
    #echo $outfile
    file_type=$(echo $outfile | cut -d'.' -f 2)
    # echo $file_type
  	outfile=$(echo $outfile | cut -d'.' -f 1)
    #echo $outfile
  	outfile="./com/as_fisher_m_hv/$outfile$i.ifs" #change here
    # echo $outfile

  	./fcom -FM3 -r $i $file $outfile  -d 2 #change here
     # echo $outfile
  	outfile1=$(echo $outfile | cut -d'/' -f 4)
    # echo $outfile1
  	outfile1=$(echo $outfile1 | cut -d'.' -f 1)
     # echho h$outfile1
    outfile1="./dcom/as_fisher_m_hv/$outfile1.dec.$file_type" #change here
     # echo $outfile1
    # echo  $outfile1
  	./fdcom -i $outfile $outfile1 -q 
    ./qm $file $outfile1 > ./result/as_fisher_m_hv/quality.txt #change here
  TEMP=$(cat ./result/as_fisher_m_hv/quality.txt | cut -d':' -f 2) # change here
  SSIM=$(echo $TEMP | cut -d' ' -f 1)
  PSNR=$(echo $TEMP | cut -d' ' -f 2)
  echo "$SSIM \t $PSNR" >> $result_file
  
 
done	

echo "" >> $result_file
done
