
main_script="/mydata/code/generic_vcf_ml_v2.sh"
test_set="test_set.txt"
train_set="training.txt"
#sen_lev="6"

###prepare test sample
#test_set : 1143 sentieon  vs edico (truthset) for examing the result
sample_test=$(cat $test_set | awk '{print $1}')
sample_test_name=$(echo $sample_test| sed 's/\..*$//g')

sample_truth=$(cat $test_set | awk '{print $2}')
sample_truth_name=$(echo $sample_truth| sed 's/\..*$//g')

echo "test: $sample_test over truth: $sample_truth"
tabix $sample_test
tabix $sample_truth
#make sure use default filter
bcftools filter -s "ECNT3" -e "ECNT=3" $sample_test -Oz > $sample_test_name.test.vcf.gz
tabix $sample_test_name.test.vcf.gz

#remove any unesccary field to reduce errors
bcftools annotate -x "INFO/FUNCOTATION,FORMAT" $sample_truth | bcftools view -f PASS -Oz > $sample_truth_name.truth.vcf.gz
tabix $sample_truth_name.truth.vcf.gz

#prepare test ready vcf
echo '##INFO=<ID=TRUTH,Number=1,Type=String,Description="Variant exist in truthset">' > ann.hdr 
bcftools annotate -a $sample_truth_name.truth.vcf.gz -h ann.hdr -m+'TRUTH=TP' $sample_test_name.test.vcf.gz -Oz > $sample_test_name.tr.vcf.gz

tot_truth_num=$(zcat $sample_truth_name.truth.vcf.gz | grep -v '^#' | wc -l)
tot_test_num=$(zcat $sample_test | grep -v '^#' | wc -l)
tot_catch_num=$(zcat $sample_test_name.tr.vcf.gz | grep "TP" | wc -l)
tot_miss_num=$(echo "$tot_truth_num - $tot_catch_num" | bc)

# training.txt training set with truthset
echo "model,iMeIeT,iMiIeT,iMiIiT,iMeIiT,eMiIeT,eMiIiT,eMeIiT,tot_test_num,tot_miss_num,M_TP,M_FP,M_TN,M_FN" > ML_summary.csv

for sen_lev in $(echo "2 3 4 5 6"); do
    echo "### sensitivity level: $sen_lev"
    i=0
    while read sample_train; do 
        i=$(expr $i + 1)
        echo $i $sample_train
        if [ $i -eq 1 ]; then 
            bash $main_script train $sample_train $sen_lev
        else
            bash $main_script train $sample_train $sen_lev $lat_model
        fi 
        pre_model=$(ls -t *.sav | head -n1)
        mv $pre_model vcf_rf_model.sen_${sen_lev}_vcf_${i}.sav
        lat_model="vcf_rf_model.sen_${sen_lev}_vcf_${i}.sav"
        #use trained model on test sample 
        bash $main_script pred $sample_test_name.tr.vcf.gz $lat_model
        
        lat_vcf=$(ls -t *.ml.vcf.gz | head -n1)
        iMeIeT=$(zcat $lat_vcf | grep -v "^#" | grep "ML_PASS" | grep -vw "PASS" | grep -v "TP" |wc -l)
        iMiIeT=$(zcat $lat_vcf | grep -v "^#" | grep "ML_PASS" | grep -w "PASS" | grep -v "TP" |wc -l)
        iMiIiT=$(zcat $lat_vcf | grep -v "^#" | grep "ML_PASS" | grep -w "PASS" | grep "TP" |wc -l)
        iMeIiT=$(zcat $lat_vcf | grep -v "^#" | grep "ML_PASS" | grep -vw "PASS" | grep "TP" |wc -l)
        eMiIeT=$(zcat $lat_vcf | grep -v "^#" | grep -v "ML_PASS" | grep -w "PASS" | grep -v "TP" |wc -l)
        eMiIiT=$(zcat $lat_vcf | grep -v "^#" | grep -v "ML_PASS" | grep -w "PASS" | grep "TP" |wc -l)
        eMeIiT=$(zcat $lat_vcf | grep -v "^#" | grep -v "ML_PASS" | grep -vw "PASS" | grep "TP" |wc -l)
        #tot_truth_num
        #tot_test_num
        M_TP=$(echo "$iMiIiT + $iMeIiT" | bc)
        M_FP=$(echo "$iMiIeT + $iMeIeT" | bc)    
        M_FN=$(echo "$eMiIiT + $eMeIiT + $tot_miss_num" | bc)
        M_TN=$(echo "$tot_test_num - $M_TP - $M_FP - $M_FN" | bc)
        
        echo "sen_${sen_lev}_vcf_${i},$iMeIeT,$iMiIeT,$iMiIiT,$iMeIiT,$eMiIeT,$eMiIiT,$eMeIiT,$tot_test_num,$tot_miss_num,$M_TP,$M_FP,$M_TN,$M_FN" >> ML_summary.csv

    done < $train_set
    rm *pred*
done

