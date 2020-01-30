
#prepare sample for machine learning, with truthset target for training (only bialleic mode), wihout for testing (split multialleic)
# training mode with a truthset vcf
# bash generic_vcf_ml.sh train input.vcf.gz target.vcf.gz sensitivity_set
# training: incremental
# bash generic_vcf_ml.sh train input.vcf.gz target.vcf.gz sensitivity_set ML_model.date.sav
# sensitivity_set = 1 - 6(higher)

# prepare sample only
# bash generic_vcf_ml.sh prep input.vcf.gz

# prediction mode with a ml model
# bash generic_vcf_ml.sh pred input.vcf.gz ML_model.date.sav

#list="18-112 18-136 18-179 18-269 18-270"
#for i in $list; do echo $i; bash ../sample_prepare_v3.sh *$i*_001_*twicefiltered.vcf.gz *$i*_R1_*.vcf.gz; done

ml_pred_py="ml_vcf_rf_pred_v2.py"
ml_train_py="ml_vcf_rf_train_v2.py"

tgrepm(){
    #select by col_name
    #tgrepm $input $col_name $sep $choice_of_line $grep_option
    #tgrepm test.tsv "a|31" , "1,2" -E
    if [ $3 == "\s" ];then 
        sep=" "
    else 
        sep=$(echo -e $3)
    fi
    
    if [ $4 == "0" ]; then 
        row_num=$(cat $1 | wc -l)
        choice_line=""
        for((i=1; i<=$row_num; i++)); do choice_line=$(echo $choice_line $i);done 
    else
        choice_line=$(echo $4 | sed 's/,/ /g')
    fi

    line_num=""
    for i in $choice_line;do
    
        col_num=$(sed -n ${i}p $1 | sed "s/${sep}/\n/g" | grep $5 $2 -n | sed 's/\:.*$//g')
        line_num=$(echo $line_num $col_num)
        
    done
    chk=$(echo $line_num |wc -c)
    if [ $chk -gt 1 ]; then
        cat $1 | cut -f $(echo $line_num | sed 's/ /,/g') -d "$sep"
    fi
}

vcf2csv(){
    ## produce table of all content:# vcf2csv abc.vcf.gz 
    ## prepare table for ML, integer/float only:# vcf2csv abc.vcf.gz ML
    vcf_file=$1
    add_req=$2
    add_com=$3
    #echo $vcf_file $add_req
    vcf_name=$(echo $vcf_file | sed 's/\.*$//g')
    rm -f *temp.txt
    #Head_COMMON="%CHROM %POS %ID %REF %ALT %QUAL %FILTER"

    Head_COMMON=$(zcat $vcf_file | grep "^#CHROM" | sed 's/\tINFO.*$//g' | sed 's/^#/\%/g' | sed 's/\t/\ \%/g')

    Head_COMMON_num=$(echo $Head_COMMON | awk -F '%' '{print NF-1}')
    #echo $Head_COMMON
    HEAD_NAME=$(bcftools query -l $vcf_file )

    #if [ $add_req == "ML" ]; then
    if [ -z $add_req ]; then
        add_req=""
        #covert all content for csv
        Head_INFO=$(echo $(zcat $vcf_file | grep "^##INFO" |sed 's/^.*\<ID=/\%INFO\//g'|sed 's/\,.*$//g'))
        #%INFO/AC %INFO/AF %INFO/AN %INFO/BaseQRankSum %INFO/ClippingRankSum %INFO/DB %INFO/DP %INFO/END %INFO/ExcessHet %INFO/FS %INFO/InbreedingCoeff %INFO/MLEAC %INFO/MLEAF %INFO/MQ %INFO/MQRankSum %INFO/QD %INFO/RAW_MQ %INFO/ReadPosRankSum %INFO/SOR

        Head_FORMAT=$(echo $(zcat $vcf_file |grep "^##FORMAT" |sed 's/^.*\<ID=/\[\ \%/g'|sed 's/\,.*$/]/g')| sed 's/\] /]/g')
        #[ %AD][ %DP][ %GQ][ %GT][ %MIN_DP][ %PGT][ %PID][ %PL][ %RGQ][ %SB]
        zcat $vcf_file|grep "^##FORMAT" |sed 's/^.*\<ID=/\%/g'|sed 's/\,.*$//g' > temp.txt

        touch all-temp.txt
        #zcat $vcf_file|grep "^##FORMAT" |sed 's/^.*\<ID=/\%/g'|sed 's/\,.*$//g' > temp.txt
        for HEAD_NAME_TEMP in $HEAD_NAME; do
            #echo $HEAD_NAME_TEMP
            sed "s/^/$HEAD_NAME_TEMP/g" temp.txt > $HEAD_NAME_TEMP-temp.txt
            cp all-temp.txt all1-temp.txt
            paste all1-temp.txt $HEAD_NAME_TEMP-temp.txt > all-temp.txt
        done
        #echo $Head_COMMON $Head_INFO$(cat all-temp.txt) | sed 's/ /\,/g' > $vcf_name.${add_req}vcf2csv.temp.csv
        echo $Head_COMMON $Head_INFO$(cat all-temp.txt) | sed 's/ /\,/g'
        rm *temp.txt
        #bcftools query -f "$Head_COMMON $Head_INFO$Head_FORMAT\n" $vcf_file | sed 's/\,/|/g' | sed 's/ /\,/g' >> $vcf_name.${add_req}vcf2csv.temp.csv
        bcftools query -f "$Head_COMMON $Head_INFO$Head_FORMAT\n" $vcf_file | sed 's/\,/|/g' | sed 's/ /\,/g' 
    else
        add_req="ML_"
        #prepare matrics table for ML 
        Head_INFO=$(echo $(zcat $vcf_file | grep -Ew "Number=1|Number=A" | grep -Ew "Type=Integer|Type=Float" | grep "^##INFO" |sed 's/^.*\<ID=/\%INFO\//g'|sed 's/\,.*$//g'))
        #Head_INFO=$(echo $(zcat $vcf_file | grep -Ew "Number=1|Number=A" | grep -Ew "Type=Integer|Type=Float" | grep "^##INFO" |sed 's/^.*\<ID=/\%INFO\//g'|sed 's/\,.*$//g'))
        #%INFO/AC %INFO/AF %INFO/AN %INFO/BaseQRankSum %INFO/ClippingRankSum %INFO/DB %INFO/DP %INFO/END %INFO/ExcessHet %INFO/FS %INFO/InbreedingCoeff %INFO/MLEAC %INFO/MLEAF %INFO/MQ %INFO/MQRankSum %INFO/QD %INFO/RAW_MQ %INFO/ReadPosRankSum %INFO/SOR

        #Head_FORMAT=$(echo $(zcat $vcf_file | grep -Ew "Type=Integer|Type=Float" | grep "^##FORMAT" |sed 's/^.*\<ID=/\[\ \%/g'|sed 's/\,.*$/]/g')| sed 's/\] /]/g')
        Head_FORMAT=$(echo $(zcat $vcf_file | grep -Ew "Number=1|Number=A" | grep -Ew "Type=Integer|Type=Float" | grep "^##FORMAT" |sed 's/^.*\<ID=/\[\ \%/g'|sed 's/\,.*$/]/g')| sed 's/\] /]/g')
        #Head_FORMAT_m=$(echo $(zcat $vcf_file | grep -Ew "Number=R|Number=G" | grep -Ew "Type=Integer|Type=Float" | grep "^##FORMAT" |sed 's/^.*\<ID=/\[\ \%/g'|sed 's/\,.*$/]/g')| sed 's/\] /]/g')
        Head_FORMAT_m=$(echo $(zcat $vcf_file | grep -Ew "Number=R" | grep -Ew "Type=Integer|Type=Float" | grep "^##FORMAT" |sed 's/^.*\<ID=/\[\ \%/g'|sed 's/\,.*$/]/g')| sed 's/\] /]/g')

        zcat $vcf_file | grep -Ew "Number=1|Number=A" | grep -Ew "Type=Integer|Type=Float" | grep "^##FORMAT" |sed 's/^.*\<ID=/\%/g'|sed 's/\,.*$//g' > temp.txt
        #zcat $vcf_file | grep -Ew "Number=R|Number=G" | grep -Ew "Type=Integer|Type=Float" | grep "^##FORMAT" |sed 's/^.*\<ID=/\%/g'|sed 's/\,.*$//g' |awk '{print $1"_R",$1"_A"}' > s_temp.txt
        zcat $vcf_file | grep -Ew "Number=R" | grep -Ew "Type=Integer|Type=Float" | grep "^##FORMAT" |sed 's/^.*\<ID=/\%/g'|sed 's/\,.*$//g' |awk '{print $1"_R",$1"_A"}' > s_temp.txt
        #paste -c " " s_temp.txt s_temp.txt > s2_temp.txt
        
        touch all-temp.txt
        touch all-2temp.txt
        #zcat $vcf_file|grep "^##FORMAT" |sed 's/^.*\<ID=/\%/g'|sed 's/\,.*$//g' > temp.txt
        i=0
        for HEAD_NAME_TEMP in $HEAD_NAME; do
            #echo $HEAD_NAME_TEMP
            i=$(expr $i + 1)
            sed "s/^/$i/g" temp.txt > $HEAD_NAME_TEMP-temp.txt
            cp all-temp.txt all1-temp.txt
            paste all1-temp.txt $HEAD_NAME_TEMP-temp.txt > all-temp.txt
            
            #sed "s/^/$HEAD_NAME_TEMP/g" s_temp.txt > $HEAD_NAME_TEMP-2temp.txt
            awk -v var=$i '{print var$1,var$2}' s_temp.txt > $HEAD_NAME_TEMP-2temp.txt
            cp all-2temp.txt all1-2temp.txt
            paste all1-2temp.txt $HEAD_NAME_TEMP-2temp.txt > all-2temp.txt

        done

        if [ -z $add_com ]; then
            #echo $Head_COMMON $Head_INFO$(cat all-temp.txt)$(cat all-2temp.txt) | sed 's/ /\,/g' > $vcf_name.${add_req}vcf2csv.temp.csv
            echo $Head_COMMON $Head_INFO$(cat all-temp.txt)$(cat all-2temp.txt) | sed 's/ /\,/g' 
            rm *temp.txt
            #bcftools query -f "$Head_COMMON $Head_INFO$Head_FORMAT$Head_FORMAT_m\n" $vcf_file | sed 's/\,/ /g' | sed 's/ /\,/g' >> $vcf_name.${add_req}vcf2csv.temp.csv
            bcftools query -f "$Head_COMMON $Head_INFO$Head_FORMAT$Head_FORMAT_m\n" $vcf_file | sed 's/\,/ /g' | sed 's/ /\,/g' 
        else
            #include special column at the right end of table
            echo $Head_COMMON $Head_INFO$(cat all-temp.txt)$(cat all-2temp.txt) $add_com | sed 's/ /\,/g' 
            rm *temp.txt
            bcftools query -f "$Head_COMMON $Head_INFO$Head_FORMAT$Head_FORMAT_m $add_com\n" $vcf_file | sed 's/\,/ /g' | sed 's/ /\,/g' 
        fi
    fi
    #cat $vcf_name.${add_req}vcf2csv.temp.csv
    #rm *temp.csv
}

mode=$1
target=$3
input=$2
sen_set=$4
model=$5
input_name=$(echo $input | sed 's/\.vcf\.gz//g')

if [ $mode == "prep" ]; then
    echo "#Preparation mode: prep input.vcf.gz"
    if [ -z $target ]; then
        echo "treat as testing set: input: $input"
        #split multialliec call into multiple line
        bcftools norm -m- $input -Oz > $input_name.norm.vcf.gz

        vcf2csv $input_name.norm.vcf.gz ML | sed 's/\,\.$/\,0/g' | sed 's/\,TP/\,1/g' > $input_name.test_matrix.csv
        #echo "CHROM,POS,ID,REF,ALT,QUAL,FILTER,DP,ECNT,NLOD,POP_AF,P_CONTAM,P_GERMLINE,TLOD,N_AD_R,N_AD_A,T_AD_R,T_AD_A,N_AF,T_AF,N_MBQ,T_MBQ,N_MFRL_R,N_MFRL_A,T_MFRL_R,T_MFRL_A,N_MMQ,T_MMQ,N_MPOS,T_MPOS" > $input_name.test_matrix.csv
        #bcftools query -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER %INFO/DP %INFO/ECNT %INFO/NLOD %INFO/POP_AF %INFO/P_CONTAM %INFO/P_GERMLINE %INFO/TLOD[ %AD][ %AF][ %MBQ][ %MFRL][ %MMQ][ %MPOS]\n'  $input_name.norm.vcf.gz \
        #| sed 's/\,/ /g' | sed 's/ /\,/g' | sed 's/\,\.$/\,0/g' | sed 's/\,TP/\,1/g' >>  $input_name.test_matrix.csv
        output=$input_name.test_matrix

        #remove unnesssary columns to reduce burden of ML training
        #selcol $output.csv .*_R , -v | sed 's/\,nan/\,0/g' > ${output}_cor.csv
        #tgrepm $output.csv .*_R , 1 -v | sed 's/\,nan/\,0/g'| sed 's/\,\./\,0/g' > ${output}_cor.csv
        cat $output.csv| sed 's/\,nan/\,0/g'| sed 's/\,\./\,0/g' > ${output}_cor.csv
    else
        echo "ERROR: multiple inputs for mode: prep"
    fi
elif [ $mode == "pred" ]; then
    echo "#Prediction mode: pred input.vcf.gz ML_model.date.sav"
    if [ -z $target ]; then
        echo " ERROR: no input model for mode: pred "
    else
        chk_target=$(echo $target | grep ".vcf.gz$" | wc -l)
        if [ $chk_target -eq 0 ]; then
            ml_model=$target
            echo "treat as testing set for prediction: input: $input; predict with ml model: $ml_model"
            #split multialliec call into multiple line
            bcftools norm -m- $input -Oz > $input_name.norm.vcf.gz

            vcf2csv $input_name.norm.vcf.gz ML | sed 's/\,\.$/\,0/g' | sed 's/\,TP/\,1/g'> $input_name.test_matrix.csv
            #echo "CHROM,POS,ID,REF,ALT,QUAL,FILTER,DP,ECNT,NLOD,POP_AF,P_CONTAM,P_GERMLINE,TLOD,N_AD_R,N_AD_A,T_AD_R,T_AD_A,N_AF,T_AF,N_MBQ,T_MBQ,N_MFRL_R,N_MFRL_A,T_MFRL_R,T_MFRL_A,N_MMQ,T_MMQ,N_MPOS,T_MPOS" > $input_name.test_matrix.csv
            #bcftools query -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER %INFO/DP %INFO/ECNT %INFO/NLOD %INFO/POP_AF %INFO/P_CONTAM %INFO/P_GERMLINE %INFO/TLOD[ %AD][ %AF][ %MBQ][ %MFRL][ %MMQ][ %MPOS]\n'  $input_name.norm.vcf.gz \
            #| sed 's/\,/ /g' | sed 's/ /\,/g' | sed 's/\,\.$/\,0/g' | sed 's/\,TP/\,1/g' >>  $input_name.test_matrix.csv
            output=$input_name.test_matrix

            #remove unnesssary columns to reduce burden of ML training
            #selcol $output.csv .*_R , -v | sed 's/\,nan/\,0/g' > ${output}_cor.csv
            #tgrepm $output.csv .*_R , 1 -v | sed 's/\,nan/\,0/g'| sed 's/\,\./\,0/g' > ${output}_cor.csv
            cat $output.csv| sed 's/\,nan/\,0/g'| sed 's/\,\./\,0/g' > ${output}_cor.csv
            #echo "predict with ml model: $ml_model"
            
            python $(dirname $0)/$ml_pred_py ${output}_cor.csv $ml_model
            
            bedfile=$(ls *$output*.pred.bed | tail -n1)
            bgzip -f $bedfile 
            tabix $bedfile.gz
            #echo '##INFO=<ID=ML_FILTER,Number=1,Type=String,Description="Pass ML filter">' > ann.hdr
            echo '##INFO=<ID=ML_PASS,Number=0,Type=Flag,Description="Pass ML filter">' > ann.hdr
            #bcftools annotate -a $bedfile.gz -h ann.hdr -c CHROM,FROM,TO,TAG -m+'ML_FILTER=PASS' $input_name.norm.vcf.gz \
            bcftools annotate -a $bedfile.gz -h ann.hdr -c CHROM,FROM,TO,ML_PASS $input_name.norm.vcf.gz \
            | bcftools norm -m+any - -Oz > $input_name.ml.vcf.gz
        else 
            echo "ERROR: no input ML model for prediction"
        fi
    fi
elif [ $mode == "train" ]; then
    echo "#Training mode: train input.vcf.gz target.vcf.gz sensitivity_set[1,2,3,4] (ML_model.date.sav)"
    if [ -z $target ]; then
        echo " ERROR: no target truthset for mode: train"
    else
        echo " treat as training set: truthset: $target; input: $input"
        # remove funcotation  and format
        bcftools annotate -x "INFO/FUNCOTATION,FORMAT" $target | bcftools view -f PASS -Oz > target.vcf.gz
    
        #keep bialleic sites
        bcftools view -m2 -M2 $input -Oz > input.vcf.gz
        bcftools index target.vcf.gz
        bcftools index input.vcf.gz

        echo '##INFO=<ID=TRUTH,Number=1,Type=String,Description="Variant exist in truthset">' > ann.hdr
        bcftools annotate -a target.vcf.gz -h ann.hdr -m+'TRUTH=TP' input.vcf.gz -Oz > $input_name.ann.vcf.gz

        vcf2csv $input_name.ann.vcf.gz ML "%INFO/TRUTH" | sed 's/\,\.$/\,0/g' | sed 's/\,TP/\,1/g' > $input_name.train_matrix.csv   
        #echo "CHROM,POS,ID,REF,ALT,QUAL,FILTER,DP,ECNT,NLOD,POP_AF,P_CONTAM,P_GERMLINE,TLOD,N_AD_R,N_AD_A,T_AD_R,T_AD_A,N_AF,T_AF,N_MBQ,T_MBQ,N_MFRL_R,N_MFRL_A,T_MFRL_R,T_MFRL_A,N_MMQ,T_MMQ,N_MPOS,T_MPOS,TRUTH" > $input_name.train_matrix.csv
        #bcftools query -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER %INFO/DP %INFO/ECNT %INFO/NLOD %INFO/POP_AF %INFO/P_CONTAM %INFO/P_GERMLINE %INFO/TLOD[ %AD][ %AF][ %MBQ][ %MFRL][ %MMQ][ %MPOS] %INFO/TRUTH\n'  $input_name.ann.vcf.gz \
        #| sed 's/\,/ /g' | sed 's/ /\,/g' | sed 's/\,\.$/\,0/g' | sed 's/\,TP/\,1/g' >>  $input_name.train_matrix.csv
        
        output=$input_name.train_matrix
        #remove unnesssary columns to reduce burden of ML training
        #selcol $output.csv .*_R , -v | sed 's/\,nan/\,0/g' > ${output}_cor.csv
        #tgrepm $output.csv .*_R , 1 -v | sed 's/\,nan/\,0/g' | sed 's/\,\./\,0/g' > ${output}_cor.csv
        cat $output.csv| sed 's/\,nan/\,0/g'| sed 's/\,\./\,0/g' > ${output}_cor.csv
        #convert sensitvity level to ml setting
        senset_t=$(echo "0 16 8 4 2 1" | awk -v var=$sen_set '{if (var>6) var=6 ;print $(var)}')
        if [ -z $model ]; then
            echo "start training mode without previous model at sensitivity level $sen_set:$senset_t"
            python $(dirname $0)/$ml_train_py ${output}_cor.csv 0.1 $senset_t            
        else
            echo "start training mode with previous model: $model at sensitivity level $sen_set:$senset_t"
            ml_model=$model
            python $(dirname $0)/$ml_train_py ${output}_cor.csv 0.1 $senset_t $ml_model
        fi
    fi
else
    echo "ERROR: wrong mode; please select from train, prep, pred"
fi

rm -f target.vcf.gz target.vcf.gz.csi input.vcf.gz input.vcf.gz.csi $output.csv ann.hdr

