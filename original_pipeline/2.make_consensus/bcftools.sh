#! /bin/sh

bcftools mpileup --threads 4 -C 50 -Ob -f $RefSeqFileName $Rajaoutname -o Raja.mpileup.bcf              #makes mipleup, adjusting MQ for reads with lots of mismatches(-C 50)
bcftools norm --threads 4 -d all -Ob -f $RefSeqFileName Raja.mpileup.bcf -o Raja.norm.mpileup.bcf       #left aligns indels (-f ref.fa) and removes duplicates (-d)
bcftools call --threads 4 -mv -P 0.9e-3 -Ob Raja.norm.mpileup.bcf -o Raja.calls.norm.mpileup.bcf        #multiallelic calling of SNPs/INDELs and more stringent calls(-P)
bcftools filter --threads 4 -i 'DP>7' -Ob Raja.calls.norm.mpileup.bcf -o Raja.DP8.calls.norm.mpileup.bcf        #filter for depth of 8+

bcftools mpileup --threads 4 -C 50 -Ou -f $RefSeqFileName $Rajaoutname |\
bcftools norm --threads 4 -d all -Ou -f $RefSeqFileName |\
bcftools call --threads 4 -mv -P 0.9e-3 -Ou |\
bcftools filter --threads 4 -g 3 -i 'DP>7' -Ou -o Raja.pipedoutput3.bcf
