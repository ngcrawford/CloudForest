Generating Genetrees from Loci:

1.) run `nexus2oneliner.py -- `

Bootstrapping Instructions:
---------------------------

1.) Install mrjob. `$ easy_install mrjob` should do the trick. If you want to run the scripts locally you'll also need to install numpy. 

2.) Convert your phylip alignments to 'oneliners.' Essentially MapReduce operates one line at a time so you need to convert your alignments to a one line per alignment format. Actually, since you want to bootstrap the entire dataset you need to convert all your alignments into a single oneline format. The following command will do this for you:

        cat /phylip_folder/*.phylip | ./prep_aligns.py > oneliners.txt

3.) To submit the job to Amazon's Elastic MapReduce the following command will spin up 15 ec2 instances and generate 5 bootstrap replicates.

        python aws_mrjob_boostrap.py -v -r emr --num-ec2-instances 15 \
        --bootreps 5 --archive=aws.phylo.tar.gz#bin \
        < oneliners.txt > trees_out

        python cloudtree.py \
        --bootreps 5 \
        --archive=osx.phylo.tar.gz#bin \
        < reptiles-extended-with-outgroup.oneliners \
        > aws.reptile.trees \

        python cloudtree.py \
        -r emr \
        --num-ec2-instances 5 \
        --jobconf mapred.map.tasks=4 \
        --jobconf mapred.reduce.tasks=1 \
        --jobconf mapred.reduce.tasks.speculative.execution=True \
        --full-analysis \
        --bootreps 2 \
        --archive=aws.phylo.tar.gz#bin \
        < practice_alignments/3.align.oneliners.txt \
        > step4.out

Species Tree Estimation:
-----------------------



