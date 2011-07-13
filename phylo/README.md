Elastic MapReduce:

Amazon webserver


Converting Alignments to 'Oneliners':

    


1.) run `nexus2oneliner.py -- `





Bootstrapping Instructions:
---------------------------

1.) Install mrjob. `$ easy_install mrjob` should do the trick. If you want to run the scripts locally you'll also need to install numpy. 


    runners:
      emr:
        aws_access_key_id: YOURIDYOURIDYOURIDYOURIDYOURID
        aws_region: us-east-1
        aws_secret_access_key: YOURSECRETKEYYOURSECRETKEYYOURSECRETKEY
        ec2_key_pair: mr-keypair
        ec2_key_pair_file: /Users/YourUserName/.ssh/mr-keypair.pem
        ec2_instance_type: m1.small
        num_ec2_instances: 1
        setup_cmds: &setup_cmds
        ssh_tunnel_is_open: true
        ssh_tunnel_to_job_tracker: true

2.) Convert your phylip alignments to 'oneliners.' Essentially MapReduce operates one line at a time so you need to convert your alignments to a one line per alignment format. Actually, since you want to bootstrap the entire dataset you need to convert all your alignments into a single oneline format. The following command will do this for you:

        cat /phylip_folder/*.phylip | ./prep_aligns.py > oneliners.txt

3.) To submit the job to Amazon's Elastic MapReduce the following command will spin up 15 ec2 instances and generate 5 bootstrap replicates.



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



