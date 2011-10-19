How it works:
-------------

blah blah blah... need to write a blerb for this. 

Basic Setup:
------------

1.) Install mrjob. `$ easy_install mrjob` should do the trick. If you want to run the scripts locally you'll also need to install numpy.

    I haven't written any unittests, but you can 'test' the 'alignments to oneliner scripts' with the following commands:

          $ python nexus2oneliner.py -i practice_alignments/nexus_primates/
          $ python phylip2oneliner.py -i practice_alignments/phylip_primates/

    Both commands should print out three long lines of DNA sequences without any error messages. 

2.) Follow the [instructions][1] for setting up mrjob. You'll need to make a mrjob.conf file. Mine looks something like this:
    
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

3.) Run the test analyses described below in the **Species Tree Estimation** and **Bootstrapping** sections

   
Species Tree Estimation:
------------------------

Here's the command to run it locally with the practice data:

        python nexus2oneliner.py -i practice_alignments/nexus_primates/ |
        python cloudtree.py \
        --gene-trees \
        --archive=osx.phylo.tar.gz#bin \
        > primate.gene.trees \

And, this will run it on Amazon:

        python nexus2oneliner.py -i practice_alignments/nexus_primates/ |
        python cloudtree.py \
        --num-ec2-instances 2 \
        --jobconf mapred.map.tasks=1 \
        --jobconf mapred.reduce.tasks=1 \
        --jobconf mapred.reduce.tasks.speculative.execution=True \
        --gene-trees \
        --archive=aws.phylo.tar.gz#bin \
        > primate.aws.gene.trees \

Real life example:

    The main things of note here are the extra map tasks and the extended timeout settings. The latter are particularly important if you have longish alignments such as those from exons and/or are using MrAIC to infer models.

        python nexus2oneliner.py -i some/real/data/nexus |
        python cloudtree.py \
        -r emr \
        --num-ec2-instances 250 \
        --jobconf mapred.map.tasks=249 \
        --jobconf mapred.reduce.tasks=1 \
        --jobconf mapred.reduce.tasks.speculative.execution=True \
        --jobconf mapred.task.timeout=18000000 \
        --full-analysis \
        --bootreps 1000 \
        --archive=aws.phylo.tar.gz#bin \
        > genecode-merged.1000-bootreps.trees
        
Bootstrapping:
-------------

Here's the command to run it locally with the practice data:

        python nexus2oneliner.py -i practice_alignments/nexus_primates/ |
        python cloudtree.py \
        --bootreps 5 \
        --archive=osx.phylo.tar.gz#bin \
        > primates.5_bootreps.trees \

Phybase:
--------

Instructions are unfinished.

pick out-group (run --taxa command)
run phybase...


[1]: http://packages.python.org/mrjob/

