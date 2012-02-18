How it works:
-------------

1. This first line reads in a directory of nexus files files and emits oneliners. 

		python nexus2oneliner.py -i test/alignments/nexus_primates/
	
2. Then step one piped directly into `cloudforest.py` with the pipe symbol (= '|').

		python nexus2oneliner.py -i test/alignments/nexus_primates/ | 
		python cloudforest.py ...

3. Cloudforest has two main functions. You can use it to infer genetrees by setting the `--gene-trees` flag or you can use it to bootstrap your data as described in [CITE] by setting the `--bootreps=#` flag. With the default setting cloudforest calculates  `--mraic`


Basic Installation:
------------------

1. Install [setuptools easy_install][2]

1. Install [Git][5]

1. Install [mrjob][3]. `$ easy_install mrjob` should do the trick. If you want to run the scripts locally you'll also need to install numpy. I recommend using [Enthought Python][4] for this, but you can also use `$ easy_install numpy`. 

1. Clone this repository and cd to phylo folder.

        git clone git@github.com:ngcrawford/BioAWS.git
        cd BioAWS/phylo

I haven't written any unittests, but you can 'test' the 'alignments to oneliner scripts' with the following commands in the phylo folder:

        python nexus2oneliner.py -i practice_alignments/nexus_primates/
        python phylip2oneliner.py -i practice_alignments/phylip_primates/

Both commands should print out three long lines of DNA sequences without any error messages. 

1. Follow the [instructions][1] for setting up mrjob. You'll need to make a mrjob.conf file. Mine looks something like this:
    
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

1. Run the test analyses described below in the **Species Tree Estimation** and **Bootstrapping** sections

   
Species Tree Estimation:
------------------------

This first command will run cloudforest locally with the practice data. Setting `--protocol=raw_value` removes 'null' keys from the final file. `--setup-cmd 'mkdir -p tmp'` creates the appropriate `tmp/` directory for storing the phyml output files. This setting is only necessary for running local jobs. `tmp/` is present by default in on Amazon Elastic MapReduce Clusters. 


		python nexus2oneliner.py -i test/alignments/nexus_primates/ | 
		python cloudforest.py \
		--protocol=raw_value \
		--setup-cmd 'mkdir -p tmp' \
		--gene-trees \
		--archive=../gzips/osx.phylo.tar.gz#bin \
		> primate.gene.trees
		
		
Setting the `--mraic` flag will makes cloudforest.py use mraic to infer the appropriate evolutionary model. 		

		python nexus2oneliner.py -i test/alignments/nexus_primates/ |
		python process.py \
		--setup-cmd 'mkdir -p tmp' \
		--gene-trees \
		--mraic \
		--archive=../gzips/osx.phylo.tar.gz#bin \
		> primate.gene.trees


And, this will run it on Amazon:

        python nexus2oneliner.py -i test/alignments/nexus_primates/ |
        python process.py \
        --num-ec2-instances 2 \
        --jobconf mapred.map.tasks=1 \
        --jobconf mapred.reduce.tasks=1 \
        --jobconf mapred.reduce.tasks.speculative.execution=True \
		--protocol=raw_value \
        --gene-trees \
        --archive=../gzips/aws.phylo.tar.gz#bin \
        > primate.aws.gene.trees \

Real life example:

    The main things of note here are the extra map tasks and the extended timeout settings. 
	The latter are particularly important if you have longish alignments such as those from
	exons and/or are using MrAIC to infer models.

		python nexus2oneliner.py -i some/real/data/nexus |
		python process.py \
		-r emr \
		--num-ec2-instances 250 \
		--jobconf mapred.map.tasks=249 \
		--jobconf mapred.reduce.tasks=1 \
		--jobconf mapred.reduce.tasks.speculative.execution=True \
		--jobconf mapred.task.timeout=18000000 \
		--gene-trees \
		--archive=../gzips/aws.phylo.tar.gz#bin \
		> lotsa.gene.trees
        
Bootstrapping:
-------------

Here's the command to run it locally with the practice data:

		python nexus2oneliner.py -i test/alignments/nexus_primates/ |
		python process.py \
		--setup-cmd 'mkdir -p tmp' \
		--full-analysis \
		--bootreps 5 \
		--archive=../gzips/osx.phylo.tar.gz#bin  \
		> primates.5bootreps.trees 

And, this will run it on Amazon:

		python nexus2oneliner.py -i test/alignments/nexus_primates/ |
		python cloudforest.py \
		-r emr \
		--num-ec2-instances 2 \
		--jobconf mapred.map.tasks=1 \
		--jobconf mapred.reduce.tasks=1 \
		--jobconf mapred.reduce.tasks.speculative.execution=True \
		--jobconf mapred.task.timeout=18000000 \
		--full-analysis \
		--bootreps 5 \
		--archive=../gzips/osx.phylo.tar.gz#bin  \
		> primates.5bootreps.aws.trees 


Phybase:
--------

Instructions are unfinished.

[1]: http://packages.python.org/mrjob/
[2]: http://pypi.python.org/pypi/setuptools
[3]: https://github.com/Yelp/mrjob
[4]: http://www.enthought.com/products/epd.php
[5]: http://help.github.com/mac-set-up-git/