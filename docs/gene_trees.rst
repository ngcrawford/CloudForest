
Basic Concept:
==============

MORE HERE.....

1. This first line reads in a directory of nexus files files and emits oneliners. 

		python nexus2oneliner.py -i test/alignments/nexus_primates/
	
2. Then step one piped directly into `cloudforest.py` with the pipe symbol (= '|').

		python nexus2oneliner.py -i test/alignments/nexus_primates/ | 
		python cloudforest.py ...

3. Cloudforest has two main functions. You can use it to infer genetrees by setting the `--gene-trees` flag or you can use it to bootstrap your data as described in [CITE] by setting the `--bootreps=#` flag. With the default setting cloudforest calculates  `--mraic`


Species Tree Estimation:
========================

This first command will run cloudforest locally with the practice data and produce a file called `genetrees.tre` in the `test` directory. ::

		python cloudforest_mpi.py \
		tests/alignments/phylip_primates/ \
		tests \
		--parallelism=multiprocessing \
		genetrees \
		binaries/PhyML3OSX

Details about all the modifier flags can be obtained by typing `python cloudforest_mpi.py -h`

You can also run essentially the same command using the emr script. Our recommendation is to use the command when you want to make sure you data looks ok prior to running a more extenstive analysis using a remote hadoop cluster. ::


		python nexus2oneliner.py -i test/alignments/nexus_primates/ | 
		python cloudforest_emr.py \
		--protocol=raw_value \
		--setup-cmd 'mkdir -p tmp' \
		--gene-trees \
		--archive=gzips/osx.phylo.tar.gz#bin \
		> primate.gene.trees

Setting `--protocol=raw_value` removes 'null' keys from the final file. `--setup-cmd 'mkdir -p tmp'` creates the appropriate `tmp/` directory for storing the phyml output files. This setting is only necessary for running local jobs. `tmp/` is present by default in on Amazon Elastic MapReduce Clusters.	
		
Setting the `--mraic` flag will makes cloudforest.py use mraic to infer the appropriate evolutionary model. :: 		

		python nexus2oneliner.py -i test/alignments/nexus_primates/ |
		python process.py \
		--setup-cmd 'mkdir -p tmp' \
		--gene-trees \
		--mraic \
		--archive=../gzips/osx.phylo.tar.gz#bin \
		> primate.gene.trees


And, this will run it on Amazon. ::

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
	exons and/or are using MrAIC to infer models. ::

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

Here's the command to run it locally with the practice data. ::

		python nexus2oneliner.py -i test/alignments/nexus_primates/ |
		python process.py \
		--setup-cmd 'mkdir -p tmp' \
		--full-analysis \
		--bootreps 5 \
		--archive=../gzips/osx.phylo.tar.gz#bin  \
		> primates.5bootreps.trees 

And, this will run it on Amazon. ::

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


