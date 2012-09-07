Basic Setup for OS X
--------------------

1. Download hadoop from one of the [mirrors](http://www.apache.org/dyn/closer.cgi/hadoop/core/)

2. Open `conf/hadoop-env.sh` and change `# export JAVA_HOME=/usr/lib/j2sdk1.6-sun ` to `export 
JAVA_HOME=/System/Library/Frameworks/JavaVM.framework/Versions/1.6.0`

3. Open `bin/hadoop-config.sh` and change `JAVA=$JAVA_HOME/bin/java` to `JAVA=$JAVA_HOME/Commands/java`

4. Open `conf/hdfs-site.xml` and add the following lines:

        <configuration>
            <property>
                <name>dfs.replication</name>
                <value>1</value>
            </property>
        </configuration>
        
5. Open `conf/core-site.xml` and add the following lines:

        <configuration>
            <property>
                <name>fs.default.name</name>
                <value>hdfs://localhost:9000</value>
            </property>
        </configuration>

6. Open `conf/mapred-site.xml` and add the following lines:

        <configuration>
            <property>
                <name>mapred.job.tracker</name>
                <value>localhost:9001</value>
            </property>
        </configuration>

7. Format the namenode with this command: `bin/hadoop namenode -format`

8. Start up Hadoop with: `bin/start-all.sh`

9. Type `jps` in the terminal and if everything worked you should see something like:

> 12543 DataNode

> 12776 Jps

> 12677 JobTracker

> 12755 TaskTracker

> 12619 SecondaryNameNode

Using Your Local Hadoop Node:
-----------------------------

1. Copy the files to the hadoop filesystem.

        $ bin/hadoop dfs -copyFromLocal /Users/nick/Desktop/Code/BioAWS/phylo/mpest.py mpest.py
        $ bin/hadoop dfs -copyFromLocal /Users/nick/Desktop/Code/BioAWS/phylo/phyml.py phyml.py
        $ bin/hadoop dfs -copyFromLocal /Users/nick/Desktop/Code/BioAWS/phylo/practice_alignments/3.align.oneliners.txt 3.align.oneliners.txt
        $ bin/hadoop dfs -copyFromLocal /Users/nick/Desktop/Code/BioAWS/phylo/osx.phylo.tar.gz osx.phylo.tar.gz 
        $ bin/hadoop dfs -copyFromLocal /Users/nick/Desktop/Code/BioAWS/phylo/setup.cfg setup.cfg
        
        $ bin/hadoop dfs -copyFromLocal /Users/nick/Desktop/Code/BioAWS/phylo/183Loci_29Species.oneliners.align 183Loci_29Species.oneliners.align
        
        
2. Run the streaming command.

        $ ./bin/hadoop jar /Users/nick/Desktop/hadoop-0.21.0/mapred/contrib/streaming/hadoop-0.21.0-streaming.jar \
        -file /Users/nick/Desktop/Code/BioAWS/phylo/phyml.py \
        -file /Users/nick/Desktop/Code/BioAWS/phylo/mpest.py \
        -file /Users/nick/Desktop/Code/BioAWS/phylo/osx.phylo.tar.gz \
        -cacheArchive osx.phylo.tar.gz#bin \
        -input 3.align.oneliners.txt \
        -output output/ \
        -mapper phyml.py \
        -reducer mpest.py \
        -verbose
        
3. Run the streaming command on AWS. The order of the commands is important with the mapper and reducer called last.

        # Test command:
        
         elastic-mapreduce --create --alive --stream \
        --cache-archive s3n://ngc-practice/exes/aws.phylo.tar.gz#bin \
        --input s3n://ngc-practice/data/3.align.oneliners.txt \
        --output s3n://ngc-practice/output_new \
        --mapper s3n://ngc-practice/exes/phyml.py \
        --reducer s3n://ngc-practice/exes/mpest.py \


        # Real command: MAKE SURE YOU DON'T HAVE UNDERSCORES IN THE FILE NAMES

        elastic-mapreduce --create --stream --num-instances 5 \
        --jobconf mapreduce.job.reduces=1 \
        --cache-archive s3n://ngc-practice/exes/aws.phylo.tar.gz#bin \
        --input s3n://ngc-practice/data/183Loci.29Species.oneliners.align \
        --output s3n://ngc-practice/output_new \
        --mapper s3n://ngc-practice/exes/phyml.py \
        --reducer s3n://ngc-practice/exes/mpest.py \
       
       

The trick is that you need to use 'cacheFile' on all the modules and binaries you want to 'call' from your mapper and reducer.  This ensures that they are accessible from the compute nodes.  The other trick is to include the following lines in your python code before the homemade module imports.

        import sys 
        sys.path.append('.')

Anyway, this works on my single node cluster. Brilliant! 


**Dumbo Instructions:**

Now you should be able to run Dumbo jobs on Elastic MapReduce. To start a cluster, you can use the Ruby client as so:

$ elastic-mapreduce --create --alive

SSH into the cluster using your EC2 keypair as user hadoop and install Dumbo with the following two commands:

$ wget http://bit.ly/ezsetup
$ sudo python ezsetup dumbo

Then you can run your Dumbo scripts. I was able to run the ipcount.py demo with the following command.

$ dumbo start ipcount.py -hadoop /home/hadoop \
-input s3://anhi-test-data/wordcount/input/ \
-output s3://anhi-test-data/output/dumbo/wc/



**Instructions derived from:**

http://www.infosci.cornell.edu/hadoop/mac.html

http://hadoop-karma.blogspot.com/2009/05/running-hadoop-020-pseudo-distributed.html

**More here:**

http://www.michael-noll.com/wiki/Running_Hadoop_On_Ubuntu_Linux_(Single-Node_Cluster)