Basic Setup for OS X
--------------------

1. Download hadoop from one of the [mirrors](http://www.apache.org/dyn/closer.cgi/hadoop/core/)

2. Open `conf/hadoop-env.sh` and change `# export JAVA_HOME=/usr/lib/j2sdk1.6-sun ` to `export 
JAVA_HOME=/System/Library/Frameworks/JavaVM.framework/Versions/1.6.0`

3. Open `bin/hadoop` and change `JAVA=$JAVA_HOME/bin/java` to `JAVA=$JAVA_HOME/Commands/java`

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

        $ bin/hadoop dfs -copyFromLocal /Users/nick/Desktop/Code/BioAWS/temp/phyml.py phyml.py
        $ bin/hadoop dfs -copyFromLocal /Users/nick/Desktop/Code/BioAWS/temp/tree.py tree.py
        $ bin/hadoop dfs -copyFromLocal /Users/nick/Desktop/Code/BioAWS/temp/newick.py newick.py
        $ bin/hadoop dfs -copyFromLocal /Users/nick/Desktop/Code/BioAWS/temp/PhyML_3.0 PhyML_3.0
        $ bin/hadoop dfs -copyFromLocal /Users/nick/Desktop/Code/BioAWS/temp/mpestBFV1 mpestBFV1

2. Run the streaming command.

        $ ./bin/hadoop jar /Users/nick/Desktop/hadoop-0.21.0/mapred/contrib/streaming/hadoop-0.21.0-streaming.jar \
        -file /Users/nick/Desktop/Code/BioAWS/temp/practice_alignments/3.align.oneliners.txt \
        -file /Users/nick/Desktop/Code/BioAWS/temp/mpest.py \
        -file /Users/nick/Desktop/Code/BioAWS/temp/phyml.py \
        -cacheFile hdfs://localhost:9000/user/nick/newick.py#newick.py \
        -cacheFile hdfs://localhost:9000/user/nick/tree.py#tree.py \
        -cacheFile hdfs://localhost:9000/user/nick/mpestBFV1#mpestBFV1 \
        -cacheFile hdfs://localhost:9000/user/nick/PhyML_3.0#PhyML_3.0 \
        -input 3.align.oneliners.txt \
        -output output/ \
        -mapper phyml.py \
        -reducer mpest.py \
        -verbose

The trick is that you need to use 'cacheFile' on all the modules and binaries you want to 'call' from your mapper and reducer.  This ensures that they are accessible from the compute nodes.  The other trick is to include the following lines in your python code before the homemade module imports.

        import sys 
        sys.path.append('.')

Anyway, this works on my single node cluster. Brilliant! 

**Instructions derived from:**

http://www.infosci.cornell.edu/hadoop/mac.html

http://hadoop-karma.blogspot.com/2009/05/running-hadoop-020-pseudo-distributed.html

**More here:**

http://www.michael-noll.com/wiki/Running_Hadoop_On_Ubuntu_Linux_(Single-Node_Cluster)