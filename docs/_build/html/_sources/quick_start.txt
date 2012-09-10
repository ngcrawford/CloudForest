Quick Start:
============

The first thing you need to do to get CloudForest running to get it installed. We've tried to make this as painless as possible by using package installers and , but you'll still need to type some commands at the command line. 


1. Open terminal.
  Prepare to cut-n-paste!


2. Install `Python`_. 
  Python should already be installed on OS X or Linux. CloudForest requires version Python2.7 so enter ``python --version`` at the commandline and make sure you're uptodate. If you're not, I recommend installing `Enthought Python`_ 2.7 which includes almost all the dependancies necessary to get CloudForest running.

3. Install `R`_. 
  R is *the* open source statistical software package. You should be able to download and install a graphical package installer from the R website. We recommend using a recent version such as R 2.15 or greater. 

4. Install `Git`_. 
  Git is a distributed version control systerm. It's open source, easy to use, and integrates with github for easy collaboration and distribution.
  
  * We recommend first installing `HomeBrew`_ and then running ``brew install git`` at the command line to install git.

5. Install `Pip`_.
  Pip is a package installer for python that makes adding and managing packages and modules a breeze.

6. Install `CloudForest`_ by running ``pip install -U cloudforest`` at the commandline. Alternately you can install the most cutting edge development version by running ``pip install git+git://github.com/ngcrawford/CloudForest.git`` at the commandline.

  * CloudForest's ``setup.py`` script should install all the dependancies you need. However, Numpy can be troublesome. If you get errors when numpy is building you may need to install numpy manually. Running ``sudo pip install numpy`` should do the trick.

7. Install `Phybase`_.   
  * First install `ape`_ Phybase's only dependancy. 

    * To do this open R. At the R command line type: ``install.packages('ape')``. Follow the instructions.

  * Then download the gziped `source code`_ for phybase.
  * ``cd`` to that directory and type ``R CMD INSTALL phybase_1.3.tar.gz``. If the version of phybase is newer, you may have to adjust the gzip filename to refect this. 


That should do it for installation.


Configure Software:
-------------------

CloudForest itsself doesn't require any configuration, but one of it's dependancies does. `MrJob`_ has a pretty simple config file you'll need to setup if you want to run analyses on EMR. A full explaination is available `here`_.

You'll need to make a mrjob.conf file. Mine looks something like this::
    
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

Test the Installation:
----------------------

#. Run the UnitTests with `nose`_. 

I'm not yet sure how to do this with an 'installed package.'

.. _CloudForest: http://github.com/ngcrawford/CloudForest
.. _ape: http://ape.mpl.ird.fr/
.. _source code: http://phybase.googlecode.com/files/phybase_1.3.tar.gz
.. _Phybase: http://code.google.com/p/phybase/
.. _R: www.r-project.org/
.. _nose: http://nose.readthedocs.org/en/latest/
.. _Python: http://www.python.org/
.. _Git: http://git-scm.com/
.. _Pip: http://pypi.python.org/pypi/pip
.. _MrJob: https://github.com/Yelp/mrjob
.. _Enthought: Python http://www.enthought.com/products/epd.php
.. _HomeBrew: http://mxcl.github.com/homebrew/
.. _here: http://packages.python.org/mrjob/writing-and-running.html