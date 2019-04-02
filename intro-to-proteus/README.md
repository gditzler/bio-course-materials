# About 

This repo presents some of the basic functionality for using `Proteus` on [Drexel University's Research Computing Facility (UCRF)](http://www.drexel.edu/research/urcf/). Proteus uses the [Univa Grid Engine](http://www.univa.com/products/grid-engine) as the batch queuing system to submit and handle jobs. Before getting started, read through the [Proteus wiki](https://proteusmaster.urcf.drexel.edu/urcfwiki) for more information on using the batch queuing system. 

# SSHing into Proteus

Proteus has two head nodes that you can use for compting. Think of the nodes as being a computer or server. The two head nodes either use: (a) AMD, or (b) Intel processors. To connect to one of the head nodes, run the follwing commands in the shell. 

```bash
  # sshing into AMD servers 
  ssh $USERNAME@proteusa01.urcf.drexel.edu
  # sshing into Intel servers 
  ssh $USERNAME@proteusi01.urcf.drexel.edu
```
where `$USERNAME` is your username for proteus. Note that in the shell, the `$` is used to denote a variable. For example, if my username was `greg`, I could use something like: 

```bash 
  USERNAME=greg
  echo "My user name is $USERNAME" 
  ssh $USERNAME@proteusi01.urcf.drexel.edu
```

# Using Modules 

Proteus is using modules to setup the environment in the shell and the software that the user has available to them. The Proteus wiki has some useful information on navigating though the [modules](https://proteusmaster.urcf.drexel.edu/urcfwiki/index.php/Environment_Modules). 

```bash
  # view available modules
  module avail 
  # view the modules that are currently loaded 
  module list
```

## Loading Modules & Some Useful Modules 

You can manully load each of the modules as you need them; however, it may make your life easier to add the following commands to your `~/.bashrc` file. 

```bash
  module load python/2.7.6
  module load R/3.0.2
  module load matlab/R2013a
  module load ncbi-blast/gcc/64/2.2.29 
```

# Simple Qsub and Scripts 

## Writing the Script
The outline for this script can be found on the [Proteus wiki](https://proteusmaster.urcf.drexel.edu/urcfwiki/index.php/Writing_Job_Scripts). As an example, consider `simple-script.sh`. The `#$` tell the scheduler that these lines are to be interpreted as flags. 

```bash
  #!/bin/bash
  # tell SGE to use bash for this script
  #$ -S /bin/bash
  # execute the job from the current working directory, i.e. the directory in which the qsub command is given
  #$ -cwd
  # set email address for sending job status
  #$ -M fixme@drexel.edu
  # project - basically, your research group name with "Grp" replaced by "Prj"
  # -P fixmePrj
  # select parallel environment, and number of job slots
  #$ -pe openmpi_ib 256
  # request 15 min of wall clock time "h_rt" = "hard real time" (format is HH:MM:SS, or integer seconds)
  #$ -l h_rt=00:15:00
  # a hard limit 8 GB of memory per slot - if the job grows beyond this, the job is killed
  #$ -l h_vmem=8G
  # want at least 6 GB of free memory
  #$ -l mem_free=6G
  # select the queue all.q, using hostgroup @intelhosts
  #$ -q all.q@@intelhosts 

  module load python/2.7.6
  python -c "print 'The task ID is $SGE_TASK_ID'""
```

## Calling the Script

We can submit the previous script to the scheduler using: 

```bash
  newgrp fixmeGrp
  qsub simple-script.sh  
```
at the shell. Note that the previous commands will produce an error since the project and group names were fictitious.

## Monitoring Progress

```bash 
  qstat
  qstat -f 
```
