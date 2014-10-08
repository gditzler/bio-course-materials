# What are we doing here?

This tutorial goes over how to use  Qiime  to process data from a high-throughput 16S rRNA sequencing studie collected from *Crawford et al., 2009*. The original tutorial can be found [here](http://qiime.org/tutorials/tutorial.html); however, due to changes in the software some components of the tuturial have been changed from its original to somethings that can be deployed on Protues's version of Qiime.  Before getting started, lets point a few things out: 

* `data/` contains all of the data that we would like to analyze
* `submitter.sh` is the script we're going to submit to SGE
* Before attempting to run any code you'll need to create a directory for the results of `submitter.sh` to be placed. create a folder called `output` in `bio-course-materials/proteus-demo/` (i.e., from this directory run `mkdir output`). Also, make sure nothing is in `output/` when you submit your job. This can be acheived by:
```bash 
  # assuming you are in bio-course-materials/proteus-demo
  rm -Rf output/ 
  mkdir output/
```


# Submitting to the Proteus cluster

You need to make sure you change from your group to the courses group before you submit your job to the queuing system. This is true for any job that you run. Once you have changed your group, run: 

```bash 
  qsub submitter.sh
```

