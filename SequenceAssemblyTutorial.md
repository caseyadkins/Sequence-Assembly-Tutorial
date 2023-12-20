# Sequence Assembly Tutorial

Authored by Casey Adkins and Ellen McMullen

### Overview

This workflow will take you through the steps of assembly for restriction digest methods (RAD, ddRAD, GBS) or related amplification-based (e.g., NextRAD, RApture) sequences, and inferring population structure of large SNP datasets in a variational Bayesian framework.

**Outline**


1. [Installing Packages](#installing)  
   A) `conda`  
   B) `ipyrad`  
   C) `structure_threader`  
   D) `fastStructure`
2. [Run `ipyrad`](#ipyrad)
3. [Run `structure_threader`](#structure_threader)
4. [Extract Information from Result Files](#results)
5. [Plot Important Results](#plots)
6. [Conclusion](#end)

## 1. Installing Packages<a name="installing"></a>

### Installing `conda`

1. Login to pronghorn with `ssh netID@pronghorn.rc.unr.edu` or via PuTTY or open your terminal.
2. Go to the [linux installation page](https://docs.anaconda.com/free/anaconda/install/linux/) and copy the prompt for Linux x86; it should be something like `curl -O https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh`.
3. `cd` to the directory you want to install conda in (your home directory should be fine - /data/gpfs/home/netID)
4. Run the `curl` command.
5. When the download is complete, run `bash Anaconda3-2023.09-0-Linux-x86_64.sh` (filled in with the up-to-date version you downloaded).
6. Press enter to view the License, and press and hold enter to scroll through the License, or press `q` to skip to the agreement. Enter `yes` to agree. Press enter again to begin the installation.
7. When the installation finishes, enter `yes` again to update. your shell profile to use conda.
8. Run the command `source ~/.bashrc` to refresh the shell settings. You should now see (base) before your prompts, indicating that you are in the base conda environment. Alternatively, you can exit Pronghorn and re-login.
9. Test that conda is working correctly by running `conda list`. You should see a list of installed packages.
10. After you have confirmed that conda is working correctly, you can delete the installer with `rm Anaconda3-2023.09-0-Linux-x86_64.sh`

**Note:** On Windows, after installing conda, if (base) or (your\_environment\_name) does not show up in the prompt, conda prompts may not work. Try running `bash` to get the environment name to show up in the prompt and test that conda is working with `conda list`. This applies to all future conda prompts; if they are not working, try to run `bash` first.

### Creating your conda environment

You may want to create a conda environment for specific projects and software. By organizing software into environments, you can work with different software versions. This can be important because some software requires other software (e.g., python) to be a specific or older version than others.

```
conda create -n myenvname
conda activate myenvname
```

#### To assemble DNA sequences you want to create three different conda environments.

One environment will house the ipyrad software, the second will house the faststructure software, and the third will house the structure_threader software. Like before, when you are in your base conda environment run the following to create and activate or enter that environment. You will need to enter `y` after the conda create command to approve the new environment.

```
conda create -n ipyrad
conda activate ipyrad
```

And make another one for fast_structure.

```
conda create -n fs
conda activate fs
```
Now you can create a new environment for structure_threader.

```
conda create -n str
conda activate str
```

Check that your new environments show up with `conda env list`.

### Install `ipyrad`

The package `ipyrad` is an assembly and analysis software for RADseq data. Additional documentation for ipyrad can be found [here](https://ipyrad.readthedocs.io/en/latest/index.html).

First, you should activate your ipyrad environment:

```
conda activate ipyrad
```

Next, we will use conda to install ipyrad. There are multiple locations or 'channels' where packages are stored. Using the -c option (for channel) will allow you to specify which of these channels conda should search first: 

```
conda install ipyrad -c conda-forge -c bioconda
```

Enter `y` when prompted to install.

There are additional packages that may also be needed. You can install these with:

```
conda install notebook -c conda-forge
conda install mpi4py -c conda-forge
```

**Note:** These may automatically install with some versions of ipyrad.

The conda-forge and bioconda channels are well vetted. You can add them as default channels to search in with:

```
conda config --add channels conda-forge
conda config --add channels bioconda
```

### Install `structure_threader`

The package `structure_threader` is a wrapper software for `structure`, `faststructure`, `MavericK`, and `ALStructure`. It can run one or all of these software in parallel decreasing the run time for DNA assembly. Additional documentation can be found [here](https://structure-threader.readthedocs.io/en/latest/).

You can install `structure_threader` using `pip3`, which is a package manager. `conda` is both a package manager and a virtual environment manager.

#### Checking for `pip3` and `python`

To use `pip3` with `conda`, you first need to check that `pip3` is installed. If you installed the full anaconda version, `pip3` is automatically installed in your base environment. However, you need to have it in your str environment as well.

1. Activate the str environment: `conda activate str`
2. Check if `pip3` is installed: `pip3 -- version`
3. If `pip3` is installed it will give you a version number with a path to the installation.
4. If `pip3` is not installed you will get an error message: "bash: pip: command not found..."

If `pip3` is not installed:

1. Check if `python` is installed in the environment: `which python`
2. If `python` is installed in the environment it will give you a path to the installation that includes your conda installation in the path. You can check which version it is with `python --version`
3. If `python` is not installed in the environment it will give you a path to the installation that does not include your conda installation path, for example, /usr/bin/python

If you need to install `python`:

1. Search which python versions are available with `conda search python`
2. Install your desired version. For example: `conda install python=3.12.0` and type `y` to install.
3. Check again for `pip3` with `pip3 --version`.

If you still need to install `pip3`

1. Run `python -m pip`
2. Confirm with `pip3 --version`

#### Using `pip3` to install `structure_threader`

1. Make sure your str environment is active.
2. Run `pip3 install structure_threader` 


### Install `faststructure`

We will be using `faststructure` within `structure_threader`, so now we need to install `faststructure` as well. Documentation for `faststructure` can be found [here](https://rajanil.github.io/fastStructure/).

`faststructure` needs python version 2.7, so we need to install python 2.7 first:

```
conda activate fs
conda install python=2.7
conda install faststructure -c bioconda
```

## 2. Run `ipyrad`<a name="ipyrad"></a>

The goal of ipyrad is to take raw sequence data and assemble it to a reference genome or denovo. For ipyrad to run you need a set of fasta files (.fasta or .fasta.gz), a parameter file (.txt), and a slurm file (.txt).

### A) Create ipyrad parameter file

Below are the contents of the ipyrad parameter file (.txt) with mostly default parameters. Documentation on the parameters can be found [here](https://ipyrad.readthedocs.io/en/latest/6-params.html).

You can create a default parameter file to edit by running `ipyrad -n test` in your ipyrad environment, in the folder that is storing your .fastq files. You can then edit the file using `nano params-test.txt`. For example, we will change the assembly method to *reference* rather than *denovo*, which requires us to also add the path to the reference genome sequence. Check out the other changes to the file:

```
------- ipyrad params file (v.0.9.24)-------------------------------------------
Project-Name  ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
/YOUR/DIRECTORY ## [1] [project_dir]: Project dir (made in curdir if not present)
## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
## [3] [barcodes_path]: Location of barcodes file
/YOUR/DIRECTORY/*.fastq.gz ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
reference ## [5] [assembly_method]: Assembly method (denovo, reference)
/YOUR/DIRECTORY/genomefile.fasta ## [6] [reference_sequence]: Location of reference sequence file
ddrad ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
5 ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
33 ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)
6 ## [11] [mindepth_statistical]: Min depth for statistical base calling
6  ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
10000  ## [13] [maxdepth]: Max cluster depth within samples
0.85  ## [14] [clust_threshold]: Clustering threshold for de novo assembly
2 ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
2 ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
35 ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2 ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
0.05 ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus
0.05 ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus
5 ## [21] [min_samples_locus]: Min # samples per locus for output
0.2 ## [22] [max_SNPs_locus]: Max # SNPs per locus
8 ## [23] [max_Indels_locus]: Max # of indels per locus
0.5 ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus
0, 130, 0, 130 ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0 ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
* ## [27] [output_formats]: Output formats (see docs)
## [28] [pop_assign_file]: Path to population assignment file
## [29] [reference_as_filter]: Reads mapped to this reference are removed in step 3
```

### B) Create Slurm Script

Slurm is a job management tool for HPCs. When using Pronghorn, we need to submit jobs using slurm. The basic command is `sbatch path/to/slurm/script.txt`. The slurm script file will contain some HPC-related parameters as well as the actual bash command(s) you want to run.

The following are the contents of slurm file for running ipyrad:

```
#!/bin/bash
#SBATCH --job-name=jobname
#SBATCH --output=outputfilename.out
#SBATCH --error=errorfilename.err
#SBATCH --account=cpu-s2-pronghornassociation-0
#SBATCH --partition=cpu-s2-core-0
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@unr.edu

## directory with params file
PARAMS=/YOUR/DIRECTORY/params.txt

## call ipyrad 
ipyrad -p $PARAMS -s 1234567 -c 32 -d
```

The first line (#!/bin/bash) denotes this file is a bash script.

The next 12 lines are slurm notation to inform the run.

1. `--job-name` name the job
2. `--output` name the output file
3. `--error` name the error file
4. `--account` what pronghorn account you are working on. Could be your personal account, or a lab group association account.
5. `--partition` which partitions (groups of nodes) are you requesting
6. `--time` how much run-time are you requesting
7. `--nodes` how many nodes are you requesting
8. `--tasks-per-node` how many tasks should be assigned to each node
9. `--cpus-per-task` how many CPUs should be assigned to each task
10. `--mem-per-cpu` how much memory (RAM) do you need per CPU
11. `--mail-type` what do you want to get email notifications about
12. `mail-user` your email address to receive email notifications

After the slurm information, you put the commands you want to run in your job.

First, create a bash parameter, `$PARAMS` that contains the path to the parameters file we created in part 4a.

Next, run the ipyrad assembly command.

* `-p` refers to the parameter file
* `-s` refers to the steps to run
* `-c` refers to the number of cores to use
* `-d` means debug mode, so if an error occurs, it will provide an informative error message.

### C) Understand `ipyrad` steps

`ipyrad` performs several steps in a row to execute the assembly:

1. Demultiplexing / Loading fastq files.
2. Filtering / Editing reads
3. Clustering / Mapping reads within Samples and alignment
4. Joint estimation of heterozygosity and error rate
5. Consensus base calling and filtering
6. Clustering / Mapping reads among Samples and alignment
7. Filtering and formatting output files

### D) Submit the slurm job

```
conda activate ipyrad

#navigate to your ipyrad directory
cd /data/gpfs/assoc/YOUR/DIRECTORY/

#submit slurm job to run ipyrad
sbatch ./ipyrad-file.txt
	#batch job ID #######
```

Record your job number so you can access the job later if you need to check on the progress or cancel it. Under these computing parameters, it takes approximately 20 hours for ipyrad to assemble 200 individual samples.

Once completed, you will have an output folder within your directory with all of the files selected in the parameters file, but importantly there is a text file with statistics of the run which provides information about read depth, SNP calls, etc. You should look through those data to look for any errors or outliers.

Most importantly ipyrad creates a Structure file (.str) for you within your output directory. You will use this file as part of an input for structure_threader.

## 3. Run `structure_threader`<a name="structure_threader"></a>

The goal of this step is to infer population structure by taking Bayesian estimates of the underlying ancestry of individuals. You will need a str file from the ipyrad output, an individual file, a main parameters file, and a slurm file.

### A) Copy your str file to your structure_threader folder

`cp file.str ../structure_step`

### B) Create Individual File

We need a file that contains the individual IDs for each individual or sample in the analysis. Copy the ipyrad output file that contains the heterozygosity estimates for each individual (.txt extension) to your working directory. Run `python extract_ind_file.py file.txt`. This will produce a file called 'individuals.indfile' that will be used by `structure_threader`.

Here is the python script we wrote for this purpose, extract_ind_file.py:

```
'''
This python script takes 1 input file, the ipyrad output file that contains the 
heterozygosity estimates for each individual. The file must have a '.txt' extension.

The script will output a .csv file with a column of individual IDs (no header)
'''

import sys
import itertools
import pandas as pd

het_file = sys.argv[1]

# extract the id column
with open(het_file, 'r') as het:
	for i, line in enumerate(het):
		if line.strip('\n') == '## Final Sample stats summary':
			#print('section found')
			break
id_list = []
with open(het_file, 'r') as het:
	for line in itertools.islice(het, i+2, None):  # start=17, stop=None
		parts = line.strip('\n').split()
		if len(parts) == 12:
			id_part = parts[0]
			id_list.append(id_part)
		else:
			pass
#print(id_list)

# convert to series
s = pd.Series(id_list)
#print(s)

# save output
s.to_csv('individuals.indfile', header=False, index=False)
```

### C) Create a Main Parameters File

You can use the param mode of `structure_threader` to get a new param file within your working directory.

`structure_threader params`

The main parameters file will look similar to this.

```
KEY PARAMETERS FOR THE PROGRAM structure. YOU WILL NEED TO SET THESE IN ORDER TO RUN THE PROGRAM. VARIOUS OPTIONS CAN BE ADJUSTED IN THE FILE extraparams.

"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean (ie insert 1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes!)

Basic Program Parameters

#define MAXPOPS   #      // (int) number of populations assumed - K
#define BURNIN    100000   // (int) length of burnin period
#define NUMREPS   10000000   // (int) number of MCMC reps after burnin

Input/Output files

#define INFILE   /YOUR/DIRECTORY/file.str   // (str) name of input data file
#define OUTFILE  /YOUR/DIRECTORY/output_directory   //(str) name of output data file

Data file format

#define NUMINDS    ###  // (int) number of diploid individuals in data file
#define NUMLOCI     ###    // (int) number of loci in data file
#define PLOIDY       2    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line

#define LABEL     0     // (B) Input file contains individual labels
#define POPDATA   0     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier

#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data before the genotype data start.

#define MARKERNAMES      0  // (B) data file contains row of marker names
#define RECESSIVEALLELES 0 // (B) data file contains dominant markers (eg AFLPs) // and a row to indicate which alleles are recessive
#define MAPDISTANCES     0  // (B) data file contains row of map distances // between loci

Advanced data file options

#define PHASED           0 // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO        0 // (B) the data for each individual contains a line indicating phase (linkage model)
#define MARKOVPHASE      0 // (B) the phase info follows a Markov model.
#define NOTAMBIGUOUS  -999 // (int) for use in some analyses of polyploid data

Command line options:

-m mainparams
-e extraparams
-s stratparams
-K MAXPOPS
-L NUMLOCI
-N NUMINDS
-i input file
-o output file
-D SEED
```

You can change these parameters from their default to get your desired results, but you have to define `MAXPOPS`, `NUMINDS`, and `NUMLOCI`.

1. `MAXPOPS` is the number of populations you expect within your individuals. It can complete a range of populations or K, for example, if you set the `MAXPOPS` to 4, then it will predict population assignment for one population, two, three, and lastly four populations. Results for all four population sizes will be created and it will estimate the best K for your data. It's important to select the best K and to move forward based on the structure_threader estimate, but also the biology and natural history of your samples.

2. `NUMINDS` is simply the number of individuals in your dataset

3. `NUMLOCI` refers to the number of loci within your dataset, it's easiest to find this number from the `ipyrad` stats output file (.txt).

### D) Create Slurm Script

The following are the contents of slurm file for running structure_threader:

```
#!/bin/bash
#SBATCH --job-name=jobname
#SBATCH --output=outputfilename.out
#SBATCH --error=errorfilename.err
#SBATCH --account=cpu-s2-pronghornassociation-0
#SBATCH --partition=cpu-s2-core-0
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@unr.edu

structure_threader run -i file.str -o outputfilename --params ./mainparams.txt -fs /YOUR/DIRECTORY/conda/envs/fs/bin/fastStructure --ind ./individuals.indfile -K # -t 5

structure_threader plot -i outputfilename -f faststructure -K # -o out_plots --ind individuals.indfile --use-ind-labels
```

See above for descriptions of SBATCH parameters.

The first line is the  `structure_threader` call and `run` mode.

* `-i` denotes the input file, which is the .str file from `ipyrad`.
* `-o` denotes the outputfile name, the same name as the SBATCH output.
* `--params` denotes the file path to the parameter file.
* `-fs` denotes the file path to your `fastStructure` software, which is in your `fs` conda environment.
* `--ind` denotes the file path to your individual file, with your individual sample id. Note: the file extension must be .indfile
* `-K` (int) is the number of populations - this value must match the `MAXPOPS` value in the parameters file.
* `-t` (int) is the number of threads to use

The second line is the `structure_threader` call and `plot`
 mode. This line will take the results from the `run` and create the classic structure barplot.
 
 * `-i` denotes the input file, which is the `structure_threader run` output.
 * `-f` denotes the external program used.
 * `-K`  values that you want to plot. Each individual K value that is provided will be plotted individually and in the end, a comparative plot will all K values will also be generated. (-K; Example: -K 2 3 4).
 * `-o` denotes the file name of output plots.
 * `--ind` denotes the file with individual IDs/names.
 * `--use` use individual sample labels from the .indfile.

### E) Submit the Slurm Job

```
conda activate str

# navigate to your working directory
cd /data/gpfs/assoc/YOUR/DIRECTORY/

# submit slurm job to run structure threader
sbatch ./str_script.txt
```

## 4. Extract Information from Result Files<a name="results"></a>

There are many output files from ipyrad and structure_threader and many analyses can be done with these outputs. For our purposes, we take 3 files and combine them into a .csv file for further analysis and plotting.

Files needed:

1. The ipyrad output file that contains the heterozygosity estimates for each individual (.txt)
2. The structure_threader output file that contains the q value estimates for each individual (.meanQ)
3. The structure_threader input file that provides individual IDs for the output q values (.indfile)

Put these three files in a new directory with the combine\_output.py script. 

Run (with the appropriate file names):

```
python combine_output.py file.txt file.meanQ file.indfile
``` 

The script will generate a new .csv file that has columns for Individual_ID, Q1, Q2, and Heterozygosity, and a row for each individual.

Here is the python script we wrote for this purpose, combine\_ouputs.py:

```
'''
This python script takes 3 input files:
1) The ipyrad output file that contains the heterozygosity estimates for each individual
2) The structure_threader output file that contains the q value estimates for each individual
3) The structure_threader input file that provides individual IDs for the output q values

Files can be provided in any order, but must have the file extensions '.txt',
'.meanQ', and '.indfile' respectively.

The script will output a .csv file with columns:
1) Individual_ID
2) Q1
3) Q2
4) Heterozygosity
'''

import sys
import itertools
import pandas as pd

# sort the input files by filetype
het_file = 0
q_file = 0
id_file = 0
for f in sys.argv[1:]:
	if f.endswith('.txt') and not het_file:
		het_file = f
		print(f'het_file: {f}')
	elif f.endswith('.txt') and het_file:
		raise ValueError('User provided multiple ipyrad output files')
	if f.endswith('.meanQ') and not q_file:
		q_file = f
		print(f'q file: {f}')
	elif f.endswith('.meanQ') and q_file:
		raise ValueError('User provided multiple structure_threader output files')
	if f.endswith('.indfile') and not id_file:
		id_file = f
		print(f'id file: {f}')
	elif f.endswith('.indfile') and id_file:
		raise ValueError('User provided multiple structure_threader input files')
		
assert het_file, 'het_file not set'
assert q_file, 'q_file not set'
assert id_file, 'id_file not set'
		
# extract the individual id column
ind_col = []
with open(id_file, 'r') as indivs:
	for line in indivs:
		l = line.strip('\n')
		ind_col.append(l)
print(f'{len(ind_col)} individuals extracted from id file')

# extract the q columns
q1_col = []
q2_col = []
with open(q_file, 'r') as qs:
	for line in qs:
		l = line.strip('\n')
		l_parts = l.split('  ')
		q1_col.append(l_parts[0])
		q2_col.append(l_parts[1])
print(f'{len(q1_col)} individuals extracted from q1 in q file')
print(f'{len(q2_col)} individuals extracted from q2 in q file')

assert len(q1_col)==len(q2_col), 'Uneven column lengths'
assert len(q1_col)==len(ind_col), 'id file and q file have unequal numbers of individuals'

# extract the het column
# make dictionary of id: het
with open(het_file, 'r') as het:
	for i, line in enumerate(het):
		if line.strip('\n') == '## Final Sample stats summary':
			#print('section found')
			break
het_dict = {}
with open(het_file, 'r') as het:
	for line in itertools.islice(het, i+2, None):  # start=17, stop=None
		parts = line.strip('\n').split()
		if len(parts) == 12:
			id_part = parts[0]
			het_part = parts[8]
			het_dict[id_part] = het_part
		else:
			pass
#print(het_dict)
print(f'{len(het_dict)} individuals extracted from het file')

assert len(het_dict)==len(ind_col), 'het file has a different number of individuals than id file or q file'
		
# make dict of id:[q1, q2, het]
final_dict = {}
for i, ind in enumerate(ind_col):
	vals = [q1_col[i], q2_col[i], het_dict[ind]]
	final_dict[ind]=vals
#print(final_dict)
		
# turn into df
df = pd.DataFrame.from_dict(final_dict, orient='index').reset_index()
df.columns = ['Individual_ID', 'Q1', 'Q2', 'Heterozygosity']
df.to_csv('combined_output.csv', index=False)
```

## 5. Plot Important Results<a name="plots"></a>

Finally, we will use the .csv file to create a common plot for this kind of data, a cluster assignment structure plot (plots q proportions for each individual).

You could also take the .csv file and import it into your program of choice for further analysis.

In the directory that contains your .csv file, run `python plotting.py combined_output.csv`

Here is the python script we wrote for this purpose, plotting.py:

```
import sys
import pandas as pd
import matplotlib.pyplot as plt

# read in the dataframe
df = pd.read_csv(sys.argv[1])

# create a structure plot
# sort by q1 value
df_sorted = df.sort_values(by=['Q1','Individual_ID'], ascending=[False, True])
# extract x axis labels
inds = df_sorted['Individual_ID']
# extract Q1 bar lengths
q1 = df_sorted['Q1']
# extract Q@ bar lengths
q2 = df_sorted['Q2']
# create a fig object
fig = plt.subplots(layout='constrained')
# add first level (q1)
p1 = plt.bar(inds, q1, color = 'blue')
# add second level (q2) starting where q1 left off
p2 = plt.bar(inds, q2, bottom = q1, color = 'red')
# add labels
plt.ylabel('Cluster Assignment')
#plt.xlabel('Individual')
plt.xticks(rotation=90)
plt.legend((p1[0], p2[0]), ('N. fuscipes', 'N. macrotis'))
plt.margins(x=0, tight=True)
# plt.show()
plt.savefig('StructurePlot.png')
```

This results in the structure plot:

![image](./StructurePlot.png)

## 6. Conclusion<a name="end"></a>

We now know the ancestry of all these individuals!
