#!/usr/bin/env python
# Author Mahmoud Al-Bassam
# Program: prober1.py
# Version: 2.3
# Date: 05/24/2022
# Update note(s): tested code in python 3.4, 3.5, 3.6, 3.7, 3.8, and 3.9. The code successfully runs.

# Input:
# 1. Sorted and indexed bam file
# 2. Bed file containing the probe IDs and coordinates

# Output:
# log file containing the following information:
# 1. No. of reads in the input or the downsized bam file
# 2. No. of on target probes
# 3. No. of dropout probes
# 4. Percentage of dropout probes (probes with no reads)
# 5. Uniformity (coefficient of variance)

# Dependencies:
# argparse, numpy, datetime, subprocess, scipy and pandas are standard libraries that come with Anaconda
# samtools is available via TSCC.
# bedtools

# example command line:
'''prober2.3.py -bam rep1_10M.bam -bed 20220323_S3196374_XTHS_A_1_Probes.txt.bed -outdir . -prefix test4 -bam_depth 500000'''
#--------------------------------------------------------------------------------------------------

########################### 1. Setting up user command line options ############################

import argparse
import sys
#using argparse to print help page which allows the user to understand the option
# use prober1.py --help to see the help page
parser = argparse.ArgumentParser(description='Takes Bam and Bed files, calculates aligned reads on target and outputs a log file. The program can reduce the depth of the bam file if needed, see -bam_depth option')
parser.add_argument("-bed", help="bed file containing probe information", type=str)
parser.add_argument("-bam", help="bam file (single file with paired-end reads)", type=str)
parser.add_argument("-prefix", help="Basename for output log file ", type=str)
parser.add_argument("-outdir", help="Basename for output directory", type=str)
parser.add_argument("-bam_depth", help="Reduces the bam file to an integer number of reads (optional). Otherwise the program takes the full-depth input bam file", type=int)
args = parser.parse_args()

# print help and error if any arguments are missing
try:
    options = parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

# importing libraries

from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats
from subprocess import check_output
import matplotlib.pyplot as plt
############################# 2. Downsampling of large bam files ##############################

bed = args.bed
prefix = args.prefix
outdir = args.outdir
bam_depth= args.bam_depth
bam = args.bam

def run_probe_analysis(df, argparsetotal_reads):
	'''A function that takes count dataframe and
	total number of reads to calculate dropout rate
	and coefficient of variance'''

	############################### 3. Performing the calculations ################################

	total_number_of_probes = check_output("wc -l {}".format(bed), shell=True)
	total_number_of_probes = total_number_of_probes.decode("utf-8")
	total_number_of_probes = float( total_number_of_probes.split()[0] )
	total_number_of_probes = int( total_number_of_probes )

	# calculate total number of probes with 1 or more reads aligned
	drop_out = len(df[df['count'] < 1])
	on_target = len(df[df['count'] >= 1])

	# calculate percent_on_target
	percent_drop_out =  drop_out / total_number_of_probes * 100
	percent_drop_out = '{:.4f}'.format(percent_drop_out)

	# calculate Coefficient of Variation
	Coefficient_of_Variation = stats.variation(list(df['count']), axis = 0)
	Coefficient_of_Variation = '{:.2f}'.format(Coefficient_of_Variation)

	if Coefficient_of_Variation == 'nan':
		Coefficient_of_Variation = 0
	# can also be calculated without scipy as follows:
	# Coefficient_of_Variation = np.std(df['{}'.format(args.bam)]) / np.mean(df['{}'.format(args.bam)])

        # calculate variance
	variance = np.var(df['count'])
	variance = '{:.2f}'.format(variance)
	if variance == 'nan':
		variance = 0

	################################## 4. writing up the results ##################################

	now = datetime.now() #today's date
	dt_string = now.strftime("%d/%m/%Y %H:%M:%S") #formatting the date

	writing_ar = ['total_number_of_reads','on_target_probes','dropout_probes', 'percent_dropout','CV', 'variance'] # this is the header
	with open('{}/{}.prober.log'.format(outdir,prefix), 'w+') as f:

	# writing the final data
		#f.write('{}\n\n'.format(dt_string))
		f.write ('\t'.join(writing_ar) + '\n')
		f.write ('{}\t{}\t{}\t{}%\t{}\t{}\n'.format(total_reads, on_target, drop_out, percent_drop_out, Coefficient_of_Variation, variance))

	################################## 5. make a KDE plot of read counts distribution #####################

	df['count'].plot.kde(linewidth=3, label = '{}'.format(prefix))
	plt.xlabel('Counts', fontsize=16)
	plt.ylabel('Density', fontsize=16)
	plt.legend()
	plt.xlim(-100,500)
	plt.title('{}\nCV={}, Var={}'.format(prefix,Coefficient_of_Variation, variance))
	plt.savefig('{}/{}_kde_plot.pdf'.format(outdir,prefix), bbox_inches='tight')


if bam_depth: #if the bam_depth option is used
	# calculating the ratio required to get the downsized depth reads from input bam
	total_reads = check_output("samtools view -c {}".format (bam), shell = True)
	total_reads = float(total_reads.decode("utf-8"))
	total_reads = int(total_reads)
	depth_ratio = (bam_depth / total_reads)
	depth_ratio = float(depth_ratio)

	if float(depth_ratio) < 1.0 : # if bam file size >downsized_depth reads then downsizing makes sense

		# writing the downsized bam file
		check_output ("samtools view -bh {} | samtools view -b -s {} >{}/{}_{}.bam".format(bam, depth_ratio, outdir, prefix, bam_depth), shell=True)

		# running bedtools to generate bed file containing the number of reads per each probe
		check_output("bedtools coverage -a {} -b {}/{}_{}.bam -counts >{}/{}.counts".format(bed, outdir, prefix, bam_depth, outdir, prefix), shell=True,)

		total_number_of_reads = check_output("samtools view -c {}/{}_{}.bam".format(outdir,prefix,bam_depth), shell =True)
		total_number_of_reads = total_number_of_reads.decode("utf-8")
		total_number_of_reads = int(float(total_number_of_reads))

		fc_table = pd.read_csv('{}/{}.counts'.format(outdir,prefix), sep='\t',
							names = ['chr','start','end','feature','count'])
	# calculate the total number of probes from the bed file

		# run the function to print the log file
		run_probe_analysis(fc_table, total_number_of_reads)
	
	elif float(depth_ratio) >= 1.0: # if input bam file depth is less than the target downsized_depth then raise an error
		raise ValueError('The input bam file reads cannot be smaller than or equal to the downsized bam file! Please reduce the downsized bam file size or run the command without the -bam_depth option')

else: # else just run bedtools coverage on the full-depth bam file
	check_output("bedtools coverage -a {} -b {} -counts >{}/{}.counts".format(bed, bam, outdir, prefix), shell=True)

	fc_table = pd.read_csv('{}/{}.counts'.format(outdir,prefix), sep='\t',
							names = ['chr','start','end','feature','count'])
	total_number_of_reads = int(check_output("samtools view -c {}".format(bam), shell=True))
	# run the analysis on the full-depth bam file
	run_probe_analysis(fc_table, total_number_of_reads)
