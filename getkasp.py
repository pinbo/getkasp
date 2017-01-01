#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  getkasp
#  
#  Copyright 2016 Junli Zhang <zhjl86@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

### Imported
from subprocess import call
import getopt, sys, os

#########################
usage="""
getkasp.py 
	-i <alignment.fa>
	-p <getkasp path> 
	-o <output file name>"
	-s <SNP site>
	-t <target_sequence_ID>
	-a <altanative allele>

Example:
	./getkasp.py -i alignment_raw.fa -p . -o primer3008.txt -s 489 -t 2AS -a T

"""

snpname = "T4-3008"
snp_site = 489
getkasp_path = "."
RawAlignFile = "alignment_raw.fa"
out = ""
target = "2AS" # target sequence name
product_min = 50
product_max = 250
alt_allele = "N"
#######################
try:
	opts, args = getopt.getopt(sys.argv[1:], "i:p:o:s:t:a:h", ["help"])
except getopt.GetoptError as err:
	# print help information and exit:
	print str(err)  # will print something like "option -a not recognized"
	print usage
	sys.exit(2)
for o, a in opts:
	if o == "-i":
		RawAlignFile = a
	elif o in ("-h", "--help"):
		print usage
		sys.exit()
	elif o in ("-p"):
		getkasp_path = a
	elif o in ("-o"):
		out = a
	elif o in ("-s"):
		snp_site = int(a)
	elif o in ("-t"):
		target = a
	elif o in ("-a"):
		alt_allele = a
	else:
		assert False, "unhandled option"

########################
getkasp_path = os.path.expanduser(getkasp_path)
#primer3_path = getkasp_path + "/bin/primer3_core"
primer3_parameter_path = getkasp_path + "/primer3web_v4_JZ.txt"

#from sys import platform
if sys.platform.startswith('linux'): # linux
	primer3_path = getkasp_path + "/bin/primer3_core"
elif sys.platform == "win32" or sys.platform == "cygwin": # Windows...
	primer3_path = getkasp_path + "/bin/primer3_core.exe"

if not out:
	out = "selected_primers.txt"
outfile = open(out, 'w') # output file

snp_site = snp_site - 1 # 0-based

########################
# read alignment file
fasta = {} # dictionary for alignment

with open(RawAlignFile) as file_one:
	for line in file_one:
		if line.startswith(">"):
			sequence_name = line.rstrip().lstrip(">")
		else:
			fasta.setdefault(sequence_name, "")
			fasta[sequence_name] += line.rstrip()


# get the variaiton site among sequences
ids = [] # all other sequence names
for kk in fasta.keys():
	if kk != target:
		ids.append(kk)

print "The other groups: ", ids
	
alignlen = len(fasta[target])
print "Alignment length: ", alignlen

variation = [] # variation sites
ngap = 0 # gap number

for i in range(alignlen):
	nd = 0 # number of difference
	if fasta[target][i] == "-":
		ngap += 1
	for k in ids:
		if fasta[k][i] != fasta[target][i]:
			nd+=1
	if nd == len(ids): # different from all other sequences
		if i - ngap not in variation: # in case multiple gaps
			variation.append(i - ngap)

#############
# loop to write primer3 input for each variation site
# primer3 inputfile
primer3input = "primer3.input"
p3input = open(primer3input, 'w')

seq_template = fasta[target].replace("-","") # remove all gaps

# because A and T give lower Tm, so use them as template
if alt_allele in "ATat":
	seq_template = seq_template[:snp_site] +  alt_allele + seq_template[snp_site + 1:]

for i in variation:
	if i < snp_site:
		left_end = i
		right_end = snp_site
	else:
		left_end = snp_site
		right_end = i
	if right_end - left_end > product_max - 35: # suppose both primers are 18 bp
		continue 
	settings = "PRIMER_TASK=generic" + "\n" + \
	"SEQUENCE_ID=" + target + "-" + str(i) + "\n" + \
	"SEQUENCE_TEMPLATE=" + seq_template + "\n" + \
	"PRIMER_PRODUCT_SIZE_RANGE=50-250" + "\n" + \
	"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + getkasp_path + "/bin/primer3_config/"  + "\n" + \
	"PRIMER_MAX_SIZE=25" + "\n" + \
	"PRIMER_PAIR_MAX_DIFF_TM=2.0" + "\n" + \
	"PRIMER_FIRST_BASE_INDEX=1" + "\n" + \
	"PRIMER_LIBERAL_BASE=1" + "\n" + \
	"PRIMER_NUM_RETURN=5"  + "\n" + \
	"PRIMER_EXPLAIN_FLAG=1"  + "\n" + \
	"SEQUENCE_FORCE_LEFT_END=" + str(left_end + 1) + "\n" + \
	"SEQUENCE_FORCE_RIGHT_END=" + str(right_end + 1) + "\n" + \
	"="
	p3input.write(settings + "\n")

p3input.close()

# primer3 output file
primer3output = "primer3.output"
#p3cmd = primer3_path + " -default_version=2 -p3_settings_file=" + primer3_parameter_path + " -output=" + primer3output + " " + primer3input
p3cmd = primer3_path + " -default_version=2 -output=" + primer3output + " " + primer3input
print "Primer3 command 1st time: ", p3cmd
call(p3cmd, shell=True)

##########################
# function to get reverse complement
def ReverseComplement(seq):
	# too lazy to construct the dictionary manually, use a dict comprehension
	seq1 = 'ATCGTAGCatcgtagc'
	seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
	return "".join([seq_dict[base] for base in reversed(seq)])
	
# classes
class Primers(object):
	"""A primer set designed by Primer3"""
	def __init__(self):
		self.start = 0
		self.end = 0
		self.length = 0
		self.tm = 0.0
		self.gc = 0.0
		self.anys = 0.0
		self.three = 0.0
		self.hairpin = 0.0
		self.end_stability = 0.0
		self.seq = ""
	def formatprimer(self):
		formatout = "\t".join(str(x) for x in [self.start, self.end, self.length, self.tm, self.gc, self.anys, self.three, self.end_stability, self.hairpin, self.seq, ReverseComplement(self.seq)])
		return(formatout)

class PrimerPair(object):
	"""A pair of primers designed by Primer3"""
	def __init__(self):
		self.left = Primers()
		self.right = Primers()
		self.compl_any = "NA"
		self.compl_end = "NA"
		self.penalty = "NA"
		self.product_size = 0



### parse primer 3 primer check output
primerpairs = {}

with open(primer3output) as infile:
	for line in infile:
		line = line.strip()
		if "SEQUENCE_ID" in line:
			seqid = line.split("=")[1]
			primerpairs[seqid] = PrimerPair()
		elif "PRIMER_PAIR_0_PENALTY" in line:
			primerpairs[seqid].penalty = line.split("=")[1]
		elif "PRIMER_PAIR_0_COMPL_ANY" in line:
			primerpairs[seqid].compl_any = line.split("=")[1]
		elif "PRIMER_PAIR_0_COMPL_END" in line:
			primerpairs[seqid].compl_end = line.split("=")[1]
		elif "PRIMER_PAIR_0_PRODUCT_SIZE" in line:
			primerpairs[seqid].product_size = int(line.split("=")[1])
		elif "PRIMER_LEFT_0_SEQUENCE" in line:
			primerpairs[seqid].left.seq = line.split("=")[1]
		elif "PRIMER_LEFT_0=" in line:
			primerpairs[seqid].left.start = int(line.split("=")[1].split(",")[0])
			primerpairs[seqid].left.length = int(line.split("=")[1].split(",")[1])
			primerpairs[seqid].left.end = primerpairs[seqid].left.start + primerpairs[seqid].left.length - 1
		elif "PRIMER_LEFT_0_TM" in line:
			primerpairs[seqid].left.tm = float(line.split("=")[1])
		elif "PRIMER_LEFT_0_GC_PERCENT" in line:
			primerpairs[seqid].left.gc = float(line.split("=")[1])
		elif "PRIMER_LEFT_0_SELF_ANY_TH" in line:
			primerpairs[seqid].left.anys = float(line.split("=")[1])
		elif "PRIMER_LEFT_0_SELF_END_TH" in line:
			primerpairs[seqid].left.three = float(line.split("=")[1])
		elif "PRIMER_LEFT_0_HAIRPIN_TH" in line:
			primerpairs[seqid].left.hairpin = float(line.split("=")[1])
		elif "PRIMER_LEFT_0_END_STABILITY" in line:
			primerpairs[seqid].left.end_stability = float(line.split("=")[1])
		elif "PRIMER_RIGHT_0_SEQUENCE" in line:
			primerpairs[seqid].right.seq = line.split("=")[1]
		elif "PRIMER_RIGHT_0=" in line:
			primerpairs[seqid].right.start = int(line.split("=")[1].split(",")[0])
			primerpairs[seqid].right.length = int(line.split("=")[1].split(",")[1])
			primerpairs[seqid].right.end = primerpairs[seqid].right.start - primerpairs[seqid].right.length + 1
		elif "PRIMER_RIGHT_0_TM" in line:
			primerpairs[seqid].right.tm = float(line.split("=")[1])
		elif "PRIMER_RIGHT_0_GC_PERCENT" in line:
			primerpairs[seqid].right.gc = float(line.split("=")[1])
		elif "PRIMER_RIGHT_0_SELF_ANY_TH" in line:
			primerpairs[seqid].right.anys = float(line.split("=")[1])
		elif "PRIMER_RIGHT_0_SELF_END_TH" in line:
			primerpairs[seqid].right.three = float(line.split("=")[1])
		elif "PRIMER_RIGHT_0_HAIRPIN_TH" in line:
			primerpairs[seqid].right.hairpin = float(line.split("=")[1])
		elif "PRIMER_RIGHT_0_END_STABILITY" in line:
			primerpairs[seqid].right.end_stability = float(line.split("=")[1])
		

# write to file
outfile = open(out, 'w')
outfile.write("index\tproduct_size\ttype\tstart\tend\tlength\tTm\tGCcontent\tany\t3'\tend_stability\thairpin\tprimer_seq\tReverseComplement\tpenalty\tcompl_any\tcompl_end\n")

#for pp in primerpairs:
for i, pp in primerpairs.items():
	if pp.product_size != 0:
		pl = pp.left
		pr = pp.right
		outfile.write("\t".join([i, str(pp.product_size), "LEFT", pl.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end]) + "\n")
		outfile.write("\t".join([i, str(pp.product_size), "RIGHT", pr.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end]) + "\n")

outfile.close()
