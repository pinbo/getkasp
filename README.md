# GetKasp
Get genomic specific KASPars from user provided sequences.

Polymarker (http://polymarker.tgac.ac.uk/) is a great tool to design genomic specific KASPars. But sometimes I would like to get extra primers or made modifications or for different genomes. That is why I wrote this python script to do that for my own purpose. I do not know whether it is useful for others. You can also directly clone and modify the polymarker scripts (https://github.com/TGAC/bioruby-polyploid-tools) to meet your own requirements.

# Dependencies

GetKasp use "Primer3" to design primers.

1. Primer3: program for designing PCR primers (http://primer3.sourceforge.net/).

# How it works
1. Find all the different sites that can differ all other sequences from the user provided alignment file;
2. Use these sites and the SNP site as SEQUENCE_FORCE_RIGHT_END in primer3 to design all possible left and right primers in the target sequence.

# Usage
```
getkasp.py 
	-i <alignment.fa>
	-p <getkasp path> 
	-o <output file name>"
	-s <SNP position in the raw target sequence (not the alignment sequence)>
	-t <target_sequence_ID>
	-a <altanative allele>
	-h help
```

**Example**
./getkasp.py -i alignment_raw.fa -p . -o primer3008.txt -s 489 -t 2AS -a T

It can be used in Linux (may need to recompile primer3_core for your specific Linux version) and windows 7 (terminal or Cygwin). The INPUT file is an multipe sequence alignment file in fasta format.

# Credits
I borrowed ideas from the polymarker scripts (https://github.com/TGAC/bioruby-polyploid-tools), a great tool for Genome Specific KASPar design in polyploid species. Thanks to the author of Polymarker.

I also borrowed some codes from biopython (https://github.com/biopython/biopython/blob/master/Bio/Emboss/Primer3.py). Thanks to them too.

Thanks to the open source software **Primer3** (http://primer3.sourceforge.net/).

