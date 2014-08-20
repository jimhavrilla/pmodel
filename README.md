pmodel
======

For modeling genetic variation based on protein domains.  Extracts 6500ESP vcf SNPs and then intersects with Pfam protein domains.  Finds top and bottom 10 occurring domains by count and by ratio of domain count to domain length.  extractcommands.sh contains all programming work.

Basically I discovered that ccds_id was the culprit in creating a variable number of fields.  Sometimes an entry had one sometimes it didn't.  I formatted everything in python (rearrange.py) to keep tabs and spaces where appropriate,  pulled the important pfamA_id to the front of the INFO field, and pushed ccds_id (if it existed) to the back.  I also added an interval length column to make calculating the total nucleic acid distance for each domain easier when I divided the number of variants for a domain by that number later (in divide.py).  I also have a version of the chromosome counts for just ratio counts that cuts out the fields after pfamA_id after the rearrange to get the data we care about.
