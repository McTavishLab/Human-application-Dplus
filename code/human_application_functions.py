import numpy as np
import pandas as pd
import math
from collections import defaultdict

def get_chimp_allele(string):
  index=string.find("CAnc")
  if index < 0:
    return("")
  substring=string[index:]
  chimp_allele=substring.split(";")[0].split("=")[1]
  return(chimp_allele)


#Function to cut introgressed region into window sizes
def get_window_cut_introgressed_region(chromosome,start,stop,window_size,pos_dict):
    results=[]
    #Greatest lower bound and Least upper bound for the windows per chromosome
    gub=pos_dict[chromosome][0];lub=pos_dict[chromosome][1]
    #If introgressed region within first window or last window
    if(stop < (gub + window_size) or start > (lub - window_size)):
        return None
    #Discard the first window of the chromosome
    if(start < (gub + window_size)):
        start = (gub + window_size)
    while start+window_size <= stop:
        #Discard the last window of the chromosome
        if(start > (lub - window_size)):
            break
        results.append([str(chromosome),str(start),str(start+window_size),str(start+window_size-start)])
        start=start+window_size + 1
    leftover = stop-(int(results[-1][2]))
    return(results, leftover)

# Function to get a bigger window with introgressed region within that window
def get_window_with_intro_region_within(chromosome,start,stop,window_size,pos_dict):
    middle_point=(stop-start)/2.
    #Greatest lower bound and Least upper bound for the windows per chromosome
    gub=pos_dict[chromosome][0];lub=pos_dict[chromosome][1]
    if(start < (gub + window_size) or stop > (lub - window_size)):
        return None
    half_window_size=window_size/2 #window size assumed to be even
    new_start = math.floor(middle_point - half_window_size)
    new_stop = math.floor(middle_point + half_window_size)
    #If the start is too close to the first window
    if(new_start < (gub+window_size)):
        new_start = (gub + window_size)
        new_stop = new_start + window_size
    elif(new_stop > (lub - window_size)):
        new_stop = (lub - window_size)
        new_start = (lub - window_size) - window_size
    return([str(int(chromosome)),str(new_start),str(new_stop),str(new_stop-new_start)])

# Function to get windows of introgressed regions, both introgressed regions shorter and larger than the window size
def get_window_from_introgressed_region(chromosome,start,stop,window_size,physical_length,pos_dict,outfile):
    #Open outfile
    file = open(outfile,"a")
    if(physical_length < window_size):
        #Make window bigger than intro region
        results = get_window_with_intro_region_within(chromosome,start,stop,window_size,pos_dict)
        leftover=0
        file.write("\t".join(results)+"\n")
    #If the introgressed region is bigger than the window size, cut into window sizes
    elif(physical_length >= window_size):
        results,leftover = get_window_cut_introgressed_region(chromosome,start,stop,window_size,pos_dict)
        for row in results:
            file.write("\t".join(row)+"\n")
    else:
        print("Error from within get_window_from_introgressed_region. Physical length being recognized")
    file.close()
    return leftover

# Function to get the first and last position of each chromosome
def first_last_pos_per_chromosome():
  start_stop_position_chromosomes_file=thousand_genomes_path+"start_stop_position_autosomal_chromosomes.csv"
  start_stop_pos_chrom_dicc = defaultdict(None)
  file = open(start_stop_position_chromosomes_file,"r")
  line=file.readline()
  for line in file:
    spline = line.split(",")
    start_stop_pos_chrom_dicc[int(spline[0])] = (int(spline[1]),int(spline[2].strip()))
  file.close()
  return(start_stop_pos_chrom_dicc)


def get_ancestral_derived(ref,alt,ancestral,alleles):
  """Get whether the alleles are ancestral or derived
  Input: ref: string with the reference allele, alt: string with the alternative allele, ancestral: string with allele from
         outgroup, alleles: string of the haplotype or haplotypes ("0","0|1","1/1" with 0s and 1s
  Output:string with the alleles denoted as A for ancestral allele and B for derived allele. ex: "A","A|B"
  """
  if(ref==ancestral):
    processed_alleles=alleles.replace("0","A")
    processed_alleles=processed_alleles.replace("1","B")
  elif(alt==ancestral):
    processed_alleles=alleles.replace("1","A")
    processed_alleles=processed_alleles.replace("0","B")
  return(processed_alleles)


def biallelic_four_populations(ref,human_alt,archaic_ref,archaic_alt,ancestral_allele):
  """Check whether the SNP at a position is biallelic across four populations
  Input: ref: string with the nucleotide of the reference allele, human_alt: string with alternative allele for humans,
         archaic_alt: string with alternative allele for archaic human, ancestral_allele: string with ancestral allele
         from outgroup
  Output: boolean, True if SNP is biallelic, False if SNP is not biallelic
  """
  nucleotides = ["A","T","C","G"]
  if(archaic_alt == "."):
    alleles=np.asarray((ref,human_alt,archaic_ref,ancestral_allele))
  else:
    alleles=np.asarray((ref,human_alt,archaic_ref,archaic_alt,ancestral_allele))
  #Are there only two unique alleles and are the alleles the nucleotides?
  return( len(np.unique(alleles))==2 and not False in np.in1d(alleles,nucleotides)   )

def get_first_last_position_per_chrom():
  ##Open the start stop position file for all chromosomes
  first_last_position_file="../start_stop_position_autosomal_chromosomes.csv"
  first_last_position_per_chrom=defaultdict()
  with open(first_last_position_file,"r") as file:
    for line in file:
      split_line=line.strip().split(",")
      try:
        first_last_position_per_chrom[split_line[0]]=(int(split_line[1]),int(split_line[2]))
      except ValueError:
        continue
  return (first_last_position_per_chrom)
