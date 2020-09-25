from collections import defaultdict
import human_application_functions as func
import numpy as np
import tabix
import sys

chromosome_number = sys.argv[1]
yri="NA18486"
gbr="HG00148"

yri_index = 1764
gbr_index = 54 #header.index(gbr)

thousand_genomes_path="/branchinecta/shared/humans/"  #thousandgenomes/"
#archaic_data_path="/branchinecta/shared/humans/archaic/"
archaic_data_path="/courses/ehuerta-sanchez/AltaiNea_Den_combined/"
human_infile=thousand_genomes_path #+"chr"+str(chromosome_number)+".biallelic.vcf.gz"
human_infile+="ALL.chr"+chromosome_number+".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
archaic_infile=archaic_data_path+"AltaiNea_Den_combined."+str(chromosome_number)+".vcf.gz"
outfile_path ="/branchinecta/llopez84/intermediate_data_yri-gbr/"

#Open files with tabix
atb = tabix.open(archaic_infile)
tb = tabix.open(human_infile)
#atb = tabix.open(archaic_infile)

#Open outfile
outfile=outfile_path+"fraw_chr-"+chromosome_number+".yri-"+yri+".gbr-"+gbr+".nean.vcf"
fout=open(outfile,"w")
header="#CHROM\tPOS\tNA18486\tHG00148\tNEAN\tCHIMP\n"
fout.write(header)

#Replace with another function?
nucleotides=["A","T","C","G"]

window_size=100000
#Get first and last position of the chromosome
first_last_position_per_chrom=func.get_first_last_position_per_chrom()
start = first_last_position_per_chrom[chromosome_number][0]
stop = first_last_position_per_chrom[chromosome_number][1]

while start < stop:
  human_iterator = tb.query(str(chromosome_number),start,(start+window_size))
  archaic_iterator = atb.query(str(chromosome_number),start,(start+window_size))

  data = defaultdict(defaultdict)
  for spline in human_iterator:
    #is the snp biallelic
    if(nucleotides.count(spline[3])==1 and nucleotides.count(spline[4])==1):
      #Does pop1 or pop2 have a derived and ancestral allele
      if("0" in spline[yri_index]+spline[gbr_index] and "1" in spline[yri_index]+spline[gbr_index]):
        #Get data
        data[spline[1]].update({"human reference":spline[3],"human alt":spline[4],"yri":spline[yri_index],
        "gbr":spline[gbr_index]})

  #Get neanderthal and denisovan data from the biallelic human SNPs
  for spline in archaic_iterator:
    #Is there data missing for neanderthals
    if "." in spline[9].split(":")[0]:
      continue
    #is it from biallelic human SNPs
    if(spline[1] in list(data.keys())):
      data[spline[1]].update({"archaic reference":spline[3],"archaic alt":spline[4],
      "chimp allele":func.get_chimp_allele(spline[7]),"neanderthal":spline[9].split(":")[0]})

  #Go through data from window:
  for key in list(data.keys()):
    #Chimp allele has one nucleotide
    try:
      if(nucleotides.count(data[key]["chimp allele"]) != 1):
        continue
    #Delete positions that were not found in archaic data
    except KeyError:
      del data[key]
      continue
    #Is this position biallelic across all four populations
    if not func.biallelic_four_populations(data[key]["human reference"],data[key]["human alt"],data[key]["archaic reference"],data[key]["archaic alt"],data[key]["chimp allele"]):
      continue
    #Assert that the human and archaic reference are the same
    assert(data[key]["human reference"]==data[key]["archaic reference"]),"Human reference different from archaic reference at position: "+key
    #Are the alleles ancestral or derived
    outstring = chromosome_number+"\t"+key+"\t"+func.get_ancestral_derived(data[key]["human reference"],data[key]["human alt"],data[key]["chimp allele"],data[key]["yri"])
    outstring+="\t"+func.get_ancestral_derived(data[key]["human reference"],data[key]["human alt"],data[key]["chimp allele"],data[key]["gbr"])
    outstring+="\t"+func.get_ancestral_derived(data[key]["human reference"],data[key]["human alt"],data[key]["chimp allele"],data[key]["neanderthal"])
    outstring+="\t"+"A" #Ancestral allele will always be A
    fout.write(outstring+"\n")

  #Move onto next window of the chromosome
  start += (window_size + 1) 

fout.close()

print("Ran through chromosome "+chromosome_number)

