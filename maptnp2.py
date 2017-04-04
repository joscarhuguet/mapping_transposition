from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import os

os.system('ls *.ab1 | sed s/.ab1//g > abi_list.txt')
#1) from a list of ab1 files extract fasta files
list_abis=open("abi_list.txt")
for abi in list_abis:
	handle = open(str(abi).rstrip()+".ab1", "rb")
	for record in SeqIO.parse(handle, "abi"):
		head =record.id
		SeqIO.convert(str(abi.rstrip())+".ab1", "abi", str(head.rstrip())+"_parsed.fasta", "fasta")
	
#concatenate all the fasta files
os.system('cat *_parsed.fasta > multi.fasta')


#blast the concateante with the end transposon
blastn_cline = NcbiblastnCommandline(query="multi.fasta" , subject="tn5_end", evalue=0.001,outfmt=6, out="positives_hits_tranposon.txt", word_size='7')
stdout, stderr = blastn_cline()

#blast the hits with tranposon IR against the CMM genome
os.system('cut -f1 positives_hits_tranposon.txt> positives.txt')
os.system('makeblastdb -in crenshaw_9_12_16_V4.fasta -dbtype nucl')
positives=open("positives.txt")
for positive in positives:
	blastn_cline2 = NcbiblastnCommandline(query=str(positive.rstrip())+"_parsed.fasta" , db="crenshaw_9_12_16_V4.fasta", evalue=0.001,outfmt=6, out=str(positive.rstrip())+"_blast.txt")
	print(blastn_cline2)
	stdout, stderr = blastn_cline2()
#os.system('rm -f *_parsed.fasta')
os.system('cat *_blast.txt > positives_hits_with_tranposonEnd_and_chromosome.txt')
os.system('rm -f *_blast.txt')