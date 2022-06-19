# Script that do something with the nucleotides sequence/s of the QBeta Bacteriophage

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper as h



print("\nLoading paths and variables")

# Remember that index starts at zero, so we substract 1 from all the positions!
# Dictionary containing the position of each coding region
ORF = {"A2": [60, 1322], "A1": [1343, 2332], "Coat": [1343, 1744], "Replicase": [2351, 4120]}

# Dictionary containing the Reference sequences of each region (provided by Pilar)
refSeqs = {"R1": "CAACAAGGTCAGCTATATCATAATATCGATATTGTAGACGGCTTTGACAGACGTGACATCCGGCTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATACGAGCAATCTTCCGTTCGCTACACTAGATCTTGATACTACCTTTAGTTCGTTTAAACACGTTCTTGATAGTATCTTTTTATTAACCCAACGCGTAAAGCGTTGAAACTTTG",
			"R2": "CGATTCATTTATCTTAAGTCGATAAATGCTTATTGCTCTCTTAGCGATATTGCGGCCTATCACGCCGATGGCGTGATAGTTGGCTTTTGGCGCGATCCATCCAGTGGTGGTGCCATACCGTTTGACTTCACTAAGTTTGATAAGACTAAATGTCCTATTCAAGCCGTGATAGTCGTTCCTCGTGCTTAGTAACTAAGGATGAAATGCATGTCTAAGACAGCATCTTCGCGTAACTCTCTCAGCGCACAATTGCGCCGAGCCGCGAACACAAGAATTGAGGTTGAAGGTAACCTCGCACTTTCCATTGCCAACGATTTACTGTTGGCCTA",
			"R3": "TTACACATTCGAGCTCGAGTCGCTTATTTTTGCTTCTCTCGCTCGTTCCGTTTGTGAGATACTGGACTTAGACTCGTCTGAGGTCACTGTTTACGGAGACGATATTATTTTACCGTCCTGTGCAGTCCCTGCCCTCCGGGAAGTTTTTAAGTATGTTGGTTTTACGACCAATACTAAAAAGACTTTTTCCGAGGGGCCGTTCAGAGAGTCGTGCGGCAAGCACTACTATTCTGGCGTAGATGTTACTCCCTTTTACATACGTCACCGTATAGTGA"}

# Dictionary containing the position for the reference sequences, obtained from NCBI
refSeqsPos = {"R1": [1059, 1329], "R2": [2144, 2472], "R3": [3328, 3602]}

# Translation table 11: The Bacterial, Archaeal and Plant Plastid Code from NCBI
transl_table = {"F": ["TTT", "TTC"], "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
				"I": ["ATT", "ATC", "ATA"], "M": ["ATG"], "V": ["GTT", "GTC", "GTA", "GTG"],
				"S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "P": ["CCT", "CCC", "CCA", "CCG"],
				"T": ["ACT", "ACC", "ACA", "ACG"], "A": ["GCT", "GCC", "GCA", "GCG"],
				"Y": ["TAT", "TAC"], "*": ["TAA", "TAG", "TGA"],
				"H": ["CAT", "CAC"], "Q": ["CAA", "CAG"],
				"N": ["AAT", "AAC"], "K": ["AAA", "AAG"],
				"D": ["GAT", "GAC"], "E": ["GAA", "GAG"],
				"C": ["TGT", "TGC"], "W": ["TGG"],
				"R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
				"G": ["GGT", "GGC", "GGA", "GGG"]}

# Reference amino acid sequences
A2_prot_NCBI = "MPKLPRGLRFGADNEILNDFQELWFPDLFIESSDTHPWYTLKGRVLNAHLDDRLPNVGGRQVRRTPHRVTVPIASSGLRPVTTVQYDPAALSFLLNARVDWDFGNGDSANLVINDFLFRTFAPKEFDFSNSLVPRYTQAFSAFNAKYGTMIGEGLETIKYLGLLLRRLREGYRAVKRGDLRALRRVIQSYHNGKWKPATAGNLWLEFRYGLMPLFYDIRDVMLDWQNRHDKIQRLLRFSVGHGEDYVVEFDNLYPAVAYFKLKGEITLERRHRHGISYANREGYAVFDNGSLRPVSDWKELATAFINPHEVAWELTPYSFVVDWFLNVGDILAQQGQLYHNIDIVDGFDRRDIRLKSFTIKGERNGRPVNVSASLSAVDLFYSRLHTSNLPFATLDLDTTFSSFKHVLDSIFLLTQRVKR*"
Coat_prot_NCBI = "MAKLETVTLGNIGKDGKQTLVLNPRGVNPTNGVASLSQAGAVPALEKRVTVSVSQPSRNRKNYKVQVKIQNPTACTANGSCDPSVTRQAYADVTFSFTQYSTDEERAFVRTELAALLASPLLIDAIDQLNPAY*"
Replicase_prot_NCBI = "MSKTASSRNSLSAQLRRAANTRIEVEGNLALSIANDLLLAYGQSPFNSEAECISFSPRFDGTPDDFRINYLKAEIMSKYDDFSLGIDTEAVAWEKFLAAEAECALTNARLYRPDYSEDFNFSLGESCIHMARRKIAKLIGDVPSVEGMLRHCRFSGGATTTNNRSYGHPSFKFALPQACTPRALKYVLALRASTHFDIRISDISPFNKAVTVPKNSKTDRCIAIEPGWNMFFQLGIGGILRDRLRCWGIDLNDQTINQRRAHEGSVTNNLATVDLSAASDSISLALCELLLPPGWFEVLMDLRSPKGRLPDGSVVTYEKISSMGNGYTFELESLIFASLARSVCEILDLDSSEVTVYGDDIILPSCAVPALREVFKYVGFTTNTKKTFSEGPFRESCGKHYYSGVDVTPFYIRHRIVSPADLILVLNNLYRWATIDGVWDPRAHSVYLKYRKLLPKQLQRNTIPDGYGDGALVGSVLINPFAKNRGWIRYVPVITDHTRDRERAELGSYLYDLFSRCLSESNDGLPLRGPSGCDSADLFAIDQLICRSNPTKISRSTGKFDIQYIACSSRVLAPYGVFQGTKVASLHEA*"
A1_prot_NCBI = "MAKLETVTLGNIGKDGKQTLVLNPRGVNPTNGVASLSQAGAVPALEKRVTVSVSQPSRNRKNYKVQVKIQNPTACTANGSCDPSVTRQAYADVTFSFTQYSTDEERAFVRTELAALLASPLLIDAIDQLNPAYWTLLIAGGGSGSKPDPVIPDPPIDPPPGTGKYTCPFAIWSLEEVYEPPTKNRPWPIYNAVELQPREFDVALKDLLGNTKWRDWDSRLSYTTFRGCRGNGYIDLDATYLATDQAMRDQKYDIREGKKPGAFGNIERFIYLKSINAYCSLSDIAAYHADGVIVGFWRDPSSGGAIPFDFTKFDKTKCPIQAVIVVPRA*"


# Define paths
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/data"
dataPathOut = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/"

print("\nPaths and variables loaded!\n")

print("\nLoading reference genome as fasta")

fIn = os.path.join(dataPathIn, "Qbeta" + ".fasta")

sequences = h.LoadFasta(fIn)

seq = list(sequences.values())[0]

print("\nReference genome loaded!")

print("\nCalculating the position of non-coding regions for R1 and R2\n")

print("\tLength of sequences")
for k, v in refSeqs.items():
	print("\t", k, ":", len(v))

# For R1 the start and end positions of the non-coding region is:
startPos_NC_R1 = ORF["A2"][1]
endPos_NC_R1 = refSeqsPos["R1"][1]
# Extract the sequence
seq_NC_R1 = h.ExtractSequence(seq, startPos_NC_R1, endPos_NC_R1)

# Get the end position for A2 from the reference sequence
endPos_A2_R1 = refSeqs["R1"].find(seq_NC_R1)
# Extract the sequence of A1 from the reference Sequence
seq_A2_R1 = refSeqs["R1"][0:endPos_A2_R1]
print("\tThe position of A2 (0-based index) for R1 is:", 0, "-", endPos_A2_R1)
print("\tThe non-coding position of A2 (0-based index) for R1 is:", endPos_A2_R1 + 1, "-", len(refSeqs["R1"]) - 1)


# For R2 the start and end positions of the replicase is:
startPos_Replicase_R2 = ORF["Replicase"][0]
endPos_Replicase_R2 = refSeqsPos["R2"][1]
# Extract the sequence
seq_R2_Replicase = h.ExtractSequence(seq, startPos_Replicase_R2, endPos_Replicase_R2)
##codons_sub, codonsPos_sub = ExtractSeqCodons(seq, startPos_Replicase_R2, endPos_Replicase_R2)

# Get the end position of seq_R2_Replicase from the reference sequence
endPos_NC_R2 = refSeqs["R2"].find(seq_R2_Replicase)
# Extract the sequence of the non-coding part of R2 from the reference sequence
seq_R2_NC = refSeqs["R2"][0:endPos_NC_R2]
print("\tThe position of Replicase (0-based index) for R2 is:", endPos_NC_R2, len(refSeqs["R2"]) - 1)

# Get the non-coding positions between Replicase and A1 
endPos_A1_Replicase = ORF["Replicase"][0] - 1
startPos_A1_Replicase = ORF["A1"][1] + 1
##print("Heere", startPos_A1_Replicase, "-", endPos_A1_Replicase)
# Extract the non-coding sequence between A1 and Replicase
seq_A1_Replicase_NC = h.ExtractSequence(seq, startPos_A1_Replicase, endPos_A1_Replicase)
##print("SEEQ", seq_A1_Replicase_NC, len(seq_A1_Replicase_NC))

# Get the start position for the non-coding sequence between A1 and Replicase
startPos_A1_Replicase_NC = refSeqs["R2"].find(seq_A1_Replicase_NC)
# Extract the sequence of the non-coding sequence between A1 and Replicase
seq_A1_Replicase_NC_refSeq = refSeqs["R2"][startPos_A1_Replicase_NC:startPos_A1_Replicase_NC + len(seq_A1_Replicase_NC)]
##print("SEEQ2", seq_A1_Replicase_NC_refSeq, len(seq_A1_Replicase_NC_refSeq))
print("\tThe non-coding position between A1 and Replicase (0-based index) for R2 is:", startPos_A1_Replicase_NC, "-", startPos_A1_Replicase_NC + len(seq_A1_Replicase_NC_refSeq) - 1)
print("\tThe position of A1 (0-based index) for R2 is:", 0, "-", startPos_A1_Replicase_NC - 1)

print("\nCalculated!\n")

print("\nBuild a dictionary for starting and ending positions of coding sequences for sequenced regions\n")

ORF_sequenced = {"A2_R1": [0, endPos_A2_R1], "A1_R2": [0, startPos_A1_Replicase_NC -1],
"Replicase_R2": [endPos_NC_R2, len(refSeqs["R2"]) - 1], "Replicase_R3": [0, len(refSeqs["R3"]) - 1]}

print("\nDictionary build!")

for k, v in ORF_sequenced.items():
	print(k, v)

###########################################################################

print("\nCalculating the codon positions of a region")

# For R1
startPos_A2 = refSeqsPos["R1"][0]
endPos_A2 = ORF["A2"][1]

# For R2
starPos_A1_Replicase = ORF["Replicase"][0]
endPos_A1_Replicase = refSeqsPos["R3"][1]

# For R3
startPos_Replicase = refSeqsPos["R3"][0]
endPos_Replicase = refSeqsPos["R3"][1]

#codons, codonsPos = SplitSeq(seq, ORF, "Replicase")
#codons_R2, codonsPos_R2 = CheckCodonPosition(seq, starPos_A1_Replicase, endPos_A1_Replicase, codonsPos)

print("\nCodon positions of a region calculated!\n")

# Translating sequences
print("\nChecking translation of reference sequences\n")

#Extract sequence for A2
A2 = h.ExtractSequence(seq, ORF["A2"][0], ORF["A2"][1])

# Get the amino acid sequence
A2_prot = h.Translation(A2, transl_table)

print("\tLength A2 is the same?", len(A2_prot) == len(A2_prot_NCBI))

print("\tAre the same A2 sequence?", A2_prot == A2_prot_NCBI)

# Extract sequence for coat protein
Coat = h.ExtractSequence(seq, ORF["Coat"][0], ORF["Coat"][1])

# Get the amino acid sequence
Coat_prot = h.Translation(Coat, transl_table)

print("\tLength Coat is the same?", len(Coat_prot) == len(Coat_prot_NCBI))

print("\tAre the same Coat sequence?", Coat_prot == Coat_prot_NCBI)

# Extract the sequence for Replicase
Replicase = h.ExtractSequence(seq, ORF["Replicase"][0], ORF["Replicase"][1])

# Get the amino acid sequence
Replicase_prot = h.Translation(Replicase, transl_table)

print("\tLength Replicase is the same?", len(Replicase_prot) == len(Replicase_prot_NCBI))

print("\tAre the same Replicase sequence?", Replicase_prot == Replicase_prot_NCBI)

# Extract the sequence for A1
A1 = h.ExtractSequence(seq, ORF["A1"][0], ORF["A1"][1])

# Get the amino acid sequence
A1_prot = h.Translation(A1, transl_table, leaky = True)

print("\tLength A1 is the same?", len(A1_prot) == len(A1_prot_NCBI))

print("\tAre the same A1 sequence?", A1_prot == A1_prot_NCBI)

##########################################################################

print("\nAll translations ok!")


print("\nTranslating Region sequences")

#codons_Replicase_r, codonsPos_Replicase_r = h.SplitSeq(seq, ORF, "Replicase")
#codons_Replicase = h.CheckCodonPosition(refSeqs["R3"], refSeqsPos["R3"][0], codonsPos_Replicase_r, ORF_sequenced["Replicase_R3"][0], ORF_sequenced["Replicase_R3"][1])
#codons_Replicase_ref, codonsPos_Replicase_ref = h.SplitSeq(refSeqs["R3"], ORF_sequenced, "Replicase_R3")
#prot_Replicase_sequenced = h.TranslationCodons(codons_Replicase, transl_table) 
#seq_Replicase_sequenced = refSeqs["R3"][ORF_sequenced["Replicase_R3"][0]:ORF_sequenced["Replicase_R3"][1] + 1]
#seq_Replicase = seq[ORF["Replicase"][0]:ORF["Replicase"][1] + 1]
#startIndex_genome = seq.find(seq_Replicase_sequenced)
#startIndex_Replicase = seq_Replicase.find(seq_Replicase_sequenced)
#startIndex_prot = Replicase_prot_NCBI.find(prot_Replicase_sequenced)
#print("Index", startIndex_prot)
#print(prot_Replicase_sequenced)
#print(Replicase_prot_NCBI)
#print(seq_Replicase)
#print(seq_Replicase_sequenced)



#codons_Replicase_r, codonsPos_Replicase_r = h.SplitSeq(seq, ORF, "Replicase")
#codons_Replicase = h.CheckCodonPosition(refSeqs["R2"], ORF["Replicase"][0], codonsPos_Replicase_r, ORF_sequenced["Replicase_R2"][0], ORF_sequenced["Replicase_R2"][1])
#codons_Replicase_ref, codonsPos_Replicase_ref = h.SplitSeq(refSeqs["R2"], ORF_sequenced, "Replicase_R2")
#prot_Replicase_sequenced = h.TranslationCodons(codons_Replicase, transl_table) 
#seq_Replicase_sequenced = refSeqs["R2"][ORF_sequenced["Replicase_R2"][0]:ORF_sequenced["Replicase_R2"][1] + 1]
#seq_Replicase = seq[ORF["Replicase"][0]:ORF["Replicase"][1] + 1]
#startIndex_genome = seq.find(seq_Replicase_sequenced)
#startIndex_Replicase = seq_Replicase.find(seq_Replicase_sequenced)
#startIndex_prot = Replicase_prot_NCBI.find(prot_Replicase_sequenced)
#print(prot_Replicase_sequenced)
#print(seq_Replicase)
#print(seq_Replicase_sequenced)


#codons_A1_r, codonsPos_A1_r = h.SplitSeq(seq, ORF, "A1")
#codons_A1 = h.CheckCodonPosition(refSeqs["R2"], refSeqsPos["R2"][0], codonsPos_A1_r, ORF_sequenced["A1_R2"][0], ORF_sequenced["A1_R2"][1])
#codons_A1_ref, codonsPos_A1_ref = h.SplitSeq(refSeqs["R2"], ORF_sequenced, "A1_R2")
#prot_A1_sequenced = h.TranslationCodons(codons_A1, transl_table) 
#seq_A1_sequenced = refSeqs["R2"][ORF_sequenced["A1_R2"][0]:ORF_sequenced["A1_R2"][1] + 1]
#seq_A1 = seq[ORF["A1"][0]:ORF["A1"][1] + 1]
#startIndex_genome = seq.find(seq_A1_sequenced)
#startIndex_A1 = seq_A1.find(seq_A1_sequenced)
#startIndex_prot = A1_prot_NCBI.find(prot_A1_sequenced)
#print(seq_A1)
#print(seq_A1_sequenced)


#print(prot_A1_sequenced)

codons_A2_r, codonsPos_A2_r = h.SplitSeq(seq, ORF, "A2")
codons_A2 = h.CheckCodonPosition(refSeqs["R1"], refSeqsPos["R1"][0], codonsPos_A2_r, ORF_sequenced["A2_R1"][0], ORF_sequenced["A2_R1"][1])
codons_A2_ref, codonsPos_A2_ref = h.SplitSeq(refSeqs["R1"], ORF_sequenced, "A2_R1")
prot_A2_sequenced = h.TranslationCodons(codons_A2, transl_table) 
seq_A2_sequenced = refSeqs["R1"][ORF_sequenced["A2_R1"][0]:ORF_sequenced["A2_R1"][1] + 1]
seq_A2 = seq[ORF["A2"][0]:ORF["A2"][1] + 1]
startIndex_genome = seq.find(seq_A2_sequenced)
startIndex_A2 = seq_A2.find(seq_A2_sequenced)
print("Len full seq protein", len(seq_A2))
print("\tAt index", seq_A2.find(seq_A2_sequenced))
print("Len sequenced", len(seq_A2_sequenced))
print("\tAt index", seq.find(seq_A2_sequenced))
startIndex_prot = A2_prot_NCBI.find(prot_A2_sequenced)
print(seq_A2)
print(seq_A2_sequenced)
print(prot_A2_sequenced)


##print("\tChecking Translation\n")
##if prot_Replicase_sequenced in Replicase_prot:
	##print("\tTranslatiton sucessful")
	##print("\tAt index", Replicase_prot.find(prot_Replicase_sequenced))
	##print("\tLenght of protein sequenced", len(prot_Replicase_sequenced))
	##print(Replicase_prot, "\n")
	##print(prot_Replicase_sequenced)
##else:
	##print("\tError")
	##print(Replicase_prot, "\n")
	##print(prot_Replicase_sequenced)
##
##print('Checking codons of Replicase from sequencing R3')
##for k, v in codons_Replicase.items():
	##print("Codon", k, v)
#
# Saving info
reg = "R1"
protName = "A2"
fOut = os.path.join(dataPathOut + "CodonInfo_" + protName + "_" + reg + ".csv")

h.saveCodonsInfo(fOut, codons_A2, prot_A2_sequenced, startIndex_genome, startIndex_A2, startIndex_prot, replicase_R2=False)

###############################################################################

