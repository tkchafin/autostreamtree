import os
import sys

# returns a decoded genotype given relevant vcf fields
def decode(gt, ref, alts, as_iupac=False, as_tuple=False, as_list=False):
	ret=[None, None]
	if gt[0] is None or gt[1] is None:
		ret = ["N","N"]
	else:
		if gt[0] == 0:
			ret[0] = ref
		else:
			ret[0] = alts[gt[0]-1]
		if gt[1] == 0:
			ret[1] = ref
		else:
			ret[1] = alts[gt[1]-1]

	if as_iupac:
		return(DNAconsensus("/".join(ret)))
	elif as_tuple:
		return(tuple(ret))
	elif as_list:
		return(ret)
	else:
		return("/".join(ret))

#function to compute nucleotide frequencies
#if ploidy = 1, ambiguities will be skipped
#if ploidy = 2, ambiguities will be resolved
def getNucFreqs(seqs, ploidy):
	freqs = list()
	for loc in range(0, len(seqs[list(seqs.keys())[0]])):
		temp = dict()
		allnucs = ""
		for samp in seqs.keys():
			allnucs = allnucs + DNAconsensus(seqs[samp][loc]).lower()
		badchars = ["?", "-", "n"]
		if ploidy == 1:
			badchars = ["?", "-", "n", "r", "y", "s", "w", "k", "m", "b", "d", "h", "v"]
		for char in badchars:
			allnucs = allnucs.replace(char, "")
		if ploidy == 2:
			iupacs = ["r", "y", "s", "w", "k", "m", "b", "d", "h", "v"]
			for ambig in iupacs:
				allnucs = allnucs.replace(ambig, "".join(get_iupac_caseless(ambig)))
			for nuc in ["a", "c", "t", "g"]:
				allnucs = allnucs.replace(char, char+char)
		total = len(allnucs)
		counts = {"a":0, "g":0,"c":0,"t":0}
		for c in allnucs:
			if c in counts:
				counts[c] += 1
		for nuc in counts.keys():
			counts[nuc] = float(counts[nuc]/total)
		freqs.append(counts)
	return(freqs)

#function to make a consensus from alleles separated by "/"
def DNAconsensus(seq):
	alleles = seq.split("/")
	consens = ""
	if len(alleles) < 1:
		return(None)
	elif not all(len(x) == len(alleles[0]) for x in alleles):
		print("ERROR: Not all alleles are the same length:",alleles)
		sys.exit(1)
	elif len(alleles) == 1:
		return(alleles[0])
	else:
		for i in range(len(alleles[0])):
			nucs = ""
			for a in alleles:
				nucs += a[i]
			temp = listToSortUniqueString(nucs.upper())
			consens += reverse_iupac_case(temp)
	#print(consens)
	return(consens)

#Function to translate a string of bases to an iupac ambiguity code, retains case
def reverse_iupac_case(char):
	iupac = {
		'A':'A',
		'N':'N',
		'-':'-',
		'C':'C',
		'G':'G',
		'T':'T',
		'AG':'R',
		'CT':'Y',
		'AC':'M',
		'GT':'K',
		'AT':'W',
		'CG':'S',
		'CGT':'B',
		'AGT':'D',
		'ACT':'H',
		'ACG':'V',
		'ACGT':'N',
		'a':'a',
		'n':'n',
		'c':'c',
		'g':'g',
		't':'t',
		'ag':'r',
		'ct':'y',
		'ac':'m',
		'gt':'k',
		'at':'w',
		'cg':'s',
		'cgt':'b',
		'agt':'d',
		'act':'h',
		'acg':'v',
		'acgt':'n'
	}
	return iupac[char]

#Function to return sorted unique string from list of chars
def listToSortUniqueString(l):
	sl = sorted(set(l))
	return(str(''.join(sl)))

#Function to split character to IUPAC codes, assuing diploidy
def get_iupac_caseless(char):
	lower = False
	if char.islower():
		lower = True
		char = char.upper()
	iupac = {
		"A"	: ["A"],
		"G"	: ["G"],
		"C"	: ["C"],
		"T"	: ["T"],
		"N"	: ["A", "C", "G", "T"],
		"?"	: ["A", "C", "G", "T"],
		"-"	: ["A", "C", "G", "T", "-"],
		"R"	: ["A","G"],
		"Y"	: ["C","T"],
		"S"	: ["G","C"],
		"W"	: ["A","T"],
		"K"	: ["G","T"],
		"M"	: ["A","C"],
		"B"	: ["C","G","T"],
		"D"	: ["A","G","T"],
		"H"	: ["A","C","T"],
		"V"	: ["A","C","G"]
	}
	ret = iupac[char]
	if lower:
		ret = [c.lower() for c in ret]
	return ret

def phaseSnp(snp):
	nucs = get_iupac_caseless(snp)
	if len(nucs) > 2 or len(nucs) < 1:
		return("n/n")
	elif len(nucs) == 1:
		s = nucs[0] + '/' + nucs[0]
		return(s)
	else:
		return("/".join(nucs))
