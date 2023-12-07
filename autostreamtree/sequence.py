from typing import List


def decode(gt: tuple, ref: str, alts: list, as_iupac: bool = False,
           as_tuple: bool = False, as_list: bool = False) -> str:
    """
    Decode a genotype from a VCF file and return the corresponding DNA sequence

    Args:
        gt (tuple): The genotype as a tuple of integers representing the
                    indices of the reference and alternate alleles.
        ref (str): The reference allele.
        alts (list): A list of alternate alleles.
        as_iupac (bool): Whether to return the consensus IUPAC symbol if the
                         genotype is heterozygous. Default is False.
        as_tuple (bool): Whether to return the result as a tuple of two strings
                         Default is False.
        as_list (bool): Whether to return the result as a list of two strings.
                        Default is False.

    Returns:
        ret (str): The decoded genotype as a string. By default, this is a
                   string in the format "ref/alt", but the format can be
                   customized using the optional arguments.
    """
    # initialize return value
    ret = [None, None]

    # if genotype is missing, set to "N"
    if gt[0] is None or gt[1] is None:
        ret = ["N", "N"]
    else:
        # decode each allele
        if gt[0] == 0:
            ret[0] = ref
        else:
            ret[0] = alts[gt[0] - 1]
        if gt[1] == 0:
            ret[1] = ref
        else:
            ret[1] = alts[gt[1] - 1]

    # apply optional formatting
    if as_iupac:
        return dna_consensus("/".join(ret))
    elif as_tuple:
        return tuple(ret)
    elif as_list:
        return ret
    else:
        return "/".join(ret)


def get_nuc_freqs(seqs: dict, ploidy: int) -> list:
    """
    Compute the nucleotide frequencies of a set of DNA sequences.

    Args:
        seqs (dict): A dictionary of DNA sequences. The keys are sample names
                     and the values are strings of nucleotides.
        ploidy (int): The ploidy of the sequences. If ploidy=1, ambiguities
                      will be skipped. If ploidy=2, ambiguities will be
                      resolved.

    Returns:
        freqs (list): A list of dictionaries, where each dictionary contains
                      the nucleotide frequencies for a single position in the
                      sequences.

    Raises:
        ValueError: If any sequence is empty or contains only invalid
                    characters, or if sequences are not all the same length.
    """
    if not all(seqs):
        raise ValueError("One or more sequences are empty.")

    sequence_length = len(seqs[list(seqs.keys())[0]])
    if not all(len(seq) == sequence_length for seq in seqs.values()):
        raise ValueError("Sequences are not all the same length.")

    freqs = []
    for loc in range(sequence_length):
        allnucs = ""
        for samp in seqs.keys():
            allnucs += dna_consensus(seqs[samp][loc]).lower()

        badchars = ["?", "-", "n"]
        if ploidy == 1:
            badchars += ["r", "y", "s", "w", "k", "m", "b", "d", "h", "v"]
        allnucs = ''.join([nuc for nuc in allnucs if nuc not in badchars])

        if ploidy == 2:
            iupacs = ["r", "y", "s", "w", "k", "m", "b", "d", "h", "v"]
            for ambig in iupacs:
                allnucs = allnucs.replace(
                    ambig, "".join(get_iupac_caseless(ambig))
                )
            for nuc in ["a", "c", "t", "g"]:
                allnucs = allnucs.replace(nuc, nuc + nuc)

        total = len(allnucs)
        counts = {"a": 0, "g": 0, "c": 0, "t": 0}
        for c in allnucs:
            if c in counts:
                counts[c] += 1
        if total <= 0:
            freqs.append({"a": 0, "g": 0, "c": 0, "t": 0})
        else:
            for nuc in counts.keys():
                counts[nuc] = float(counts[nuc] / total)
            freqs.append(counts)
    return freqs


def dna_consensus(seq: str) -> str:
    """
    Make a consensus DNA sequence from alleles separated by a "/" character.

    Args:
        seq (str): A string of DNA sequences separated by a "/" character.

    Returns:
        consens (str): A consensus DNA sequence derived from the input
                       sequences.
    """
    if not seq:
        return None
    if len(seq) <= 0:
        return None
    alleles = seq.split("/")
    consens = ""
    if len(alleles) < 1:
        return None
    elif not all(len(x) == len(alleles[0]) for x in alleles):
        raise ValueError(
            "Not all alleles are the same length: " + str(alleles)
        )
    elif len(alleles) == 1:
        return alleles[0]
    else:
        for i in range(len(alleles[0])):
            nucs = ""
            for a in alleles:
                nucs += a[i]
            temp = list_to_sort_unique_string(nucs.upper())
            consens += reverse_iupac_case(temp)
    return consens


def reverse_iupac_case(char: str) -> str:
    """
    Translate a string of DNA bases to an IUPAC ambiguity code, retaining case.

    Args:
        char (str): A string of DNA bases.

    Returns:
        iupac (str): An IUPAC ambiguity code.
    """
    iupac = {
        'A': 'A',
        'N': 'N',
        '-': '-',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'AG': 'R',
        'CT': 'Y',
        'AC': 'M',
        'GT': 'K',
        'AT': 'W',
        'CG': 'S',
        'CGT': 'B',
        'AGT': 'D',
        'ACT': 'H',
        'ACG': 'V',
        'ACGT': 'N',
        'a': 'a',
        'n': 'n',
        'c': 'c',
        'g': 'g',
        't': 't',
        'ag': 'r',
        'ct': 'y',
        'ac': 'm',
        'gt': 'k',
        'at': 'w',
        'cg': 's',
        'cgt': 'b',
        'agt': 'd',
        'act': 'h',
        'acg': 'v',
        'acgt': 'n'
    }
    return iupac[char]


def list_to_sort_unique_string(input: list) -> str:
    """
    Convert a list of characters to a sorted, unique string.

    Args:
        l (list): A list of characters.

    Returns:
        s (str): A string containing the characters from the input list, sorted
                 and with duplicates removed.
    """
    s = ''.join(sorted(set(input)))
    return s


def get_iupac_caseless(char: str) -> List[str]:
    """
    Split a character to IUPAC codes, assuming diploidy.

    Args:
        char (str): A character to be split.

    Returns:
        codes (list): A list of IUPAC codes corresponding to the input
                      character.
    """
    lower = False
    if char.islower():
        lower = True
        char = char.upper()
    iupac = {
        "A": ["A"],
        "G": ["G"],
        "C": ["C"],
        "T": ["T"],
        "N": ["A", "C", "G", "T"],
        "?": ["A", "C", "G", "T"],
        "-": ["A", "C", "G", "T", "-"],
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"]
    }
    codes = iupac[char]
    if lower:
        codes = [c.lower() for c in codes]
    return codes


def phase_snp(snp: str) -> str:
    """
    Phase a single nucleotide polymorphism (SNP) by splitting its IUPAC code
    into two alleles.

    Args:
        snp (str): The IUPAC code for the SNP.

    Returns:
        phase (str): The SNP alleles separated by a "/" character.
    """
    nucs = get_iupac_caseless(snp)
    if len(nucs) > 2 or len(nucs) < 1:
        return "n/n"
    elif len(nucs) == 1:
        phase = nucs[0] + '/' + nucs[0]
        return phase
    else:
        return "/".join(nucs)
