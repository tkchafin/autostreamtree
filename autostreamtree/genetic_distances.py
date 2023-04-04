import sys
import os
import itertools
import math
import scipy
import numpy as np
import pandas as pd
from sortedcontainers import SortedDict

from typing import List, Tuple, Dict, Any, Union, Optional

import autostreamtree.aggregators as agg
import autostreamtree.sequence as seq

def get_pop_genmat(dist, indmat, popmap, dat, seqs, pop_agg="ARITH", loc_agg="ARITH", ploidy=2, global_het=False):
    #make matrix
    genmat = np.zeros((len(popmap),len(popmap)))
    #establish as nan
    genmat[:] = np.nan
    #for each combination, either average ind distances or calc freq dist
    for ia,ib in itertools.combinations(range(0,len(popmap)),2):
        if dist in ["JC69", "K2P", "PDIST", "TN84", "TN93"]:
            #print(popmap.keys()[ia], popmap.values()[ia])
            inds1 = [dat.index(x) for x in popmap.values()[ia]]
            inds2 = [dat.index(x) for x in popmap.values()[ib]]
            #print(inds1)
            genmat[ia,ib] = (agg.aggregate_dist(pop_agg, ([indmat[i, j] for i in inds1 for j in inds2])))
            genmat[ib,ia] = genmat[ia,ib]
        elif dist == "JOST":
            results=list()
            for loc in range(0, len(seqs[popmap.values()[ia][-1]])):
                seqs1 = get_alleles([seqs[x][loc] for x in popmap.values()[ia]])
                seqs2 = get_alleles([seqs[x][loc] for x in popmap.values()[ib]])
                if not clean_list(set(seqs1), ["n", "N", "-", "?"]) or not clean_list(set(seqs2), ["n", "N", "-", "?"]):
                    continue
                results.append(two_pop_jost_d(seqs1, seqs2, ploidy, global_het))
            if len(results) > 1:
                genmat[ia,ib] = genmat[ib,ia] = agg.aggregate_dist(loc_agg, results)
            elif len(results) < 1:
                #print("ERROR: population",popmap.values()[ia],"or",popmap.values()[ib],"lacks any data")
                raise ValueError
            else:
                genmat[ia,ib] = genmat[ib,ia] = results[0]
        elif dist == "GST" or dist == "GSTPRIME":
            HT=list()
            HS=list()
            for loc in range(0, len(seqs[popmap.values()[ia][-1]])):
                seqs1 = get_alleles([seqs[x][loc] for x in popmap.values()[ia]])
                seqs2 = get_alleles([seqs[x][loc] for x in popmap.values()[ib]])
                if not clean_list(set(seqs1), ["n", "N", "-", "?"]) or not clean_list(set(seqs2), ["n", "N", "-", "?"]):
                    continue
                if dist == "GST" or "GSTPRIME":
                    (ht, hs) = two_pop_ht_hs(seqs1, seqs2, ploidy, global_het)
                    HT.append(ht)
                    HS.append(hs)
            Ht_global = np.mean(HT)
            Hs_global = np.mean(HS)
            if dist == "GST":
                if Ht_global <= 0.0:
                    genmat[ia,ib] = genmat[ib,ia] = 0.0
                Gst = ((Ht_global - Hs_global) / Ht_global )
                GprimeST = ((Gst * (1.0 + Hs_global)) / (1.0 - Hs_global))
                genmat[ia,ib] = genmat[ib,ia] = GprimeST
            elif dist == "GSTPRIME":
                Ghedrick = ((2.0*(Ht_global - Hs_global)) / (((2.0*Ht_global) - Hs_global) * (1.0 - Hs_global)))
                genmat[ia,ib] = genmat[ib,ia] = Ghedrick
        elif dist == "FST" or dist == "LINFST":
            num = list() #numerator; a
            denom = list() #denominator; a*b*c
            for loc in range(0, len(seqs[popmap.values()[ia][-1]])):
                seqs1 = clean_inds([seqs[x][loc] for x in popmap.values()[ia]])
                seqs2 = clean_inds([seqs[x][loc] for x in popmap.values()[ib]])
                if not all("/" in x for x in seqs1) or not all("/" in x for x in seqs2):
                    print("ERROR: FST estimates require phased data.")
                    sys.exit(1)
                if len(seqs1) == 0 or len(seqs2) == 0:
                    #print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
                    continue
                (n, d) = two_pop_weir_cockerham_fst(seqs1, seqs2)
                num.append(n)
                denom.append(d)

            #if either population lacking data, set value to nan
            if len(num) <= 0 or len(denom) <= 0:
                #print("ERROR (two_pop_weir_cockerham_fst): No data for pops "+ia+" and "+ib+".")
                np.nan
            #if denominator is 0, set Fst to 0
            elif np.sum(denom) == 0.0:
                theta = 0.0
            #otherwise, calculate as normal
            else:
                theta = np.sum(num) / np.sum(denom)

            if dist == "FST":
                genmat[ia,ib] = genmat[ib,ia] = theta
            elif dist == "LINFST":
                if theta != 1.0:
                    genmat[ia,ib] = genmat[ib,ia] = (theta / (1-theta))
                else:
                    genmat[ia,ib] = genmat[ib,ia] = theta
        elif dist == "NEI83":
            results = list()
            loci = 0.0
            for loc in range(0, len(seqs[popmap.values()[ia][-1]])):
                seqs1 = get_alleles([seqs[x][loc] for x in popmap.values()[ia]])
                seqs2 = get_alleles([seqs[x][loc] for x in popmap.values()[ib]])
                if not clean_list(set(seqs1), ["n", "N", "-", "?"]) or not clean_list(set(seqs2), ["n", "N", "-", "?"]):
                    #print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
                    continue
                loci += 1.0
                results.append(two_pop_nei_da(seqs1, seqs2))
            Da = (1.0 - (np.sum(results) / loci))
            genmat[ia,ib] = genmat[ib,ia] = Da
        elif dist == "EUCLID":
            results = list()
            for loc in range(0, len(seqs[popmap.values()[ia][-1]])):
                seqs1 = get_alleles([seqs[x][loc] for x in popmap.values()[ia]])
                seqs2 = get_alleles([seqs[x][loc] for x in popmap.values()[ib]])
                if not clean_list(set(seqs1), ["n", "N", "-", "?"]) or not clean_list(set(seqs2), ["n", "N", "-", "?"]):
                    #print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
                    continue
                results.append(two_pop_euclid_dist(seqs1, seqs2))
            euclid = np.sum(results)
            genmat[ia,ib] = genmat[ib,ia] = euclid
        elif dist == "CHORD":
            num = list()
            denom = list()
            for loc in range(0, len(seqs[popmap.values()[ia][-1]])):
                seqs1 = get_alleles([seqs[x][loc] for x in popmap.values()[ia]])
                seqs2 = get_alleles([seqs[x][loc] for x in popmap.values()[ib]])
                if not clean_list(set(seqs1), ["n", "N", "-", "?"]) or not clean_list(set(seqs2), ["n", "N", "-", "?"]):
                    #print("WARNING: Skipping locus "+str(loc)+" in comparison of populations "+str(ia)+" and "+str(ib)+": Not enough data.")
                    continue
                (n, d) = two_pop_chord_dist(seqs1, seqs2)
                num.append(n)
                denom.append(d)
            Dch = np.sqrt((np.sum(num))/(np.sum(denom)))
            genmat[ia,ib] = genmat[ib,ia] = Dch
    np.fill_diagonal(genmat, 0.0)
    if 0.0 in genmat:
        print("WARNING: Coercing negative distances to 0.0")
        genmat[genmat<0.0] = 0.0
    return(genmat)

#function computes pairwise JC69-corrected genetic distances
def get_genmat(dist, points, seqs, ploidy, het, loc_agg):
    #make matrix
    genmat = np.zeros((len(points),len(points)))
    #establish as nan
    genmat[:] = np.nan

    #for models which relax equal nuc frequencies, get global frequencies for each locus
    #freqs will be a list of loci, with each locus as a dist of freqs
    if dist in ["TN84", "TN93"]:
        freqs = seq.get_nuc_freqs(seqs, ploidy)
        index = 1
        for f in freqs:
            print("Empirical base frequencies for locus",index, end=": [ ")
            for n in f:
                print(f'{n}={f[n]:.3f} ', end="")
            print("]")
            index = index + 1

    #for each combination, calc jukes-cantor corrected distance
    for ia, ib in itertools.combinations(range(0,len(points)),2):
        results=list()
        for loc in range(0, len(seqs[points.keys()[ia]])):
            seq1 = seqs[points.keys()[ia]][loc]
            seq2 = seqs[points.keys()[ib]][loc]
            if "/" in seq1:
                seq1 = seq.dna_consensus(seq1)
            if "/" in seq2:
                seq2 = seq.dna_consensus(seq2)
            if dist == "JC69":
                results.append(jukes_cantor_distance(seq1, seq2, het))
            elif dist == "K2P":
                results.append(k2p_distance(seq1, seq2, het))
            elif dist == "PDIST":
                if het:
                    results.append(p_distance(seq1, seq2))
                else:
                    results.append(hamming_distance(seq1, seq2))
            elif dist == "TN84":
                results.append(tn84_distance(seq1, seq2, freqs[loc], het))
            elif dist == "TN93":
                results.append(tn93_distance(seq1, seq2, freqs[loc], het))
        #aggregate results across loci
        genmat[ia,ib] = agg.aggregate_dist(loc_agg, results)
        genmat[ib,ia] = genmat[ia,ib]
    #fill diagonals
    np.fill_diagonal(genmat, 0.0)
    return(genmat)

def jukes_cantor_distance(seq1: str, seq2: str, het: bool = False) -> float:
    """
    Description: This function calculates the JC69-corrected p-distance between two DNA sequences.

    Args:
    - seq1 (str): The first DNA sequence.
    - seq2 (str): The second DNA sequence.
    - het (bool): Whether to use the Jukes-Cantor correction for heterozygous sites (True) or for all sites (False). Default is False.

    Returns:
    - dist (float): The JC69-corrected p-distance between the two sequences.
    """
    obs = 0.0
    # calculate the observed distance
    if het:
        obs = p_distance(seq1, seq2, trans=False)
    else:
        obs = hamming_distance(seq1, seq2, trans=False)

    # apply the JC69 correction
    if obs >= 0.75:
        obs = 0.74999
    dist = -0.75 * np.log(1.0 - ((4.0/3.0) * obs))

    # ensure distance is not negative and return
    if not dist > 0.0:
        return 0.0
    return dist

#function to return Kimura 2-parameter distances
def k2p_distance(seq1, seq2, het=False):
    """
    Computes the Kimura 2-parameter distance between two DNA sequences.

    Args:
    - seq1 (str): The first DNA sequence to compare.
    - seq2 (str): The second DNA sequence to compare.
    - het (bool): Whether to use the p-distance or Jukes-Cantor distance.
                  Defaults to False (Jukes-Cantor distance).

    Returns:
    - dist (float): The Kimura 2-parameter distance between the two sequences.
    """
    P=0.0
    Q=0.0
    if het:
        (P,Q)=p_distance(seq1, seq2, trans=True)
    else:
        (P,Q)=hamming_distance(seq1, seq2, trans=True)

    dist=-0.5*(np.log((1.0-(2.0*P)-Q) * math.sqrt(1.0-(2.0*Q))))
    #print(dist)
    if dist <= 0.0:
        return(0.0)
    return(dist)

def tn84_distance(seq1: str, seq2: str, freqs: Dict[str, float], het: bool = False) -> float:
    """
    Compute the TN84 distance between two DNA sequences.

    Args:
    seq1: str
        A DNA sequence.
    seq2: str
        Another DNA sequence.
    freqs: Dict[str, float]
        A dictionary containing the frequency of each nucleotide.
        For example, {'A': 0.3, 'C': 0.2, 'G': 0.2, 'T': 0.3}.
    het: bool, optional (default=False)
        Indicates whether heterozygosity is used.

    Returns:
    dist: float
        TN84 distance between seq1 and seq2.

    """
    D=0.0
    if het:
        D=p_distance(seq1, seq2, trans=False)
    else:
        D=hamming_distance(seq1, seq2, trans=False)

    ss=0.0
    for n in freqs:
        ss = ss + np.square(freqs[n])
    b=float(1.0-ss)
    dist=-b*np.log(1.0-((1.0/b)*D))
    #print(dist)
    if dist <= 0.0:
        return(0.0)
    return(dist)

def tn93_distance(seq1: str, seq2: str, freqs: Dict[str, float], het: bool = False) -> float:
    """
    Compute TN93 distances between two sequences.

    Args:
        seq1: A DNA sequence.
        seq2: A DNA sequence.
        freqs: A dictionary of nucleotide frequencies, where the keys are nucleotides (A, T, C, G)
            and the values are the frequency of the corresponding nucleotide in the alignment.
        het: Whether the nucleotide differences are considered heterozygous. Default is False.

    Returns:
        The TN93 distance between the two input sequences.

    """

    P1 = 0.0
    P2 = 0.0
    Q = 0.0
    if het:
        P1, P2, Q = p_distance(seq1, seq2, trans=False, transSplit=True)
    else:
        P1, P2, Q = hamming_distance(seq1, seq2, trans=False, transSplit=True)

    # Calculate nucleotide frequencies
    gR = freqs['g'] + freqs['a']
    gY = freqs['c'] + freqs['t']

    # Calculate k values
    k1 = (2.0 * (freqs['a'] * freqs['g'])) / gR
    k2 = (2.0 * (freqs['t'] * freqs['c'])) / gY
    k3 = 2.0 * ((gR * gY) - ((freqs['a'] * freqs['g'] * gY) / gR) - ((freqs['t'] * freqs['c'] * gR) / gY))

    # Calculate weight values
    w1 = 1.0 - (P1 / k1) - (Q / (2.0 * gR))
    w2 = 1.0 - (P2 / k3) - (Q / (2.0 * gY))
    w3 = 1.0 - (Q / (2.0 * gR * gY))

    # Calculate the distance
    dist = -(k1 * np.log(w1)) - (k2 * np.log(w2)) - (k2 * np.log(w3))

    # Check for negative or NaN distances
    if dist <= 0.0 or np.isnan(dist):
        return 0.0

    return dist



#p distance = D / L (differences / length)
#L is the UNGAPPED distance #TODO: Maybe make this optional later
#ambigs are expanded
#when trans = true, returns two values: P (transitions/L) and Q (transversions/L)
def p_distance(seq1, seq2, trans=False, transSplit=False):
    L = min(len(seq1), len(seq2))
    D=0.0
    P=0.0
    Q=0.0
    P1=0.0
    P2=0.0
    for n1, n2 in zip(seq1.lower(), seq2.lower()):
        if n1 in ["?", "-", "n"] or n2 in ["?", "-", "n"] :
            L = L-1
            continue
        elif n1 == n2:
            continue
        else: #if n1 and n2 not equal and not gaps
            if n1 in ["a", "c", "g", "t"] and n2 in ["a", "c", "g", "t"]:
                if trans or transSplit:
                    if n1 in ["a","g"] and n2 in ["a", "g"]:
                        P = P+1
                        P1 = P1+1
                    elif n1 in ["c","t"] and n2 in ["c","t"]:
                        P=P+1
                        P2=P2+1
                    else:
                        Q=Q+1
                else:
                    D = D + 1.0
                continue
            else:
                ex1 = seq.get_iupac_caseless(n1)
                ex2 = seq.get_iupac_caseless(n2)
                val=1.0
                if len(ex1) > 1 and len(ex2) > 1:
                    val=0.25
                elif len(ex1) == 1 and len(ex2)== 1:
                    val=1.0
                else:
                    val=0.50
                for nuc1, nuc2 in zip(ex1, ex2):
                    if nuc2 == nuc1:
                        continue
                    else:
                        if trans or transSplit:
                            if n1 in ["a","g"] and n2 in ["a", "g"]:
                                P = P+val
                                P1 = P1+val
                            elif n1 in ["c","t"] and n2 in ["c","t"]:
                                P = P+val
                                P2 = P2+val
                            else:
                                Q=Q+val
                        else:
                            D=D+val
    if trans:
        transitions=0.0
        transversions=0.0
        if P > 0.0:
            transitions = float(P/L)
        if Q > 0.0:
            transversions = float(Q/L)
        return(transitions,transversions)
    elif transSplit:
        transition1=0.0
        transition2=0.0
        transversions=0.0
        if P1 > 0.0:
            transition1 = float(P1/L)
        if P2 > 0.0:
            transition2 = float(P2/L)
        if Q > 0.0:
            transversions = float(Q/L)
        return(transition1, transition2, transversions)
    else:
        if D <= 0.0:
            return(0.0)
        else:
            return(float(D/L))

#p distance = D / L (differences / length)
#gaps ignored
#ambigs are treated as alleles
def hamming_distance(seq1, seq2, trans=False, transSplit=False):
    L = min(len(seq1), len(seq2))
    D=0.0
    P=0.0
    Q=0.0
    P1=0.0
    P2=0.0
    for n1, n2 in zip(seq1.lower(), seq2.lower()):
        if n1 in ["?", "-", "n"] or n2 in ["?", "-", "n"]:
            L = L-1
            continue
        elif n1 == n2:
            continue
        else:
            if n1 != n2:
                if trans or transSplit:
                    ex1 = seq.get_iupac_caseless(n1)
                    ex2 = seq.get_iupac_caseless(n2)
                    for nuc1, nuc2 in zip(ex1, ex2):
                        if nuc2 == nuc1:
                            continue
                        else:
                            if n1 in ["a","g"] and n2 in ["a", "g"]:
                                P = P+1.0
                                P1 = P1 + 1.0
                            elif n1 in ["c","t"] and n2 in ["c","t"]:
                                P=P+1.0
                                P2 = P2 + 1.0
                            else:
                                Q=Q+1.0
                else:
                    D = D + 1.0
                continue
    if trans:
        transitions=0.0
        transversions=0.0
        if P > 0.0:
            transitions = float(P/L)
        if Q > 0.0:
            transversions = float(Q/L)
        return(transitions,transversions)
    elif transSplit:
        transition1=0.0
        transition2=0.0
        transversions=0.0
        if P1 > 0.0:
            transition1 = float(P1/L)
        if P2 > 0.0:
            transition2 = float(P2/L)
        if Q > 0.0:
            transversions = float(Q/L)
        return(transition1, transition2, transversions)
    else:
        if D <= 0.0:
            return(0.0)
        else:
            return(float(D/L))

def two_pop_nei_da(s1: List[str], s2: List[str]) -> float:
    """
    Computes Nei's 1983 Da estimator for two populations.

    Args:
        s1 (List[str]): A list of phased genotypes from population 1, e.g. ['A/A', 'A/B', 'B/B', ...].
        s2 (List[str]): A list of phased genotypes from population 2, e.g. ['A/A', 'A/C', 'C/C', ...].

    Returns:
        float: The Nei's 1983 Da estimator.

    """
    # Clean the input lists by removing individuals with unknown or gap alleles.
    s1 = clean_list(s1, ["n", "?", "-", "N"])
    s2 = clean_list(s2, ["n", "?", "-", "N"])
    
    # Get the list of unique alleles from both populations.
    uniques = uniq_alleles(s1+s2)
    
    # Compute the sum of squared roots of frequencies for each allele.
    sumSqRt = 0.0
    for allele in uniques:
        if allele in ["-", "?", "n", "N"]:
            continue
        else:
            Xu = float(s1.count(allele) / len(s1))
            Yu = float(s2.count(allele) / len(s2))
            sumSqRt += np.sqrt(Xu*Yu)
    return sumSqRt


def two_pop_euclid_dist(s1: List[str], s2: List[str]) -> float:
    """
    Computes Euclidean distance between two populations for a single locus.

    Args:
    - s1 (List[str]): List of allele calls for population 1 at a given locus
    - s2 (List[str]): List of allele calls for population 2 at a given locus

    Returns:
    - Euclidean distance between the two populations for the given locus

    """
    # Get unique alleles in the two populations
    uniques = seq.uniq_alleles(s1+s2)
    
    # Clean the sequences by removing any unknown or gap alleles
    s1 = clean_list(s1, ["n", "?", "-", "N"])
    s2 = clean_list(s2, ["n", "?", "-", "N"])
    
    # Calculate Euclidean distance
    sumSq = 0.0
    for allele in uniques:
        if allele in ["-", "?", "n", "N"]:
            continue
        else:
            Xu = float(s1.count(allele) / len(s1))
            Yu = float(s2.count(allele) / len(s2))
            sumSq += np.square(Xu - Yu)
    return sumSq


#Cavalli-Sforza and Edwards 1967 chord distance
#non-nucleotide alleles are deleted
def two_pop_chord_dist(s1, s2):
    s1 = clean_list(s1, ["n", "?", "-", "N"])
    s2 = clean_list(s2, ["n", "?", "-", "N"])
    uniques = uniq_alleles(s1+s2)
    sumSqRt = 0.0
    for allele in uniques:
        if allele in ["-", "?", "n", "N"]:
            continue
        else:
            Xu = float(s1.count(allele) / len(s1))
            Yu = float(s2.count(allele) / len(s2))
            sumSqRt += np.sqrt(np.square(Xu)*np.square(Yu))
    return((1.0-sumSqRt), (len(uniques)-1.0))

def two_pop_weir_cockerham_fst(s1, s2):
    """
    Computes Weir and Cockerham's THETAst Fst approximation for two populations.

    Args:
        s1 (list): A list of phased genotypes for population 1.
        s2 (list): A list of phased genotypes for population 2.

    Returns:
        tuple: A tuple containing two floats. The first float is the numerator
        of the THETAst estimator for the locus. The second float is the denominator
        of the THETAst estimator for the locus.

    Raises:
        ValueError: If the inputs are not valid.

    """
    # Check inputs
    if not isinstance(s1, list) or not isinstance(s2, list):
        raise ValueError("Inputs must be lists.")
    if not s1 or not s2:
        raise ValueError("Inputs must not be empty.")
    if not all(isinstance(x, str) for x in s1) or not all(isinstance(x, str) for x in s2):
        raise ValueError("Inputs must only contain strings.")
    if not all("/" in x for x in s1) or not all("/" in x for x in s2):
        raise ValueError("Inputs must be phased genotypes.")
    
    # Initialize variables
    num = 0.0
    denom = 0.0

    #mean sample size
    alleles1 = get_alleles(s1) #split alleles s1
    alleles2 = get_alleles(s2) #split alleles s2
    uniques = uniq_alleles(s1+s2) #list of unique alleles only
    r = 2.0 #number of pops
    n1 = float(len(s1)) #pop size of pop 1
    n2 = float(len(s2)) #pop size of pop 2
    csd = np.std([n1, n2])
    cm = np.mean([n1, n2])
    nbar = cm
    csquare = (csd*csd) / (cm*cm)
    nC   = nbar * (1.0 - (csquare/r)) #coeff of pop size variance
    for allele in uniques:
        ac1 = float(alleles1.count(allele))
        ac2 = float(alleles2.count(allele))
        p1 = ac1 / float(len(alleles1))
        p2 = ac2 / float(len(alleles2))
        h1 = get_het_from_phased(allele, s1, count=True)
        h2 = get_het_from_phased(allele, s2, count=True)
        pbar = (ac1+ac2) / (float(len(alleles1)) + float(len(alleles2)))
        ssquare = ((np.sum( [ (n1* (np.square(p1 - pbar)) ), (n2* (np.square(p2 - pbar))) ])) / ((r-1.0)*nbar))
        hbar = ((h1+h2) / (r * nbar))
        #print(nbar)
        #print(hbar)
        if nbar != 1.0:
            a = ((nbar/nC) *
                (ssquare -
                ((1.0 / (nbar-1.0)) *
                ((pbar * (1.0-pbar)) -
                ((r - 1.0) * ssquare / r) -
                (hbar / 4.0)))))
            b = ((nbar / (nbar-1.0)) *
                ((pbar * (1.0 - pbar)) -
                ((r - 1.0) * ssquare / r) -
                (((2.0 * nbar) - 1.0) * hbar / (4.0 * nbar))))
            c = hbar/2.0
            d = a+b+c
            num += a
            denom += d
    return(num, denom)

def clean_inds(inds):
    """
    Removes individuals with unknown or gap alleles.

    Args:
        inds (list): A list of individuals.

    Returns:
        list: A list of individuals without unknown or gap alleles.
    """
    ret = []
    for ind in inds:
        if "-" not in ind and "?" not in ind and "n" not in ind and "N" not in ind:
            ret.append(ind)
    return ret


def two_pop_jost_d(seqs1, seqs2, ploidy, global_het=False):
    """
    Computes Jost's D using Nei and Chesser's Hs and Ht estimators.

    Args:
    - seqs1: list of sequences for population 1
    - seqs2: list of sequences for population 2
    - ploidy: the ploidy of the sequences
    - global_het (optional): whether to compute global heterozygosity instead of population-specific heterozygosity
    Returns:
    - The Jost's D estimate for the two populations.
    """
    # Compute Hs and Ht estimates
    if global_het:
        Ht = get_global_het(seqs1 + seqs2)
    else:
        Ht = get_average_het(seqs1, seqs2)
    Hs = np.mean([get_global_het(seqs1), get_global_het(seqs2)])
    harmN = scipy.stats.hmean([(len(seqs1) / ploidy), (len(seqs1) / ploidy)])
    Hs_est = Hs * ((2.0 * harmN) / ((2.0 * harmN) - 1.0))
    Ht_est = Ht + (Hs / harmN * 2.0 * 2.0)
    
    # Compute Jost's D
    if Ht_est == 0.0:
        return(0.0)
    D = (Ht_est - Hs_est) * 2.0 #b/c in pw estimate, N/N-1 is always 2
    if D == 2:
        return(1.0)
    elif D <= 0.0:
        return(0.0)
    return(D)


def two_pop_ht_hs(seqs1, seqs2, ploidy, global_het=False):
    """
    Computes Nei's Fst estimator (Gst) using Nei and Chessers Hs and Ht estimators
    also applies Hedrick's (2005) sample size correction, thus returning G'st.

    Args:
    - seqs1, seqs2 (list): lists of sequences for populations 1 and 2.
    - ploidy (int): ploidy of the sequences.
    - global_het (bool, optional): If True, computes global heterozygosity. Defaults to False.

    Returns:
    - A tuple containing Ht_est and Hs_est.

    """
    if global_het:
        Ht = get_global_het(seqs1 + seqs2)
    else:
        Ht = get_average_het(seqs1, seqs2)

    Hs = np.mean([get_global_het(seqs1), get_global_het(seqs2)])
    harmN = scipy.stats.hmean([(len(seqs1) / ploidy), (len(seqs1) / ploidy)])

    # Hs correction based on Hedrick (2005)
    Hs_est = Hs * ((2.0 * harmN) / ((2.0 * harmN) - 1.0))

    # Ht correction based on Hedrick (2005)
    Ht_est = Ht + (Hs / (harmN * 2.0 * 2.0))

    # Gst = ((Ht_est - Hs_est) / Ht_est )

    # GprimeST = ((Gst * (1.0 + Hs_est)) / (1.0 - Hs_est))
    return (Ht_est, Hs_est)

def get_het_from_phased(allele, phasedList, count=False):
    """
    Returns observed heterozygosity of an allele given a list of phased genotypes (e.g. allele1/allele2 for each individual). Assumes diploid.

    Args:
    allele (str): The allele for which to compute heterozygosity.
    phasedList (list): A list of phased genotypes, where each element is a string of the form 'allele1/allele2'.
    count (bool): Whether to return the number of heterozygotes instead of the proportion.

    Returns:
    float: The proportion of heterozygotes for the given allele, unless `count` is True, in which case the count of heterozygotes is returned.
    """
    hets = 0.0
    twoN = (len(phasedList)) * 2.0
    for genotype in phasedList:
        if "/" not in genotype:
            print("ERROR (get_het_from_phased): Phased genotypes are required.")
        gens = genotype.split("/")
        if gens[0] == allele and gens[1] != allele:
            hets += 1.0
            continue
        elif gens[1] == allele and gens[0] != allele:
            hets += 1.0
            continue
        else:
            continue
    if count:
        return hets
    else:
        return hets / twoN

def get_global_het(seqs):
    """
    Computes global expected heterozygosity (Ht) from a set of sequences.

    Args:
    seqs (list): A list of sequences.

    Returns:
    float: The Ht estimate for the given sequences.
    """
    hom = 0.0
    uniq_alleles = clean_list(set(seqs), ["n", "?", "-", "N"])
    freqs = [np.square(float(seqs.count(x) / len(seqs))) for x in uniq_alleles]
    hom = np.sum(freqs)
    return 1.0 - hom

def get_average_het(s1, s2):
    """Computes the mean expected heterozygosity (Ht) from two populations.

    Args:
        s1 (list): A list of alleles from population 1.
        s2 (list): A list of alleles from population 2.

    Returns:
        float: The Ht estimate per locus.
    """
    hom = 0.0
    # Get unique alleles
    uniq_alleles = clean_list(set(s1+s2), ["n", "?", "-", "N"])
    # Compute frequency for each allele
    freqs = [np.square(np.mean([float(s1.count(x)/len(s1)), float(s2.count(x)/len(s2))])) for x in uniq_alleles]
    hom = np.sum(freqs)
    return(1.0-hom)

def clean_list(l, bads):
    """Removes bad items from a list.

    Args:
        l (list): A list to clean.
        bads (list): A list of items to remove.

    Returns:
        list: A cleaned list.
    """
    if not any(item not in bads for item in set(l)):
        return(False)
    for b in bads:
        if b in l:
            l.remove(b)
    return(l)

def get_alleles(s):
    """Splits a string of alleles separated by "/" into a list.

    Args:
        s (str): A string of alleles separated by "/".

    Returns:
        list: A list of alleles.
    """
    return(sum([x.split("/") for x in s], []))

def uniq_alleles(s):
    """Returns the unique alleles in a list of alleles.

    Args:
        s (list): A list of alleles.

    Returns:
        set: A set of unique alleles.
    """
    return(set(sum([x.split("/") for x in s], [])))

