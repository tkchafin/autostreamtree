import pytest
import autostreamtree.sequence as seq_tools


# Fixture for common reference and alternate alleles
@pytest.fixture
def ref_and_alts():
    return "A", ["G", "C", "T"]


# Parametrized test for basic decoding
@pytest.mark.parametrize("gt, expected", [
    ((1, 2), "G/C"),
    ((0, 1), "A/G"),
    ((2, 2), "C/C")
])
def test_basic_decode(gt, expected, ref_and_alts):
    ref, alts = ref_and_alts
    assert seq_tools.decode(gt, ref, alts) == expected


# Test for handling missing genotypes
def test_missing_genotype(ref_and_alts):
    ref, alts = ref_and_alts
    assert seq_tools.decode((None, None), ref, alts) == "N/N"


# Parametrized test for IUPAC encoding
@pytest.mark.parametrize("gt, expected", [
    ((1, 2), seq_tools.dna_consensus("G/C")),
    ((0, 1), seq_tools.dna_consensus("A/G")),
    ((1, 1), seq_tools.dna_consensus("G/G")),
])
def test_iupac_encoding(gt, expected, ref_and_alts):
    ref, alts = ref_and_alts
    assert seq_tools.decode(gt, ref, alts, as_iupac=True) == expected


# Test for tuple output
def test_tuple_output(ref_and_alts):
    ref, alts = ref_and_alts
    gt = (1, 2)
    assert seq_tools.decode(gt, ref, alts, as_tuple=True) == ("G", "C")


# Test for list output
def test_list_output(ref_and_alts):
    ref, alts = ref_and_alts
    gt = (1, 2)
    assert seq_tools.decode(gt, ref, alts, as_list=True) == ["G", "C"]


# Parametrized test for get_nuc_freqs function
@pytest.mark.parametrize("seqs, ploidy, expected, raises_exception", [
    # Scenario 1: Basic Frequency Calculation with No Ambiguities
    (
        {"sample1": "ACGT", "sample2": "TGCA", "sample3": "GTAC"},
        1,
        [{'a': 0.3333333333333333, 'g': 0.3333333333333333, 'c': 0.0,
          't': 0.3333333333333333},
         {'a': 0.0, 'g': 0.3333333333333333, 'c': 0.3333333333333333,
          't': 0.3333333333333333},
         {'a': 0.3333333333333333, 'g': 0.3333333333333333,
          'c': 0.3333333333333333, 't': 0.0},
         {'a': 0.3333333333333333, 'g': 0.0, 'c': 0.3333333333333333,
          't': 0.3333333333333333}],
        False
    ),
    # Scenario 2: Frequency Calculation with Ambiguities (Ploidy 1)
    (
        {"sample1": "AGR", "sample2": "TCY", "sample3": "GTS"},
        1,
        [{'a': 0.3333333333333333, 'g': 0.3333333333333333, 'c': 0.0,
          't': 0.3333333333333333},
         {'a': 0.0, 'g': 0.3333333333333333, 'c': 0.3333333333333333,
          't': 0.3333333333333333},
         {'a': 0, 'g': 0, 'c': 0, 't': 0}],
        False
    ),
    # Scenario 3: Frequency Calculation with Ambiguities (Ploidy 2)
    (
        {"sample1": "AGR", "sample2": "TCY", "sample3": "GTS"},
        2,
        [{'a': 0.3333333333333333, 'g': 0.3333333333333333, 'c': 0.0,
          't': 0.3333333333333333},
         {'a': 0.0, 'g': 0.3333333333333333, 'c': 0.3333333333333333,
          't': 0.3333333333333333},
         {'a': 0.16666666666666666, 'g': 0.3333333333333333,
          'c': 0.3333333333333333, 't': 0.16666666666666666}],
        False
    ),
    # Edge Case: One or more sequences are entirely empty or invalid
    (
        {"sample1": "ACGT", "sample2": "", "sample3": "ACGT"},
        1,
        1,
        True
    ),

    # Edge Case: Sequences of mixed lengths
    (
        {"sample1": "ACG", "sample2": "T", "sample3": "GTAC"},
        1,
        1,
        True
    ),
])
def test_get_nuc_freqs(seqs, ploidy, expected, raises_exception):
    if raises_exception:
        with pytest.raises(ValueError):
            seq_tools.get_nuc_freqs(seqs, ploidy)
    else:
        assert seq_tools.get_nuc_freqs(seqs, ploidy) == expected


@pytest.mark.parametrize("test_input, expected, raises_exception", [
    ("A/C", "M", False),  # Simple heterozygous case
    ("A/A", "A", False),  # Homozygous case
    ("", None, False),    # Empty string
    ("A/C/G", "V", False),  # Multiple alleles
    ("A/AC", None, True),  # Alleles of different lengths
    ("A//T", None, True),  # Missing allele
])
def test_dna_consensus(test_input, expected, raises_exception):
    if raises_exception:
        with pytest.raises(ValueError):
            seq_tools.dna_consensus(test_input)
    else:
        assert seq_tools.dna_consensus(test_input) == expected


@pytest.mark.parametrize("char, expected", [
    ("AG", "R"),  # Basic case
    ("ct", "y"),  # Lowercase
    ("ACGT", "N"),  # All bases
])
def test_reverse_iupac_case(char, expected):
    assert seq_tools.reverse_iupac_case(char) == expected


@pytest.mark.parametrize("input, expected", [
    (["A", "C", "G", "T"], "ACGT"),  # Basic sorted string
    (["G", "G", "A", "C"], "ACG"),  # Duplicates and unsorted
    ([], ""),  # Empty list
])
def test_list_to_sort_unique_string(input, expected):
    assert seq_tools.list_to_sort_unique_string(input) == expected


@pytest.mark.parametrize("snp_input, expected, raises_exception", [
    ("R", "A/G", False),  # Valid IUPAC code for heterozygous SNP
    ("Y", "C/T", False),  # Another valid IUPAC code
    ("A", "A/A", False),  # Homozygous SNP
    ("N", "n/n", False)

])
def test_phase_snp(snp_input, expected, raises_exception):
    if raises_exception:
        with pytest.raises(ValueError):
            seq_tools.phase_snp(snp_input)
    else:
        assert seq_tools.phase_snp(snp_input) == expected
