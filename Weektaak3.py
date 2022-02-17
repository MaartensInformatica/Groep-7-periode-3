# Weektaak 3 Course 3
# Groep 7

from translation import code
import matplotlib.pyplot as plt
import numpy as np


def read_file(file):
    """" Read a fasta file and creates a dictionary with the headers as
    key and the sequences as value

    input:
    file - str - name of the file

    output:
    dict_prot - dict - dictionary with headers as key
                       and sequences as value
                       dict_prot[header]=seq
    """

    # Create empty dictionary and variables
    dict_prot = {}
    header = ""
    seq = ""

    open_file = open(file, "r")

    for line in open_file:
        # Check if line is the header
        if line.startswith(">"):
            # Set line as header if seq is empty
            if seq == "":
                header = line.strip()

            # Add header and seq to dict when new header is found
            else:
                dict_prot[header] = seq
                header = line.strip()
                seq = ""

        # Add line to sequence
        else:
            seq += line.strip()

    # No new header can be found at end of file
    # Add last header and seq to dict
    dict_prot[header] = seq

    open_file.close()

    return dict_prot


def freq_codons_organism(dict_prot):
    """" Goes through a dictionary with sequences and counts the
    frequency of used codons in the sequence

    input:
    dict_prot - dict - dictionary with headers as key
                       and sequences as value
                       dict_prot[header]=seq

    output:
    freq_codon - dict - dictionary with codons as key and their
                        frequency as value
    """

    # Create empty dictionary for frequency
    freq_codon = {}

    # Loop through the dictionary with sequences
    for header, seq in dict_prot.items():
        # Loop to go through the seq with steps of 3
        for i in range(0, len(seq), 3):
            # Create codon variable
            codon = seq[i:i+3]
            # Add 1 to frequency if codon found in dict
            if codon in freq_codon:
                freq_codon[codon] += 1
            # Add codon to dict if not found in dict
            else:
                freq_codon[codon] = 1

    return freq_codon


def freq_codons_virus(dict_prot):
    """"Goes through a dictionary with sequences and counts the
    frequency of used codons in the sequence. Separate dictionaries are
    created for the envelop protein and all the other proteins

    input:
    dict_prot - dict - dictionary with headers as key
                       and sequences as value
                       dict_prot[header]=seq

    output:
    freq_codon_int - dict - dictionary with codons as key and their
                            frequency as value
    freq_codon_env - dict - dictionary with codons as key and their
                            frequency as value
    """

    # Create empty dictionaries for frequency
    freq_codon_int = {}
    freq_codon_env = {}

    # Loop to go through the sequences
    for header, seq in dict_prot.items():
        # Check if the header is for envelop protein
        if "gene=env" in header:
            # Loop to go through the seq with steps of 3
            for i in range(0, len(seq), 3):
                # Create codon variable
                codon = seq[i:i+3]
                # Add 1 to frequency if codon found in dict
                if codon in freq_codon_env:
                    freq_codon_env[codon] += 1
                # Add codon to dict if not found in dict
                else:
                    freq_codon_env[codon] = 1

        # For non-envelop proteins
        else:
            # Loop to go through the seq with steps of 3
            for i in range(0, len(seq), 3):
                # Create codon variable
                codon = seq[i:i+3]
                # Add 1 to frequency if codon found in dict
                if codon in freq_codon_int:
                    freq_codon_int[codon] += 1
                # Add codon to dict if not found in dict
                else:
                    freq_codon_int[codon] = 1

    return freq_codon_int, freq_codon_env


def freq_amino_acid(freq_codon):
    """" Goes through a dictionary with codons and their frequency and
    counts the frequency of amino acids based on the codons

    input:
    freq_codon - dict - dictionary with codons as key and their
                        frequency as value

    output:
    dict_freq_aa - dict - dictionary with amino acid as key and their
                          frequency as value
    """

    # Create empty dictionary
    dict_freq_aa = {}

    # Loop to go through the dictionary
    for codon, frequency_codon in freq_codon.items():
        # Create the amino acid variable using the codon
        amino_acid = code[codon.lower()]
        # Add 1 to frequency if amino acid found in dict
        if amino_acid in dict_freq_aa:
            dict_freq_aa[amino_acid] += frequency_codon

        # Add amino acid to dict if not found in dict
        else:
            dict_freq_aa[amino_acid] = frequency_codon

    return dict_freq_aa


def fraction_count_codon(freq_codon, dict_freq_aa):
    """"Goes through dictionaries with sequences and and calculates the
    fraction of a codon used per amino acid

    input:
    freq_codon - dict - dictionary with codons as key and their
                        frequency as value
    dict_freq_aa - dict - dictionary with amino acid as key and their
                          frequency as value

    output:
    dict_fraction_codon - dict - dictionary with codon as key and their
                                 fraction as value
    """

    # Create empty dictionary
    dict_fraction_codon = {}

    # Loop to go through the dictionary
    for codon, frequency_codon in freq_codon.items():
        # Create the amino acid variable using the codon
        amino_acid = code[codon.lower()]

        # Calculate fraction of codon
        fraction = frequency_codon / dict_freq_aa[amino_acid]

        # Add fraction to new dictionary
        dict_fraction_codon[codon] = fraction

    return dict_fraction_codon


def main():
    CDS_HIV1 = "HIV 1 CDS.txt"
    CDS_HIV2 = "HIV 2 CDS.txt"
    CDS_SIV = "SIV CDS.txt"
    CDS_SIVmnd2 = "SIVmnd2 CDS.txt"
    CDS_H_sapiens = "Homo sapiens fumarate hydratase.fasta"
    CDS_H_annuus = "Helianthus annuus fumarate hydratase.fasta"
    CDS_P_aeruginosa = "Pseudomonas aeruginosa fumarate hydratase.fasta"
    CDS_S_aureus = "Staphylococcus aureus fumarate hydratase.fasta"

    dict_HIV1 = read_file(CDS_HIV1)
    freq_int_HIV1, freq_env_HIV1 = freq_codons_virus(dict_HIV1)
    freq_aa_int_HIV1 = freq_amino_acid(freq_int_HIV1)
    freq_aa_env_HIV1 = freq_amino_acid(freq_env_HIV1)
    fraction_codon_int_HIV1 = fraction_count_codon(freq_int_HIV1,
                                                   freq_aa_int_HIV1)
    fraction_codon_env_HIV1 = fraction_count_codon(freq_env_HIV1,
                                                   freq_aa_env_HIV1)

    dict_HIV2 = read_file(CDS_HIV2)
    freq_int_HIV2, freq_env_HIV2 = freq_codons_virus(dict_HIV2)
    freq_aa_int_HIV2 = freq_amino_acid(freq_int_HIV2)
    freq_aa_env_HIV2 = freq_amino_acid(freq_env_HIV2)
    fraction_codon_int_HIV2 = fraction_count_codon(freq_int_HIV2,
                                                   freq_aa_int_HIV2)
    fraction_codon_env_HIV2 = fraction_count_codon(freq_env_HIV1,
                                                   freq_aa_env_HIV2)

    dict_SIV = read_file(CDS_SIV)
    freq_int_SIV, freq_env_SIV = freq_codons_virus(dict_SIV)
    freq_aa_int_SIV = freq_amino_acid(freq_int_SIV)
    freq_aa_env_SIV = freq_amino_acid(freq_env_SIV)
    fraction_codon_int_SIV = fraction_count_codon(freq_int_SIV,
                                                  freq_aa_int_SIV)
    fraction_codon_env_SIV = fraction_count_codon(freq_env_SIV,
                                                  freq_aa_env_SIV)

    dict_SIVmnd2 = read_file(CDS_SIVmnd2)
    freq_int_SIVmnd2, freq_env_SIVmnd2 = freq_codons_virus(dict_SIVmnd2)
    freq_aa_int_SIVmnd2 = freq_amino_acid(freq_int_SIVmnd2)
    freq_aa_env_SIVmnd2 = freq_amino_acid(freq_env_SIVmnd2)
    fraction_codon_int_SIVmnd2 = fraction_count_codon(freq_int_SIVmnd2,
                                                      freq_aa_int_SIVmnd2)
    fraction_codon_env_SIVmnd2 = fraction_count_codon(freq_env_SIVmnd2,
                                                      freq_aa_env_SIVmnd2)

    print("--------")
    for k, v in freq_env_HIV1.items():
        print(k, v)

    print("--------")

    for k, v in freq_aa_env_HIV1.items():
        print(k, v)

    print("--------")

    for k, v in fraction_codon_env_HIV1.items():
        print(k, v)


main()
