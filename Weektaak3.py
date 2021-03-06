# Weektaak 3 Course 3
# Groep 7

from translation import code
from translation import aa3
import matplotlib.pyplot as plt
from colour import Color


class Amino_acid:

    def __init__(self, amino_acid, codons, frequencies):
        self.__amino_acid = amino_acid
        self.__codons = codons
        self.__frequencies = frequencies
        self.codon_frequencies = self.generate_freq_dict()

    def get_amino_acid(self):
        return self.__amino_acid

    def get_codons(self):
        return self.__codons

    def get_frequencies(self):
        return self.__frequencies

    def get_codon_frequencies(self):
        return self.codon_frequencies

    def generate_freq_dict(self):

        _dict = {}
        for i in range(len(self.__codons)):
            _dict[self.__codons[i]] = self.__frequencies[i]
        return _dict

class amino_acid:

    def __init__(self, amino_acid, codons):
        self.__amino_acid = amino_acid
        self.__codons = codons




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
            seq += line.strip().lower()

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
        if "env" in header:
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


def fill_amino_acids(dict_freq_aa, codon_freq_dict):

    amino_acids = []
    for aa in dict_freq_aa.keys():
        aa_freq = dict_freq_aa[aa]

        codons = []
        frequencies = []

        for codon in codon_freq_dict.keys():
            if code[codon.lower()] == aa:
                codons.append(codon)
                frequencies.append(codon_freq_dict[codon]/aa_freq)

        obj = Amino_acid(aa, codons, frequencies)
        amino_acids.append(obj)

    return amino_acids


def amino_graph(aa_codon_frequencies, title):
    """"

    """
    fig, ax = plt.subplots(1, figsize=(12, 10))
    # Plot stacked bar

    red = Color("red")
    colors = list(red.range_to(Color("yellow"), 16))
    yellow = Color("yellow")
    colors += list(yellow.range_to(Color("green"), 16))
    green = Color("green")
    colors += list(green.range_to(Color("blue"), 16))
    blue = Color("blue")
    colors += list(blue.range_to(Color("purple"), 16))

    bar_colors = []
    for color in colors:
        bar_colors.append(color.rgb)

    aa_names = []
    for aa in aa_codon_frequencies:
        aa_names.append(aa.get_amino_acid())

    all_codons = code.keys()

    left = len(aa_names) * [0]
    for c in range(len(all_codons)):
        freqs = []
        for i in range(len(aa_codon_frequencies)):
            try:
                freq = aa_codon_frequencies[i].get_codon_frequencies()[
                    list(all_codons)[c]]
                freqs.append(freq)
                plt.text(left[i] + freq * 0.5, i, round(freq, 2))
            except KeyError:
                freqs.append(0)

        plt.barh(aa_names, freqs, left=left, color=bar_colors[c])

        for i in range(len(left)):
            left[i] += freqs[i]

    plt.suptitle(title)
    plt.xlabel("Codonfrequentie")
    plt.ylabel("Aminozuren")
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    # adjust limits and draw grid lines
    plt.ylim(-0.5, ax.get_yticks()[-1] + 0.5)
    ax.set_axisbelow(True)
    ax.xaxis.grid(color='gray', linestyle='dashed')
    plt.legend(all_codons, bbox_to_anchor=(1.04, 0.5), ncol=4, frameon=False)
    plt.tight_layout()
    plt.show()

def list_of_all_data(dict_fraction_codon):
    amino_codon_freqs = []

    for amino_acids, codons in dict_fraction_codon.items():
        if dict_fraction_codon[amino_acids] in aa3


def amino_graph(dict_fraction_codon, dict_freq_aa):
    amino_labels = []
    codon_labels = []
    codon_frequence = []

    for amino_acids, codons in aa3.items():
        amino_labels.append(amino_acids)
        codon_labels.append(codons)

    for codons, fractions in dict_fraction_codon.items():
        codon_frequence.append(fractions)






def main():
    list_virus = ["HIV 1 CDS.txt", "HIV 2 CDS.txt", "SIV CDS.txt",
                  "SIVmnd2 CDS.txt"]

    for i in range(len(list_virus)):
        file_name = list_virus[i]
        dict_prot = read_file(file_name)
        freq_codon_int, freq_codon_env = freq_codons_virus(dict_prot)
        dict_freq_aa_int = freq_amino_acid(freq_codon_int)
        dict_freq_aa_env = freq_amino_acid(freq_codon_env)
        aa_codon_frequencies_int = fill_amino_acids(dict_freq_aa_int,
                                                    freq_codon_int)
        aa_codon_frequencies_env = fill_amino_acids(dict_freq_aa_env,
                                                    freq_codon_env)
        title_int = file_name.replace(".txt", " int")
        amino_graph(aa_codon_frequencies_int, title_int)
        title_env = file_name.replace(".txt", " env")
        amino_graph(aa_codon_frequencies_env, title_env)

    list_organisms = ["Homo sapiens fumarate hydratase.fasta",
                      "Helianthus annuus fumarate hydratase.fasta",
                      "Pseudomonas aeruginosa fumarate hydratase.fasta",
                      "Staphylococcus aureus fumarate hydratase.fasta"]

    for i in range(len(list_organisms)):
        file_name = list_organisms[i]
        dict_prot = read_file(file_name)
        freq_codon = freq_codons_organism(dict_prot)
        dict_freq_aa = freq_amino_acid(freq_codon)
        aa_codon_frequencies = fill_amino_acids(dict_freq_aa, freq_codon)

        title = file_name.replace(".fasta", "")
        amino_graph(aa_codon_frequencies, title)


main()
