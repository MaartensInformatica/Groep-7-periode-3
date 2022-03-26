# Weektaak 4 Course 3
# Groep 7

import matplotlib.pyplot as plt
import numpy as np
from AAclasses import aa_classes
from AAclasses import aa_colours
from AAclasses import aa_abr


class AminoAcid:

    def __init__(self, amino_acid, abbreviation, frequency, percentage,
                 colour):
        self.__amino_acid = amino_acid
        self.__abbreviation = abbreviation
        self.__frequency = frequency
        self.__percentage = percentage
        self.__colour = colour

    def get_amino_acid(self):
        return self.__amino_acid

    def get_abbreviation(self):
        return self.__abbreviation

    def get_frequency(self):
        return self.__frequency

    def get_percentage(self):
        return self.__percentage

    def get_colour(self):
        return self.__colour

    def get_hydro_interaction(self):
        if "hydrophobic" in aa_classes[self.__amino_acid]:
            return "hydrophobic"
        if "hydrophylic" in aa_classes[self.__amino_acid]:
            return "hydrophylic"
        if "inbetween" in aa_classes[self.__amino_acid]:
            return "inbetween"


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
    dict_prot_seqs = {}
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
                dict_prot_seqs[header] = seq
                header = line.strip()
                seq = ""

        # Add line to sequence
        else:
            seq += line.strip()

    # No new header can be found at end of file
    # Add last header and seq to dict
    dict_prot_seqs[header] = seq

    open_file.close()

    return dict_prot_seqs


def separate_env_int_virus(dict_prot_seqs):
    """"Goes through a dictionary with sequence(s) and separates
    sequences for envelop protein(s) from other sequences

    input:
    dict_prot - dict - dictionary with header as key
                       and sequence as value
                       dict_prot[header]=seq

    output:
    dict_prot_seq_int - dict - dictionary with header as key and their
                               sequence as value
    dict_prot_seq_env - dict - dictionary with header as key and their
                               sequence as value
    """

    # Create empty dictionaries for frequency
    dict_prot_seq_int = {}
    dict_prot_seq_env = {}

    # Loop to go through the sequences
    for header, seq in dict_prot_seqs.items():
        # Check if the header is for envelop protein
        if "env" in header:
            # Add envelop protein seq to envelop dictionary
            dict_prot_seq_env[header] = seq

        else:
            # Add other proteins to intern dictionary
            dict_prot_seq_int[header] = seq

    return dict_prot_seq_int, dict_prot_seq_env


def freq_amino_acid(dict_prot_seqs):
    """" Goes through a dictionary with sequences and
    counts the frequency of amino acids

    input:
    dict_prot_seqs - dict - dictionary with header as key and the
                            sequence as value

    output:
    dict_freq_aa - dict - dictionary with amino acid as key and their
                          frequency as value
    """

    # Create empty dictionary
    dict_freq_aa = {}

    # Empty variable for total amino acids in sequence(s)
    total_aa = 0

    # Loop to go through the dictionary
    for sequences in dict_prot_seqs.values():
        for amino_acid in sequences:
            total_aa += 1

            # Add 1 to frequency if amino acid found in dict
            if amino_acid in dict_freq_aa:
                dict_freq_aa[amino_acid] += 1

            # Add amino acid to dict if not found in dict
            else:
                dict_freq_aa[amino_acid] = 1

    return dict_freq_aa, total_aa


def make_aa_class(dict_freq_aa, total_aa):
    """"Create class objectives from amino acids and add the objects to
    a list

    input:
    dict_freq_aa - dict - dictionary with amino acid as key and their
                          frequency as value
    total_aa - int - the total amount of amino acids that where in the
                     sequence(s)

    output:
    class_amino_acids - list - list with the amino acids saved as a class
    """

    # Empty list for the amino acids
    class_amino_acids = []

    # Create object for each amino acid in AminoAcid class
    for amino_acid in dict_freq_aa.keys():
        # Set the variables
        abbreviation = aa_abr[amino_acid]
        frequencies = dict_freq_aa[amino_acid]
        percentage = frequencies / total_aa * 100
        colour = aa_colours[amino_acid]

        # Set object, add object to list
        obj = AminoAcid(amino_acid, abbreviation, frequencies, percentage,
                        colour)
        class_amino_acids.append(obj)

    return class_amino_acids


def calculate_percentages_hydro(class_amino_acids, total_aa):
    """"Calculate the percentages of hydrophobic, hydrophylic and
    inbetween amino acids and add the results to a dictionary

    input:
    class_amino_acids - list - list with the amino acids saved as a class
    total_aa - int - the total amount of amino acids that where in the
                     sequence(s)

    output:
    dict_hydro_percentage - dict - dictionary with either hydrophobic,
                                   hydrophylic or inbetween as key and
                                   float as value, rounded with 2 decimals
    """

    # Create empty dictionary
    dict_hydro_percentage = {}

    # Set start value at zero
    total_hydrophobic = 0
    total_hydrophylic = 0
    total_inbetween = 0

    # Loop to go through the list with amino acid objects
    for aa in class_amino_acids:
        # Check what what the interaction is with water
        # Add the frequency to the correct total variable
        if "hydrophobic" in aa.get_hydro_interaction():
            total_hydrophobic += aa.get_frequency()
        if "hydrophylic" in aa.get_hydro_interaction():
            total_hydrophylic += aa.get_frequency()
        if "inbetween" in aa.get_hydro_interaction():
            total_inbetween += aa.get_frequency()

    # Calculate the percentages of the water interactions
    # Add the percentages to the dictionary
    per_hydrophobic = total_hydrophobic / total_aa * 100
    dict_hydro_percentage["hydrophobic"] = per_hydrophobic

    per_hydrophylic = total_hydrophylic / total_aa * 100
    dict_hydro_percentage["hydrophylic"] = per_hydrophylic

    per_inbetween = total_inbetween / total_aa * 100
    dict_hydro_percentage["inbetween"] = per_inbetween

    return dict_hydro_percentage


def highest_lowest_aa_freq(dict_freq_aa, total_aa):
    """"Determine which three amino acids are the most and the least frequent
    in a protein and create a dictionary for both the highest and lowest
    frequencies

    input:
    dict_freq_aa - dict - dictionary with amino acid as key and their
                          frequency as value
    total_aa - int - the total amount of amino acids that where in the
                     sequence(s)

    output:
    dict_highest_three - dict - dictionary with the three most frequent amino
                                acids, the amino acids are the keys and their
                                frequency the values
    dict_lowest_three - dict -  dictionary with the three least frequent amino
                                acids, the amino acids are the keys and their
                                frequency the values
    """

    # Create empty dictionaries
    dict_highest_three = {}
    dict_lowest_three = {}

    # Create list with the three highest values
    highest_aa = sorted(dict_freq_aa, key=dict_freq_aa.get, reverse=True)[:3]

    # Create list with the three lowest values
    lowest_aa = sorted(dict_freq_aa, key=dict_freq_aa.get, reverse=False)[:3]

    # Add the keys and values to the dictionaries from the lists
    for i in highest_aa:
        dict_highest_three[i] = round(dict_freq_aa[i] / total_aa * 100, 2)

    for i in lowest_aa:
        dict_lowest_three[i] = round(dict_freq_aa[i] / total_aa * 100, 2)

    return dict_highest_three, dict_lowest_three


def amino_graph_cys_trp(class_amino_acids, title):
    """"Create pie chart which will show how much cysteine and tyrosine is present
    in the protein(s) in percentages

    input:
    class_amino_acids - list - list with the amino acids saved as a class
    title - str - the title that will be added to the graph
    """

    # Define the variables
    c_freq = 0
    w_freq = 0
    others_freq = 0

    c_label = ""
    w_label = ""
    others_label = ""

    c_colour = ""
    w_colour = ""
    other_colour = ""

    # Define the values for the variables
    for aa in class_amino_acids:
        if aa.get_amino_acid() == "C":
            c_freq += aa.get_percentage()
            c_label = aa.get_abbreviation()
            c_colour = aa.get_colour()

        if aa.get_amino_acid() == "W":
            w_freq += aa.get_percentage()
            w_label = aa.get_abbreviation()
            w_colour = aa.get_colour()

        else:
            others_freq += aa.get_percentage()
            others_label = "Others"
            other_colour = "silver"

    # Add the values to lists
    values = [c_freq, w_freq, others_freq]
    labels = [c_label, w_label, others_label]
    colours = [c_colour, w_colour, other_colour]

    # Create the pie chart
    fig = plt.figure()
    ax11 = fig.add_subplot(111)

    w, l, p = ax11.pie(values, labels=labels, colors=colours,
                       autopct='%1.1f%%', startangle=140, pctdistance=1,
                       radius=1)

    # Change the label positions
    l[1]._y = l[1]._y - 0.1

    # Set position for the values in chart
    pctdists = [.8, .5, .4]

    for t, d in zip(p, pctdists):
        xi, yi = t.get_position()
        ri = np.sqrt(xi ** 2 + yi ** 2)
        phi = np.arctan2(yi, xi)
        x = d * ri * np.cos(phi)
        y = d * ri * np.sin(phi)
        t.set_position((x, y))

    plt.title(title)
    plt.tight_layout()
    plt.show()


def amino_graph_hydro(dict_hydro_percentage, title):
    """"Create a pie chart with the amount of hydrophobic, hydrophylic and
    inbetween amino acids shown in percentages

    input:
    dict_hydro_percentage - dict - dictionary with either hydrophobic,
                                   hydrophylic or inbetween as key and
                                   float as value, rounded with 2 decimals
    title - str - the title that will be added to the graph
    """

    # Define the 'empty' variables
    hydrophobic_perc = 0
    hydrophylic_perc = 0
    inbetween_perc = 0

    # Add values to the variables
    for k, v in dict_hydro_percentage.items():
        if k == "hydrophobic":
            hydrophobic_perc += v
        if k == "hydrophylic":
            hydrophylic_perc += v
        if k == "inbetween":
            inbetween_perc += v

    # Create values and labels list
    values = [hydrophobic_perc, hydrophylic_perc, inbetween_perc]
    labels = ["Hydrophobic", "Hydrophylic", "Inbetween"]
    colours = ["darkorange", "darkmagenta", "g"]

    # Create pie chart
    fig1, ax1 = plt.subplots()
    ax1.pie(values, labels=labels, autopct='%1.1f%%', colors=colours,
            pctdistance=0.6, startangle=90, radius=1)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def amino_graph_highest_lowest(dict_lowest_three, dict_highest_three, title):
    """"Create a bar chart which shows the three most and least frequent
    amino acids in percentages

    input:
    dict_highest_three - dict - dictionary with the three most frequent amino
                                acids, the amino acids are the keys and their
                                frequency the values
    dict_lowest_three - dict -  dictionary with the three least frequent amino
                                acids, the amino acids are the keys and their
                                frequency the values
    title - str - the title that will be added to the graph
    """

    # Create empty lists
    x_labels = []
    y_values = []
    colours = []

    # Add values to the lists for the least frequent amino acids
    for amino_acid, percentage in dict_lowest_three.items():
        x_labels.append(aa_abr[amino_acid])
        y_values.append(percentage)
        colours.append(aa_colours[amino_acid])

    # Add values to the lists for the most frequent amino acids
    for amino_acid, percentage in dict_highest_three.items():
        x_labels.append(aa_abr[amino_acid])
        y_values.append(percentage)
        colours.append(aa_colours[amino_acid])

    # Create the bar chart
    plt.bar(x_labels, y_values, color=colours)

    # Add values in the graph
    for i in range(len(x_labels)):
        if y_values[i] < max(y_values) / 12:
            plt.text(i - 0.25, y_values[i] + 0.1, f"{y_values[i]}%")
        elif y_values[i] >= 10.0:
            plt.text(i - 0.35, y_values[i] - max(y_values) / 20, f"{y_values[i]}%")
        else:
            plt.text(i - 0.25, y_values[i] - max(y_values) / 20, f"{y_values[i]}%")

    plt.ylabel('Percentage (%)')
    plt.xlabel('Amino acid')
    plt.title(title)

    plt.show()


def main():
    list_virus = ["HIV 1 Protein.txt", "HIV 2 Protein.txt", "SIV Protein.txt",
                  "SIVmnd2 Protein.txt"]

    for i in range(len(list_virus)):
        file_name_virus = list_virus[i]
        dict_prot_seqs = read_file(file_name_virus)
        dict_prot_seq_int, dict_prot_seq_env = \
            separate_env_int_virus(dict_prot_seqs)

        # Call everything for internal viral protein
        dict_freq_aa_int, total_aa_int = freq_amino_acid(dict_prot_seq_int)
        class_amino_acids_int = make_aa_class(dict_freq_aa_int, total_aa_int)
        dict_hydro_percentage_int = calculate_percentages_hydro(
            class_amino_acids_int, total_aa_int)
        dict_highest_three_int, dict_lowest_three_int = \
            highest_lowest_aa_freq(dict_freq_aa_int, total_aa_int)
        title_int = file_name_virus.replace(".txt", " intern")
        # amino_graph_cys_trp(class_amino_acids_int, title_int)
        # amino_graph_hydro(dict_hydro_percentage_int, title_int)
        # amino_graph_highest_lowest(dict_lowest_three_int, dict_highest_three_int,
        #                           title_int)

        # Call everything for envelop viral protein
        dict_freq_aa_env, total_aa_env = freq_amino_acid(dict_prot_seq_env)
        class_amino_acids_env = make_aa_class(dict_freq_aa_env, total_aa_env)
        dict_hydro_percentage_env = calculate_percentages_hydro(
            class_amino_acids_env, total_aa_env)
        dict_highest_three_env, dict_lowest_three_env = \
            highest_lowest_aa_freq(dict_freq_aa_env, total_aa_env)
        title_env = file_name_virus.replace(".txt", " envelop")
        # amino_graph_cys_trp(class_amino_acids_env, title_env)
        # amino_graph_hydro(dict_hydro_percentage_env, title_env)
        # amino_graph_highest_lowest(dict_lowest_three_env, dict_highest_three_env,
        #                           title_env)

    list_proteins = ["Homo sapiens fumarate hydratase.txt",
                      "Homo sapiens glucose transporter.txt",
                      "Homo sapiens GPCR.txt",
                      "Homo sapiens receptor kinase.txt"]

    for i in range(len(list_proteins)):
        file_name_proteins = list_proteins[i]
        dict_prot_seq = read_file(file_name_proteins)
        dict_freq_aa, total_aa = freq_amino_acid(dict_prot_seq)
        class_amino_acids = make_aa_class(dict_freq_aa, total_aa)
        dict_hydro_percentage = calculate_percentages_hydro(class_amino_acids,
                                                            total_aa)
        dict_highest_three, dict_lowest_three = \
            highest_lowest_aa_freq(dict_freq_aa, total_aa)
        title = file_name_proteins.replace(".txt", "")
        amino_graph_cys_trp(class_amino_acids, title)
        # amino_graph_hydro(dict_hydro_percentage, title)
        # amino_graph_highest_lowest(dict_lowest_three, dict_highest_three,
        #                           title)


main()
