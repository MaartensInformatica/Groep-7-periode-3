# Weektaak 2
import matplotlib.pyplot as plt
import numpy as np


def strip_header(file):
    """"
    """

    list_genome = []

    open_file = open(file, "r")

    for line in open_file:
        if line.startswith(">"):
            continue
        else:
            line = line.strip()
            list_genome.append(line)

    open_file.close()

    return list_genome


def count_per_100_bp(list_genome):
    """"
    """

    list_genome_100 = []

    count = 0
    total_count = 0
    hundred_bp = ""
    for line in list_genome:
        total_count += 1
        for nucleotide in line:
            hundred_bp += nucleotide
            count += 1
            if count == 100:
                list_genome_100.append(hundred_bp)
                hundred_bp = ""
                count = 0
        if total_count == len(list_genome):
            list_genome_100.append(hundred_bp)

    return list_genome_100


def count_per_10000_bp(list_genome):
    """"
    """

    list_genome_10000 = []

    count = 0
    total_count = 0
    ten_thousand_bp = ""
    for line in list_genome:
        total_count += 1
        for nucleotide in line:
            ten_thousand_bp += nucleotide
            count += 1
            if count == 10000:
                list_genome_10000.append(ten_thousand_bp)
                ten_thousand_bp = ""
                count = 0
        if total_count == len(list_genome):
            list_genome_10000.append(ten_thousand_bp)

    return list_genome_10000


def count_per_100000_bp(list_genome):
    """"
    """

    list_genome_100000 = []

    count = 0
    total_count = 0
    hundred_thousand_bp = ""
    for line in list_genome:
        total_count += 1
        for nucleotide in line:
            hundred_thousand_bp += nucleotide
            count += 1
            if count == 100000:
                list_genome_100000.append(hundred_thousand_bp)
                hundred_thousand_bp = ""
                count = 0
        if total_count == len(list_genome):
            list_genome_100000.append(hundred_thousand_bp)

    return list_genome_100000


def count_gc_percentage(list_genome_per_bp):
    """"
    """

    list_genome_gc_percentage = []
    list_genome_n_percentage = []

    gc_count = 0
    n_count = 0
    at_count = 0

    for line in list_genome_per_bp:
        for nucleotide in line:
            if nucleotide == "G" or nucleotide == "C":
                gc_count += 1

            elif nucleotide == "A" or nucleotide == "T":
                at_count += 1

            elif nucleotide == "N":
                n_count += 1

        total_known = gc_count + at_count
        if total_known != 0 or gc_count != 0:
            gc_percentage = gc_count / total_known * 100
        else:
            gc_percentage = 0
        list_genome_gc_percentage.append(gc_percentage)
        if n_count != 0:
            n_percentage = n_count / (len(line)) * 100
            list_genome_n_percentage.append(n_percentage)
        else:
            n_percentage = 0
            list_genome_n_percentage.append(n_percentage)
        gc_count = 0
        at_count = 0
        n_count = 0

    return list_genome_gc_percentage, list_genome_n_percentage


def create_graph(list_genome_bp, list_genome_gc_percentage,
                 list_genome_n_percentage, title):
    """"
    """

    x_values = []
    x_value = 0

    for i in range(len(list_genome_bp)):
        x_len = len(list_genome_bp[i])
        x_value += x_len
        x_values.append(x_value)

    mean = []
    for i in range(len(x_values)):
        mean.append(np.mean(list_genome_gc_percentage))

    median = []
    for i in range(len(x_values)):
        median.append(np.median(list_genome_gc_percentage))

    print(title)
    print("mean: ", np.mean(list_genome_gc_percentage))
    print("median: ", np.median(list_genome_gc_percentage))
    print("variance: ", np.var(list_genome_gc_percentage))

    plt.plot(x_values, list_genome_gc_percentage, color="darkgoldenrod",
             label="GC percentage")
    plt.plot(x_values, list_genome_n_percentage, color="gold",
             label="N percentage", alpha=0.5)
    plt.plot(x_values, mean, color="forestgreen", label="Mean")
    plt.plot(x_values, median, color="maroon", label="Median")

    plt.title(title)
    plt.xlabel("Basepairs (bp)")
    plt.ylabel("Percentage (%)")

    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
    plt.tight_layout()
    plt.show()


def main():
    genome_HIV1 = "HIV 1 Genome.fasta"
    genome_HIV2 = "HIV 2 Genome.fasta"
    genome_SIV = "SIV Genome.fasta"
    genome_SIVmnd2 = "SIVmnd2 Genome.fasta"
    genome_H_sapiens = "Homo sapiens Chromosome X.fasta"
    genome_H_annuus = "Helianthus annuus Chromosome 7.fasta"
    genome_P_aeruginosa = "Pseudomonas aeruginosa PAO1 genome.fasta"
    genome_S_aureus = "Staphylococcus aureus genome.fasta"

    list_HIV1_genome = strip_header(genome_HIV1)
    list_HIV1_genome_100 = count_per_100_bp(list_HIV1_genome)
    list_HIV1_gc_percentage, list_HIV1_n_percentage = \
        count_gc_percentage(list_HIV1_genome_100)
    title_HIV1 = "HIV 1"
    create_graph(list_HIV1_genome_100, list_HIV1_gc_percentage,
                 list_HIV1_n_percentage, title_HIV1)

    list_HIV2_genome = strip_header(genome_HIV2)
    list_HIV2_genome_100 = count_per_100_bp(list_HIV2_genome)
    list_HIV2_gc_percentage, list_HIV2_n_percentage = \
        count_gc_percentage(list_HIV2_genome_100)
    title_HIV2 = "HIV 2"
    create_graph(list_HIV2_genome_100, list_HIV2_gc_percentage,
                 list_HIV2_n_percentage, title_HIV2)

    list_SIV_genome = strip_header(genome_SIV)
    list_SIV_genome_100 = count_per_100_bp(list_SIV_genome)
    list_SIV_gc_percentage, list_SIV_n_percentage = \
        count_gc_percentage(list_SIV_genome_100)
    title_SIV = "SIV"
    create_graph(list_SIV_genome_100, list_SIV_gc_percentage,
                 list_SIV_n_percentage, title_SIV)

    list_SIVmnd2_genome = strip_header(genome_SIVmnd2)
    list_SIVmnd2_genome_100 = count_per_100_bp(list_SIVmnd2_genome)
    list_SIVmnd2_gc_percentage, list_SIVmnd2_n_percentage = \
        count_gc_percentage(list_SIVmnd2_genome_100)
    title_SIVmnd2 = "SIVmnd2"
    create_graph(list_SIVmnd2_genome_100, list_SIVmnd2_gc_percentage,
                 list_SIVmnd2_n_percentage, title_SIVmnd2)

    list_H_sapiens = strip_header(genome_H_sapiens)
    list_H_sapiens_100000 = count_per_100000_bp(list_H_sapiens)
    list_H_sapiens_gc_percentage, list_H_sapiens_n_percentage = \
        count_gc_percentage(list_H_sapiens_100000)
    title_H_sapiens = "Homo sapiens chromosoom X"
    create_graph(list_H_sapiens_100000, list_H_sapiens_gc_percentage,
                 list_H_sapiens_n_percentage, title_H_sapiens)

    list_H_annuus = strip_header(genome_H_annuus)
    list_H_annuus_100000 = count_per_100000_bp(list_H_annuus)
    list_H_annuus_gc_percentage, list_H_annuus_n_percentage = \
        count_gc_percentage(list_H_annuus_100000)
    title_H_annuus = "Helianthus annuus chromosoom 7"
    create_graph(list_H_annuus_100000, list_H_annuus_gc_percentage,
                 list_H_annuus_n_percentage, title_H_annuus)

    list_P_aeruginosa = strip_header(genome_P_aeruginosa)
    list_P_aeruginosa_10000 = count_per_10000_bp(list_P_aeruginosa)
    list_P_aeruginosa_gc_percentage, list_P_aeruginosa_n_percentage = \
        count_gc_percentage(list_P_aeruginosa_10000)
    title_P_aeruginosa = "Pseudomonas aeruginosa"
    create_graph(list_P_aeruginosa_10000, list_P_aeruginosa_gc_percentage,
                 list_P_aeruginosa_n_percentage, title_P_aeruginosa)

    list_S_aureus = strip_header(genome_S_aureus)
    list_S_aureus_10000 = count_per_10000_bp(list_S_aureus)
    list_S_aureus_gc_percentage, list_S_aureus_n_percentage = \
        count_gc_percentage(list_S_aureus_10000)
    title_S_aureus = "Staphylococcus aureus"
    create_graph(list_S_aureus_10000, list_S_aureus_gc_percentage,
                 list_S_aureus_n_percentage, title_S_aureus)


main()
