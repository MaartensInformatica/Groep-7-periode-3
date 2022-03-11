import numpy
import matplotlib.pyplot as plt


def read_file(file_name):
    """
    """

    lines = []
    with open(file_name, "r") as file:
        for line in file:
            line = line.strip()
            if "#" not in line and line != "":
                lines.append(line)

    file.close()

    return lines


def get_percentages(lines):
    """

    """

    list_percentages = []
    for line in lines:
        start_index = line.index("100.00")

        value_string = ""

        index_counter = -1
        for i in line:
            index_counter += 1
            if index_counter < start_index:
                continue

            else:
                value_string += i

        list_percentages.append(value_string.split())

    return list_percentages


def count_percentages(list_percentages):
    """"
    """

    list_values = []

    for i in range(len(list_percentages)):
        for j in range(len(list_percentages[i])):
            value = list_percentages[i][j]

            if value != "100.00":
                list_values.append(float(value))

    return list_values


def calculate_mean(list_values):
    """"
    """

    mean_percentage = numpy.mean(list_values)

    return mean_percentage


def calculate_standard_deviation(list_values):
    """"
    """

    stdev_percentage = numpy.std(list_values)

    return stdev_percentage


def create_histogram(list_values, mean_percentage, stdev_percentage,
                     mean_blast):
    """"
    """

    plt.hist(list_values, density=True, bins=300, color='green')
    plt.axvline(mean_percentage, color='r', linewidth=1)
    plt.axvline(mean_blast, color='yellow', linewidth=1)

    plt.axvline(mean_percentage - stdev_percentage, color='b', linewidth=1)
    plt.axvline(mean_percentage + stdev_percentage, color='b', linewidth=1)

    plt.xlabel("Percentages (%)")
    plt.ylabel("Fraction")
    plt.title("Identities 404 sequences")

    plt.show()


if __name__ == '__main__':
    identities_blast = [34, 36, 31, 45, 34, 34]
    file_name = "results clustalO.txt"
    lines = read_file(file_name)
    list_percentages = get_percentages(lines)
    list_values = count_percentages(list_percentages)
    mean_percentage = calculate_mean(list_values)
    stdev_percentage = calculate_standard_deviation(list_values)
    mean_blast = calculate_mean(identities_blast)
    create_histogram(list_values, mean_percentage, stdev_percentage,
                     mean_blast)
