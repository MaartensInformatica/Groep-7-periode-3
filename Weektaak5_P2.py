import matplotlib.pyplot as plt


def append_identity_percentages(percentage_file):
    """
    :param percentage_file:
    :return:
    """

    lines = []
    with open(percentage_file, "r") as file:
        for line in file:
            line = line.strip()
            if "#" not in line and line != "":
                lines.append(line)

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


def create_histogram(list_values):
    """"
    """

    plt.hist(list_values, density=True, bins=300, color='green')

    plt.xlabel("Percentages (%)")
    plt.ylabel("Fraction")
    plt.title("Identities 404 sequences")

    plt.show()


if __name__ == '__main__':
    test_file = "results clustalO.txt"
    list_percentages = append_identity_percentages(test_file)
    list_values = count_percentages(list_percentages)
    create_histogram(list_values)
