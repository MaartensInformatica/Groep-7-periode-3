import matplotlib.pyplot as mpl


def create_group_dict(step):
    """
    Creates a dictionary and adds a key for every percentage group.
    (0-2%, 2-4%, ...)
    :return:
        percentage_groups -> Dict -> Dictionary containing the
                                    group percentage count.
    """
    percentage_groups = {}

    counter = 0
    while counter != 100:
        group = str(counter) + "-" + str(counter + step)
        percentage_groups[group] = 0
        counter += step

    return percentage_groups


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

    percentages_list = []
    for line in lines:
        start_index = line.index('100.00')

        value_string = ""

        index_counter = -1
        for i in line:
            index_counter += 1
            if index_counter < start_index:
                continue

            else:
                value_string += i

        percentages_list.append(value_string.split())

    return percentages_list


def count_percentages(percentage_list):
    values_dict = {}

    for i in range(len(percentage_list)):
        for j in range(len(percentage_list[i])):
            value = percentage_list[i][j]

            if value != "100.00":
                if value in values_dict.keys():
                    values_dict[value] += 1

                else:
                    values_dict[value] = 1

    return values_dict


def add_percentages_to_dict(percentage_groups, value_count_dict):
    for i in value_count_dict.keys():
        for key in percentage_groups.keys():
            key_split = key.split("-")

            if int(key_split[0]) < float(i) <= int(key_split[1]):
                percentage_groups[key] += 1

    return percentage_groups


def create_histogram(percentage_groups):
    # keys
    x_as = []
    # frequencies
    y_as = []

    for keys, values in percentage_groups.items():
        x_as.append(keys)
        y_as.append(values)

    mpl.bar(x_as, y_as)
    mpl.xlabel("Percentage groups")
    mpl.ylabel("Frequency")
    mpl.tight_layout()
    mpl.xlim(0, 24)
    mpl.show()


def calculate_values(percentage_groups):
    totaal = 0
    gemiddelde = 0

    value_counter = 0
    for key in percentage_groups:
        value = percentage_groups[key]
        totaal += value
        value_counter += 1

    gemiddelde = totaal / value_counter
    print(gemiddelde)

if __name__ == '__main__':
    test_file = "44 test.txt"
    # bepaalt de grootte van de percentage groups.
    percentage_group_step = 2
    percentage_groups = create_group_dict(percentage_group_step)
    percentage_list = append_identity_percentages(test_file)
    value_count_dict = count_percentages(percentage_list)
    percentage_groups = add_percentages_to_dict(percentage_groups,
                                                value_count_dict)
    create_histogram(percentage_groups)
    calculate_values(percentage_groups)