# Weektaak 7

import random


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

    # Create empty variables
    env_prot_seq = ""
    header_env = ""

    open_file = open(file, "r")

    for line in open_file:
        # Check if header is for envelop
        if line.startswith(">"):
            header_env += line.strip()
        else:
            env_prot_seq += line.strip()

    open_file.close()

    return env_prot_seq, header_env


def assign_position_to_aa(env_prot_seq):
    """"
    """

    dict_n_post = {}

    for i in range(len(env_prot_seq)):
        aa = env_prot_seq[i]
        dict_n_post[i] = aa

    return dict_n_post


def randomize_100_seqs(dict_n_post):
    """"
    """

    random_seqs = []

    for i in range(100):
        new_order = [0]
        while len(new_order) != len(dict_n_post):
            n = random.randint(1, len(dict_n_post)-1)
            if n not in new_order:
                new_order.append(n)

        new_seq = ""
        for j in range(len(new_order)):
            new_seq += dict_n_post[new_order[j]]

        random_seqs.append(new_seq)

    return random_seqs


def create_headers_list(header_env, header_virus, number_header, random_seqs):
    """"
    """

    list_headers = [header_env]
    header_counter = number_header

    for i in range(len(random_seqs)):
        list_headers.append(">" + str(header_counter) + "|Random|" + header_virus)
        header_counter += 1

    return list_headers, header_counter


def create_sequences_list(env_prot_seq, random_seqs):
    """"
    """

    list_sequences = [env_prot_seq]

    for i in range(len(random_seqs)):
        list_sequences.append(random_seqs[i])

    return list_sequences


def create_file_from_list(list_404_headers, list_404_seqs):
    """"
    """

    with open("random sequences.txt", "w") as new_file:
        for i in range(len(list_404_headers)):
            for j in range(len(list_404_headers[i])):
                new_file.write(str(list_404_headers[i][j]))
                new_file.write("\n")
                new_file.write(str(list_404_seqs[i][j]))
                new_file.write("\n")
                new_file.write("\n")

    new_file.close()


def main():
    list_virus = ["HIV 1 Envelop Protein.txt", "HIV 2 Envelop Protein.txt",
                  "SIV Envelop Protein.txt", "SIVmnd2 Envelop Protein.txt"]

    list_404_headers = []
    list_404_seqs = []

    number_header = 1

    for i in range(len(list_virus)):
        file_name = list_virus[i]

        env_prot_seq, header_env = read_file(file_name)
        dict_n_post = assign_position_to_aa(env_prot_seq)
        random_seqs = randomize_100_seqs(dict_n_post)
        header_virus = file_name.replace(" Envelop Protein", "")
        list_headers, header_counter = create_headers_list(header_env, header_virus, number_header, random_seqs)
        number_header = header_counter
        list_sequences = create_sequences_list(env_prot_seq, random_seqs)

        list_404_headers.append(list_headers)
        list_404_seqs.append(list_sequences)

    create_file_from_list(list_404_headers, list_404_seqs)


main()
