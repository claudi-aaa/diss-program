
import itertools
import pandas as pd
import random 
from Bio import SeqIO
import numpy as np 


BASES = ['A', 'T', 'C', 'G']

ORGAADICT = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W',
}


def reassign_aa(kmer_size, bases, seed_input):
    """
    returns random new amino acid dictionary - can use seed to make reproducible 
    """
    permutations = [i for i in itertools.product(bases, repeat=kmer_size)]
    aa_list = [list(x) for x in permutations]
    aa_perms = [''.join(item) for item in aa_list]
    
    aa_arr = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    random.seed(seed_input)

    new_aa_dict = {}

    while len(aa_perms) > 0:
        current_perm = aa_perms[0]
        random_num = random.randint(0,19)
        random_acid = aa_arr[random_num]
        new_aa_dict[current_perm] = random_acid
        aa_perms.remove(current_perm)

    return new_aa_dict


def reassign_aa_by_chart(kmer_size, bases):
    """
    return new amino acid dictionary based on percentage of amount seen in original amino acid chart
    """

    permutations = [i for i in itertools.product(bases, repeat=kmer_size)]
    aa_list = [list(x) for x in permutations]
    aa_perms = [''.join(item) for item in aa_list]

    # using the skewed percentages based on original chart 
    group_1 = 0.09375
    group_2 = 0.0625
    group_3 = 0.046875
    group_4 = 0.03125
    group_5 = 0.015625

    aa_arr = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X']

    new_aa_arr = []
    total_perms = len(aa_perms)

    for aa in aa_arr:
        # percentage = 0
        if aa == 'S' or aa == 'R' or aa == 'L':
            # percentage = group_1 
            new_aa_arr.extend([aa for i in range(int(total_perms * 0.09375))])
        if aa == 'V' or aa == 'A' or aa == 'G' or aa == 'P' or aa == 'T':
            # percentage = group_2
            new_aa_arr.extend([aa for i in range(int(total_perms * 0.0625))])
        if aa == 'I' or aa == 'X':
            # percentage = group_3
            new_aa_arr.extend([aa for i in range(int(total_perms * 0.046875))])
        if aa == 'F' or aa == 'C' or aa == 'Y' or aa == 'Q' or aa == 'N' or aa == 'H' or aa == 'E' or aa == 'D' or aa == 'K':
            # percentage = group_4
            new_aa_arr.extend([aa for i in range(int(total_perms * 0.03125))])
        if aa == 'M' or aa == 'W':
            # percentage = group_5
            new_aa_arr.extend([aa for i in range(int(total_perms * 0.015625))])


    random.shuffle(new_aa_arr)

    new_aa_dict = {}

    

    while len(aa_perms) > 0:
        current_perm = aa_perms[0]
        new_aa_dict[current_perm] = new_aa_arr[-1]
        aa_perms.remove(current_perm)
        new_aa_arr.pop()
    
    return new_aa_dict



# new_dict = reassign_aa_by_chart(4, BASES)
# print(new_dict)


def reassign_aa_by_weight(kmer_size, bases):
    """
    return new amino acid dictionary based on percentage of amount seen in original amino acid chart
    """

    keepSearching = True

    permutations = [i for i in itertools.product(bases, repeat=kmer_size)]
    aa_list = [list(x) for x in permutations]
    aa_perms = [''.join(item) for item in aa_list]

    aa_arr = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X']

    while keepSearching:
        weights = np.random.randint(1, 100, size=(21))
        sum_rand = np.sum(weights)
        weights = weights/sum_rand 

        new_aa_arr = []
        total_perms = len(aa_perms)

        pos = 0
        for aa in aa_arr:
            if pos < len(weights):
                new_aa_arr.extend([[aa for i in range(int(round(total_perms * (weights[pos])))) ]])
                pos += 1

        new_aa_arr = list(itertools.chain.from_iterable(new_aa_arr))
        # print(len(new_aa_arr))

        # escape while loop when assignments equals total permutations 
        if len(new_aa_arr) == total_perms:
            keepSearching = False

    random.shuffle(new_aa_arr)
    
    new_aa_dict = {}

    while len(aa_perms) > 0:
        current_perm = aa_perms[0]
        new_aa_dict[current_perm] = new_aa_arr[-1]
        aa_perms.remove(current_perm)
        new_aa_arr.pop()
    
    return new_aa_dict







def make_partition(arr, partition_length):
        """
        generates array split into partition of set size
        """
        print(len(arr))

        for i in range(0, len(arr), partition_length):
             yield arr[i:i + partition_length]



def gen_aa_chain(partitions, kmer_length, aa_dict):
    """
    generates amino acid chain based on kmer length, amino acid assignment, and partitions
    """
    arr_to_convert = []
    
    for partition in partitions: 
        temp_arr = []
        for i in range(len(partition) - (kmer_length - 1)):
            temp_arr.append(partition[i: i + (kmer_length)])
        arr_to_convert.append(temp_arr)


    aa_chains = []

    for arr in arr_to_convert:
        aa_chain = []
        for k_mer in arr:
            current_aa = aa_dict[k_mer]
            aa_chain.append(current_aa)
        aa_chains.append(aa_chain)

    return aa_chains

def generate_aa_chain_list(seq_list, kmer_length, aa_dict):
    """
    returns amino acid chain list based on dna seq given, kmer length and amino acid assignment
    """
    arr_to_convert = []

    for item in seq_list:
        temp_arr = []
        
        for i in range(len(item) - (kmer_length - 1)):
            temp_arr.append(item[i: i + (kmer_length)])

        arr_to_convert.append(temp_arr)
    

    aa_chains = []

    for arr in arr_to_convert:
        aa_chain = []
        for k_mer in arr:
            current_aa = aa_dict[k_mer]
            aa_chain.append(current_aa)
        aa_chains.append(aa_chain)

    return aa_chains


def record_aa_assignment(kmer_length, seed_value, aa_dict, folder_name):
    """
    writes to the file to record the amino acid reassignements
    """
    filename = "rand" + str(kmer_length) + "_seed" + str(seed_value) + "_aadict.csv" 
    final_filename = folder_name + "/" + filename
    aa_df = pd.DataFrame(aa_dict, index=[1])
    aa_df.to_csv(final_filename)


def record_aa(kmer_length, seed_value, aa_chains, folder_name):
    """
    writes to file to record amino acid chains
    """
    filename = "rand" + str(kmer_length) + "_seed" + str(seed_value) + ".phy"
    # filepath = "/aa_chains"
    # final_filename = folder_name + filepath + filename
    # file = open(final_filename, "w")


    final_filename = folder_name + "/" + filename
    file = open(final_filename, "w")
    counter = 1 

    num_of_seqs = len(aa_chains)
    len_of_seqs = len(aa_chains[0])
    first_line = str(num_of_seqs) + " " + str(len_of_seqs) + "\n"

    file.write(first_line)
    
    for chain in aa_chains:
        seq_name = "s" + str(counter)
        spaces = "        "
        curr_chain = "".join(chain)
        file_data = seq_name + spaces + curr_chain + "\n"
        file.write(file_data)
        counter += 1
    
    file.close()



def record_org_aa(kmer_length, seed_value, aa_chains, folder_name):
    """
    writes to file to record amino acid chains
    """
    filename = "org_aa.phy" 

    final_filename = folder_name + "/" + filename
    file = open(final_filename, "w")
    counter = 1 

    num_of_seqs = len(aa_chains)
    len_of_seqs = len(aa_chains[0])
    first_line = str(num_of_seqs) + " " + str(len_of_seqs) + "\n"

    file.write(first_line)
    
    for chain in aa_chains:
        seq_name = "s" + str(counter)
        spaces = "        "
        curr_chain = "".join(chain)
        file_data = seq_name + spaces + curr_chain + "\n"
        file.write(file_data)
        counter += 1
    
    file.close()

def record_ts_org_aa(s_number, kmer_length, seed_value, aa_chains, folder_name):
    """
    writes to file to record amino acid chains
    """
    filename = "s" + str(s_number + 1) + "org_aa.phy" 

    final_filename = folder_name + "/" + filename
    file = open(final_filename, "w")
    counter = 1 

    num_of_seqs = len(aa_chains)
    len_of_seqs = len(aa_chains[0])
    first_line = str(num_of_seqs) + " " + str(len_of_seqs) + "\n"

    file.write(first_line)
    
    for chain in aa_chains:
        seq_name = "s" + str(counter)
        spaces = "        "
        curr_chain = "".join(chain)
        file_data = seq_name + spaces + curr_chain + "\n"
        file.write(file_data)
        counter += 1
    
    file.close()


def record_ts_aa(s_number, kmer_length, seed_value, aa_chains, folder_name):
    """
    writes to file to record amino acid chains 
    """
    
    filename = "s" + str(s_number + 1) + "_rand" + str(kmer_length) + "_seed" + str(seed_value) + ".phy"

    final_filename = folder_name + "/" + filename
    file = open(final_filename, "w")

    counter = 1 

    num_of_seqs = len(aa_chains)
    len_of_seqs = len(aa_chains[0])
    first_line = str(num_of_seqs) + " " + str(len_of_seqs) + "\n"

    file.write(first_line)
    
    for chain in aa_chains:
        seq_name = "s" + str(counter)
        spaces = "        "
        curr_chain = "".join(chain)
        file_data = seq_name + spaces + curr_chain + "\n"
        file.write(file_data)
        counter += 1
    
    file.close()

def make_org_aa_phy(seed_value, aa_chains, folder_name):
    """
    writes to file to record amino acid chains in .phy format
    """
    
    filename = "org_seed" + str(seed_value) + ".phy"

    final_filename = folder_name + "/" + filename
    file = open(final_filename, "w")

    counter = 1 

    num_of_seqs = len(aa_chains)
    len_of_seqs = len(aa_chains[0])
    first_line = str(num_of_seqs) + " " + str(len_of_seqs) + "\n"

    file.write(first_line)
    
    for chain in aa_chains:
        seq_name = "s" + str(counter)
        spaces = "        "
        curr_chain = "".join(chain)
        file_data = seq_name + spaces + curr_chain + "\n"
        file.write(file_data)
        counter += 1
    
    file.close()


def make_new_aa_phy(kmer_length, seed_value, aa_chains, folder_name):
    """
    writes to file to record amino acid chains in .phy format
    """
    
    filename = "rand" + str(kmer_length) + "_seed" + str(seed_value) + ".phy"

    final_filename = folder_name + "/" + filename
    file = open(final_filename, "w")

    counter = 1 

    num_of_seqs = len(aa_chains)
    len_of_seqs = len(aa_chains[0])
    first_line = str(num_of_seqs) + " " + str(len_of_seqs) + "\n"

    file.write(first_line)
    
    for chain in aa_chains:
        seq_name = "s" + str(counter)
        spaces = "        "
        curr_chain = "".join(chain)
        file_data = seq_name + spaces + curr_chain + "\n"
        file.write(file_data)
        counter += 1
    
    file.close()

def make_org_aa_ts_phy(cycle, seed_value, aa_chains, folder_name):
    """
    writes to file to record amino acid chains in .phy format
    """

    filename = f"org_cycle{str(cycle)}_seed{str(seed_value)}.phy"


    final_filename = folder_name + "/" + filename
    file = open(final_filename, "w")

    counter = 1 

    num_of_seqs = len(aa_chains)
    len_of_seqs = len(aa_chains[0])
    first_line = str(num_of_seqs) + " " + str(len_of_seqs) + "\n"

    file.write(first_line)
    
    for chain in aa_chains:
        seq_name = "s" + str(counter)
        spaces = "        "
        curr_chain = "".join(chain)
        file_data = seq_name + spaces + curr_chain + "\n"
        file.write(file_data)
        counter += 1
    
    file.close()


def make_new_aa_ts_phy(cycle, kmer_length, seed_value, aa_chains, folder_name):
    """
    writes to file to record amino acid chains in .phy format
    """
    filename = f"rand{str(kmer_length)}_cycle{str(cycle)}_seed{str(seed_value)}.phy"
    final_filename = folder_name + "/" + filename
    file = open(final_filename, "w")

    counter = 1 

    num_of_seqs = len(aa_chains)
    len_of_seqs = len(aa_chains[0])
    first_line = str(num_of_seqs) + " " + str(len_of_seqs) + "\n"

    file.write(first_line)
    
    for chain in aa_chains:
        seq_name = "s" + str(counter)
        spaces = "        "
        curr_chain = "".join(chain)
        file_data = seq_name + spaces + curr_chain + "\n"
        file.write(file_data)
        counter += 1
    
    file.close()














