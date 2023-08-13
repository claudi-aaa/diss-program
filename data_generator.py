import pyvolve
import random

def generate_sequences(input_tree, seq_size, folder_name):
    """
    generates genetically evolved DNA data and outputs FASTA file with number of sequences based on number of tree branches
    k-mer length must be a factor in sequence size
    """
    model_type = "nucleotide"

    try:
        phylogeny = pyvolve.read_tree(tree=input_tree)
        my_model = pyvolve.Model(model_type)
        my_partition = pyvolve.Partition(models=my_model, size=seq_size)
        my_evolver = pyvolve.Evolver(partitions=my_partition, tree=phylogeny)

        # need to edit so that the file is being put in the same folder as the run instance + function needs to accept the folder path
        folder_path = f"{folder_name}/gendata.fasta"
        my_evolver(seqfile=folder_path, ratefile=False, infofile=False)
        print('File has been created')

    except:
        print("Sorry, k-mer length must be a multiple of sequence length  e.g. kmer length of 3 or 4 for sequence length of 144")

# Adapted from: https://stackoverflow.com/questions/69576054/random-mutation-in-a-dna-sequence-mutation-rate-in-python
# Last accessed July 17th 2023
def mutate_sequence(dna_seq, mutation_rate):
    """
    accepts string of dna sequence and evolves it and returns new modified dna string
    """
    dna_list = list(dna_seq)
    for i in range(len(dna_seq)):
        rand_no = random.random()
        # print(dna_list)
        if rand_no < mutation_rate:
            mutation_position = random.randint(0, len(dna_list) - 1)
            dna_list[mutation_position] = random.choice(list('ATCG'))
        # print(dna_list)
        return ''.join(dna_list)



def mutate_sequences(arr_of_sequences, mutation_rate):
    """
    takes list of sequences and takes each and evovles each to return list of modified sequences 
    """
    mutations = []
    for seq in arr_of_sequences:
        mutation = mutate_sequence(seq, mutation_rate)
        mutations.append(str(mutation))
    # print(mutations)
    return mutations


def evolve_seq_arr(dna_data_arr, mutation_rate):
    new_dna_arr = []
    for arr in dna_data_arr:
        result = mutate_sequences(arr, mutation_rate)
        new_dna_arr.append(result)
    return new_dna_arr


# test_seq = [['AAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAA'], ['AAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAA'], ['AAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAA']]

# result = mutate_sequences(test_seq, 0.50)
# print(result)

# test = ['AAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAA']
# test = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
# result = mutate_sequences(test, 0.75)
# print(result)

# new_sequences = []
# print(test_seq)

# new_results = evolve_seq_arr(test_seq, 0.80)
# print(new_results)

# for arr in test_seq:
#     result = mutate_sequences(arr, 0.75)
#     new_sequences.append(result)

# print(new_sequences)



# tree_1 = "(L1:0.5, (L2:0.5, (L3:0.5, (L4:0.5, (L5:0.5, (L6:0.5, (L7:0.5, (L8:0.5):0.5):0.5):0.5):0.5):0.5):0.5):0.5);"
# seq_size = 1728
# folder_test = "test_21"

# generate_sequences(tree_1, seq_size, folder_test)