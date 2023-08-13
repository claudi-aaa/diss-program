from data_generator import generate_sequences, evolve_seq_arr
from amino_acid import reassign_aa,record_aa_assignment, make_new_aa_ts_phy, make_org_aa_ts_phy, generate_aa_chain_list
from runsplits import run_workflow
from tree import parse_output
from ete3 import Tree
from Bio import SeqIO
import os
import shutil
import csv



for run_no in range(1, 11):
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

    # ------SET PROGRAM SETTINGS VIA USER --------

    # get folder name for where to store files
    # new_dir_name = input('Enter foldername to save files: ')
    # parent_dir = "C:\\Users\\rothc\\Documents\\GitHub\\diss"
    # path = os.path.join(parent_dir, new_dir_name)
    new_dir_name = f"CY9-{run_no}"
    os.mkdir(str(new_dir_name))

    PROJECT_PATH = f'/users/claudiaroth/Documents/GitHub/diss/{new_dir_name}'

    WORKFLOW_FILE = 'workflow5.wflow6'

    # get settings for programs
    # kmer_input = input('Enter the k-mer size e.g. 3: ')
    # kmer_size = int(kmer_input)
    kmer_size = 9
    # try:
    #     kmer_size = int(kmer_input)
    # except ValueError:
    #     print("Kmer input must be a number")

    # seed_input = input("Enter random seed value: ")
    # seed_value = int(seed_input)
    seed_value = run_no
    # try:
    #     seed_value = int(seed_input)
    # except ValueError:
    #     print("Seed value input must be a number")

    # -----GENERATE THE GENETIC DATA -------------

    tree_1 = "(L1:0.5, (L2:0.5, (L3:0.5, (L4:0.5, (L5:0.5, (L6:0.5, (L7:0.5, (L8:0.5):0.5):0.5):0.5):0.5):0.5):0.5):0.5);"
    seq_size = 2520
    no_of_branches = 8

    generate_sequences(tree_1, seq_size, new_dir_name)

    # -------REASSIGN AMINO ACIDS BASED ON K-MER *RANDOMLY----

    BASES = ['A', 'T', 'C', 'G']

    # reassign the amino acids and record in csv file in  user defined folder
    new_aa_dict = reassign_aa(kmer_size, BASES, seed_value)
    record_aa_assignment(kmer_size, seed_value, new_aa_dict, new_dir_name)



    # ------GENERATE AMINO ACID CHAINS FROM SEQ DATA ---------

    fasta_file = f"{new_dir_name}/gendata.fasta"
    records = list(SeqIO.parse(fasta_file, "fasta"))


    # original seq data 
    seq_data_list = [str(records[0].seq), str(records[1].seq), str(records[2].seq), str(records[3].seq), str(records[4].seq), str(records[5].seq), str(records[6].seq), str(records[7].seq)]
   

    # generate subsequent "mutated sequences"
    # !NOTE here is where user input could determine how many genetic cycles to incorporate
    NO_OF_CYCLES = 10
    MUTATION_RATE = 0.05
    seq_ref = []
    seq_ref.append(seq_data_list)


    for i in range(NO_OF_CYCLES - 1):
        starting_data = seq_ref[-1]
        evolved_data = evolve_seq_arr(starting_data, MUTATION_RATE)
        seq_ref.append(evolved_data)
        # print(evolved_data)
        # print(f"type of evovled_data{type(evolved_data)}")
    
    # print(len(seq_ref))


    for i in range(len(seq_ref)):

        # translate RNA into amino acid chains  
        current_branch = []
    
        for branch in seq_ref[i]:
            if type(branch) == list:
                new_branch = ''.join(map(str, branch))
                current_branch.append(new_branch)
            else:
                current_branch.append(branch)
                
        # print(f"current branch{str(i + 1)} = {current_branch}")


        list_of_new_aa_chains = generate_aa_chain_list(current_branch, kmer_size, new_aa_dict)
        list_of_aa_chains = generate_aa_chain_list(current_branch, 3, ORGAADICT)

        # record amino acid chains of each branch in phylip format 
        make_new_aa_ts_phy(i+1, kmer_size, seed_value, list_of_new_aa_chains, new_dir_name)
        make_org_aa_ts_phy(i+1, seed_value, list_of_aa_chains, new_dir_name)
        




    # ------EXECUTE SPLITSTREE FROM PROGRAM ----------------

    # NOTE! must have some default setting if allowing the user to determine measurement type

    # setting up workflow files to run splits tree 
    WORKFLOW_PATH = '/Applications/SplitsTreeCE/Tools'
    try:
        shutil.copy(WORKFLOW_FILE, PROJECT_PATH)
    except:
        print('error: unable to record workflow for new directory')



    for i in range(NO_OF_CYCLES):

        # running for original run 
        input_true_fn = f"org_cycle{str(i + 1)}_seed{str(seed_value)}.phy"
        output_true_fn = f"splits_cycle{str(i + 1)}_seed{str(seed_value)}_org.nex"
        run_workflow(WORKFLOW_FILE, input_true_fn, output_true_fn, WORKFLOW_PATH, PROJECT_PATH)


        # running for experimental run 
        input_test_fn = f"rand{str(kmer_size)}_cycle{str(i + 1)}_seed{str(seed_value)}.phy"
        output_test_fn = f"splits_rand{str(kmer_size)}_cycle{str(i + 1)}_seed{str(seed_value)}.nex"
        run_workflow(WORKFLOW_FILE, input_test_fn, output_test_fn, WORKFLOW_PATH, PROJECT_PATH)





    # ------MAKE TREE COMPARISONS --------------------------

    # redirect to project path (from splitstree)
    os.chdir(PROJECT_PATH)


    for i in range(len(seq_ref)):

        output_test_file = f"splits_rand{str(kmer_size)}_cycle{str(i + 1)}_seed{str(seed_value)}.nex"
        output_true_file = f"splits_cycle{str(i + 1)}_seed{str(seed_value)}_org.nex"

        # isolate only tree info from nexus output
        nwk_tree_1 = parse_output(output_test_file)
        nwk_tree_2 = parse_output(output_true_file)

        # compare the two trees using ete3 Tree to compute RF
        t1 = Tree(nwk_tree_1)
        t2 = Tree(nwk_tree_2)

        rf, max_rf, common_leaves, parts_t1, parts_t2, discard_t1, discard_t2 = t1.robinson_foulds(t2, unrooted_trees=True)

        t1_minus_t2 = parts_t1 - parts_t2
        t2_minus_t1 = parts_t2 - parts_t1

        headers = ['t1_minus_t2', 't2_minus_t1', 'rf_value', 'max_rf']
        rf_data = [t1_minus_t2, t2_minus_t1, rf, max_rf]


        rf_filename = f"cycle{str(i+1)}_rf_report.csv"
        f = open(f'{PROJECT_PATH}/{rf_filename}', 'w')
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerow(rf_data)
        f.close()
    
    os.chdir('/users/claudiaroth/Documents/GitHub/diss/')




    # for i in range(no_of_branches):
    #     branch_no = i + 1
    #     output_test_fn = f"s{str(branch_no)}splits_rand{str(kmer_size)}_seed{str(seed_value)}.nex"
    #     output_true_fn = f"s{str(seed_value)}org_aa.nex"

    #     nwk_tree_1 = parse_output(output_test_fn)
    #     nwk_tree_2 = parse_output(output_true_fn)

    #     # compare the two trees using ete3 Tree to compute RF
    #     t1 = Tree(nwk_tree_1)
    #     t2 = Tree(nwk_tree_2)

    #     rf, max_rf, common_leaves, parts_t1, parts_t2, discard_t1, discard_t2 = t1.robinson_foulds(t2, unrooted_trees=True)

    #     t1_minus_t2 = parts_t1 - parts_t2
    #     t2_minus_t1 = parts_t2 - parts_t1

    #     headers = ['t1_minus_t2', 't2_minus_t1', 'rf_value', 'max_rf']
    #     rf_data = [t1_minus_t2, t2_minus_t1, rf, max_rf]


    #     rf_filename = "rf_report.csv"
    #     f = open(f'{PROJECT_PATH}/{rf_filename}', 'w')
    #     writer = csv.writer(f)
    #     writer.writerow(headers)
    #     writer.writerow(rf_data)
    #     f.close()

