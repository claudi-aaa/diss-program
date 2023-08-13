from data_generator import generate_sequences
from amino_acid import reassign_aa, make_partition, gen_aa_chain, record_aa, record_org_aa, record_aa_assignment, make_new_aa_phy, make_org_aa_phy, generate_aa_chain_list
from runsplits import run_workflow
from tree import parse_output
from ete3 import Tree
from Bio import SeqIO
import os
import shutil
import csv
import subprocess



for run_no in range(1, 5):
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
    # new_dir_name = f"TS3-{run_no}"
    new_dir_name = f"S4-{run_no}"
    os.mkdir(str(new_dir_name))

    PROJECT_PATH = f'/users/claudiaroth/Documents/GitHub/diss/{new_dir_name}'

    WORKFLOW_FILE = 'workflow5.wflow6'

    # get settings for programs
    # kmer_input = input('Enter the k-mer size e.g. 3: ')
    # kmer_size = int(kmer_input)
    kmer_size = 4
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

    # incorrect line below as this only takes 1 of the branches seq data 
    # start_seq_data = str(records[0].seq)

    seq_data_list = [str(records[0].seq), str(records[1].seq), str(records[2].seq), str(records[3].seq), str(records[4].seq), str(records[5].seq), str(records[6].seq), str(records[7].seq)]
   
    # print(seq_data_list)

    # translate RNA into amino acid chains 
    list_of_new_aa_chains = generate_aa_chain_list(seq_data_list, kmer_size, new_aa_dict)
    list_of_aa_chains = generate_aa_chain_list(seq_data_list, 3, ORGAADICT)

    # print(len(list_of_new_aa_chains))
    # print(list_of_new_aa_chains)
    # print(len(list_of_aa_chains))
    # print(list_of_aa_chains)

    # record amino acid chains of each branch in phylip format 
    make_new_aa_phy(kmer_size, seed_value, list_of_new_aa_chains, new_dir_name)
    make_org_aa_phy(seed_value, list_of_aa_chains, new_dir_name)



    # PARTITIONLENGTH = 18
    # list_of_seq_partitions = []
    # list_of_aa_chains = []
    # list_of_new_aa_chains = []

    # make original tree to make subsequent comparisions to: this will require a loop still * for time series
    # org_partitions = list(make_partition(start_seq_data, PARTITIONLENGTH))
    # org_aa_chain = gen_aa_chain(org_partitions, 3, ORGAADICT)
    # record_org_aa(3, seed_value, org_aa_chain, new_dir_name)



    # for seq in seq_data_list:
    #     partitions = list(make_partition(seq, PARTITIONLENGTH))
    #     list_of_seq_partitions.append(partitions)
        

    #     org_aa_chain = gen_aa_chain(partitions, 3, ORGAADICT)
    #     new_aa_chain = gen_aa_chain(partitions, kmer_size, new_aa_dict)

    #     list_of_aa_chains.append(org_aa_chain)
    #     list_of_new_aa_chains.append(new_aa_chain)










    # ------EXECUTE SPLITSTREE FROM PROGRAM ----------------

    # NOTE! must have some default setting if allowing the user to determine measurement type

    # setting up workflow files to run splits tree 
    WORKFLOW_PATH = '/Applications/SplitsTreeCE/Tools'
    try:
        shutil.copy(WORKFLOW_FILE, PROJECT_PATH)
    except:
        print('error: unable to record workflow for new directory')


    # running for original run 
    input_true_fn = f"org_seed{str(seed_value)}.phy"
    output_true_fn = f"splits_seed{str(seed_value)}_org.nex"
    run_workflow(WORKFLOW_FILE, input_true_fn, output_true_fn, WORKFLOW_PATH, PROJECT_PATH)


    # running for experimental run 
    input_test_fn = f"rand{str(kmer_size)}_seed{str(seed_value)}.phy"
    output_test_fn = f"splits_rand{str(kmer_size)}_seed{str(seed_value)}.nex"
    run_workflow(WORKFLOW_FILE, input_test_fn, output_test_fn, WORKFLOW_PATH, PROJECT_PATH)




    # ------MAKE TREE COMPARISONS --------------------------

    # redirect to project path (from splitstree)
    os.chdir(PROJECT_PATH)

    # isolate only tree info from nexus output
    nwk_tree_1 = parse_output(output_test_fn)
    nwk_tree_2 = parse_output(output_true_fn)

    print(nwk_tree_1)
    print(nwk_tree_2)

    # compare the two trees using ete3 Tree to compute RF
    t1 = Tree(nwk_tree_1)
    t2 = Tree(nwk_tree_2)

    # visualizing the two trees 
    viz_tree_1 = open("viztree1.txt", "w")
    viz_tree_1.write(str(t1))
    viz_tree_1.close()

    viz_tree_2 = open("viztree2.txt", "w")
    viz_tree_2.write(str(t2))
    viz_tree_2.close()

    # saving newick tree files for comparison by quartet distance 
    tree_file_1 = open("tree1.new", "w")
    tree_file_1.write(nwk_tree_1)
    tree_file_1.close()

    tree_file_2 = open("tree2.new", "w")
    tree_file_2.write(nwk_tree_2)
    tree_file_2.close()

    # running the qDist from within program 
    command = "quartet_dist -v tree1.new tree2.new"
    res = subprocess.run(command, capture_output=True, shell=True)
    finalres = res.stdout.decode()
 
    
    thing = " ".join(finalres.split())
    finalres = thing.split()
    print(finalres)


    qheaders = ['leaves', 'total_quartets', 'qd', 'qdnorm', 'resolved_1', 'norm_resolved_q', 'unresolved_q', 'norm_unresolved_q']
    qdata = [finalres[0], finalres[1], finalres[2], finalres[3], finalres[4], finalres[5], finalres[6], finalres[7]]


    quartet_filename = "quartet_report.csv"
    qfile = open(f'{PROJECT_PATH}/{quartet_filename}', 'w')
    qwriter = csv.writer(qfile)
    qwriter.writerow(qheaders)
    qwriter.writerow(qdata)
    qfile.close()


    # calculating RF distance 

    rf, max_rf, common_leaves, parts_t1, parts_t2, discard_t1, discard_t2 = t1.robinson_foulds(t2, unrooted_trees=True)
    t1_minus_t2 = parts_t1 - parts_t2
    t2_minus_t1 = parts_t2 - parts_t1

    headers = ['t1_minus_t2', 't2_minus_t1', 'rf_value', 'max_rf']
    rf_data = [t1_minus_t2, t2_minus_t1, rf, max_rf]


    rf_filename = "rf_report.csv"
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

