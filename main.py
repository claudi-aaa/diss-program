import random_program
import charted_program
import weighted_program



# PROGRAM SETUP FUNCTIONS

def get_ksize():
    """Takes user input and returns an integer greater than 3 representing kmer size being tested for run(s)"""

    print("---------------------------------------------------")
    k_input = input('Enter the k-mer size as integer: ')
    try: 
        if k_input == 'exit':
            print('Exiting program')
            return 
        else:
            k_size = int(k_input)
            if k_size > 2:
                return k_size
            else:
                print('Sorry run must be greater than 2, try again')
                print('To exit program enter "exit"')
                get_ksize()
    except:
        print('Sorry that is not a valid number, try again')
        print('To exit program enter "exit"')
        get_ksize()


def get_seq_size(k):
    """Takes user input and checks if genetic sequence is valid"""

    print("---------------------------------------------------")
    seq_input = input('Enter genetic sequence length as number of bases: ')
    try: 
        if seq_input == 'exit':
            print('Exiting program')
            return 
        else:
            seq_size = int(seq_input)
            if seq_size > 0:
                if (seq_size % k) == 0:
                    print(seq_size)
                    return seq_size
                else:
                    print('Sorry, the k-mer size is not divisble by genetic sequence length, please try again')
                    get_seq_size(k)
            else:
                print('Sorry run must be greater than 0, try again')
                get_seq_size(k)
    except:
        print('Sorry that is not a valid number, try again')
        get_seq_size(k)




def get_runs():
    """Takes user input and returns an integer greater than 0 representing total number of times to run the program"""

    print("---------------------------------------------------")
    run_input = input('Enter the number of runs as integer: ')
    try: 
        if run_input == 'exit':
            print('Exiting program')
            return 
        else:
            total_runs = int(run_input)
            if total_runs > 0:
                return total_runs
            else:
                print('Sorry run must be greater than 0, try again')
                get_runs()
    except:
        print('Sorry that is not a valid number, try again')
        get_runs()


def get_assignment_type():
    """Takes user input and determines amino acid assignment"""

    print("---------------------------------------------------")
    print("Enter the number of the respective assignment you would like to use\n")

    print("1 for random amino acid assignment")
    print("2 for weighted-random assignment based on original 3-mer chart")
    print("3 for weighted-random assignment based on random walk\n")

    aa_type = input('Enter selection as integer: ')
    try: 
        if aa_type == 'exit':
            print('Exiting program')
            return 
        else:
            aa_type = int(aa_type)

            # adapt based on programs available 
            if aa_type > 0 and aa_type < 4: 
                return aa_type
            else:
                print('Sorry not a valid selection, try again')
                get_runs()
    except:
        print('Sorry that is not a valid selection, try again')
        get_runs()



# INITIAL SETUP BY USER PROGRAM


# program intro here
# name, credits, and description intro here by launching program 
# NOTE! --- disclaimer about having splits tree CE installed in certain file path (/aaplications/splitstreeCE/Tools)



# standard settings - could be changed to accept user input at a later date 
WORKFLOW_FILE = 'workflow5.wflow6'
SPLITS_PATH = '/Applications/SplitsTreeCE/Tools'
PROJECT_PATH = '/users/claudiaroth/Desktop/diss-program'
no_of_branches = 8
tree_1 = "(L1:0.5, (L2:0.5, (L3:0.5, (L4:0.5, (L5:0.5, (L6:0.5, (L7:0.5, (L8:0.5):0.5):0.5):0.5):0.5):0.5):0.5):0.5);"


# user defined settings 
kmer_size = get_ksize()
seq_size = get_seq_size(kmer_size)
runs = get_runs()
assignment_type = get_assignment_type()



# run subsequent program based on the assignment type 
if assignment_type == 1:
    random_program.program_1(kmer_size, seq_size, runs, WORKFLOW_FILE, no_of_branches, tree_1, PROJECT_PATH, SPLITS_PATH)
elif assignment_type == 2:    
    charted_program.program_2(kmer_size, seq_size, runs, WORKFLOW_FILE, no_of_branches, tree_1, PROJECT_PATH, SPLITS_PATH)
elif assignment_type == 3:
    weighted_program.program_3(kmer_size, seq_size, runs, WORKFLOW_FILE, no_of_branches, tree_1, PROJECT_PATH, SPLITS_PATH)
else:
    print("sorry that is not a valid program, please try again")







# new_dir_name = f"S{str(kmer_size)}-{run_no}"


# PROJECT_PATH = f'/users/claudiaroth/Desktop/diss-program/{new_dir_name}'
