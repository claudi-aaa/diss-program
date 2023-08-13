import os
import shutil

WORKFLOW_PATH = '/Applications/SplitsTreeCE/Tools'
PROJECT_PATH = '/users/claudiaroth/Documents/GitHub/diss'


def move_input(input_file, project_path, wf_path):
    """
    takes file from project directory and moves it into to where splitstree workflow terminal is located
    """
    try:
        os.chdir(project_path)
        shutil.move(input_file, wf_path)
        print('Input file moved!')
    except:
        print('error: unable to move input files')


def move_output(output_file, project_path, wf_path):
    """
    moves the output file from the splitstree tool folder and moves it into the defined directory
    """
    try:
        # os.chdir(wf_path)
        shutil.move(output_file, project_path)
        print('Ouput file moved!')
    except:
        print('error: unable to move ouput file')


def run_workflow(workflow_file, input_file, output_file, wf_path, project_path):
    wf_command = f'./workflow-run -w {workflow_file} -i {input_file} -f Phylip -o {output_file}'
    try:
        os.chdir(project_path)
        move_input(workflow_file, project_path, wf_path)
        move_input(input_file, project_path, wf_path)
        os.chdir(wf_path)
        os.system(wf_command)
        move_output(output_file, project_path, wf_path)

        # move the files that were analysed back to project folder
        move_output(workflow_file, project_path, wf_path)
        move_output(input_file, project_path, wf_path)
    except:
        print('A problem occured')


# workflow_file = 'workflow5.wflow6'
# input_file = 'aligndata2.phy'
# output_file = 'output1.nex'
# run_workflow(workflow_file, input_file, output_file, WORKFLOW_PATH, PROJECT_PATH)
# project_dir_path = '/users/claudiaroth/Documents/GitHub/diss/test5'

# move_input("rand3_seed7.phy", project_dir_path, WORKFLOW_PATH)
# move_input("workflow5.wflow6", project_dir_path, WORKFLOW_PATH)


# run_workflow('workflow5.wflow6', 'rand3_seed7.phy',
#              'newoutput.nex', WORKFLOW_PATH, project_dir_path)
