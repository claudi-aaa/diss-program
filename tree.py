import re


regex1 = r"""
(?<=\[TREES\])(.*)(?=\[TREES\])
"""

regex2 = r"""
[0-9:.,()]+
"""

def parse_output(filepath):
    """
    takes nexus output and returns only tree block within ()
    """
    with open(filepath, 'r') as nexus_file:
        data = nexus_file.read()
        data = data.strip()

    matches = re.findall(regex1, data, re.VERBOSE | re.MULTILINE | re.IGNORECASE | re.DOTALL)
        
    new_matches = re.findall(regex2, str(matches), re.VERBOSE | re.MULTILINE | re.IGNORECASE | re.DOTALL)

    tree_data = str(new_matches[1])
    tree_string = tree_data + ';'
    return tree_string

# print(parse_output('./test10/splits_rand3_seed5.nex'))


# parse_output('./test7/splits_rand3_seed6.nex')

# with open('./test7/splits_rand3_seed6.nex', 'r') as nexus_file:
#     data = nexus_file.read()
#     data = data.strip()




# regex = r"""
# (?<=\[TREES\])(.*)(?=\[TREES\])
# """

# # regex2 = r"""
# # (?<=\()(.*)(?=[:punct:])
# # """

# regex3 = r"""
# [0-9:.,()]+
# """

# matches = re.findall(regex, data, re.VERBOSE |
#                      re.MULTILINE | re.IGNORECASE | re.DOTALL)
# # print(matches)

# new_matches = re.findall(regex3, str(matches), re.VERBOSE |
#                          re.MULTILINE | re.IGNORECASE | re.DOTALL)


# tree_data = str(new_matches[1])
# print(tree_data)
# print(str(new_matches[1]))



