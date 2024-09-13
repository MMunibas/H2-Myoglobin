import os
import numpy as np
from glob import glob

# MDAnalysis
import MDAnalysis

import matplotlib.pyplot as plt

from ase import Atoms
from ase import io
from ase.io import read


def build_Exec_string(input_file_path_specific, output_file_path_specific):
    returnval = "$" + "my_orca" + " " + input_file_path_specific + " > " + output_file_path_specific
    return returnval


def build_JOBNAME_line(file_index):
    return "#SBATCH --job-name=" + "MyoScan" + f"{file_index:03d}"


def build_input_output_line(prefix, file_index):
    return prefix + f"{file_index:03d}"


template_file_path = "orca_tmpl.sh"

input_file_path_prefix = "input/input_orca_scan_"

output_file_path_prefix_on_beethoven = "output/output_orca_scan_"

exec_file_path_prefix = "execfiles/exec_file_"

with open(template_file_path, 'r') as file:
    temp_lines = file.read()

tag1 = "header to replace"
tag2 = "run orca file to replace"

list_files_in_inp_folder = glob(input_file_path_prefix + "*")
num_files = len(list_files_in_inp_folder)
print(num_files)

for i in range(num_files):
    exec_file_path = exec_file_path_prefix + f"{i:03d}"

    input_file_path_beethoven = input_file_path_prefix + str(i)
    output_file_path_beethoven = build_input_output_line(output_file_path_prefix_on_beethoven, i)

    lines = temp_lines[:]

    replacement = build_JOBNAME_line(i)
    lines = lines.replace(tag1, replacement)

    replacement = build_Exec_string(input_file_path_beethoven, output_file_path_beethoven)
    lines = lines.replace(tag2, replacement)

    with open(exec_file_path, 'w') as file:
        file.write(lines)

# #SBATCH --job-name=JOBNAME
# srun $my_orca input_file_path > output_file_path
