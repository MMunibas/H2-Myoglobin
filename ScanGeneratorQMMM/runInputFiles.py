import subprocess
from glob import glob

input_file_path_prefix = "/WeirdMyoScan/input/input_orca_scan_"
output_file_path_prefix = "/WeirdMyoScan/output/output_orca_scan_"

input_file_list = glob(input_file_path_prefix + "*")
print(input_file_list)
print(len(input_file_list))
num_files = len(input_file_list)

for i in range(num_files):
    input_file_path = input_file_path_prefix + str(i)
    result = subprocess.run(['/home/kaeserj/program/orca5/orca', input_file_path], capture_output=True, text=True)
    text = result.stdout

    with open(output_file_path_prefix + str(i), "w") as file:
        file.write(text)