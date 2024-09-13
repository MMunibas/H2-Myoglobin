import subprocess
from subprocess import PIPE
from glob import glob

exec_file_path_prefix = "/data/kaeserj/GaussianOpt/OrcaOpt/MyoScan/execfiles/exec_file_"
output_file_path_prefix = "/data/kaeserj/GaussianOpt/OrcaOpt/MyoScan/output/output_file_orca_scan_"
#
# exec_file_path_prefix = "/home/kaeserj/PycharmProjects/CurveFitMorse/ScanGenerator/execfiles/exec_file_"
# output_file_path_prefix = "/home/kaeserj/PycharmProjects/CurveFitMorse/ScanGenerator/output/output_orca_scan_"

exec_file_list = glob(exec_file_path_prefix + "*")
print(exec_file_list)
print(len(exec_file_list))
num_files = len(exec_file_list)

# exec_file_path = exec_file_path_prefix + f"{i:03d}"

for i in range(21, 25):
    exec_file_path = exec_file_path_prefix + f"{i:03d}"
    output_file_path = output_file_path_prefix + f"{i:03d}"
    with open(output_file_path, "w") as file:
        result = subprocess.run(['sbatch', exec_file_path], stdout=file, stderr=file)



