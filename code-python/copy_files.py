import os
import shutil

# Define the list of prefixes
file_prefixes = ['R117-2007-06-07-TT08_1','R117-2007-06-11-TT03_1','R117-2007-06-12-TT03_7','R117-2007-06-13-TT03_1','R117-2007-06-15-TT03_3','R117-2007-06-15-TT07_1','R117-2007-06-16-TT07_1','R117-2007-06-19-TT07_3','R117-2007-06-20-TT07_1','R119-2007-06-29-TT01_1','R119-2007-07-03-TT11_1','R119-2007-07-07-TT01_3','R131-2007-09-07-TT12_5','R132-2007-10-21-TT05_3','R132-2007-10-21-TT05_7','R132-2007-10-23-thr30-TT05_6','R132-2007-10-23-thr30-TT12_3','R132-2007-10-29-TT05_3']
# Define input and output folder paths
input_folder = 'D:\\RandomVstrAnalysis\\CellFingerprint'
output_folder = 'D:\\RandomVstrAnalysis\\lfr_fsi_only'

# Ensure output directory exists
os.makedirs(output_folder, exist_ok=True)

# Loop through files in the input folder
for filename in os.listdir(input_folder):
    for prefix in file_prefixes:
        if filename.startswith(prefix):
            src_path = os.path.join(input_folder, filename)
            dst_path = os.path.join(output_folder, filename)
            shutil.copy2(src_path, dst_path)
            break  # Avoid checking other prefixes once matched