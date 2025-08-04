import os
import shutil

# Define the list of prefixes
file_prefixes = ['R117-2007-06-02-TT03_1','R117-2007-06-02-TT04_1','R117-2007-06-16-TT05_1','R117-2007-06-21-TT03_5','R117-2007-06-21-TT05_1','R119-2007-06-25-TT12_1','R119-2007-06-26-TT08_5','R119-2007-06-27-TT06_4','R119-2007-06-27-TT08_4','R119-2007-06-29-TT11_4','R119-2007-06-30-TT10_4','R119-2007-07-01-TT06_6','R119-2007-07-03-TT05_2','R119-2007-07-03-TT06_3','R119-2007-07-10-TT11_3','R119-2007-07-14-TT05_7','R119-2007-07-14-TT07_8','R119-2007-07-14-TT11_1','R131-2007-09-04-TT02_8','R131-2007-09-04-TT02_9','R131-2007-09-05-TT02_3','R131-2007-09-06-TT02_3','R131-2007-09-07-TT02_2','R131-2007-09-10-TT06_2','R131-2007-09-13-TT12_3','R131-2007-09-14-TT12_5','R131-2007-09-15-TT02_2','R131-2007-09-15-TT06_6','R131-2007-09-17-TT06_1','R131-2007-09-18-TT02_2','R132-2007-10-09-TT06_5','R132-2007-10-09-TT07_2','R132-2007-10-09-TT08_2','R132-2007-10-09-TT09_2','R132-2007-10-11-TT02_1','R132-2007-10-11-TT02_7','R132-2007-10-11-TT03_4','R132-2007-10-11-TT07_1','R132-2007-10-11-TT07_2','R132-2007-10-11-TT09_2','R132-2007-10-16-TT06_3','R132-2007-10-16-TT09_1','R132-2007-10-17-TT09_1','R132-2007-10-18-TT02_6','R132-2007-10-18-TT06_2','R132-2007-10-18-TT08_2','R132-2007-10-19-TT01_3','R132-2007-10-19-TT06_2','R132-2007-10-20-thr30-TT02_7','R132-2007-10-20-thr30-TT08_4','R132-2007-10-21-TT08_1','R132-2007-10-21-TT09_4','R132-2007-10-21-TT10_4','R132-2007-10-23-thr30-TT01_6','R132-2007-10-23-thr30-TT02_3','R132-2007-10-23-thr30-TT02_4','R132-2007-10-23-thr30-TT05_2','R132-2007-10-23-thr30-TT06_5','R132-2007-10-23-thr30-TT07_2','R132-2007-10-29-TT01_9','R132-2007-10-29-TT02_1',] 
# Define input and output folder paths
input_folder = 'D:\\RandomVstrAnalysis\\CellFingerPrint2'
output_folder = 'D:\\RandomVstrAnalysis\\CellFingerPrint2\\lfr_msn_only'

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