# optional imports:
import os
import time
# necessary import:
import difflib

f1path = r'C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Code\python\decoding\detect_activations.py'
f2path = r'C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Code\python\decoding\detect_activations_old.py'

with open(f1path, 'rU') as f1:
    with open(f2path, 'rU') as f2:
        readable_last_modified_time1 = time.ctime(os.path.getmtime(f1path)) # not required
        readable_last_modified_time2 = time.ctime(os.path.getmtime(f2path)) # not required
        print(''.join(difflib.unified_diff(
          f1.readlines(), f2.readlines(), fromfile=f1path, tofile=f2path,
          fromfiledate=readable_last_modified_time1, # not required
          tofiledate=readable_last_modified_time2, # not required
          )))
        difftext = ''.join(difflib.unified_diff(
                f1.readlines(), f2.readlines(), fromfile=f1path, tofile=f2path,
                fromfiledate=readable_last_modified_time1,  # not required
                tofiledate=readable_last_modified_time2,  # not required
                ))
        with open('diffon1and2', 'w') as diff_file:
            diff_file.write(difftext)
print()
