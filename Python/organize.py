import glob
import shutil
import os
filelist = glob.glob('Output/*.pkl')

for file in filelist:
    filename = file.split('/')[-1]
    if filename[7]=='_':
        print('Moving '+file)
        if not os.path.exists('Output/'+filename[:7]+'/'):
            os.mkdir('Output/'+filename[:7])
        shutil.move(file, 'Output/'+filename[:7]+'/'+filename)

print('Organizing output files finished.')