import os
import re
from glob import glob

path = '/data/user/alldata/sequencing/20250106_huangc_hd/J23002058_2'
sample_name = path.split('/')[-1]
patterm = fr'({sample_name})-(\d+)_S1_L001_R(\d)_001\.fastq\.gz'

# 定义正则表达式模式
pattern = re.compile(patterm)

# 遍历当前目录下的所有文件
for filename in sorted(glob(path + '/*.fastq.gz')):
    _filname = os.path.basename(filename)
    match = pattern.match(_filname)
    if match:
        # 提取匹配的组
        sample_name = match.group(1)
        lane_number = match.group(2)
        read_number = match.group(3)
        
        # 构建新的文件名
        new_lane = f'L00{lane_number}'
        new_filename = f'{sample_name}_S1_{new_lane}_R{read_number}_001.fastq.gz'
        
        print(f'ln -s {filename} {new_filename}')
