#!/software/miniconda/bin/python
import os
import re
import sys
import time
import hashlib
import datetime
import subprocess
import pandas as pd
from pathlib import Path
from argparse import ArgumentParser
from ruamel.yaml import YAML

v_type = ['ALL', 'DNA', 'RNA']
pwd = Path(os.path.split(os.path.realpath(__file__))[0])
old_yml = 'input/.old_files.yaml'
old_dic = YAML().load(open(old_yml))
today = str(datetime.date.today())

def getInputPath():
	"""获取输入文件
	"""
	input_dir = pwd.joinpath('input')
	pn_path = ''
	re_path = ''
	atime1 = 0
	atime2 = 0
	for file_path in list(input_dir.glob('*.xlsx')):
		f_atime = os.path.getatime(file_path)
		if re.search('positive', str(file_path), re.I):
			md5hash = hashlib.md5(str(pd.read_excel(file_path)))
			md5 = md5hash.hexdigest()
			if md5 in old_dic:
				continue
			if f_atime > atime1:
				pn_path = file_path
				md5_1  = md5
				atime1 = f_atime
		elif re.search('requantification', str(file_path), re.I):
			md5hash = hashlib.md5(str(pd.read_excel(file_path)))
			md5 = md5hash.hexdigest()
			if md5 in old_dic:
				continue
			if f_atime > atime2:
				re_path = file_path
				md5_2  = md5
				atime2 = f_atime
	
	if pn_path == '' or re_path == '':
		sys.stderr.write('ERROR: Check input!\n')
		raise SystemExit(1)
	else:
		old_dic.update({md5_1: pn_path, md5_2: re_path})
		with open(old_yml, 'w') as fw:
			YAML().dump(old_dic, fw)
		return pn_path, re_path

def checkInputFile(pn_path, re_path):
	"""检查输入文件表头
	"""
	pn_df = pd.read_excel(pn_path)
	re_df = pd.read_excel(re_path)
	if 'sample' not in pn_df.columns:
		sys.stderr.write(f'ERROR: Check input file {pn_path} column name!\n')
		raise SystemExit(1)
	if 'sample' not in re_df or 'vsname' not in re_df.columns or 'vtype' not in re_df:
		sys.stderr.write(f'ERROR: Check input file {re_path} column name!\n')
		raise SystemExit(1)
	return

def callWorker():
	# pn_path, re_path = getInputPath()
	pn_path = 'input/batfeces_positive_0724_rerun_v3_filtered.xlsx'
	re_path = 'input/EastAfrica_Requantification_0724_rerun_v1_batfecesfiltedspecies_filtered.xlsx'
	checkInputFile(pn_path, re_path)
	work_r = pwd.joinpath('src/work.R')
	output_dir = pwd.joinpath('output')

	# print('>>>> 准备输入文件...')
	# cmd = f'{work_r} -p {pn_path} -r {re_path}'
	# try:
	# 	subprocess.run(cmd, check=True, shell=True)
	# except subprocess.CalledProcessError as e:
	# 	sys.stderr.write(f'ERROR: in trimming data.\n')
	# 	raise SystemExit(1)
	
	print('>>>> 投递任务拟合模型...')
	shell_fun = """
 	qsubst () {qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:$1 -l vf=$2,num_proc=$1 $3}
	qsubsuper () {qsub -cwd -l vf=$2,num_proc=$1 -P P20Z10200N0206_super -binding linear:$1 -q st_supermem.q $3}
	"""
	shell_root = output_dir.joinpath(today)
	shell_root.mkdir(parents=True, exist_ok=True)
	shell_path = str(output_dir.joinpath(today)) + f'/GAM_submit.sh'
	for v in v_type:
		script_path = str(output_dir.joinpath(today)) + f'/fit_{v}.sh'
		script_cmd  = f'{work_r} -v {v}\n'
		with open(script_path, 'w') as f:
			f.write(script_cmd)
		shell_fun += f'qsubst 10 20G {script_path}\n'
		shell_fun += f'qsubsuper 5 100G {script_path}\n'
	
	with open(shell_path, 'w') as f:
		f.write(shell_fun)
	cmd = f'sh {shell_path}'
	try:
		subprocess.run(cmd, check=True, shell=True)
	except subprocess.CalledProcessError as e:
		sys.stderr.write(f'ERROR: in submit GAM fitting.\n')
		raise SystemExit(1)
	return

def gitlabMonitor():
	"""实现CICD
	"""
	return

if __name__=='__main__':
	print('>>>> Running...')
	callWorker()
	print('>>>> All GAM submitted.')




