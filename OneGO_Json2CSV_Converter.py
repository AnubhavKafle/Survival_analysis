import subprocess
import os
from glob import glob
count=1
for i in subprocess.check_output(["find","-iname","*.json"]).splitlines():
	#print i 
	c=i.strip().split("/")[1]
	s="/media/guest-vxtpfp/KINGSTON/Gene-Xpression-Quant_TCGA/"+c
	os.chdir("%s"%s)
	#print os.getcwd()
	subprocess.call(['python', '/media/guest-vxtpfp/KINGSTON/Gene-Xpression-Quant_TCGA/json_to_csv.py', '%s'%c, '%s.csv'%c]+glob('*.json'))
	print count,"Done converting for >>>>>", os.getcwd()
	os.chdir("/media/guest-vxtpfp/KINGSTON/Gene-Xpression-Quant_TCGA")
	count=count+1
	 

print "Done" 	
