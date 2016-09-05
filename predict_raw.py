#Author : Lewis Mervin lhm30@cam.ac.uk
#Supervisor : Dr. A. Bender
#All rights reserved 2016
#Protein Target Prediction Tool trained on SARs from PubChem (Mined 21/06/16) and ChEMBL21
#Molecular Descriptors : 2048bit Morgan Binary Fingerprints (Rdkit) - ECFP4
#Dependencies : rdkit, sklearn, numpy

#libraries
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.naive_bayes import BernoulliNB
import cPickle
import glob
import os
import sys
import math
import numpy as np
from multiprocessing import Pool
import multiprocessing
multiprocessing.freeze_support()

def introMessage():
	print '=============================================================================================='
	print ' Author: Lewis Mervin\n Email:  lhm30@cam.ac.uk\n Supervisor: Dr. A. Bender'
	print ' Address: Centre For Molecular Informatics, Dept. Chemistry, Lensfield Road, Cambridge CB2 1EW'
	print '==============================================================================================\n'
	return

	
#calculate 2048bit morgan fingerprints, radius 2
def calcFingerprints(smiles):
	m1 = Chem.MolFromSmiles(smiles)
	fp = AllChem.GetMorganFingerprintAsBitVect(m1,2, nBits=2048)
	binary = fp.ToBitString()
	return list(binary) 

#calculate fingerprints for chunked array of smiles
def arrayFP(input):
	outfp = []
	for i in input:
		outfp.append(calcFingerprints(i))
	return outfp
	
#import user query
def importQuery(query):
	matrix = []
	smiles_per_core = int(math.ceil(len(query) / N_cores)+1)
	chunked_smiles = [query[x:x+smiles_per_core] for x in xrange(0, len(query), smiles_per_core)]
	pool = Pool(processes=N_cores)  # set up resources
	jobs = pool.imap(arrayFP, chunked_smiles)
	for i, result in enumerate(jobs):
		matrix += result
	pool.close()
	pool.join()
	return np.array(matrix,dtype=np.uint8)

#get info for uniprots
def getUniprotInfo():
	model_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + '/classes_in_model.txt').read().splitlines()]
	return_dict = {l[0] : l[0:5] for l in model_info}
	return return_dict

#prediction worker	
def doTargetPrediction(pickled_model_name):
	with open(pickled_model_name, 'rb') as fid:
		clf = cPickle.load(fid)
	preds = clf.predict(querymatrix)
	if sum(preds) > 0:
		return pickled_model_name.split('/')[-1][:-4],preds
	else: return None

#prediction runner
def performTargetPrediction(models):
	prediction_results = []
	pool = Pool(processes=N_cores)  # set up resources
	jobs = pool.imap_unordered(doTargetPrediction, models)
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(models)))*100 + 1
		sys.stdout.write(' Performing Classification on Query Molecules: %3d%%\r' % percent)
		sys.stdout.flush()
		if result is not None: prediction_results.append(result)
	pool.close()
	pool.join()
	return prediction_results

#main
input_name = sys.argv[1]
N_cores = int(sys.argv[2])
introMessage()
print ' Predicting Targets for ' + input_name
print ' Using ' + str(N_cores) + ' Cores'
models = [modelfile for modelfile in glob.glob(os.path.dirname(os.path.abspath(__file__)) + '/models/*.pkl')]
model_info = getUniprotInfo()
print ' Total Number of Classes : ' + str(len(models))
output_name = input_name + '_out_raw.txt'
out_file = open(output_name, 'w')
smiles = open(input_name).read().splitlines()
querymatrix = importQuery(smiles)
print ' Total Number of Query Molecules : ' + str(len(querymatrix))
prediction_results = performTargetPrediction(models)
out_file.write('Uniprot\tPref_Name\tGene ID\tTarget_Class\tOrganism\t' + '\t'.join(map(str,smiles)) + '\n')
for row in sorted(prediction_results):
	out_file.write('\t'.join(map(str,model_info[row[0]])) + '\t' + '\t'.join(map(str,row[1])) + '\n')
print '\n Wrote Results to: ' + output_name
out_file.close()