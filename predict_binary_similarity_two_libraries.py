#Author : Lewis Mervin lhm30@cam.ac.uk
#Supervisor : Dr. A. Bender
#All rights reserved 2016
#Protein Target Prediction Tool trained on SARs from PubChem (Mined 21/06/16) and ChEMBL21
#Molecular Descriptors : 2048bit Morgan Binary Fingerprints (Rdkit) - ECFP4
#Dependencies : rdkit, sklearn, numpy

#libraries
from rdkit import Chem
from rdkit.Chem import AllChem
import cPickle
import glob
import zipfile
import os
import sys
import math
import numpy as np
from multiprocessing import Pool
import multiprocessing
from scipy.spatial.distance import rogerstanimoto
from scipy.spatial.distance import jaccard

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
def arrayFP(inp):
	outfp = []
	outsmi = []
	for i in inp:
		try:
			outfp.append(calcFingerprints(i))
			outsmi.append(i)
		except:
			print 'SMILES Parse Error: ' + i
	return outfp,outsmi

#import user query
def importQuery(in_file):
	query = open(in_file).read().splitlines()
	matrix = np.empty((len(query), 2048), dtype=np.uint8)
	smiles_per_core = int(math.ceil(len(query) / N_cores)+1)
	chunked_smiles = [query[x:x+smiles_per_core] for x in xrange(0, len(query), smiles_per_core)]
	pool = Pool(processes=N_cores)  # set up resources
	jobs = pool.imap(arrayFP, chunked_smiles)
	current_end = 0
	processed_smi = []
	for i, result in enumerate(jobs):
		matrix[current_end:current_end+len(result[0]), :] = result[0]
		current_end += len(result[0])
		processed_smi += result[1]
	pool.close()
	pool.join()
	return matrix[:current_end], processed_smi

#get info for uniprots
def getUniprotInfo():
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	model_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + sep + 'classes_in_model.txt').read().splitlines()]
	return_dict = {l[0] : l[0:8] for l in model_info}
	return return_dict

#unzip a pkl model
def open_Model(mod):
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	with zipfile.ZipFile(os.path.dirname(os.path.abspath(__file__)) + sep + 'models' + sep + mod + '.pkl.zip', 'r') as zfile:
		with zfile.open(mod + '.pkl', 'r') as fid:
			clf = cPickle.load(fid)
	return clf

#prediction worker	
def doTargetPrediction(pickled_model_name):
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	mod = pickled_model_name.split(sep)[-1].split('.')[0]
	clf = open_Model(mod)
	probs = clf.predict_proba(querymatrix)[:,1]
	preds = map(int,probs > threshold)
	return preds

#prediction runner
def performTargetPrediction(models):
	prediction_results = []
	pool = Pool(processes=N_cores, initializer=initPool, initargs=(querymatrix,threshold,))  # set up resources
	jobs = pool.imap(doTargetPrediction, sorted(models))
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(models)))*100 + 1
		sys.stdout.write(' Performing Classification on Query Molecules: %3d%%\r' % percent)
		sys.stdout.flush()
		if result is not None: prediction_results.append(result)
	pool.close()
	pool.join()
	return np.array(prediction_results)

#initializer for the pool
def initPool(querymatrix_, threshold_):
	global querymatrix, threshold
	querymatrix = querymatrix_
	threshold = threshold_

#main
if __name__ == '__main__':
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	input_name = sys.argv[1]
	input_name2 = sys.argv[2]
	N_cores = int(sys.argv[3])
	introMessage()
	print ' Predicting Targets for ' + input_name
	print ' Using ' + str(N_cores) + ' Cores'
	try:
		threshold = float(sys.argv[3])
	except ValueError:
		print 'ERROR: Enter a valid float (max 2 decimal places) for threshold'
		quit()
	models = [modelfile for modelfile in glob.glob(os.path.dirname(os.path.abspath(__file__)) + sep + 'models' + sep + '*.zip')]
	model_info = getUniprotInfo()
	print ' Total Number of Classes : ' + str(len(models))
	print ' Using TPR threshold of : ' + str(threshold)
	output_name = input_name + '_out_binary_' + str(threshold) + '.txt'
	out_file = open(output_name, 'w')
	querymatrix,smiles = importQuery(input_name)
	querymatrix2,smiles2 = importQuery(input_name2)
	print ' Total Number of Query Molecules file 1 : ' + str(len(querymatrix))
	print ' Total Number of Query Molecules file 2 : ' + str(len(querymatrix2))
	prediction_results = performTargetPrediction(models)
	prediction_results2 = performTargetPrediction(models)
	sim_output = []
	sim_output2 = []
	for i in len(prediction_results):
		sim_output.append(rogerstanimoto(prediction_results[:,i],prediction_results2[:,i]))
		sim_output2.append(jaccard(prediction_results[:,i],prediction_results2[:,i]))
	out_file.write('Compound Pair No.\tSmiles 1\tSmiles 2\tJaccard Sim\n')
	for idx, comp1 in enumerate(smiles):
		comp2 = smiles2[idx]
		s = sim_output[idx]
		s2 = sim_output[idx]
		out_file.write('\t'.join(map(str,[idx,comp1,comp2,s,s2])) + '\n')
	print '\n Wrote Results to: ' + output_name
	out_file.close()