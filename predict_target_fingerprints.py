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
import zipfile
import glob
import os
import sys
import math
import numpy as np
from multiprocessing import Pool
import multiprocessing
from collections import Counter
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
	model_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + '/classes_in_model.txt').read().splitlines()]
	return_dict = {l[0] : l[0:7] for l in model_info}
	return return_dict
	
#unzip a pkl model
def open_Model(mod):
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
	probs = probs.round(2)
	if threshold is not None: probs = np.array(probs >threshold,dtype=int)
	return mod,probs

#prediction runner
def performTargetPrediction(models):
	prediction_results = dict()
	pool = Pool(processes=N_cores, initializer=initPool, initargs=(querymatrix,))  # set up resources
	jobs = pool.imap_unordered(doTargetPrediction, models)
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(models)))*100 + 1
		sys.stdout.write(' Performing Classification on Query Molecules: %3d%%\r' % percent)
		sys.stdout.flush()
		prediction_results[result[0]] = result[1]
	pool.close()
	pool.join()
	prediction_matrix = np.array([j for i,j in sorted(prediction_results.items())]).transpose()
	return sorted(prediction_results.keys()), prediction_matrix

#initializer for the pool
def initPool(querymatrix_):
	global querymatrix
	querymatrix = querymatrix_

#main
if __name__ == '__main__':
	input_name, N_cores,  = sys.argv[1], int(sys.argv[2])
	try:
		threshold = float(sys.argv[3])
	except KeyError:
		if sys.argv[3] == 'None': threshold = 'None'
		else: print 'Enter valid threshold or "None"'
	introMessage()
	print ' Using ' + str(N_cores) + ' Cores'
	try:
		desired_organism = sys.argv[4]
	except IndexError:
		desired_organism = None

	model_info = getUniprotInfo()
	if desired_organism is not None:
		if os.name == 'nt': sep = '\\'
		else: sep = '/'
	models = [modelfile for modelfile in glob.glob(os.path.dirname(os.path.abspath(__file__)) + sep + 'models' + sep + '*.zip')]
		models = [mod for mod in models if model_info[mod.split(sep)[-1].split('.')[0]][4] == desired_organism]
	print ' Total Number of Classes : ' + str(len(models))
	if desired_organism is not None:
		print ' Predicting for organism : ' + desired_organism
		out_name = input_name + '_out_target_fingerprints_' + desired_organism[:3] + '_' + str(threshold) + '.txt'
		out_file = open(out_name, 'w')
	else:
		out_name = input_name + '_out_target_fingerprints_' + str(threshold) + '.txt'
		out_file = open(out_name, 'w')

	#perform target predictions and tp fingerprints to file 
	querymatrix,smiles = importQuery(input_name)
	print ' Total Number of Query Molecules : ' + str(len(querymatrix))
	print ' Using threshold : ' + str(threshold)
	sorted_targets,prediction_matrix = performTargetPrediction(models)
	out_file.write('Compound\t' + '\t'.join(map(str,[i for i in sorted_targets])) + '\n')
	for i, row in enumerate(prediction_matrix):
		#write target prediction fp
		out_file.write(smiles[i] + '\t' + '\t'.join(map(str,row)) + '\n')
	print '\n Wrote Results to: ' + out_name
	out_file.close()
