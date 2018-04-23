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
	#collect IDs, if present
	if len(query[0].split()) > 1:
		ids = [line.split()[1] for line in query]
		query = [line.split()[0] for line in query]
	else:
		ids = None
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
	#remove IDs of SMILES parsing errors
	if ids:
		processed_ids = []
		for idx, smi in enumerate(query):
			if smi in processed_smi:
				processed_ids.append(ids[idx])
		ids = processed_ids
	#if IDs weren't present, use SMILES as IDs
	else:
		ids = processed_smi
	return matrix[:current_end], processed_smi, ids

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
	preds = clf.predict_proba(querymatrix)[:,1]
	if sum(preds) > 0:
		return mod,preds
	else: return None

#prediction runner
def performTargetPrediction(models):
	prediction_results = []
	pool = Pool(processes=N_cores, initializer=initPool, initargs=(querymatrix,))  # set up resources
	jobs = pool.imap_unordered(doTargetPrediction, models)
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(models)))*100 + 1
		sys.stdout.write(' Performing Classification on Query Molecules: %3d%%\r' % percent)
		sys.stdout.flush()
		if result is not None: prediction_results.append(result)
	pool.close()
	pool.join()
	return prediction_results

#initializer for the pool
def initPool(querymatrix_):
	global querymatrix
	querymatrix = querymatrix_

#main
if __name__ == '__main__':
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	multiprocessing.freeze_support()
	input_name = sys.argv[1]
	N_cores = int(sys.argv[2])
	introMessage()
	print ' Predicting Targets for ' + input_name
	print ' Using ' + str(N_cores) + ' Cores'
	try:
		desired_organism = sys.argv[3]
	except IndexError:
		desired_organism = None
	models = [modelfile for modelfile in glob.glob(os.path.dirname(os.path.abspath(__file__)) + sep + 'models' + sep + '*.zip')]
	model_info = getUniprotInfo()
	if desired_organism is not None:
		models = [mod for mod in models if model_info[mod.split(sep)[-1].split('.')[0]][4] == desired_organism]
		print ' Predicting for organism : ' + desired_organism
		out_name = input_name + '_out_raw_' + desired_organism[:3] + '.txt'
		out_file = open(out_name, 'w')
	else:
		out_name = input_name + '_out_raw.txt'
		out_file = open(out_name, 'w')
	print ' Total Number of Classes : ' + str(len(models))
	querymatrix,smiles,ids = importQuery(input_name)
	print ' Total Number of Query Molecules : ' + str(len(querymatrix))
	prediction_results = performTargetPrediction(models)
	out_file.write('Uniprot\tPref_Name\tGene ID\tTarget_Class\tOrganism\tPDB_ID\tDisGeNET_Diseases_0.06\tChEMBL_First_Published\t' + '\t'.join(map(str,ids)) + '\n')
	for row in sorted(prediction_results):
		out_file.write('\t'.join(map(str,model_info[row[0]])) + '\t' + '\t'.join(map(str,row[1])) + '\n')
	print '\n Wrote Results to: ' + out_name
	out_file.close()
