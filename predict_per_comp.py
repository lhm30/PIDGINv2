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
	model_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + sep + 'classes_in_model.txt').read().splitlines()]
	return_dict = {l[0] : l[0:7] for l in model_info}
	return return_dict

#get info for diseases
def getDisgenetInfo():
	return_dict1 = dict()
	return_dict2 = dict()
	disease_file = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + sep + 'DisGeNET_diseases.txt').read().splitlines()]
	for l in disease_file:
		try:
			return_dict1[l[0]].append(l[1])
		except KeyError:
			return_dict1[l[0]] = [l[1]]
		try:
			return_dict2[(l[1],l[0])] = float(l[2])
		except ValueError: pass
	return return_dict1, return_dict2

#get info for biosystems pathways
def getPathwayInfo():
	return_dict1 = dict()
	return_dict2 = dict()
	pathway_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + sep + 'biosystems.txt').read().splitlines()]
	for l in pathway_info:
		try:
			return_dict1[l[0]].append(l[1])
		except KeyError:
			return_dict1[l[0]] = [l[1]]
		return_dict2[l[1]] = l[2:]
	return return_dict1, return_dict2 

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
	preds = map(int,probs > threshold)
	if sum(preds) > 0:
		return mod,preds
	else: return None

#prediction runner
def performTargetPrediction(models):
	total_pw = set()
	total_disease = set()
	prediction_results = dict()
	pool = Pool(processes=N_cores, initializer=initPool, initargs=(querymatrix,threshold,))  # set up resources
	jobs = pool.imap_unordered(doTargetPrediction, models)
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(models)))*100 + 1
		sys.stdout.write(' Performing Classification on Query Molecules: %3d%%\r' % percent)
		sys.stdout.flush()
		if result is not None:
			prediction_results[result[0]] = result[1]
			total_pw = total_pw | set(pathway_links.get(result[0],[]))
			try:
				total_disease = total_disease | set([dis for dis in disease_links[result[0]] if disease_score[(dis,result[0])] > dgn_threshold])
			except KeyError: pass
	pool.close()
	pool.join()
	prediction_matrix = np.array([j for i,j in sorted(prediction_results.items())]).transpose()
	return prediction_results, prediction_matrix, sorted(list(total_pw)), sorted(list(total_disease))
	
#initializer for the pool
def initPool(querymatrix_, threshold_):
	global querymatrix, threshold
	querymatrix = querymatrix_
	threshold = threshold_

#main
if __name__ == '__main__':
	input_name, N_cores,  = sys.argv[1], int(sys.argv[2])
	introMessage()
	print ' Using ' + str(N_cores) + ' Cores'
	try:
		threshold = float(sys.argv[3])
	except ValueError:
		print 'ERROR: Enter a valid float (2DP) for threshold'
		quit()
	try:
		dgn_threshold = float(sys.argv[4])
	except IndexError:
		dgn_threshold = 0
	try:
		desired_organism = sys.argv[5]
	except IndexError:
		desired_organism = None
	model_info = getUniprotInfo()
	models = [modelfile for modelfile in glob.glob(os.path.dirname(os.path.abspath(__file__)) + sep + 'models' + sep + '*.zip')]
	if desired_organism is not None:
		if os.name == 'nt': sep = '\\'
		else: sep = '/'
		models = [mod for mod in models if model_info[mod.split(sep)[-1].split('.')[0]][4] == desired_organism]
	disease_links, disease_score = getDisgenetInfo()
	pathway_links, pathway_info = getPathwayInfo()
	print ' Total Number of Classes : ' + str(len(models))
	print ' Using TPR threshold of : ' + str(threshold)
	print ' Using DisGeNET score threshold of : ' + str(dgn_threshold)
	if desired_organism is not None:
		print ' Predicting for organism : ' + desired_organism
		out_file = open(input_name + '_out_percomp_target_' + str(threshold) + '_' + desired_organism[:3] + '.txt', 'w')
		out_file2 = open(input_name + '_out_percomp_pathway_' + str(threshold) + '_' + desired_organism[:3] + '.txt', 'w')
		out_file3 = open(input_name + '_out_percomp_disease_' + str(threshold) + '_' + str(dgn_threshold) + '_' + desired_organism[:3] + '.txt', 'w')
	else:
		out_file = open(input_name + '_out_percomp_target_' + str(threshold) + '.txt', 'w')
		out_file2 = open(input_name + '_out_percomp_pathway_' + str(threshold) + '.txt', 'w')
		out_file3 = open(input_name + '_out_percomp_disease_' + str(threshold) + '_' + str(dgn_threshold) + '.txt', 'w')

	#perform target predictions and tp fingerprints to file 
	querymatrix,smiles = importQuery(input_name)
	print ' Total Number of Query Molecules : ' + str(len(querymatrix))
	prediction_results,prediction_matrix,sorted_pws,sorted_diseases = performTargetPrediction(models)
	sorted_targets = sorted(prediction_results.keys())
	out_file.write('Compound\t' + '\t'.join(map(str,[i for i in sorted(prediction_results.keys())])) + '\n')
	out_file.write('-\t' + '\t'.join(map(str,[model_info[i][4] for i in sorted(prediction_results.keys())])) + '\n')
	out_file2.write('Compound\t' + '\t'.join(map(str,sorted_pws)) + '\n')
	out_file2.write('Compound\t' + '\t'.join(map(str,[pathway_info[i][0] for i in sorted_pws])) + '\n')
	out_file3.write('Compound\t' + '\t'.join(map(str,sorted_diseases)) + '\n')
	for i, row in enumerate(prediction_matrix):
		#write target prediction fp
		out_file.write(smiles[i] + '\t' + '\t'.join(map(str,row)) + '\n')
		pred_targets = [sorted_targets[i2] for i2, pred in enumerate(row) if pred == 1]
		#write pathway fp
		comp_pws = Counter(sum([[pw for pw in pathway_links.get(t,[])] for t in pred_targets],[]))
		pw_line = [smiles[i]] + [comp_pws.get(pwid,0) for pwid in sorted_pws]
		out_file2.write('\t'.join(map(str,pw_line)) + '\n')
		#write disease fp
		comp_disease = Counter(sum([[dis for dis in disease_links.get(t,[]) if disease_score[(dis,t)] > dgn_threshold] for t in pred_targets],[]))
		disease_line = [smiles[i]] + [comp_disease.get(dis,0) for dis in sorted_diseases]
		out_file3.write('\t'.join(map(str,disease_line)) + '\n')
	print '\n Wrote Results'
	out_file.close()
	out_file2.close()
	out_file3.close()