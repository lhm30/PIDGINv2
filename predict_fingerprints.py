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
import operator
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
	for i in inp:
		try:
			outfp.append(calcFingerprints(i))
		except:
			print 'SMILES Parse Error: ' + i
	return outfp
	
#import user query
def importQuery(in_file):
	query = open(in_file).read().splitlines()
	matrix = np.empty((len(query), 2048), dtype=np.uint8)
	smiles_per_core = int(math.ceil(len(query) / N_cores)+1)
	chunked_smiles = [query[x:x+smiles_per_core] for x in xrange(0, len(query), smiles_per_core)]
	pool = Pool(processes=N_cores)  # set up resources
	jobs = pool.imap(arrayFP, chunked_smiles)
	current_end = 0
	for i, result in enumerate(jobs):
		matrix[current_end:current_end+len(result), :] = result
		current_end += len(result)
	pool.close()
	pool.join()
	return matrix[:current_end]

#get info for uniprots
def getUniprotInfo():
	model_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + '/classes_in_model.txt').read().splitlines()]
	return_dict = {l[0] : l[0:5] for l in model_info}
	return return_dict

#get info for diseases
def getDisgenetInfo():
	return_dict1 = dict()
	return_dict2 = dict()
	disease_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + '/DisGeNET_diseases.txt').read().splitlines()]
	for l in disease_info:
		try:
			return_dict1[l[0]].append(l[1])
		except KeyError:
			return_dict1[l[0]] = [l[1]]
		return_dict2[l[1]] = l[2]
	return return_dict1, return_dict2 
	
#get info for biosystems pathways
def getPathwayInfo():
	return_dict1 = dict()
	return_dict2 = dict()
	pathway_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + '/biosystems.txt').read().splitlines()]
	for l in pathway_info:
		try:
			return_dict1[l[0]].append(l[1])
		except KeyError:
			return_dict1[l[0]] = [l[1]]
		return_dict2[l[1]] = l[2:]
	return return_dict1, return_dict2 

#calculate prediction ratio for two sets of predictions
def calcPredictionRatio(preds1,preds2):
	preds1_percentage = float(preds1)/float(len(querymatrix1))
	preds2_percentage = float(preds2)/float(len(querymatrix2))
	if preds1 == 0: return 999.0, round(preds1_percentage,3), round(preds2_percentage,3) 
	if preds1 == 0 and preds2 == 0: return None
	if preds2 == 0: return 0.0, round(preds1_percentage,3), round(preds2_percentage,3)
	return round(preds2_percentage/preds1_percentage,3), round(preds1_percentage,3), round(preds2_percentage,3)

#prediction worker	
def doTargetPrediction(pickled_model_name):
	with open(pickled_model_name, 'rb') as fid:
		clf = cPickle.load(fid)
	uniprot = pickled_model_name.split('/')[-1][:-4]
	preds = clf.predict_proba(querymatrix1)[:,1] > threshold
	return uniprot, list(map(int,preds))

#prediction runner
def performTargetPrediction(models):
	prediction_results = []
	pool = Pool(processes=N_cores)  # set up resources
	jobs = pool.imap_unordered(doTargetPrediction, models)
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(models)))*100 + 1
		sys.stdout.write(' Performing Classification on Query Molecules: %3d%%\r' % percent)
		sys.stdout.flush()
		prediction_results.append(result[0])
		for p in pathway_links[result[0][0]]
	pool.close()
	pool.join()
	prediction_results = np.array(prediction_results).transpose()
	return prediction_results

#main
#set up environment
input_name1, N_cores, threshold  = 'out3.txt', int(10), 0.5
introMessage()
print ' Using ' + str(N_cores) + ' Cores'
models = [modelfile for modelfile in glob.glob('models/*.pkl')][:3]
model_info = getUniprotInfo()
disease_links, disease_info = getDisgenetInfo()
pathway_links, pathway_info = getPathwayInfo()
print ' Total Number of Classes : ' + str(len(models))
print ' Using TPR threshold of : ' + str(threshold)
output_name = input_name1 + '_out_target_fingerprint_' + str(threshold) + '.txt'
out_file = open(output_name, 'w')

#perform target predictions and write to file
querymatrix1 = importQuery(input_name1)
disease_hits = dict()
pathway_hits = dict()
print ' Total Number of Molecules in ' +input_name1+ ' : ' + str(len(querymatrix1))
prediction_results = performTargetPrediction(models)
for row in prediction_results:
	out_file.write('\t'.join(map(str,row)) + '\n')
print '\n Wrote Results to: ' + output_name
out_file.close()

#write disease results to file
processed_diseases, inp1_total, inp2_total = processDiseaseHits()
output_name = input_name1 + '_vs_' + input_name2 + '_out_enriched_diseases_' + str(threshold) + '.txt'
out_file = open(output_name, 'w')
out_file.write('Disease_Name\t'+input_name1+'_Hits\t'+input_name1+'_Precent_Hits\t'+input_name2+'_Hits\t'+input_name2+'Precent_Hits\tRatio\tDisease_Confidence_Score\n')
for disease, ratio in sorted(processed_diseases.items(), key=operator.itemgetter(1)):
	inp1_stats = [disease_hits[disease][0],round(float(disease_hits[disease][0])/float(inp1_total),3)]
	inp2_stats = [disease_hits[disease][1],round(float(disease_hits[disease][1])/float(inp2_total),3)]
	out_file.write(disease + '\t' + '\t'.join(map(str,inp1_stats)) + '\t' + '\t'.join(map(str,inp2_stats)) + '\t' + str(ratio) + '\t' + disease_info[disease] + '\n')
out_file.close()

#write pathway results to file
processed_pathways, inp1_total, inp2_total = processPathwayHits()
output_name = input_name1 + '_vs_' + input_name2 + '_out_enriched_pathways_' + str(threshold) + '.txt'
out_file = open(output_name, 'w')
out_file.write('Pathway_Name\tPathway_Name\tSource\tClass\t'+input_name1+'_Hits\t'+input_name1+'_Precent_Hits\t'+input_name2+'_Hits\t'+input_name2+'Precent_Hits\tRatio\n')
for pathway, ratio in sorted(processed_pathways.items(), key=operator.itemgetter(1)):
	inp1_stats = [pathway_hits[pathway][0],round(float(pathway_hits[pathway][0])/float(inp1_total),3)]
	inp2_stats = [pathway_hits[pathway][1],round(float(pathway_hits[pathway][1])/float(inp2_total),3)]
	out_file.write(pathway + '\t' + '\t'.join(map(str,pathway_info[pathway])) + '\t' + '\t'.join(map(str,inp1_stats)) + '\t' + '\t'.join(map(str,inp2_stats)) + '\t' + str(ratio) + '\n')
out_file.close()