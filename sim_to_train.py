#Author : Lewis Mervin lhm30@cam.ac.uk
#Supervisor : Dr. A. Bender
#All rights reserved 2016
#Protein Target Prediction Tool trained on SARs from PubChem (Mined 21/06/16) and ChEMBL21
#Molecular Descriptors : 2048bit Morgan Binary Fingerprints (Rdkit) - ECFP4
#Dependencies : rdkit, sklearn, numpy

#libraries
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import cPickle
import zipfile
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
	return fp

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
	smiles_per_core = int(math.ceil(len(query) / N_cores)+1)
	chunked_smiles = [query[x:x+smiles_per_core] for x in xrange(0, len(query), smiles_per_core)]
	pool = Pool(processes=N_cores)  # set up resources
	jobs = pool.imap(arrayFP, chunked_smiles)
	processed_fp = []
	processed_smi = []
	for i, result in enumerate(jobs):
		processed_fp += result[0]
		processed_smi += result[1]
	pool.close()
	pool.join()
	#if IDs aren't present, use SMILES as IDs
	if not ids:
		ids = processed_smi
	return processed_fp, processed_smi, ids

#get info for uniprots
def getUniprotInfo():
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	model_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + sep + 'classes_in_model.txt').read().splitlines()]
	return_dict = {l[0] : l[0:7] for l in model_info}
	return return_dict

#sim worker
def doSimSearch(model_name):
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	mod = model_name.split(sep)[-1].split('.')[0]
	try:
		with zipfile.ZipFile(os.path.dirname(os.path.abspath(__file__)) + sep + 'actives' + sep + mod + '.smi.zip', 'r') as zfile:
			comps = [i.split('\t') for i in zfile.open(mod + '.smi', 'r').read().splitlines()]
	except IOError: return
	comps2 = []
	afp = []
	for comp in comps:
		try:
			afp.append(calcFingerprints(comp[1]))
			comps2.append(comp)
		except: pass
	ret = []
	for i,fp in enumerate(querymatrix):
		sims = DataStructs.BulkTanimotoSimilarity(fp,afp)
		idx = sims.index(max(sims))
		ret.append([sims[idx], mod] + comps2[idx] + [smiles[i]])
	return ret

#prediction runner
def performSimSearch(models):
	sims_results = []
	pool = Pool(processes=N_cores, initializer=initPool, initargs=(querymatrix,smiles))  # set up resources
	jobs = pool.imap_unordered(doSimSearch, models)
	out_file2.write('Uniprot\tPref_Name\tGene ID\tTarget_Class\tOrganism\tPDB_ID\tDisGeNET_Diseases_0.06\t' + '\t'.join(map(str,ids)) + '\n')
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(models)))*100 + 1
		sys.stdout.write(' Calculating Sims for Query Molecules: %3d%%\r' % percent)
		sys.stdout.flush()
		if result is not None:
			sims_results += result
			out_file2.write('\t'.join(map(str,model_info[result[0][1]])))
			for sim in result:
				out_file2.write('\t' + str(round(sim[0],3)))
			out_file2.write('\n')
	pool.close()
	pool.join()
	return sims_results

#initializer for the pool
def initPool(querymatrix_,smiles_):
	global querymatrix, smiles
	querymatrix = querymatrix_
	smiles = smiles_

#main
if __name__ == '__main__':
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	input_name = sys.argv[1]
	N_cores = int(sys.argv[2])
	try:
		desired_organism = sys.argv[3]
	except IndexError:
		desired_organism = None
	introMessage()
	print ' Calculating Near-Neighbors for ' + input_name
	print ' Using ' + str(N_cores) + ' Cores'
	models = [modelfile for modelfile in glob.glob(os.path.dirname(os.path.abspath(__file__)) + sep + 'models' + sep + '*.zip')]
	model_info = getUniprotInfo()

	if desired_organism is not None:
		models = [mod for mod in models if model_info[mod.split(sep)[-1].split('.')[0]][4] == desired_organism]
		output_name = input_name + '_out_similarity_details' + desired_organism[:3] + '.txt'
		output_name2 = input_name + '_out_similarity_matrix' + desired_organism[:3] + '.txt'
		print ' Predicting for organism : ' + desired_organism
	else:
		output_name = input_name + '_out_similarity_details.txt'
		output_name2 = input_name + '_out_similarity_matrix.txt'
	print ' Total Number of Classes : ' + str(len(models))
	output_name = input_name + '_out_similarity_details.txt'
	output_name2 = input_name + '_out_similarity_matrix.txt'
	out_file = open(output_name, 'w')
	out_file2 = open(output_name2, 'w')
	querymatrix,smiles,ids = importQuery(input_name)
	print ' Total Number of Query Molecules : ' + str(len(querymatrix))
	sims_results = performSimSearch(models)
	out_file.write('Uniprot\tPref_Name\tGene ID\tTarget_Class\tOrganism\tNear_Neighbor_ChEMBLID\tNear_Neighbor_Smiles\tNear_Neighbor_Bioactive_organism\tNear_Neighbor_conf_score\tNN_activity\tNN_Units\tInput_Compound\tSimilarity\n')
	for row in sorted(sims_results,reverse=True):
		out_file.write('\t'.join(map(str,model_info[row[1]][:5])) + '\t' + '\t'.join(map(str,row[2:])) + '\t' + str(row[0]) + '\n')
	print '\n Wrote Results to: ' + output_name
	print ' Wrote Results to: ' + output_name2
	out_file.close()
