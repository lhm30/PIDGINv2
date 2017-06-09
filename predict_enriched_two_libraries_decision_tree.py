#Author : Lewis Mervin lhm30@cam.ac.uk
#Supervisor : Dr. A. Bender
#All rights reserved 2016
#Protein Target Prediction Tool trained on SARs from PubChem (Mined 21/06/16) and ChEMBL21
#Molecular Descriptors : 2048bit Morgan Binary Fingerprints (Rdkit) - ECFP4
#Dependencies : rdkit, sklearn, numpy

#libraries
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn import tree
from sklearn.externals.six import StringIO
import cPickle
import zipfile
import glob
import os
import sys
import math
import numpy as np
import scipy.stats as stats
from multiprocessing import Pool
import multiprocessing
import operator
import pydot
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
	model_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + sep + 'classes_in_model.txt').read().splitlines()]
	return_dict = {l[0] : l[0:8] for l in model_info}
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

#calculate prediction ratio for two sets of predictions
def calcPredictionRatio(preds1,preds2):
	preds1_percentage = float(preds1)/float(len(querymatrix1))
	preds2_percentage = float(preds2)/float(len(querymatrix2))
	if preds1 == 0 and preds2 == 0: return None
	if preds1 == 0: return 999.0, round(preds1_percentage,3), round(preds2_percentage,3) 
	if preds2 == 0: return 0.0, round(preds1_percentage,3), round(preds2_percentage,3)
	return round(preds2_percentage/preds1_percentage,3), round(preds1_percentage,3), round(preds2_percentage,3)

#unzip a pkl model
def open_Model(mod):
	with zipfile.ZipFile(os.path.dirname(os.path.abspath(__file__)) + sep + 'models' + sep + mod + '.pkl.zip', 'r') as zfile:
		with zfile.open(mod + '.pkl', 'r') as fid:
			clf = cPickle.load(fid)
	return clf

#prediction worker to predict targets and calculate Fisher's test in parallel
def doTargetPrediction(pickled_model_name):
	mod = pickled_model_name.split(sep)[-1].split('.')[0]
	clf = open_Model(mod)
	probs1 = map(int, clf.predict_proba(querymatrix1)[:,1] > threshold)
	preds1 = sum(probs1)
	probs2 = map(int, clf.predict_proba(querymatrix2)[:,1] > threshold)
	preds2 = sum(probs2)
	oddsratio, pvalue = stats.fisher_exact([[preds2,len(querymatrix2)-preds2],[preds1,len(querymatrix1)-preds1]])
	try:
		ratio, preds1_percentage, preds2_percentage = calcPredictionRatio(preds1,preds2)
		return ratio, mod, preds1, preds1_percentage, preds2, preds2_percentage, oddsratio, pvalue, probs1 + probs2
	except TypeError: return None

#prediction runner
def performTargetPrediction(models):
	prediction_results = []
	decision_tree_matrix = []
	decision_tree_node_label = []
	pool = Pool(processes=N_cores, initializer=initPool, initargs=(querymatrix1,querymatrix2,threshold))  # set up resources
	jobs = pool.imap_unordered(doTargetPrediction, models)
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(models)))*100 + 1
		sys.stdout.write(' Performing Classification on Query Molecules: %3d%%\r' % percent)
		sys.stdout.flush()
		if result is not None: 
			prediction_results.append(result[:8])
			updateHits(disease_links,disease_hits,result[1],result[2],result[4])
			updateHits(pathway_links,pathway_hits,result[1],result[2],result[4])
			decision_tree_matrix.append(result[8])
			decision_tree_node_label.append(model_info[result[1]][2])
	pool.close()
	pool.join()
	decision_tree_matrix = np.array(decision_tree_matrix,dtype=np.uint8).transpose()
	return prediction_results, decision_tree_matrix, decision_tree_node_label
	
#update counts for each pathway/disease that is hit by predictions	
def updateHits(links,hits,uniprot,hit1,hit2):
	try:
		for idx in links[uniprot]:
			#try checks if pw or dnet
			try:
				if disease_score[(idx,uniprot)] < dgn_threshold: continue
			except KeyError: pass
			try:
				hits[idx] = hits[idx] + np.array([hit1,hit2])
			except KeyError:
				hits[idx] = np.array([hit1,hit2])
	except KeyError: return
	return
	
#worker for the processHits to calculate the prediction ratio, Chi-square test in parallel
def doHitProcess(inp):
	idx, hits, n_f1_hits, n_f2_hits = inp
	if hits[0] == 0 and hits[1] == 0: return
	if hits[0] == 0: return idx, 999.0, 0, 0, hits[1], float(hits[1])/float(n_f2_hits), 'NA', 'NA'
	if hits[1] == 0: return idx, 0.0, hits[0], float(hits[0])/float(n_f1_hits), 0, 0, 'NA', 'NA'
	h1_p = float(hits[0])/float(n_f1_hits)
	h2_p = float(hits[1])/float(n_f2_hits)
	chi, pvalue, _, _ = stats.chi2_contingency([[hits[1],n_f2_hits-hits[1]],[hits[0],n_f1_hits-hits[0]]])
	return idx, round(h2_p/h1_p,3), hits[0], h1_p, hits[1], h2_p, chi, pvalue

#calculate the enrichment ratio between predictions
def processHits(inp_dict):
	out_dict = dict()
	total_hits = np.array(inp_dict.values()).sum(axis=0)
	if total_hits.shape is (): return out_dict, 0, 0
	n_f1_hits = total_hits[0]
	n_f2_hits = total_hits[1]
	tasks = [[idx,hits,n_f1_hits,n_f2_hits] for idx, hits in inp_dict.iteritems()]
	pool = Pool(processes=N_cores)  # set up resources
	jobs = pool.imap_unordered(doHitProcess, tasks)
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(tasks)))*100 + 1
		sys.stdout.write(" Calculating Fisher's test: %3d%%\r" % percent)
		sys.stdout.flush()
		if result is None: continue
		out_dict[result[0]] = result[1:]
	return out_dict, n_f1_hits, n_f2_hits

#train decision tree on predictions and output graph for pdf
def createTree(matrix,label):	
	vector = [1] * len(querymatrix1) + [0] * len(querymatrix2)
	ratio = float(len(vector)-sum(vector))/float(sum(vector))
	sw = np.array([ratio if i == 1 else 1 for i in vector])
	pc_10 = int(len(querymatrix1)*0.01)
	clf = tree.DecisionTreeClassifier(min_samples_split=min_sampsplit,min_samples_leaf=min_leafsplit,max_depth=max_d)
	clf.fit(matrix,vector)
	dot_data = StringIO()  
	tree.export_graphviz(clf, out_file=dot_data,  
							feature_names=label,  
							class_names=['File2','File1'],  
							filled=True, rounded=True,  
							special_characters=True,
							proportion=False,
							impurity=True)
	out_tree = dot_data.getvalue()
	out_tree = out_tree.replace('True','Inactive').replace('False','Active').replace(' &le; 0.5', '')
	graph = pydot.graph_from_dot_data(str(out_tree))
	try:
		graph.write_jpg(output_name_tree)
	except AttributeError:
		graph = pydot.graph_from_dot_data(str(out_tree))[0]
		graph.write_jpg(output_name_tree)
	return

#initializer for the pool
def initPool(querymatrix1_, querymatrix2_, threshold_):
	global querymatrix1, querymatrix2, threshold
	querymatrix1 = querymatrix1_
	querymatrix2 = querymatrix2_
	threshold = threshold_

#main
#set up environment
if __name__ == '__main__':
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	input_name1, input_name2, N_cores  = sys.argv[1], sys.argv[2], int(sys.argv[3])
	introMessage()
	print ' Using ' + str(N_cores) + ' Cores'
	try:
		threshold = float(sys.argv[4])
	except ValueError:
		print 'ERROR: Enter a valid float (2DP) for threshold'
		quit()
	try:
		dgn_threshold = float(sys.argv[5])
	except IndexError:
		dgn_threshold = 0

	min_sampsplit = int(sys.argv[6])
	min_leafsplit = int(sys.argv[7])
	max_d = int(sys.argv[8])
	try:
		desired_organism = sys.argv[9]
	except IndexError:
		desired_organism = None
	model_info = getUniprotInfo()
	models = [modelfile for modelfile in glob.glob(os.path.dirname(os.path.abspath(__file__)) + sep + 'models' + sep + '*.zip')]
	if desired_organism is not None:
		models = [mod for mod in models if model_info[mod.split(sep)[-1].split('.')[0]][4] == desired_organism]
	disease_links, disease_score = getDisgenetInfo()
	pathway_links, pathway_info = getPathwayInfo()
	print ' Total Number of Classes : ' + str(len(models))
	print ' Using TPR threshold of : ' + str(threshold)
	print ' Using DisGeNET score threshold of : ' + str(dgn_threshold)
	if desired_organism is not None:
		print ' Predicting for organism : ' + desired_organism
		output_name = input_name1 + '_vs_' + input_name2 + '_out_enriched_targets_' + str(threshold) + '_' + desired_organism[:3] + '.txt'
		output_name_tree = input_name1 + '_vs_' + input_name2 + '_decision_tree_' + str(threshold) + '_' + desired_organism[:3] + '.jpg'
		output_name2 = input_name1 + '_vs_' + input_name2 + '_out_enriched_diseases_' + str(threshold) + '_' + str(dgn_threshold) + '_' + desired_organism[:3] + '.txt'
		output_name3 = input_name1 + '_vs_' + input_name2 + '_out_enriched_pathways_' + str(threshold) + '_' + desired_organism[:3] + '.txt'
	else:
		output_name = input_name1 + '_vs_' + input_name2 + '_out_enriched_targets_' + str(threshold) + '.txt'
		output_name_tree = input_name1 + '_vs_' + input_name2 + '_decision_tree_' + str(threshold) + '.jpg'
		output_name2 = input_name1 + '_vs_' + input_name2 + '_out_enriched_diseases_' + str(threshold) + '_' + str(dgn_threshold) + '.txt'
		output_name3 = input_name1 + '_vs_' + input_name2 + '_out_enriched_pathways_' + str(threshold) + '.txt'
	print ' Using max sample split, max leaves and max depth of : ' + ', '.join(map(str,[min_sampsplit,min_leafsplit,max_d]))

	#perform target predictions and write to file
	querymatrix1 = importQuery(input_name1)
	querymatrix2 = importQuery(input_name2)
	disease_hits, pathway_hits = dict(), dict()
	print ' Total Number of Molecules in ' +input_name1+ ' : ' + str(len(querymatrix1))
	print ' Total Number of Molecules in ' +input_name2+ ' : ' + str(len(querymatrix2))
	prediction_results, decision_tree_matrix, decision_tree_node_label = performTargetPrediction(models)
	out_file = open(output_name, 'w')
	out_file.write('Uniprot\tPref_Name\tGene ID\tTarget_Class\tOrganism\tPDB_ID\tDisGeNET_Diseases_0.06\tChEMBL_First_Published\t'+input_name1+'_Hits\t'+input_name1+'_Precent_Hits\t'+input_name2+'_Hits\t'+input_name2+'_Precent_Hits\tOdds_Ratio\tFishers_Test_pvalue\tPrediction_Ratio\n')
	for row in sorted(prediction_results):
		out_file.write('\t'.join(map(str,model_info[row[1]])) + '\t' + '\t'.join(map(str, row[2:])) + '\t' + str(row[0]) + '\n')
	print '\n Wrote Results to: ' + output_name
	out_file.close()

	#perform decision tree function and write to file
	createTree(decision_tree_matrix,decision_tree_node_label)
	print 'Wrote Results to: ' + output_name_tree

	#write disease results to file
	processed_diseases, inp1_total, inp2_total = processHits(disease_hits)
	out_file = open(output_name2, 'w')
	out_file.write('Disease_Name\t'+input_name1+'_Hits\t'+input_name1+'_Precent_Hits\t'+input_name2+'_Hits\t'+input_name2+'Precent_Hits\tchi2_test_statistic\tchi2_pvalue\tPrediction_Ratio\n')
	for disease, ratio in sorted(processed_diseases.items(), key=operator.itemgetter(1)):
		out_file.write(disease + '\t' + '\t'.join(map(str,ratio[1:])) + '\t' + str(ratio[0]) + '\n')
	print '\n Wrote Results to: ' + output_name2
	out_file.close()

	#write pathway results to file
	processed_pathways, inp1_total, inp2_total = processHits(pathway_hits)
	out_file = open(output_name3, 'w')
	out_file.write('Pathway_Name\tPathway_Name\tSource\tClass\t'+input_name1+'_Hits\t'+input_name1+'_Precent_Hits\t'+input_name2+'_Hits\t'+input_name2+'Precent_Hits\tchi2_test_statistic\tchi2_pvalue\tPrediction_Ratio\n')
	for pathway, ratio in sorted(processed_pathways.items(), key=operator.itemgetter(1)):
		out_file.write(pathway + '\t' + '\t'.join(map(str,pathway_info[pathway])) + '\t' + '\t'.join(map(str,ratio[1:])) + '\t' + str(ratio[0]) + '\n')
	print '\n Wrote Results to: ' + output_name3
	out_file.close()