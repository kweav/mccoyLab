#!/usr/bin/env python3

import sys
import numpy as np
from sklearn import linear_model, metrics
from scipy import stats
from itertools import chain, combinations
import os
import argparse as ap
import logging
import datetime

'''add a chr_sizes file to the parser and the initialization rather than the dictionary DONE'''
'''let's do one with the correlation cutoff and one without'''
'''--> therefore, add a store_true for correlation cutoff to the parser'''
'''--> Also consider randomly generating the betas; add a lim argument for this to the parser'''
'''--> Also consider naively setting -1 +1 or something similar based on knowledge of promoting or repressing'''
'''need a generator for the powerset'''
'''add a store_true log_props value to the parser'''
'''add a desired base for log_props to the parser'''
'''add a desired add value for log_props to the parser'''
'''add a store_true for filter all 0 cres to the parser'''
'''add a store_true for cutting state 0 out of props to the parser'''
'''add a store_true for excluding CREs within the promoter region of a gene'''
'''add a desired length of region when excluding CREs within the promoter region of a gene'''


'''How do I want to balance the gene by gene pairing with the finding of beta values and the final total MSE calculation?'''

def main():
	parser = generate_parser()
	args = parser.parse_args()
	setup_logging(args.logfile)
	setup_file_locs(args, args.where_run, args.other_path)
	check_pc(args, args.pc_only)
	setup_threads(args.threads)
	model = regress_gene_cre_random(args.train_cre_state_file, args.train_tss_state_file, args.train_exp_file,
									args.test_cre_state_file, args.test_tss_state_file, args.test_exp_file,
									args.output_pair_file, args.output_beta_file,
									args.chr_sizes,
									args.dist_gamma,
									args.use_corr, args.abs_bool,
									abc_dist=args.abc_dist, element_dist=args.element_dist, abc_thresh=args.abc_thresh, exclude_wp=args.exclude_wp, exclude_wp_dist=args.exclude_wp_dist, dn_cross_bound=args.dn_cross_bound,
									log_props=args.log_props, log_props_base=args.log_props_base, log_props_add_val=args.log_props_add_val,
									filter0=args.filter0, cutState0=args.cutState0, zero_val=args.zero_val,
									iter_val=args.iter_val, seed=args.seed,
									run_Gibbs=args.run_Gibbs, burn_in=args.burn_in, init_sigma_sqr=args.init_sigma_sqr, vary_sigma_sqr=args.vary_sigma_sqr)

	if args.chroms == "all":
		chrom_list = ['chr{}'.format(i) for i in np.hstack((np.arange(1,20), 'X'))]
	else:
		chrom_list = args.chroms

	for chrom in chrom_list:
		model.subset_based_on_chrom(chrom)
		model.subset_based_ongroup(args.exp_type, args.m_thresh, args.s_thresh)
		model.set_initial_betas(args.init_beta_method, args.init_beta_rlim, args.init_beta_file, args.init_beta_PrRe_file)
		model.find_initial_weighted_sum()
		model.drive_random_iter(3)

def setup_logging(logfile_oi):
	logfile=logfile_oi
	logging.basicConfig(filename=logfile, level=logging.INFO, filemode='w', format='%(name)s - %(levelname)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
	logging.info("command used to run analysis:\n {}".format(' '.join(sys.argv)))

def setup_file_locs(args, where_run, other_path):
	argumentToAdd = {'mine' : '/Users/kateweaver/taylorLab/VISION/regress_dist_gibbs/inputs/',
					 'marcc' : 'NA',
					 'comp' : '/project/vision/Target_Genes/npz_inputs/',
					 'other' : other_path}

	args.train_cre_state_file = argumentToAdd[where_run] + args.train_cre_state_file
	args.test_cre_state_file = argumentToAdd[where_run] + args.test_cre_state_file
	args.train_tss_state_file = argumentToAdd[where_run] + args.train_tss_state_file
	args.test_tss_state_file = argumentToAdd[where_run] + args.test_tss_state_file
	args.train_exp_file = argumentToAdd[where_run] + args.train_exp_file
	args.test_exp_file = argumentToAdd[where_run] + args.test_exp_file

	args.chr_sizes = argumentToAdd[where_run] + args.chr_sizes
	args.init_beta_file = argumentToAdd[where_run] + args.init_beta_file

def check_pc(args, pc_only):
	if pc_only:
		args.train_exp_file = args.train_exp_file.replace('.npz', '_pc.npz')
		args.test_exp_file = args.test_exp_file.replace('.npz', '_pc.npz')
		args.train_tss_state_file = args.train_tss_state_file.replace('.npz', '_pc.npz')
		args.test_tss_state_file = args.test_tss_state_file.replace('.npz', '_pc.npz')

def setup_threads(threads):
	os.environ["OMP_NUM_THREADS"] = threads
	os.environ["OPENBLAS_NUM_THREADS"] = threads
	os.environ["MKL_NUM_THREADS"] = threads
	os.environ["VECLIB_MAXIMUM_THREADS"] = threads
	os.environ["NUMEXPR_NUM_THREADS"] = threads

def generate_parser():
	parser = ap.ArgumentParser(description = 'VISION regression of state and distance to assign CREs to genes based on ability to predict gene expression, using an IC approach to penalize larger feature sets and randomly selecting from a powerset for the CRE feature set used each iteration')
	'''arguments related to running logistics'''
	parser.add_argument('--where', action='store', dest='where_run', type=str, default='comp', help='{comp, mine, marcc, other}; adds path to files; if "other" used, make sure to provide the other path in the --otherpath argument')
	parser.add_argument('--otherpath', action='store', dest='other_path', type=str, default='NA', help='use to give other path if not one of the 3 previously specified. MUST be the same for ALL input files')
	parser.add_argument('--threads', action='store', dest='threads', type=str, default="1")
	parser.add_argument('--logfile', action='store', dest='logfile', type=str, default='base_random_linear.log', help='name of logfile')
	parser.add_argument('--seed_value', action='store', dest='seed', type=int, default=42, help='random seed to be used for reproducibility')

	'''arguments related to input data files'''
	parser.add_argument('-c', '--train_cre_state_file', action='store', dest='train_cre_state_file', type=str, default='train_cre_state_prop.npz', help='training file; default stores the file name, though a different one can be provided, path will be added based on --where argument')
	parser.add_argument('-t', '--train_tss_state_file', action='store', dest='train_tss_state_file', type=str, default='trainTSS_window_state_prop.npz', help='training file; default stores the file name, though a different one can be provided, path will be added based on --where argument')
	parser.add_argument('-e', '--train_exp_file', action='store', dest='train_exp_file', type=str, default='trainTPM.npz', help='training file; default stores the file name, though a different one can be provided, path will be added based on --where argument')
	parser.add_argument('--test_cre_state_file', action='store', dest='test_cre_state_file', type=str, default='test_cre_state_prop.npz', help='testing file; default stores the file name, though a different one can be provided, path will be added based on --where argument')
	parser.add_argument('--test_tss_state_file', action='store', dest='test_tss_state_file', type=str, default='testTSS_window_state_prop.npz', help='testing file; default stores the file name, though a different one can be provided, path will be added based on --where argument')
	parser.add_argument('--test_exp_file', action='store', dest='test_exp_file', type=str, default='testTPM.npz', help='testing file; default stores the file name, though a different one can be provided, path will be added based on --where argument')
	parser.add_argument('--init_beta_file', action='store', dest='init_beta_file', type=str, default='NA', help='the file containing the initial beta parameters saved as a numpy object readable with np.load, only necessary if --init_beta_method=file')
	parser.add_argument('--init_beta_PrRe_file', action='store', dest='init_beta_PrRe_file', type=str, default='NA', help='the file containing the initial beta parameters saved as numpy object readable with np.load, only necessary if --init_beta_method=PrRe')
	parser.add_argument('--chr_sizes_file', action='store', dest='chr_sizes', type=str, default='chromosomes.txt', help='this tab-delimited file contains the sizes of the chromosomes and chromosomes stored ')
	parser.add_argument('--pc_only', action='store_true', dest='pc_only', help='provide if you want to assess protein coding genes only. WARNING: Will overwrite the input file names, replacing .npz with _pc.npz')

	'''arguments related to what to run'''
	parser.add_argument('--chroms', action='store', nargs = '+', dest='chroms', type=str, default='all')
	parser.add_argument('--iter_val', action='store', dest='iter_val', type=int, default=10**5, help='the number of iterations to run')
	parser.add_argument('--burn_in', action='store', dest='burn_in', type=int, default=500, help='if using Gibbs sampling, the number of iterations to use for the burn-in period, but not analysis')
	parser.add_argument('-g', '--group', action='store', dest='exp_type', type=int, default=0, help='{0, 1, 2, 3, 4} Group subselection. 0 is for all and is default')
	parser.add_argument('-m', '--mean_thresh', action='store', dest='m_thresh', type=int, default=-4)
	parser.add_argument('-s', '--stdev_thresh', action='store', dest='s_thresh', type=int, default=2)


	'''arguments related to variables or distance thresholds in the model'''
	parser.add_argument('--init_beta_method', action='store', dest='init_beta_method', type=str, default='regress', help='{regress, random, PrRe, file} use of regress will use regression of TSS and exp to set, random will randomly set, limiting generation using --input_beta_rlim, and file will require a file input using the --init_beta_file argument, and PrRe will require a file input using the --init_beta_PrRe_file argument')
	parser.add_argument('--init_beta_rlim', action='store', dest='init_beta_rlim', type=float, default=40, help='absolute value of upper limit to use when randomly generating initial betas')
	parser.add_argument('--dist_gamma', action='store', dest='dist_gamma', type=float, default=0.5, help='gamma value for distance correction')
	parser.add_argument('--abc_thresh', action='store', dest='abc_thresh', type=float, default=0.2, help='the threshold to use for cutoff for CRE - gene pairing when using the ABC relationship as a measurement')
	parser.add_argument('--abc_dist', action='store', dest='abc_dist', type=int, default=5000000, help='the number of base pairs to use in a two sided window for the abc measurement')
	parser.add_argument('--cre_dist', action='store', dest='element_dist', type=int, default=1000000, help='the number of base pairs to use in a two sided window for finding CREs within a certain distance of a promoter')
	parser.add_argument('--dist_zero_val', action='store', dest='zero_val', type=float, default=1, help='value to in ABC calculations when distance is 0')
	parser.add_argument('--init_sigma_sqr', action='store', dest='init_sigma_sqr', type=float, default=1.0, help='initial sigma square value, if using Gibbs sampling')

	'''store_true or store_false arguments that affect behavior'''
	parser.add_argument('--use_corr', action='store_true', dest='use_corr', help='provide this flag if and only if you want to use a correlation threshold before constructing the powerset')
	parser.add_argument('--abs_bool', action='store_false', dest='abs_bool', help='include this argument if and only if you want to have signed distance calculations rather than using the absolute value of all distance calculations')
	parser.add_argument('--exclude_wp', action='store_true', dest='exclude_wp', help='include this argument if and only if you want to exclude CREs within a certain distance from the promoter region. Tune the region length using --exclude_wp_dist argument')
	parser.add_argument('--dn_cross_bound', action='store_false', dest='dn_cross_bound', help='include this argument if and only if you do not want to include cre elements that cross the within promoter boundaries')
	parser.add_argument('--log_props', action='store_true', dest='log_props', help='include this argument if and only if you want to log the proportions (of each CRE/promoter region of the inpute for each epigenetic state). Use --log_props_base and --log_props_add_val to tune the behavior')
	parser.add_argument('--filter0', action='store_true', dest='filter0', help='include this argument if and only if you want to filter out CREs that are completely quiescent state (epigenetic IDEAS state 0)')
	parser.add_argument('--cutState0', action='store_true', dest='cutState0', help='include this argument if and only if you want to exclude the quiescent state (epigenetic state 0) from consideration in the model. Recommend that you also include --filter0 if you use this argument')
	parser.add_argument('--runGibbs', action='store_true', dest='run_Gibbs', help='include this argument if and only if you want to run Gibbs sampling')
	parser.add_argument('--vary_sigma_sqr', action='store_true', dest='vary_sigma_sqr', help='include this argument if and only if you are performing Gibbs sampling and want to update sigma square from the init_sigma_sqr value')

	'''arguments directly related to store_true arguments'''
	parser.add_argument('--exclude_wp_dist', action='store', dest='exclude_wp_dist', type=int, default=1000, help='the number of base pairs to use in a two sided window when excluding CREs within a certain distance from the promoter region')
	parser.add_argument('--log_props_base', action='store', dest='log_props_base', type=int, default=2, help='the base to use if you log the proportions of the input data')
	parser.add_argument('--log_props_add_val', action='store', dest='log_props_add_val', type=float, default=0.001, help='the pseudocount value to add to all proportions (because of 0 values) of the input data before taking the log')


	'''arguments related to the output files'''
	parser.add_argument('--ouput_pair_file', action='store', dest='output_pair_file', type=str, default='base_random_linear_pairing_{}.txt', help='name for output file storing predicted pairings for each chromosome. Name should include {} for later addition of chromosome')
	parser.add_argument('--output_beta_file', action='store', dest='output_beta_file', type=str, default='base_random_linear_betas_{}.txt', help='name for output file storing beta values for each chromosome. Name should include {} for later addition of chromosome')

	return parser

class regress_gene_cre_random():
	def __init__(self, train_cre_state_file, train_tss_state_file, train_exp_file,
					   test_cre_state_file, test_tss_state_file, test_exp_file,
					   output_pair_file, output_beta_file,
					   chr_sizes,
					   dist_gamma,
					   use_corr, abs_bool,
					   abc_dist=5000000, element_dist=1000000, abc_thresh=0.2, exclude_wp=True, exclude_wp_dist=1000, dn_cross_bound=True,
					   log_props=False, log_props_base=2, log_props_add_val = 0.001,
					   filter0 = True, cutState0=False, zero_val=1,
					   iter_val=400, seed=42,
					   run_Gibbs=False, burn_in=500, init_sigma_sqr=1, vary_sigma_sqr=False):

		self.dist_gamma = dist_gamma
		self.use_corr = use_corr
		self.abs_bool = abs_bool
		self.abc_dist = abc_dist
		self.element_dist = element_dist
		self.abc_thresh = abc_thresh
		self.exclude_wp = exclude_wp
		self.exclude_wp_dist = exclude_wp_dist
		self.dn_cross_bound = dn_cross_bound
		self.log_props = log_props
		self.log_props_base = log_props_base
		self.log_props_add_val = log_props_add_val
		self.filter0 = filter0
		self.cutState0 = cutState0
		self.zero_val = zero_val
		self.iter_val = iter_val
		self.seed = seed
		self.run_Gibbs = run_Gibbs
		self.burn_in = burn_in
		self.init_sigma_sqr = init_sigma_sqr
		self.vary_sigma_sqr = vary_sigma_sqr

		self.output_pair_file = output_pair_file
		self.output_beta_file = output_beta_file

		self.exp_values_all, self.cellIndex, self.cell_to_index, self.TSS_chr_all, self.TSSs_all = self.load_expression(train_exp_file)
		self.cellN = self.cellIndex.shape[0]
		print(self.cellN)
		self.TSS_window_props_all, self.TSS_window_chr_all = self.load_TSS_window_states(train_tss_state_file)
		self.cre_props_all, self.cre_chr_all, self.cre_coords_all = self.load_cre_states(train_cre_state_file)
		self.stateN = self.cre_props_all.shape[2] #Note this value will be indicitive of whether state 0 was cut or not
		print(self.stateN)
		self.chrSizes = {x: int(y) for x,y in (line.strip('\r\n').split() for line in open(chr_sizes))}

	def find_valid_cellTypes(self, cellIndexOI):
		valid = np.array([np.where(cellIndexOI == x) for x in self.cellIndex if np.isin(x, cellIndexOI)]).reshape(-1)
		return valid

	def load_expression(self, exp_file):
		npzfile = np.load(exp_file, allow_pickle=True)
		exp_values = npzfile['exp'].astype(np.float64)
		exp_cellIndex = npzfile['cellIndex'] #index to celltype
		exp_cell_to_index = {exp_cellIndex[i]: i for i in range(exp_cellIndex.shape[0])} #celltype to index
		TSS_info = npzfile['TSS']
		TSS_chr = TSS_info[:,0]
		TSSs = np.around(TSS_info[:,1].astype(np.float64)).astype(np.int32)
		return exp_values, exp_cellIndex, exp_cell_to_index, TSS_chr, TSSs

	def load_TSS_window_states(self, tss_state_file):
		npzfile = np.load(tss_state_file, allow_pickle=True)
		valid = self.find_valid_cellTypes(npzfile['cellIndex'])
		TSS_window_props_valid = npzfile['props'].astype(np.float64)[:,valid,:]
		if self.log_props:
			TSS_window_props_valid = np.log(TSS_window_props_valid + self.log_props_add_val)/np.log(self.log_props_base)
		if self.cutState0:
			TSS_window_props_valid = TSS_window_props_valid[:,:,1:]
		TSS_window_index = npzfile['ccREIndex']
		TSS_window_chr = TSS_window_index[:,0]
		return TSS_window_props_valid, TSS_window_chr

	def load_cre_states(self, cre_state_file):
		npzfile = np.load(cre_state_file, allow_pickle=True)
		valid = self.find_valid_cellTypes(npzfile['cellIndex'])
		cre_props_valid = npzfile['props'].astype(np.float64)[:,valid,:]
		cre_index = npzfile['ccREIndex']
		cre_chr = cre_index[:,0]
		cre_coords = cre_index[:,1:].astype(np.int32)
		if self.filter0: #filter ouut any CREs which are fully state 0 in all cell types
			mask = np.array(np.hstack(([0], np.tile([1], cre_props_valid.shape[2]-1))))
			cre_props_adjusted = np.sum(cre_props_valid * mask, axis=2)
			where_row_not_zero = np.sum(cre_props_adjusted, axis=1) != 0
			cre_props_valid = cre_props_valid[where_row_not_zero]
			cre_chr = cre_cre[where_row_not_zero]
			cre_coords = cre_coords[where_row_not_zero]
		if self.log_props:
			cre_props_valid = np.log(cre_props_valid + self.log_props_add_val)/np.log(self.log_props_base)
		if self.cutState0:
			cre_props_valid = cre_props_valid[:,:,1:]

		return cre_props_valid, cre_chr, cre_coords

	def subset_based_on_chrom(self,chrom):
		self.chrom = chrom
		where_chr = np.where(self.TSS_chr_all == self.chrom)[0]
		self.exp_values = self.exp_values_all[where_chr]
		self.TSSs = self.TSSs_all[where_chr]

		self.TSS_window_props = self.TSS_window_props_all[where_chr]

		where_chr = np.where(self.cre_chr_all == self.chrom)[0]
		self.cre_props = self.cre_props_all[where_chr]
		self.cre_coords = self.cre_coords_all[where_chr]
		self.creM = self.cre_props.shape[0]
		self.creIndex_range = np.arange(self.creM)

	def subset_based_ongroup(self, group, m_thresh, s_thresh):
		self.m_thresh = m_thresh
		self.s_thresh = s_thresh
		self.group = group
		if self.group != 0:
			m = np.mean(self.exp_values, axis=1)
			s = np.std(self.exp_values, axis=1, ddof=1)
			group_dict= { 1: np.where((m <= self.m_thresh) & (s <= self.s_thresh))[0],
						  2: np.where((m <= self.m_thresh) & (s > self.s_thresh))[0],
						  3: np.where((m > self.m_thresh) & (s > self.s_thresh))[0],
						  4: np.where((m > self.m_thresh) & (s <= self.s_thresh))[0]}
			self.exp_values = self.exp_values[group_dict[self.group]]
			self.TSSs = self.TSSs[group_dict[self.group]]
			self.TSS_window_props = self.TSS_window_props[group_dict[self.group]]
		self.tssN = self.TSSs.shape[0]

	def linear_regression(self, X, Y, intercept=False):
		X = X.reshape(-1, self.stateN)
		Y = Y.reshape(-1,)
		fit_model = linear_model.LinearRegression(fit_intercept=intercept).fit(X,Y)
		model_coeffs = fit_model.coef_
		r_squared = fit_model.score(X,Y)
		return {'coeffs': model_coeffs, 'rsquare': r_squared, 'fit_model': fit_model}

	def set_initial_betas(self, method, lim, file_classic, file_PrRe):
		if method == "regress":
			self.set_initial_betas_byRegression()
		elif method == "random":
			self.set_initial_betas_byRandom(lim)
		elif method == "PrRe":
			self.set_initial_betas_byReadIn(file_PrRe)
		elif method == "file":
			self.set_initial_betas_byReadIn(file_classic)
		else:
			logging.error("not an accepted method to set the initial beta coefficients. Make sure to use one from the following set: {regress, random, PrRe, file}")

	def set_initial_betas_byRegression(self):
		self.initial_coeffs = self.linear_regression(self.TSS_window_props, self.exp_values)['coeffs']

	def set_initial_betas_byRandom(self, lim):
		loc = -np.abs(lim)
		scale = np.abs(lim) - loc
		self.initial_coeffs = stats.uniform.rvs(loc=loc, scale=scale, size=self.stateN, random_state=self.seed)

	def set_initial_betas_byReadIn(self, read_in_file):
		self.initial_coeffs = np.load(read_in_file, allow_pickle=True)

	def find_initial_weighted_sum(self):
		'''find weighted sum of ccRE props with initial coeffs
		given a 3 dimensional array ccreN * cellN * stateN with p_ijk equal to the proportion of ccre_i for celltype_j in state_k
				1 dimensional array stateN with c_i equal to the initial coefficient for state_i
		return a 2 dimensional array ccreN * cellN with w_ij equal to the sum of l = 0 to 26 of c_l*p_ijl'''
		self.cre_weighted_sum = np.sum(self.cre_props * self.initial_coeffs, axis=2)

	def build_pairing_array(self):
		self.pairing_array = np.zeros((self.creM, self.cellN, 3), dtype=np.bool) #3D array of creM * cellN * 3. Last dimension is for initial, best, and current

	def build_overall_pairing(self):
		self.pairing_overall = np.zeros((self.tssN, self.creM, self.cellN))

	def get_dim_of_pairing(self, cellN_in, last_dim_in):
		'''if calling for initial pairing, last_dim_in of 0,
		   if calling for best so far, last_dim_in of 1,
		   if calling for current, last_dim_in of 2'''
		return np.sum(self.pairing_array[:, cellN_in, last_dim_in] == 1) #will return number of CREs paired in that slice

	# def get_dim_of_overall_pairing(self):
	# 	return np.sum(self.pairing_overall[something] == 1)

	def build_beta_array(self):
		self.beta_array = np.zeros((self.stateN, self.iter_val+1), dtype=np.float32) #2D array to store betas for initial pairing 0 and then for every random iteration after that

	def compute_adj_distance(self, starts, stops, TSS):
		locsTA = (starts + stops) // 2
		distance = np.abs(locsTA - TSS) #using abs because don't care if upsteam or downstream
		adjusted_distance = distance**self.dist_gamma
		adjusted_distance[np.where(adjusted_distance == 0)[0]] = self.zero_val #anywhere that distance is 0, set to zero_val
		return adjusted_distance

	def compute_distance(self, starts, stops, TSS):
		midLocs = (starts + stops) // 2
		if self.abs_bool:
			distance = np.abs(midLocs - TSS)
		else:
			distance = midLocs - TSS
		return distance

	def adjust_by_distance(self, to_adjust, TSS, starts, stops):
		adj_dist = self.compute_adj_distance(starts, stops, TSS)
		Y = np.ones(to_adjust.shape[0])
		Y /= adj_dist
		if to_adjust.ndim == 2:
			adjusted = np.multiply(to_adjust, Y.reshape(-1, 1))
		elif to_adjust.ndim == 3:
			adjusted = np.multiply(to_adjust, Y.reshape(-1,1,1))
		elif to_adjust.ndim == 1:
			adjusted = np.multiply(to_adjust, Y.reshape(-1))
		return adjusted

	def find_within(self, within_dist, TSS):
		none_within=False
		windowMin = max(0, TSS - within_dist)
		e_windowMin = max(0, TSS - self.exclude_wp_dist)
		e_windowMax = min(TSS + self.exclude_wp_dist, self.chrSizes[self.chrom])
		windowMax = min(TSS + within_dist, self.chrSizes[self.chrom])
		if self.exclude_wp:
			within_bool1 = (self.cre_coords[:,1]>=windowMin) & (self.cre_coords[:,0]<=windowMax)
			if self.dn_cross_bound:
				not_within_bool2 = not ((self.cre_coords[:,1] >= e_windowMin) & (self.cre_coords[:,0] <= e_windowMax))
			else:
				not_within_bool2 = not ((self.cre_coords[:,1] <= e_windowMax) & (self.cre_coords[:,0] >= e_windowMin))
			within_bool = (within_bool1 & not_within_bool2)
		else:
			within_bool = (self.cre_coords[:,1]>=windowMin) & (self.cre_coords[:,0]<=windowMax)
		if np.sum(within_bool) == 0:
			none_within = True
		return within_bool, none_within

	def adjust_abc(self, to_adjust, to_adjust_parent, TSS, starts, starts_parent, stops, stops_parent):
		print(to_adjust.shape)
		'''First we multiply the activity by the contact'''
		adjusted_pt1 = self.adjust_by_distance(to_adjust, TSS, starts, stops)
		print(adjusted_pt1.shape)
		'''Next we need to divide by the summed activity*contact for everything within 5Mbp'''
		'''find everything within 5Mbp'''
		within_bool, none_within_flag = self.find_within(self.abc_dist, TSS)
		for_denom = np.sum(self.adjust_by_distance(to_adjust_parent[within_bool], TSS, starts_parent[within_bool], stops_parent[within_bool]))
		print(for_denom.shape)
		adjusted = adjusted_pt1/for_denom
		return adjusted


	def test_abc_adjust(self):
		for i, TSS in enumerate(self.TSSs):
			if i < 20:
				CREs_within_bool, none_within_bool = self.find_within(self.element_dist, TSS)
				if none_within_bool:
					continue
				else:
					CREs_within_weighted_sum = self.cre_weighted_sum[CREs_within_bool]
					CREs_within_starts = self.cre_coords[CREs_within_bool, 0]
					CREs_within_stops = self.cre_coords[CREs_within_bool, 1]
					self.adjust_abc(CREs_within_weighted_sum, self.cre_weighted_sum, TSS, CREs_within_starts, self.cre_coords[:,0], CREs_within_stops, self.cre_coords[:,1])
			else:
				quit()

	def powerset_generator(self, num_passing):
		'''create an iterator object which makes a powerset, excluding the first empty set and the last full set
		By definition, powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
		In practice, powerset_generator([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3)
		Further, because it's a generator and I'll use next() to grab subsets (while I want to randomly pick subsets) I've got some fancy shuffling.
		Sadly, that means all of the subsets of a certain length will not be shuffled with respect to each other'''
		np.random.seed(num_passing//10)
		s = np.random.permutation(num_passing) #a shuffled np.arange(num_passing), these are the elements that will be combined
		t = np.arange(1, num_passing) #these are the elements that will be used for the length of the combinations, we start with 1 to avoid the empty set
		np.random.shuffle(t) #we randomly shuffle t so that we don't always get all the 1 element subsets and then all the 2 element subsets etc.
		powerset_iterator = chain.from_iterable(combinations(s, r) for r in t)
		return powerset_iterator


	def find_num_in_powerset(self, num_passing):
		return ((2**num_passing)-2)

	# def check_and_get_i(self, checked, current_i, powerset_size):
	# 	'''Should be called before self.move_iterator()
	# 	Inputs: checked, list, which indices (i's) have been sampled/used in analysis already
	# 	current_i, integer, where the powerset iterator '''
	# 	new_i = False
	# 	next_i = 0
	# 	if (current_i in checked) or (current_i >= powerset_size):
	# 		new_i = True
	# 		np.random.seed(current_i)
	# 		next_i = np.random.choice(np.setdiff1d(np.arange(powerset_size), checked, assume_unique=True))
	# 	return new_i, next_i
	#
	# def move_iterator(self, current_i, new_i, next_i, checked, iterator_obj, num_passing):
	# 	'''This function is super cool cause it'll move through the iterator to get to a specific subset, even resetting the iterator if it needs to
	# 		Should be called after self.check_and_get_i()
	# 	Inputs: current_i, integer, where the powerset iterator is at following a single next(iterator_obj) call
	# 		   new_i, boolean, whether or not we need to move the powerset iterator to be at next_i instead of current_i, from self.check_and_get_i()
	# 		   next_i, integer, where the powerset iterator needs to be at if new_i is True, from self.check_and_get_i()
	# 		   checked, list, which indices (i's) have been sampled/used in analysis already
	# 		   iterator_obj, an iterator object, for the powerset that will be at current_i with a single next(iterator_obj) call
	# 		   num_passing, an integer, the number of CREs passing the within boundaries requirements
	# 	Returns: new_subset, a numpy array, the new subset randomly selected from the powerset iterator
	# 	  		 checked, list, updated with the new indice (i) that has been sampled/used in analysis already
	# 			 iterator_obj, an iterator object, for the powerset that will be at to_return_current with a single next(iterator_obj) call.
	# 			 to_return_current, an integer, the indice (i) where the powerset iterator will be following a single next(iterator_obj) call.
	# 		   '''
	# 	if new_i:
	# 		while True:
	# 			try:
	# 				new_subset = next(iterator_obj)
	# 				current_i += 1
	# 			except StopIteration:
	# 				'''reset the iterator'''
	# 				current_i = 0
	# 				iterator_obj = self.powerset_generator(num_passing)
	# 			if current_i == next_i:
	# 				new_subset = next(iterator_obj)
	# 				break
	# 		to_return_current = current_i + 1
	# 		checked.append(current_i-1)
	# 	else:
	# 		new_subset = next(iterator_obj)
	# 		to_return_current = current_i + 1
	# 		checked.append(current_i)
	# 	return np.array(new_subset), checked, iterator_obj, to_return_current

	def move_iterator_for_sorted(self, iterator_obj, current_i, next_i):
		while True:
			try:
				new_subset = next(iterator_obj)
				current_i += 1
			except StopIteration:
				'''shouldn't actually need this'''
				logging.error("Something bad happened with iteration")
				break
			if current_i == next_i:
				new_subset = next(iterator_obj)
				break
			elif current_i > next_i:
				break
		yield new_subset

	def indice_iterator(self, powerset_size, gene_index):
		np.random.seed(int((gene_index+1)*self.seed))
		indice_iterator = iter(sorted(np.random.choice(np.arange(powerset_size), size=self.iter_val, replace=False)))
		return indice_iterator

	# def get_from_powerset(self, iterator_obj, checked, current_i, powerset_size, num_passing):
	# 	new_i, next_i = self.check_and_get_i(checked, current_i, powerset_size)
	# 	new_subset, checked, iterator_obj, returned_current = self.move_iterator(current_i, new_i, next_i, checked, iterator_obj, num_passing)
	# 	return new_subset, checked, iterator_obj, returned_current

	def evalute_pairing(self):
		return 0

	def drive_random_iter(self, i):
		'''
		for each gene, get the num passing using whatever cutoffs <-- can I parallelize this part of it?
		initialize empty list, checked
		create the iterator object by calling self.powerset_generator()
		find the powerset size by calling self.find_num_in_powerset()
		make an iterator of all of the random indices we'll check, sorted.

		for each iter
		call move_iterator_for_sorted
		evaluate pairing
		'''
		#get num_passing
		TSS = self.TSSs[i]
		CREs_within_bool, none_within_bool = self.find_within(self.element_dist, TSS)
		if none_within_bool: #no pairing possible
			return
		else:
			CREs_within_weighted_sum = self.cre_weighted_sum[CREs_within_bool]
			CREs_within_starts = self.cre_coords[CREs_within_bool, 0]
			CREs_within_stops = self.cre_coords[CREs_within_bool, 1]

			'''adjusted_by_abc should be np.sum(CREs_within_bool) X self.cellN -- and it is now!'''
			adjusted_by_abc = self.adjust_abc(CREs_within_weighted_sum, self.cre_weighted_sum, TSS, CREs_within_starts, self.cre_coords[:,0], CREs_within_stops, self.cre_coords[:,1])
			#Find num_passing for each cell type using the abc_thresh
			#add code for using correlation if use_corr

			#need gene_index either from a for loop or by calling this function in a parallel manner
			powerset_it_obj = self.powerset_generator(num_passing)
			powerset_size = self.find_num_in_powerset(num_passing)
			indice_it_obj = self.indice_iterator(powerset_size, gene_index)
			current_i = 0
			while True:
				try:
					next_indice = next(indice_it_obj)
				except StopIteration:
					break
				new_subset = np.array(list(self.move_iterator_for_sorted(powerset_it_obj, current_i, next_indice)))
				current_i = next_indice + 1
				#evaluate pairing for each cell type
main()
