#mem usage 20GB
#carries out igeGWAS conditioning on direct genotypes at each SNP

import sys
import os
import h5py
import numpy as np
sys.path.insert(0,'/homes/abaud/CFW/code/reproduce/LIMIX')
from social_data_CFW import SocialData
# noNoise_qtl runs GWAS without estimating the covariance structure: uses the total covariance matrix provided (built from DGE IGE DEE IEE and cage effects)
import noNoise_qtl as qtl
import h5py
import pdb
import csv

task = 'VD_CFW'

kinship_type = sys.argv[2]
subset = sys.argv[3]

phenos_to_run_file = open('/homes/abaud/CFW/data/reproduce/phenotypes/phenos_to_run.txt')
phenos_to_run = csv.reader(phenos_to_run_file)
for i, row in enumerate(phenos_to_run):
	if i==(int(sys.argv[1])-1):
		col = int(row[0])-1
		break

data = SocialData(task,kinship_type,subset)
doto = data.get_data(col)
trait = doto['trait']
print trait

DGE = "DGE"
IGE = "IGE"
IEE = "IEE"
cageEffect = "cageEffect"

#try opening pvalues file early so that lmm doesnt run if file opening is going to fail...
pvalues_file_dir = "".join(['/homes/abaud/CFW/output/reproduce/pvalues_pointwise_conditional/IGE/',kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect]))])

if not os.path.exists(pvalues_file_dir):
	os.makedirs(pvalues_file_dir)
pvalues_file_name = "".join([pvalues_file_dir,'/',trait,'.h5'])
#if os.path.exists(pvalues_file_name):
#	sys.exit(0)

covar_outfile_dir = "".join(["/homes/abaud/CFW/output/reproduce/null_covars/",kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect]))])
covar_outfile_name = "".join([covar_outfile_dir,"/",trait,".h5"])
covar_outfile = h5py.File(covar_outfile_name,'r')
saved = {}
saved['sampleID'] = covar_outfile['sampleID'][:]
saved['pheno'] = covar_outfile['pheno'][:]
saved['covs'] = covar_outfile['covs'][:]
saved['covar_mat'] = covar_outfile['covar_mat'][:]
covar_outfile.close()

#import genotypes. direct and social geno already matched on rows and cols
input_file_name = '/homes/abaud/CFW/data/reproduce/CFWmice.h5'
input_file = h5py.File(input_file_name,'r')
geno  = input_file['direct_pruned_dosages']
geno_matrix = geno['matrix'][:].T
geno_sample_ID = geno['row_header']['sample_ID'][:]
social_geno  = input_file['social_pruned_dosages']
social_geno_matrix = social_geno['matrix'][:].T
position = {
	"chr" : geno['col_header']['chr'][:],
	"pos"   : geno['col_header']['pos'][:],
	"cumpos" : geno['col_header']['cumpos'][:]
}
input_file.close()

#match genotypes with the rest
Imatch = np.nonzero(saved['sampleID'][:,np.newaxis]==geno_sample_ID)
print "Imatch"+str(len(Imatch[0]))
saved['sampleID'] = saved['sampleID'][Imatch[0]]
saved['pheno'] = saved['pheno'][Imatch[0]]
saved['covs'] = saved['covs'][Imatch[0],:]
saved['covar_mat'] = saved['covar_mat'][Imatch[0],:][:,Imatch[0]]
geno_matrix = geno_matrix[Imatch[1],:]
geno_sample_ID=geno_sample_ID[Imatch[1]]
social_geno_matrix = social_geno_matrix[Imatch[1],:]

pvalues_file = h5py.File(pvalues_file_name,'w')

pvalues = [-1] * social_geno_matrix.shape[1]
betas = [-1] * social_geno_matrix.shape[1]
lmm = qtl.test_lmm(snps=social_geno_matrix[:,0:1],pheno=saved['pheno'],covs = np.concatenate([saved['covs'],geno_matrix[:,0:1]], axis=1),K=saved['covar_mat'],NumIntervalsDelta0=2,searchDelta=False)
for i in xrange(social_geno_matrix.shape[1]):
	covs = np.concatenate([saved['covs'],geno_matrix[:,i:i+1]], axis=1)
	snp = social_geno_matrix[:,i:i+1]
	lmm.setCovs(covs)
	lmm.setSNPs(snp)
	lmm.process()
	pvalues[i] = lmm.getPv()[0,0]
	betas[i] = lmm.getBetaSNP()[0,0]

pvalues_file.create_dataset(name = 'chr',data = position['chr'])
pvalues_file.create_dataset(name = 'pos',data = position['pos'])
pvalues_file.create_dataset(name = 'cumpos',data = position['cumpos'])
pvalues_file.create_dataset(name = 'pvalues',data = pvalues)
pvalues_file.create_dataset(name = 'betas',data = betas)
pvalues_file.close()

