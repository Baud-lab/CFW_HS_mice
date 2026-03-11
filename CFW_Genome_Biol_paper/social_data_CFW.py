import scipy as sp
import h5py
#re for regular expressions
import re
import pdb

class SocialData():
    
    def __init__(self, task = None, kinship_type = "", subset = None, effect = None, chr = None):
        assert task is not None, 'Specify task!'
        self.task=task
        self.kinship_type=kinship_type
        self.subset = subset
        self.effect = effect
        self.chr = chr
        self.load()
    
    
    def load(self):
        
        if 'VD_CFW' in self.task:
            in_file = '/homes/abaud/CFW/data/reproduce/CFWmice.h5'
            print(in_file)

            f = h5py.File(in_file,'r')
    
            self.measures = f['phenotypes']['col_header']['phenotype_ID'].asstr()[:]
            self.all_pheno = f['phenotypes']['matrix'][:].T
            self.pheno_ID = f['phenotypes']['row_header']['sample_ID'][:]
            self.all_covs = f['covariates2']['matrix'][:].T
            self.covs_ID = f['covariates2']['row_header']['sample_ID'][:]
            self.covariates = f['covariates2']['col_header']['covariate_ID'].asstr()[:]
            self.cage_full = f['phenotypes']['row_header']['cage'].asstr()[:]
            self.cage_full_ID = f['phenotypes']['row_header']['sample_ID'][:]
            self.all_covs2use = f['phenotypes']['col_header']['covariatesUsed']

            if self.chr is not None:
                self.kinship_full = f['GRM'][self.kinship_type][''.join(['chr',str(chr)])],['matrix'][:]
            else:    
                self.kinship_full = f['GRM'][self.kinship_type]['matrix'][:]
            self.kinship_full_ID = f['GRM'][self.kinship_type]['row_header']['sample_ID'][:]
 
            if self.subset is None:
                self.subset_IDs = None
            else:
                self.subset_IDs = f['subsets'][self.subset][:]

        else:
            print("Nothing done: task unknown!")


    def get_data(self,col):
        
        self.trait = self.measures[col]
        self.pheno = self.all_pheno[:,col]

        #that's if no covs in entire study
        if self.all_covs2use is None:
            self.covs = None
            self.covs_ID = None
            covariates_names = None
        else:
            covs2use = self.all_covs2use.asstr()[col].split(',')
            Ic = sp.zeros(self.covariates.shape[0],dtype=bool)
            for cov in covs2use:
                Ic = sp.logical_or(Ic,self.covariates==cov)
            # covariates_names will be empty list rather than None if no cov for that phenotype (col)
            covariates_names = self.covariates[Ic]
            print('Initial covs in social_data are ' + str(covariates_names))
            if len(self.all_covs.shape)==1:
                if Ic:
                    self.covs = self.all_covs
                else:
                    self.covs = None
            else:
                self.covs = self.all_covs[:,Ic]
        
        return {'trait' : self.trait,
                'pheno' : self.pheno,
                'pheno_ID' : self.pheno_ID,
                'covs' : self.covs,
                'covs_ID' : self.covs_ID,
                'covariates_names' : covariates_names,
                'kinship_type' : self.kinship_type,
                'kinship_full' : self.kinship_full,
                'kinship_full_ID' : self.kinship_full_ID,
                'cage_full' : self.cage_full,
                'cage_full_ID' : self.cage_full_ID,
                'subset_IDs' : self.subset_IDs}






