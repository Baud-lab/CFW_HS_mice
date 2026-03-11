#runs the VD with the correlation rho between DGE and IGE constrained to be 0

import os
import sys
sys.path.insert(0,'/homes/abaud/CFW/code/reproduce/LIMIX')
from social_data_CFW import SocialData
from dirIndirVD_wSS import DirIndirVD
import pdb
import csv

if __name__=='__main__':
    
    phenos_to_run_file = open('/homes/abaud/CFW/data/reproduce/phenotypes/phenos_to_run.txt')
    phenos_to_run = csv.reader(phenos_to_run_file)
    for i, row in enumerate(phenos_to_run):
        if i==(int(sys.argv[1])-1):
            col = int(row[0])-1
            break

    task = 'VD_CFW'
   
    kinship_type = sys.argv[2]
    
    subset = sys.argv[3]

    DGE = "DGE"
    IGE = "IGE"
    IEE = "IEE"
    cageEffect = "cageEffect"

    data = SocialData(task,kinship_type,subset)
    doto = data.get_data(col)
    
    trait=doto['trait']
    print(trait)
    
    VD_outfile_dir = "".join(['/homes/abaud/CFW/output/reproduce/VD/',kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect])),'/'])
    if not os.path.exists(VD_outfile_dir):
        os.makedirs(VD_outfile_dir)
    constrained_VD_outfile_name="".join([VD_outfile_dir,'/',trait,'_constrained1.txt'])
   
    if os.path.exists(constrained_VD_outfile_name):
        sys.exit(0)

    #ff is free form
    vc_ff = DirIndirVD(vc_init_type = None, vc_init = None, pheno = doto['pheno'],pheno_ID = doto['pheno_ID'],covs = doto['covs'], covs_ID = doto['covs_ID'], covariates_names = doto['covariates_names'], kinship_all = doto['kinship_full'], kinship_all_ID = doto['kinship_full_ID'], cage_all = doto['cage_full'], cage_all_ID = doto['cage_full_ID'], DGE = (DGE is not None), IGE = (IGE is not None), IEE = (IEE is not None), cageEffect = (cageEffect is not None), calc_ste=False, subset_IDs = doto['subset_IDs'], SimplifNonIdableEnvs = False)
    
    rv_ff=vc_ff.getOutput()
    toWrite_ff=(trait,rv_ff['sample_size'],rv_ff['sample_size_cm'],rv_ff['var_Ad'],rv_ff['var_As'],rv_ff['total_var'],rv_ff['corr_Ads'],rv_ff['conv'],rv_ff['LML'])
    constrained_VD_outfile=open(constrained_VD_outfile_name,'w')
    constrained_VD_outfile.write("\t".join(str(e) for e in toWrite_ff)+'\n')

    #diag is diagonal. equivalent to imposing correlation 0 between DGE and IGE.
    vc_diag = DirIndirVD(vc_init_type = 'diagonal', vc_init = vc_ff, pheno = doto['pheno'],pheno_ID = doto['pheno_ID'],covs = doto['covs'], covs_ID = doto['covs_ID'], covariates_names = doto['covariates_names'], kinship_all = doto['kinship_full'], kinship_all_ID = doto['kinship_full_ID'], cage_all = doto['cage_full'], cage_all_ID = doto['cage_full_ID'], DGE = (DGE is not None), IGE = (IGE is not None), IEE = (IEE is not None), cageEffect = (cageEffect is not None), calc_ste=False, subset_IDs = doto['subset_IDs'], SimplifNonIdableEnvs = False)
    rv=vc_diag.getOutput()
    toWrite_diag=(trait,rv['sample_size'],rv['sample_size_cm'],rv['var_Ad'],rv['var_As'],rv['total_var'],rv['corr_Ads'],rv['conv'],rv['LML'])
    constrained_VD_outfile.write("\t".join(str(e) for e in toWrite_diag)+'\n')

    constrained_VD_outfile.close()
