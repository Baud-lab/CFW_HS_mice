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
    
    if len(sys.argv) == 4:
        subset = sys.argv[3]
    else:
        subset = None

    DGE = "DGE"
    IGE = "IGE"
    IEE = "IEE"
    cageEffect = "cageEffect"
    
    data = SocialData(task,kinship_type,subset)
    doto = data.get_data(col)
    
    trait=doto['trait']
    print('trait is ' + trait)
    
    VD_outfile_dir = "".join(['/homes/abaud/CFW/output/reproduce/VD/',kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect]))])
    if not os.path.exists(VD_outfile_dir):
        os.makedirs(VD_outfile_dir)
    VD_outfile_name="".join([VD_outfile_dir,'/',trait,'.txt'])
    #if os.path.exists(VD_outfile_name):
    #    sys.exit(0)
   
    vc = DirIndirVD(vc_init_type = None, vc_init = None, pheno = doto['pheno'],pheno_ID = doto['pheno_ID'],covs = doto['covs'], covs_ID = doto['covs_ID'], covariates_names = doto['covariates_names'], kinship_all = doto['kinship_full'], kinship_all_ID = doto['kinship_full_ID'], cage_all = doto['cage_full'], cage_all_ID = doto['cage_full_ID'], DGE = (DGE is not None), IGE = (IGE is not None), IEE = (IEE is not None), cageEffect = (cageEffect is not None), calc_ste=True, subset_IDs = doto['subset_IDs'], SimplifNonIdableEnvs = False)

    rv=vc.getOutput()
    
    toWrite=(trait,rv['sample_size'],rv['sample_size_all'],rv['var_Ad'],rv['var_As'],rv['STE_Ad'],rv['STE_As'],rv['total_var'],rv['corr_Ads'],rv['STE_corr_Ads'],rv['corr_params'],rv['conv'],rv['LML'],rv['var_Ed'],rv['var_Es'],rv['corr_Eds'],rv['var_C'],rv['env_rho'],rv['env_sigma2'])
    VD_outfile=open(VD_outfile_name,'w')
    VD_outfile.write("\t".join(str(e) for e in toWrite)+'\n')
    VD_outfile.close()



