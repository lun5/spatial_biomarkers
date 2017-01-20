import sqlite3, csv, os, math, time
from datetime import datetime
import numpy as np
from scipy.spatial.distance import cdist, pdist
from sklearn.neighbors import NearestNeighbors, KDTree
from sklearn import linear_model
import scipy.io
from itertools import product
import math
import multiprocessing
from functools import partial
import sys
import argparse
os.system("taskset -p 0xff %d" % os.getpid())

def parse_args(argv):
    print 'in parse args'
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_fname", default='Filtered.QCed.NoTcorr.SegExclude.HighDAPI.log2.Median.Norm.Cells.csv',
                       help='biomarker csv file')
    parser.add_argument("--data_path", default="/home/lun5/Documents/multiplex/data/",
                       help='directory of spots.db')
    parser.add_argument("--input_path", default="/home/lun5/Documents/multiplex/input/",
                       help='directory of biomarker csv file')
    parser.add_argument("--output_path", default="/home/lun5/Documents/output/",
                       help='directory of output file')
    return parser.parse_args(argv)


def read_data(input_path, data_fname, data_path):
    f = open(os.path.join(input_path,data_fname),'rb')
    f_csv = csv.DictReader(f)
    # load data
    biomarkers = np.load(os.path.join(data_path,'biomarkers.npy'))
    # get biomarker names from csv file
    fields = f_csv.fieldnames[3:15] + f_csv.fieldnames[16:]
    # only take nuclear biomarkers
    nuc_bio_indx = ['.Nuc.' in fields[i] for i in xrange(len(fields))]
    split_name = [str.split(fields[i],'.')[-1] for i in xrange(len(fields))]
    nuc_bio_indx = np.nonzero(np.array(nuc_bio_indx))[0]
    nuc_biomarkers_names = [split_name[nuc_bio_indx[i]] for i in xrange(len(nuc_bio_indx))]
    # only take high quality biomarkers
    f_quality = open(os.path.join(input_path,'biomarkers_quality.csv'),'rb')
    f_csv_quality = csv.DictReader(f_quality)
    marker = []
    quality = []
    for row in f_csv_quality:
        marker.append(row['marker'])
        quality.append(row['quality'])
    # compare lower case names
    nuc_biomarkers_names = [nuc_biomarkers_names[i].lower() for i in xrange(len(nuc_biomarkers_names))]
    marker = [marker[i].lower() for i in xrange(len(marker))]
    nuc_biomarker_qual = np.array([int(quality[marker.index(nuc_biomarkers_names[i])]) 
                               for i in xrange(len(nuc_bio_indx))])
    indx_good_nuc_biomarker = nuc_bio_indx[nuc_biomarker_qual == 1]
    
    # delete biomarkers with too high variance
    delete_bm_indx = [12,19,24,27]
    delete_bm = [str.split(fields[nuc_bio_indx[i]],'.')[-1] for i in delete_bm_indx]
    good_nuc_biomarkers = [split_name[indx_good_nuc_biomarker[i]] for i in xrange(len(indx_good_nuc_biomarker))]
    indx_delete = [good_nuc_biomarkers.index(delete_bm[i]) for i in xrange(1,len(delete_bm))]
    indx_good_nuc_biomarker = np.delete(indx_good_nuc_biomarker,indx_delete)
    good_nuc_biomarkers = [split_name[indx_good_nuc_biomarker[i]] for i in xrange(len(indx_good_nuc_biomarker))]
    bm_data = biomarkers[:,indx_good_nuc_biomarker]
    meta_data = biomarkers[:,[0,1,2,7,8]]
    
    # read ingested spot
    sqlite_file = os.path.join(data_path, 'spots.db')
    conn = sqlite3.connect(sqlite_file)
    conn.text_factory = str
    c = conn.cursor()
    sel_cmd = ('SELECT {tn1}.{coi1} '
           +' FROM {tn1} JOIN {tn2} ON {tn1}.{coi1} = {tn2}.{coi1}'
           + ' WHERE {tn2}.{coi2} == {0}').format(1, 
           coi1='spot_id',coi2='ingested', tn1='spots',tn2 = 'ingestion')
    c.execute(sel_cmd)
    results = c.fetchall()
    ingested_spot_id =  np.asarray([results[i][0] for i in xrange(len(results))])
    
    return (bm_data, meta_data, ingested_spot_id, conn)


def calculate_spot_features(spot_id, all_input, all_params):
    bm_data, meta_data, db_path, NN_OUTPUT, scramble_seed = all_input
    alpha, radius, l1_ratio = all_params
    
    conn = sqlite3.connect(db_path)
    conn.text_factory = str

    c = conn.cursor()
    sel_cmd = ('SELECT {coi1} FROM {tn1} WHERE {coi2}=={0}').format(spot_id, 
           coi1='cell_id',coi2='spot_id', tn1='cells')
    
    c.execute(sel_cmd)
    results = c.fetchall()
    cell_id =  np.asarray([results[i][0] for i in xrange(len(results))])
    spot_bm = bm_data[cell_id - 1, :] # biomarker values of the spot            
    spot_meta = meta_data[cell_id -1, :] # x,y coordinates of each cell
    
    sel_cmd = ('SELECT {coi1} FROM {tn1} WHERE {coi2}=={0}').format(spot_id, 
           coi1='spot_name',coi2='spot_id', tn1='spots')
    
    c.execute(sel_cmd)
    results = c.fetchall()
    spot_name =  results[0][0]
    
    # get nearest neighbors
    spot_xy = spot_meta[:,:2]
    tree = KDTree(spot_xy)
    ind, dist = tree.query_radius(spot_xy, r = radius, return_distance = True)
    
    # if scramble
    if scramble_seed:
        np.random.seed(scramble_seed)
        spot_bm = spot_bm[np.random.permutation(np.arange(len(cell_id))),:]
        output_name = os.path.join(NN_OUTPUT,spot_name +  '_seed_' + str(scramble_seed) + '.mat')
    else:
        output_name = os.path.join(NN_OUTPUT,spot_name + '.mat')
    if os.path.isfile(output_name):
        return 1 
    #num_non_zeros_nn = [] # number of contributing neighbors
    entropies = [] # entropy
    angle_std = [] # std
    neighbor_id = [] # id of neighbor
    residuals = [] # residuals at each cells
    coeff_neighbors = [] # coefficient of neighbor each cell
    pmf_neighbors = []
    
    for j in xrange(spot_bm.shape[0]):
        clf = linear_model.ElasticNet(alpha=alpha, l1_ratio=l1_ratio)
        
        cell_bm = spot_bm[j,:].reshape(1,-1) # cell of question
        cell_xy = spot_xy[j,:]
        
        nb_bm = spot_bm[ind[j][dist[j]!=0], :] # eliminate itself
        nb_xy = spot_xy[ind[j][dist[j]!=0], :]
        
        curr_neigh_ind =   ind[j][dist[j]!=0]
        neighbor_id.append(curr_neigh_ind)      

        if np.min(nb_bm.shape) > 0:
            clf.fit(nb_bm.T,cell_bm.T)
            pred = clf.predict(nb_bm.T).reshape(1,-1)
            resid_biomarkers = np.square(cell_bm - pred)
            residuals.append(np.sqrt(resid_biomarkers.sum())/cell_bm.shape[1])
            coeff_neighbors.append(clf.coef_)
            
            nb_diff_xy = nb_xy - cell_xy
            nb_diff_xy = nb_diff_xy/np.linalg.norm(nb_diff_xy,axis = 1).reshape(-1,1) # normalize the vector
            
            if len(clf.coef_) >= 2: # and (np.linalg.norm(clf.coef_) > 0):
                nb_angles = [0]
                for k in xrange(len(clf.coef_)-1):
                    a = np.dot(nb_diff_xy[0,:],nb_diff_xy[k+1,:])
                    a = np.sign(a)*np.minimum(1,np.abs(a))
                    nb_angles.append(math.acos(a))
                nb_angles = [(nb_angles[i] + np.pi/100) % 2*np.pi for i in xrange(len(nb_angles))]
                bin_vec = np.linspace(0, 2*np.pi, num=36)
                nb_bin = [(bin_vec -nb_angles[i] >=0).nonzero()[0][0] for i in range(len(nb_angles))]
                pk = np.zeros(36)
                for i in xrange(len(nb_angles)):
                    pk[nb_bin[i]] = pk[nb_bin[i]] + np.abs(clf.coef_[i])
                
                pmf_neighbors.append(np.array(pk))
                cell_entropy = scipy.stats.entropy(pk + 1e-6)
                if np.isinf(cell_entropy):
                    mat_dict={'passed': 0, 'spot_meta':spot_meta,'entropies':entropies,
                            'residuals':residuals,'angle_std':angle_std,'pmf':pmf_neighbors,
                            'coeff_neigh':coeff_neighbors, 'neighbor_id':neighbor_id
                           }
                    scipy.io.savemat(os.path.join(NN_OUTPUT,spot_name + '.mat'), mat_dict)                    
                    return 0
                entropies.append(cell_entropy)
                coef_ = np.abs(clf.coef_)/np.sum(np.abs(clf.coef_))
                cos_mean_resultant = np.sum([coef_[i]*math.cos(nb_angles[i]) for i in range(len(coef_))])
                sin_mean_resultant = np.sum([coef_[i]*math.sin(nb_angles[i]) for i in range(len(coef_))])
                mean_resultant = np.sqrt(cos_mean_resultant**2+sin_mean_resultant**2)
                cell_angle_std = np.sqrt(2*(1-mean_resultant))
                if np.isnan(cell_angle_std):
                    mat_dict={'passed': 0, 'spot_meta':spot_meta,'entropies':entropies,
                            'residuals':residuals,'angle_std':angle_std,'pmf':pmf_neighbors,
                            'coeff_neigh':coeff_neighbors, 'neighbor_id':neighbor_id
                           }
                    scipy.io.savemat(output_name, mat_dict)                    
                    return 0
                angle_std.append(cell_angle_std)
            else:
                pmf_neighbors.append([0])
                entropies.append(0)
                angle_std.append(0)
        else:
            #print spot_name, j, nb_bm.shape                   
            residuals.append(0)
            entropies.append(0)
            angle_std.append(0)
            pmf_neighbors.append([])
            coeff_neighbors.append([])
    
    mat_dict={'passed': 1, 'spot_meta':spot_meta,'entropies':entropies,
              'residuals':residuals,'angle_std':angle_std,'pmf':pmf_neighbors,
              'coeff_neigh':coeff_neighbors, 'neighbor_id':neighbor_id}
    scipy.io.savemat(output_name, mat_dict) 
    
    return 1

def run_calculation():
    args = {}
    args['output_path'] = "/home/lun5/Documents/multiplex/output/"
    args['data_fname'] = 'Filtered.QCed.NoTcorr.SegExclude.HighDAPI.log2.Median.Norm.Cells.csv'
    args['input_path'] = "/home/lun5/Documents/multiplex/input/"
    args['data_path'] = "/home/lun5/Documents/multiplex/data/"
    
    
    if not os.path.isdir(args['output_path']):
        os.mkdir(args['output_path'])
    
    # COMBINATIONS OF PARAMETER (WILL BE IN ARGS)
    alphas = [0.1, 1]
    radii = [50,100] # Values: 50, 100
    l1_ratios = [0,0.5,1] # need to add 0 in there
    
    # READ IN DATA
    bm_data, meta_data, ingested_spot_id, conn = read_data(
        args['input_path'], args['data_fname'], args['data_path'])
    
    db_path = os.path.join(args['data_path'], 'spots.db')
    for alpha, radius, l1_ratio in product(alphas, radii,l1_ratios):
        NN_OUTPUT = os.path.join(args['output_path'], ('Entropy_radius_' + str(radius)) + '_alpha_' +str(alpha) 
                          + '_L1ratio_' + str(l1_ratio))
        if os.path.isdir(NN_OUTPUT):
            print(('Elastic: Already started calculation for alpha = %0.2f radius = %d, l1_ratio = %0.2f') %(
                alpha, radius, l1_ratio))
        else:
            print(('Elastic: Have not calculated for alpha = %0.2f radius = %d, l1_ratio = %0.2f') %(
                alpha, radius, l1_ratio))
            os.makedirs(NN_OUTPUT)

        start_time = time.time()
        
        partial_spot_features = partial(calculate_spot_features, all_input=(
                    bm_data, meta_data, db_path, NN_OUTPUT, []), all_params = (alpha, radius, l1_ratio))
        #feature_results = partial_spot_features(ingested_spot_id[0])
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-1)
        os.system("taskset -p 0xff %d" % os.getpid())
        feature_results = pool.map(partial_spot_features, ingested_spot_id, 500)
        pool.close()
        pool.join()
        non_passed_cases = ingested_spot_id[np.nonzero(np.array(feature_results) == 0)[0]]
        print 'Non scramble: id of cases not passed are', non_passed_cases
        
        
        # scramble
        
        for seed in [1,50,200]:
            partial_spot_features = partial(calculate_spot_features, all_input=(
                        bm_data, meta_data, db_path, NN_OUTPUT, seed), 
                                            all_params = (alpha, radius, l1_ratio))
            #feature_results = partial_spot_features(ingested_spot_id[0])
            pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-1)
            os.system("taskset -p 0xff %d" % os.getpid())
            feature_results = pool.map(partial_spot_features,ingested_spot_id, 50)
            pool.close()
            pool.join()
            non_passed_cases = ingested_spot_id[np.nonzero(np.array(feature_results) == 0)[0]]
            print 'Scramble with seed ', seed, ', id of cases not passed are', non_passed_cases
        
        print('Done with alpha=%0.1f, radius = %d in %0.2f seconds' %(
                alpha,radius, time.time() - start_time))    

def main():
    run_calculation()

if __name__ == "__main__":
    sys.exit(main())

