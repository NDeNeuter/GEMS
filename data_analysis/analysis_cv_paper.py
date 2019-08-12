import os

from joblib import Parallel, delayed

import GeneExpressionAnalysis as GEA


def deseq_iteration(readcount_path, col_data_path):
    
    run = GEA.DESeq2run(readcount_path,
                    col_data_path,
                    sample_column='tags')
    
    result_paths = run.calc_DEGs(design = 'group2',
                                 to_test={'group2': ['EV', 'BM']},
                                 R_script_dir='/'.join(readcount_path.split('/')[:-1]))
    
    for path in result_paths:
        deg_result = GEA.DESeq2results(path)
        # write out significant DEG results with different log fold thresholds
        for fold_thres in [0, 1]:
            deg_result.get_sig_DEGs(path.replace('.txt', '_fold{}_sig.txt'.format(fold_thres)),
                                                 fold_threshold=fold_thres)
        # write out all results
        deg_result.get_sig_DEGs(path.replace('.txt', '_all.txt'), padj_threshold=1)
        
        
if __name__ == '__main__':
    
    maindir = '/Users/nicolasdeneuter/Dropbox/PhD/Projects/GOA/GEMS'
    
    readcount_df = GEA.collect_data(maindir+'/original_data', file_pattern = '/readcounts.txt')
    readcount_df = readcount_df.reindex(sorted(readcount_df.columns), axis=1)
    
    preproc = GEA.Preprocessor()
    readcount_df = preproc.clean(readcount_df)
    readcount_df = preproc.combine_lane_counts(readcount_df)
    
    ''' Viral EV1: meningitis vs bacterial BM1: meningitis - leave-one-out cross validation setup'''
    
    analysis_dir = '{}/data_paper/plan_classification_samples_without_BM2'.format(maindir)
    try:
        os.mkdir(analysis_dir)
    except FileExistsError:
        pass
    
    lp = GEA.LabelProcessor(readcount_df)
    lp.add_samples('group', 'EV1')
    lp.add_samples('group', 'BM1')
    lp.add_samples('group', 'EV2')
    #lp.add_samples('group', 'BM2')
    lp.regroup_meta_tags('sample_type', 'M_vs_C', {'M': ['M', ''], 'C': ['C']})
    lp.add_samples('group', 'VM')
    lp.add_samples('group', 'PN')
    lp.add_samples('group', 'REU1')
    lp.add_samples('group', 'REU2')
    
    lp.remove_samples('M_vs_C', 'C')
    lp.regroup_meta_tags('group', 'group2', {'EV': ['EV1', 'EV2'], 'BM': ['BM1']})
    
    # generate whole dataset and write to disk
    lp.generate_DESeq_data('{}/read_counts.txt'.format(analysis_dir),
                           '{}/col_data.txt'.format(analysis_dir))
    
    # run DESeq2 analysis on all the data
    run_1 = GEA.DESeq2run('{}/read_counts.txt'.format(analysis_dir),
                      '{}/col_data.txt'.format(analysis_dir),
                      sample_column='tags')
    
    result_paths = run_1.calc_DEGs(design = 'run + group2',
                                   to_test={'group2': ['EV', 'BM']},
                                   R_script_dir=analysis_dir)
    
    for path in result_paths:
        deg_result = GEA.DESeq2results(path)
        for fold_thres in [0, 1]:
            deg_result.get_sig_DEGs(path.replace('.txt', '_fold{}_sig.txt'.format(fold_thres)),
                                                 fold_threshold=fold_thres)
            deg_result.PCA('group', filename='PCA_normalized_sig_genes_foldthres{}'.format(fold_thres))
            deg_result.heatmap(filename='heatmap_normalized_sig_genes_foldthres{}'.format(fold_thres))
    
    
    ###############################
    ## do cross validation setup ##
    ###############################
    
    analysis_dir = '{}/data_paper/plan_classification_samples_without_BM2'.format(maindir)
    try:
        os.mkdir(analysis_dir)
    except FileExistsError:
        pass
    
    # includes writing results to disk
    outputpaths = lp.generate_CV_DESeq_data(outputdir=analysis_dir)
    Parallel(n_jobs=7)(delayed(deseq_iteration)(readcount_path, col_data_path) for readcount_path, col_data_path in outputpaths)