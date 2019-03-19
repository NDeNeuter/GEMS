#!/usr/bin/env python3

from GeneExpressionAnalysis import *

if __name__ == '__main__':
    
    maindir = '..'
    
    readcount_df = collect_data(maindir+'/original_data', file_pattern = '/readcounts.txt')
    readcount_df = readcount_df.reindex(sorted(readcount_df.columns), axis=1)
    
    preproc = Preprocessor()
    readcount_df = preproc.clean(readcount_df)
    readcount_df = preproc.combine_lane_counts(readcount_df)
    
   
    ''' Viral EV1EV2: meningitis vs bacterial BM1BM2: meningitis'''
    
    analysis_dir = '{}/data_paper/plan_EV1EV2_M_vs_B1B2_M'.format(maindir)
    try:
        os.mkdir(analysis_dir)
    except FileExistsError:
        pass
    
    lp = LabelProcessor(readcount_df)
    lp.add_samples('group', 'EV1')
    lp.add_samples('group', 'BM1')
    lp.add_samples('group', 'EV2')
    lp.add_samples('group', 'BM2')
    lp.regroup_meta_tags('sample_type', 'M_vs_C', {'M': ['M', ''], 'C': ['C']})
    lp.remove_samples('M_vs_C', 'C')
    lp.regroup_meta_tags('group', 'EV_vs_B', {'EV': ['EV1', 'EV2'], 'B': ['BM1', 'BM2']})
    
    lp.generate_DESeq_data('{}/read_counts.txt'.format(analysis_dir),
                           '{}/col_data.txt'.format(analysis_dir))
    
    R_script_dir = '{}/code/generated_R_scripts'.format(maindir)
    
    run_1 = DESeq2run('{}/read_counts.txt'.format(analysis_dir),
                      '{}/col_data.txt'.format(analysis_dir),
                      sample_column='tags')
    
    result_paths = run_1.calc_DEGs(design = 'run + EV_vs_B',
                                   to_test={'EV_vs_B': ['EV', 'B']},
                                   R_script_dir=analysis_dir)    
    run_1.PCA('EV_vs_B')
    for path in result_paths:
        deg_result = DESeq2results(path)
        for fold_thres in [0, 1]:
            print(path.replace('.txt', '_fold{}_sig.txt'.format(fold_thres)))
            deg_result.get_sig_DEGs(path.replace('.txt', '_fold{}_sig.txt'.format(fold_thres)),
                                                 fold_threshold=fold_thres)
            deg_result.heatmap(filename='heatmap_normalized_sig_genes_foldthres{}'.format(fold_thres), color_label='EV_vs_B')
            deg_result.PCA('EV_vs_B', filename='PCA_normalized_sig_genes_foldthres{}'.format(fold_thres))
            
            
    ''' Viral EV1: control vs bacterial BM1: control'''
    
    analysis_dir = '{}/data_paper/plan_EV1_C_vs_BM1_C'.format(maindir)
    try:
        os.mkdir(analysis_dir)
    except FileExistsError:
        pass
    
    lp = LabelProcessor(readcount_df)
    lp.add_samples('group', 'EV1')
    lp.add_samples('group', 'BM1')
    lp.regroup_meta_tags('sample_type', 'M_vs_C', {'M': ['M', ''], 'C': ['C']})
    lp.remove_samples('M_vs_C', 'M')
    
    lp.generate_DESeq_data('{}/read_counts.txt'.format(analysis_dir),
                           '{}/col_data.txt'.format(analysis_dir))
    
    R_script_dir = '{}/code/generated_R_scripts'.format(maindir)
    
    run_1 = DESeq2run('{}/read_counts.txt'.format(analysis_dir),
                      '{}/col_data.txt'.format(analysis_dir),
                      sample_column='tags')
    
    result_paths = run_1.calc_DEGs(design = 'run + group',
                                   to_test={'group': ['EV1', 'BM1']},
                                   R_script_dir=analysis_dir)    
    run_1.PCA('group')
    
    for path in result_paths:
        deg_result = DESeq2results(path)
        for fold_thres in [0, 1]:
            deg_result.get_sig_DEGs(path.replace('.txt', '_fold{}_sig.txt'.format(fold_thres)),
                                                 fold_threshold=fold_thres)
             
    
    ''' Viral EV1: matched meningitis vs control.
        doesn't include meningitis samples without a control sample '''
    
    analysis_dir = '{}/data_paper/plan_EV1_M_vs_C'.format(maindir)
    try:
        os.mkdir(analysis_dir)
    except FileExistsError:
        pass
    
    lp = LabelProcessor(readcount_df)
    lp.add_samples('group', 'EV1')
    lp.regroup_meta_tags('sample_type', 'M_vs_C', {'M': ['M', ''], 'C': ['C']})
    
    lp.generate_DESeq_data('{}/read_counts.txt'.format(analysis_dir),
                           '{}/col_data.txt'.format(analysis_dir))
    
    R_script_dir = '{}/code/generated_R_scripts'.format(maindir)
    
    run_1 = DESeq2run('{}/read_counts.txt'.format(analysis_dir),
                      '{}/col_data.txt'.format(analysis_dir),
                      sample_column='tags')
    
    result_paths = run_1.calc_DEGs(design = 'tags + M_vs_C',
                                   to_test={'M_vs_C': ['M', 'C']},
                                   R_script_dir=analysis_dir)    
    run_1.PCA('M_vs_C')
    
    for path in result_paths:
        deg_result = DESeq2results(path)
        for fold_thres in [0, 1]:
            deg_result.get_sig_DEGs(path.replace('.txt', '_fold{}_sig.txt'.format(fold_thres)),
                                                 fold_threshold=fold_thres)
            deg_result.heatmap(filename='heatmap_normalized_sig_genes_foldthres{}'.format(fold_thres),
                               label='M_vs_C', color_label='M_vs_C')
            deg_result.PCA('M_vs_C', filename='PCA_normalized_sig_genes_foldthres{}'.format(fold_thres))
    
    
    ''' Bacterial BM1: BM1 meningitis vs BM1 control.
        doesn't include meningitis samples without a control sample '''
    
    analysis_dir = '{}/data_paper/plan_BM1_M_vs_C'.format(maindir)
    try:
        os.mkdir(analysis_dir)
    except FileExistsError:
        pass
    
    lp = LabelProcessor(readcount_df)
    lp.add_samples('group', 'BM1')
    lp.remove_samples('sample', 'STA2')
    lp.remove_samples('sample', 'UZB11')
    lp.remove_samples('sample', 'STLUC1')
    lp.remove_samples('sample', 'UZB12')
    lp.regroup_meta_tags('sample_type', 'M_vs_C', {'M': ['M', ''], 'C': ['C']})
    
    lp.generate_DESeq_data('{}/read_counts.txt'.format(analysis_dir),
                           '{}/col_data.txt'.format(analysis_dir))
    
    R_script_dir = '{}/code/generated_R_scripts'.format(maindir)
    
    run_1 = DESeq2run('{}/read_counts.txt'.format(analysis_dir),
                      '{}/col_data.txt'.format(analysis_dir),
                      sample_column='tags')
    
    result_paths = run_1.calc_DEGs(design = 'tags + M_vs_C',
                                   to_test={'M_vs_C': ['M', 'C']},
                                   R_script_dir=analysis_dir)    
    run_1.PCA('M_vs_C')
    
    for path in result_paths:
        deg_result = DESeq2results(path)
        for fold_thres in [0, 1]:
            deg_result.get_sig_DEGs(path.replace('.txt', '_fold{}_sig.txt'.format(fold_thres)),
                                                 fold_threshold=fold_thres)
            deg_result.heatmap(filename='heatmap_normalized_sig_genes_foldthres{}'.format(fold_thres),
                               label='M_vs_C', color_label='M_vs_C')
            deg_result.PCA('M_vs_C', filename='PCA_normalized_sig_genes_foldthres{}'.format(fold_thres))
            
    
    ''' Bacterial BM1: BM1 meningitis vs BM1/EV1 control '''
    
    analysis_dir = '{}/data_paper/plan_BM1_M_vs_all_C'.format(maindir)
    try:
        os.mkdir(analysis_dir)
    except FileExistsError:
        pass
    
    lp = LabelProcessor(readcount_df)
    lp.add_samples('group', 'EV1')
    lp.remove_samples('sample_type', 'M')
    lp.remove_samples('sample_type', '')
    lp.add_samples('group', 'BM1')
    lp.regroup_meta_tags('sample_type', 'M_vs_C', {'M': ['M', ''], 'C': ['C']})
    
    lp.generate_DESeq_data('{}/read_counts.txt'.format(analysis_dir),
                           '{}/col_data.txt'.format(analysis_dir))
    
    R_script_dir = '{}/code/generated_R_scripts'.format(maindir)
    
    run_1 = DESeq2run('{}/read_counts.txt'.format(analysis_dir),
                      '{}/col_data.txt'.format(analysis_dir),
                      sample_column='tags')
    
    result_paths = run_1.calc_DEGs(design = 'run + M_vs_C',
                                   to_test={'M_vs_C': ['M', 'C']},
                                   R_script_dir=analysis_dir)    
    run_1.PCA('M_vs_C')
    
    for path in result_paths:
        deg_result = DESeq2results(path)
        for fold_thres in [0, 1]:
            deg_result.get_sig_DEGs(path.replace('.txt', '_fold{}_sig.txt'.format(fold_thres)),
                                                 fold_threshold=fold_thres)
            deg_result.heatmap(filename='heatmap_normalized_sig_genes_foldthres{}'.format(fold_thres),
                               label='group', color_label='group')
            deg_result.PCA('M_vs_C', filename='PCA_normalized_sig_genes_foldthres{}'.format(fold_thres))
            
            
    ''' Enterovirus EV: EVM versus everything else '''
    
    #EV M vs EV C, BM, BM C, VM, REU1, PN
    
    analysis_dir = '{}/data_paper/plan_EVM_vs_rest'.format(maindir)
    try:
        os.mkdir(analysis_dir)
    except FileExistsError:
        pass
    
    lp = LabelProcessor(readcount_df)
    lp.add_samples('group', 'EV1')
    lp.add_samples('group', 'EV2')
    lp.regroup_meta_tags('sample_type', 'EVM', {'EVM': ['M', ''], 'other': ['C']})
    for group_type in ['BM1', 'BM2', 'PN', 'VM', 'REU1']:
        lp.add_samples('group', group_type)
    lp.regroup_meta_tags('EVM', 'EVM2', {'EVM': ['EVM'], 'other': ['other', '']})
        
    lp.generate_DESeq_data('{}/read_counts.txt'.format(analysis_dir),
                           '{}/col_data.txt'.format(analysis_dir))
    
    R_script_dir = '{}/code/generated_R_scripts'.format(maindir)
    
    run = DESeq2run('{}/read_counts.txt'.format(analysis_dir),
                    '{}/col_data.txt'.format(analysis_dir),
                    sample_column='tags')
    
    result_paths = run.calc_DEGs(design = 'run + EVM2',
                                 to_test={'EVM2': ['EVM', 'other']},
                                 R_script_dir=analysis_dir)
    
    run.PCA('EVM2')
    
    for path in result_paths:
        deg_result = DESeq2results(path)
        for fold_thres in [0, 1]:
            deg_result.get_sig_DEGs(path.replace('.txt', '_fold{}_sig.txt'.format(fold_thres)),
                                                 fold_threshold=fold_thres)
            deg_result.heatmap(filename='heatmap_normalized_sig_genes_foldthres{}'.format(fold_thres),
                               label='group', color_label='EVM2')
            deg_result.PCA('EVM2', filename='PCA_normalized_sig_genes_foldthres{}'.format(fold_thres))