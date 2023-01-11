# -*- coding: utf-8 -*-
"""Conform functional facename task behavior data during scan session."""

from pathlib import Path
# from pathsjson.automagic import PATHS as project
import numpy as np
import pandas as pd
import scipy.io as sio

# Get directories
bids_dir = Path('/brain/guixue/Jintao/geneproject/facename')

root_dir = Path('/brain/guixue/geneproject/fMRIPrep/Beijing')
src_dir = root_dir.joinpath('sourcedata', 'behav', 'all')

# Get subjects list
subj_info = pd.read_csv(bids_dir.joinpath('facename_subj_info_478.csv'))
subj_list = subj_info['subID'].tolist()
#subj_list = [i.replace('sub-', '') for i in subj_info['labelID'].tolist()]
#subj_list = ['2004'] # for test

# In[1]
# Convert beh data
for subj_id in subj_list:
    
    out_dir = bids_dir.joinpath(f'sub-{subj_id}', 'func')
    out_dir.mkdir(parents=True, exist_ok=True)

    # facename learning
    # Check file existence
    fid = list(src_dir.glob(f'fMRI_FN_learning_sub{subj_id}_*.mat'))
    if len(fid) != 1:
        print(f'Mismatching number of perceiving mat file for sub-{subj_id}')
        break

    # Generate dataframe
    raw_dat = sio.loadmat(fid[0])
    
    src = pd.DataFrame(
        raw_dat['FL'],
        columns=[
            'trial_id', 'face_id1', 'name_id', 'num_pres', 'lag', 'gender', 'response',
            'response_time', 'Aonset', 'Donset', 'name_id1','name_id2', 'name_id3', 
            'face_id'
        ])
    
    beh = pd.DataFrame(
        {
            'onset': src['Aonset'],
            'duration': [2.5] * src.shape[0],
            'trial_type': ['n/a'] * src.shape[0],
            'response_time': src['response_time'],
            'response_type': ['n/a'] * src.shape[0],
            'task_type': ['encoding'] * src.shape[0],
            'face_id': src['face_id'].astype('int'),
            'name_id': src['name_id'].astype('int'),
            'num_pres': src['num_pres'].astype('int'),
            'lag': src['lag'].astype('int'),
            'response': src['response'].astype('int'),
            'trial_id': src['trial_id'].astype('int'),
            'design_onset': src['Donset'],
            'grate_orientation': ['n/a'] * src.shape[0]
            },
            columns=[
                'onset', 'duration', 'trial_type', 'response_time', 'response_type', 
                'task_type', 'face_id','name_id', 'num_pres', 'lag', 'response', 
                'trial_id', 'design_onset', 'grate_orientation'
            ])
    beh.loc[beh['response'] == 1, 'trial_type'] = "fit"
    beh.loc[beh['response'] == 2, 'trial_type'] = "unfit"
    beh.loc[beh['response'] == 0, 'trial_type'] = "no response"
    beh.loc[beh['response'] == 0, 'response_time'] = np.NaN
    beh.loc[beh['response'] == 1, 'response_type'] = 'fit'
    beh.loc[beh['response'] == 2, 'response_type'] = 'unfit'
    beh.loc[beh['response'] == 0, 'response_type'] = 'no response'
    
    # red box
    red_box = pd.DataFrame(
        {
            'onset': src['Aonset'],
            'duration': src['response_time'],
            'trial_type': ['n/a'] * src.shape[0],
            'response_time': src['response_time'],
            'response_type': ['n/a'] * src.shape[0],
            'task_type': ['judgement'] * src.shape[0],
            'face_id': src['face_id'].astype('int'),
            'name_id': src['name_id'].astype('int'),
            'num_pres': src['num_pres'].astype('int'),
            'lag': src['lag'].astype('int'),
            'response': src['response'].astype('int'),
            'trial_id': src['trial_id'].astype('int'),
            'design_onset': ['n/a'] * src.shape[0],
            'grate_orientation': ['n/a'] * src.shape[0]
            },
            columns=[
                'onset', 'duration', 'trial_type', 'response_time', 'response_type', 
                'task_type', 'face_id','name_id', 'num_pres', 'lag', 'response', 
                'trial_id', 'design_onset', 'grate_orientation'
            ])
    red_box.loc[:, 'onset'] += 2.5 
    red_box.loc[red_box['response'] == 1, 'trial_type'] = "fit"
    red_box.loc[red_box['response'] == 2, 'trial_type'] = "unfit"
    red_box.loc[red_box['response'] == 0, 'trial_type'] = "no response"
    red_box.loc[red_box['response'] == 0, 'duration'] = 1.5
    red_box.loc[red_box['response'] == 0, 'response_time'] = np.NaN
    red_box.loc[red_box['response'] == 1, 'response_type'] = 'fit'
    red_box.loc[red_box['response'] == 2, 'response_type'] = 'unfit'
    red_box.loc[red_box['response'] == 0, 'response_type'] = 'no response'
    
    # Merge dataframes
    beh = beh.append(red_box)
    beh = beh.sort_values(by=['trial_id', 'onset'])
    
    # grate
    src = pd.DataFrame(
        raw_dat['grate'],
        columns=[
            'trial_id', 'orientation', 'response', 'response_time', 'onset2response', 
            'contrast'
        ])
    # delete first row with trial_id == 0
    src = src[src.trial_id != 0]
    grate = pd.DataFrame(
        {
            'onset': (src['onset2response'] - src['response_time']),
            'duration': src['response_time'],
            'trial_type': ['n/a'] * src.shape[0],
            'response_time': src['response_time'],
            'task_type': ['grate'] * src.shape[0],
            'response_type': ['n/a'] * src.shape[0],
            'response': src['response'].astype('int'),
            'trial_id': src['trial_id'].astype('int'),
            'face_id': np.nan,
            'name_id': np.nan,
            'design_onset': np.nan,
            'grate_orientation': src['orientation'].astype('int')
            },
            columns=[
                'onset', 'duration', 'trial_type', 'response_time', 'task_type', 
                'response_type', 'response', 'trial_id', 'face_id', 'name_id', 
                'design_onset', 'grate_orientation'
            ])
    grate.loc[:, 'duration'] -= 0.5  # fix a bug in grating code
    grate.loc[grate['duration'] < 0, 'duration'] = 0.1
    grate.loc[:, 'response_time'] = grate.loc[:, 'duration']
    grate.loc[grate['response'] == 1, 'trial_type'] = 'correct'
    grate.loc[grate['response'] == -1, 'trial_type'] = 'incorrect'
    grate.loc[grate['grate_orientation'] == 1, 'grate_orientation'] = 'left'
    grate.loc[grate['grate_orientation'] == 2, 'grate_orientation'] = 'right'
    
    # Merge dataframes
    beh = beh.append(grate)
    beh = beh.sort_values(by=['trial_id', 'onset'])
    # In[2]
    # Save file
    out_file = out_dir.joinpath(
        f'sub-{subj_id}_task-facename_events.tsv')
    beh.to_csv(out_file, sep='\t', header=True, index=False, na_rep='n/a', float_format='%.3f')
    # In[3]
    # Store memory status for retrieval
    tmp = []
    tmp.append(beh.loc[beh.task_type == 'encoding', ['face_id', 'num_pres', 'response_type']])
    # Merge subsequent memroy table
    tmp = pd.concat(tmp)
    # Reserving subjests' twice response 
    sub_mem1 = tmp[tmp.num_pres == 1]
    sub_mem1 = sub_mem1.drop('num_pres', axis = 1)
    sub_mem1.rename(columns={'response_type': 'subs_response_first'}, inplace=True)
    sub_mem2 = tmp[tmp.num_pres == 2]
    sub_mem2 = sub_mem2.drop('num_pres', axis = 1)
    sub_mem2.rename(columns={'response_type': 'subs_response_second'}, inplace=True)
    sub_mem = pd.merge(sub_mem1, sub_mem2, on='face_id')
    
    # Clear variables
    del fid, src, beh, out_file, grate, tmp, sub_mem1, sub_mem2
    
    # In[4]
    # Retrieval
    # Check file existence
    fid = list(src_dir.glob(f'fMRI_FN_test_sub{subj_id}_*.mat'))
#    if len(fid) != 1:
#        print(f'Mismatching number of perceiving tsv file for sub-{subj_id}')
#        break

    # Generate dataframe
    raw_dat = sio.loadmat(fid[0])
    
    src = pd.DataFrame(
        raw_dat['FN'],
        columns=[
            'trial_id', 'face_id', 'old_new', 'right_name', 'right_response', 'response',
            'response_time', 'Donset', 'Aonset', 'name_id1','name_id2', 'name_id3'
        ])
    src = pd.merge(src, sub_mem, on='face_id', how='outer')
    beh = pd.DataFrame(
        {
            'onset': src['Aonset'].values,
            'duration': src['response_time'],
            'trial_type': ['n/a'] * src.shape[0],
            'response_time': src['response_time'],
            'response_type': (src['response'] - src['right_response']).astype('int'),
            'subs_response_first': src['subs_response_first'],
            'subs_response_second': src['subs_response_second'],
            'task_type': ['retrieval'] * src.shape[0],
            'response': src['response'],
            'right_response': src['right_response'],
            'trial_id': src['trial_id'],
            'design_onset': src['Donset'],
            'face_id': (src['face_id']).astype('int'),
            'old_new': src['old_new']
        },
        columns=[
                'onset', 'duration', 'trial_type', 'response_time', 'response_type', 
                'subs_response_first', 'subs_response_second', 'task_type', 'response', 
                'right_response', 'trial_id', 'design_onset', 'face_id', 'old_new'
            ])
    
    beh.loc[beh['response_type'] == 0, 'trial_type'] = "correct"
    beh.loc[beh['response_type'] != 0, 'trial_type'] = "incorrect"
    beh.loc[beh['response'] == 0, 'trial_type'] = "no response"
    beh.loc[beh['response'] == 0, 'response_time'] = np.NaN
    beh.loc[beh['response'] == 0, 'duration'] = 4
    beh.loc[beh['response_type'] != 0, 'response_type'] = -1
    beh.loc[beh['response_type'] == 0, 'response_type'] = 1
    beh.loc[beh['response_type'] == -1, 'response_type'] = 0
    beh.loc[beh['old_new'] == 2, 'subs_response_first'] = 'new'
    beh.loc[beh['old_new'] == 2, 'subs_response_second'] = 'new'
    beh.loc[beh['old_new'] == 1, 'old_new'] = 'old'
    beh.loc[beh['old_new'] == 2, 'old_new'] = 'new'
    
    # In[4]
    # Fixation
    # Generate dataframe
    fixation = pd.DataFrame(
        {
            'onset': (src['Aonset'] + src['response_time']),
            'duration': ['n/a'] * src.shape[0],
            'trial_type': ['fixation'] * src.shape[0],
            'response_time': ['n/a'] * src.shape[0],
            'response_type': ['n/a'] * src.shape[0],
            'subs_response_first': ['n/a'] * src.shape[0],
            'subs_response_second': ['n/a'] * src.shape[0],
            'task_type': ['fixation'] * src.shape[0],
            'response': ['n/a'] * src.shape[0],
            'right_response': ['n/a'] * src.shape[0],
            'trial_id': src['trial_id'],
            'design_onset': ['n/a'] * src.shape[0],
            'face_id': ['n/a'] * src.shape[0],
            'old_new': ['n/a'] * src.shape[0],
            'Tonset': abs(beh['onset'].diff(periods=-1))
        },
        columns=[
                'onset', 'duration', 'trial_type', 'response_time', 'response_type', 
                'subs_response_first', 'subs_response_second', 'task_type', 'response', 
                'right_response', 'trial_id', 'design_onset', 'face_id', 'old_new', 'Tonset'
            ])
    lc = 153*2 - beh.onset[-1:] # compute duration for the last row, 153 volumes total
    # replace the last row
    fixation.loc[(src.shape[0] - 1), 'Tonset'] = lc[beh.shape[0] - 1]
    # In[5]
    fixation.loc[:, 'duration'] = (fixation.loc[:, 'Tonset'] - beh.loc[:, 'duration'])
    fixation = fixation.drop('Tonset', axis = 1)
    
    # Merge dataframes
    beh = beh.append(fixation)
    beh = beh.sort_values(by=['trial_id', 'onset'])
    
    # In[6]
    # Save file
    out_file = out_dir.joinpath(
        f'sub-{subj_id}_task-facetest_events.tsv')
    beh.to_csv(out_file, sep='\t', header=True, index=False, float_format='%.3f')

    # Clear variables
    #del fid, src, beh, out_file
