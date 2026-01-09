"""
Compute distance score between original and SCRaMbLEd melodies.
"""

import os
import pretty_midi
import numpy as np
from scipy.stats import wasserstein_distance

def extract_features(
    midi_data
):
    """
    Extract infos from MIDI data
    
    Parameter:
        midi_data: audio data parsed by pretty_midi.PrettyMIDI
    
    Return:
        notes: list containing note infos, list
    """
    notes = []
    for instrument in midi_data.instruments:
        for note in instrument.notes:
            notes.append({
                'pitch': note.pitch,
                'start': note.start,
                'end': note.end,
                'duration': note.end - note.start,
                'velocity': note.velocity
            })
    return notes

def edit_distance(list1, list2):
    """
    Calculate Levenshtein distance between two lists
    
    Parameters:
        list1: input list 1 containing infos
        list2: input list 2 containing infos
    
    Return:
        dp[m][n]: edit distance score between list1[:m] and list2[:n]
    """
    m, n = len(list1), len(list2)
    
    #  Initialize edit distance list
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j
    
    #  Fill edit distance list
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if list1[i-1] == list2[j-1]:
                # Identical element -> unchange
                dp[i][j] = dp[i-1][j-1]
            else:
                # Different element -> get minimal value of deletion, insertion or replacement
                dp[i][j] = min(
                    dp[i-1][j] + 1,    #  Delete list1[i-1]
                    dp[i][j-1] + 1,    #  Insert list2[j-1]
                    dp[i-1][j-1] + 1   #  Replace list1[i-1] with list2[j-1]
                )
    
    return dp[m][n]

#  Default parameters
original_midi = 'Canon.synIXR.150_bit_per_note.mid'
original_note_info = 'Canon.synIXR.150_bit_per_note.txt'
SCRaMbLEd_midi_dir = 'synIXR.SCRaMbLE.melody'
SCRaMbLEd_note_info_dir = 'synIXR.SCRaMbLE.note_info'


### Edit distance
#----------------------------------------------------------------------------------------------------------------
#  Parse MIDI audio
midi_orig = pretty_midi.PrettyMIDI(original_midi)
notes_orig = extract_features(midi_orig)
pitch_seq_orig = [(n['pitch'], n['duration']) for n in notes_orig]

with open('SCRaMbLEd_Edit_distance.txt', 'w') as hd:
    for file in os.listdir(SCRaMbLEd_midi_dir):
        #  Parse MIDI audios of SCRaMbLEd strains
        midi_shuff = pretty_midi.PrettyMIDI(SCRaMbLEd_midi_dir + '/' + file)
        notes_shuff = extract_features(midi_shuff)
        pitch_seq_shuff = [(n['pitch'], n['duration']) for n in notes_shuff]
        
        #  Calculate distance score
        distance = edit_distance(pitch_seq_orig, pitch_seq_shuff)
        hd.write(file.split('.')[0] + '\t' + str(distance) + '\n')
#----------------------------------------------------------------------------------------------------------------


### New melodic fragments
#----------------------------------------------------------------------------------------------------------------
#  Read note infos and split to 16-note melodic fragments
with open(original_note_info, 'r') as f:
    num_lists_ori = []
    num_list_ful = f.read().splitlines()
    for i in range(len(num_list_ful) - 16 + 1):
        num_lists_ori.append(num_list_ful[i:i+16])

#  Calculate number of new 16-note melodic fragments
with open('SCRaMbLEd_new_melody_fragment_count.txt', 'w') as hd:
    for file in os.listdir(SCRaMbLEd_note_info_dir):
        prefix = file.split('.')[0]
        new_frag = 0
        with open(SCRaMbLEd_note_info_dir + '/' + file, 'r') as f:
            num_lists = []
            num_list_ful = f.read().splitlines()
            for i in range(len(num_list_ful) - 16 + 1):
                num_lists.append(num_list_ful[i:i+16])
        for num_list in num_lists:
            if not num_list in num_lists_ori:
                new_frag += 1
        hd.write(prefix + '\t' + str(new_frag) + '\n')
#----------------------------------------------------------------------------------------------------------------


### Pitch difference
#----------------------------------------------------------------------------------------------------------------
#  Parse MIDI audio
midi_orig = pretty_midi.PrettyMIDI(original_midi)
notes_orig = extract_features(midi_orig)
pitches_orig = sorted(np.array([n['pitch'] for n in notes_orig]))

with open('SCRaMbLEd_pitch_dist.txt', 'w') as hd:
    for file in os.listdir(SCRaMbLEd_midi_dir):
        #  Parse MIDI audios of SCRaMbLEd strains
        midi_shuff = pretty_midi.PrettyMIDI(SCRaMbLEd_midi_dir + '/' + file)
        notes_shuff = extract_features(midi_shuff)
        pitches_shuff = sorted(np.array([n['pitch'] for n in notes_shuff]))
        
        #  Calculate distance score
        dist = wasserstein_distance(pitches_orig, pitches_shuff)
        hd.write(file.split('.')[0] + '\t' + str(dist) + '\n')
#----------------------------------------------------------------------------------------------------------------


### Duration difference
#----------------------------------------------------------------------------------------------------------------
#  Parse MIDI audio
midi_orig = pretty_midi.PrettyMIDI(original_midi)
notes_orig = extract_features(midi_orig)
duration_orig = sorted(np.array([n['duration'] for n in notes_orig]))

with open('SCRaMbLEd_duration_dist.txt', 'w') as hd:
    for file in os.listdir(SCRaMbLEd_midi_dir):
        #  Parse MIDI audios of SCRaMbLEd strains
        midi_shuff = pretty_midi.PrettyMIDI(SCRaMbLEd_midi_dir + '/' + file)
        notes_shuff = extract_features(midi_shuff)
        duration_orig_shuff = sorted(np.array([n['duration'] for n in notes_shuff]))
        
        #  Calculate distance score
        dist = wasserstein_distance(duration_orig, duration_orig_shuff)
        hd.write(file.split('.')[0] + '\t' + str(dist) + '\n')
#----------------------------------------------------------------------------------------------------------------
