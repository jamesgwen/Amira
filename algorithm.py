import pandas as pd
import json
import re
import numpy as np
import ast
import itertools
from enum import IntEnum

#path = "~/amira/amira_rr_mock_data/"
#data = pd.read_table("/users/bag/clme1992/amira/slassify-error.tsv", sep = "\t")

def remove_val_element(lst, val):
    result = lst
    if val in lst:
        result.remove(val)
    return result

def remove_by_indices(lst, indices):
    for ind in sorted(indices, reverse=True):
        del lst[ind]
    return(lst)

def classify_error(txt_align, tran_align):
    if tran_align == txt_align:
        return "correct", None
    else: # there are errors.. may have multiple
        error_record = []
        ind = 0
        while ind <= len(tran_align) - 1:
            tran_word = tran_align[ind]
            txt_word = txt_align[ind]
            if tran_word == "-":
                print("tran gap...")
                error_record.append(["error - skip", ind])
                ind = ind + 1
                continue
            elif txt_word == "-": # either repeat or self correction
                print("txt gap...")
                # find out where the gap ends
                ind_start = ind
                txt_gap_inds = []
                while txt_word == "-":
                    txt_gap_inds.append(ind)
                    ind = ind + 1
                    if (ind > len(tran_align) - 1):
                        break
                    txt_word = txt_align[ind]
                insert = [tran_align[i] for i in txt_gap_inds]
                control_word = ["." for item in range(len(insert))]
                tran_wo_insert = remove_by_indices(tran_align, txt_gap_inds)
                txt_score = needleman_wunsch(insert, remove_val_element(txt_align, "-"))[2]
                tran_score = needleman_wunsch(insert, remove_val_element(tran_align, "-"))[2]
                control_txt_score = needleman_wunsch(control_word, remove_val_element(txt_align, "-"))[2]
                control_tran_score = needleman_wunsch(control_word, remove_val_element(tran_wo_insert, "-"))[2]
                if (txt_score <= control_txt_score and tran_score <= control_tran_score):
                    error_record.append(["error - random words", ind_start])
                elif (txt_score > tran_score):
                    error_record.append(["error - self correct", ind_start])
                else:
                    error_record.append(["error - repeat", ind_start])
                continue
            elif txt_word != tran_word: # either substitution or reverse
                print("mismatch...")
                # find out where the mismatch end
                ind_start = ind
                mismatch_inds = []
                while txt_word != tran_word and txt_word != "-" and tran_word != "-":
                    mismatch_inds.append(ind)
                    ind = ind + 1
                    if (ind > len(tran_align) - 1):
                        break
                    print(ind)
                    txt_word = txt_align[ind]
                    tran_word = tran_align[ind]
                mismatch_tran = [tran_align[i] for i in mismatch_inds]
                mismatch_txt = [txt_align[i] for i in mismatch_inds]
                if len(mismatch_tran) == 1:
                    error_record.append(["error - substitute", ind_start])
                else:
                    is_rearrange = 0
                    for new_tran_mismatch in itertools.permutations(mismatch_tran):
                        new_tran_mismatch = list(new_tran_mismatch)
                        if (new_tran_mismatch == mismatch_txt):
                            error_record.append(["error - rearrange", ind_start])
                            is_rearrange = 1
                            break
                    if is_rearrange == 0:
                        error_record.append(["error - substitute", ind_start])
                        ind = ind_start
            ind = ind + 1
        return(error_record)