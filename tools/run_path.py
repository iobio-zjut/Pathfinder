"""
!/usr/bin/env python
 -*- coding: utf-8 -*-
 @Time    : 2021/12/10 16:35
 @Author  : hzh
 @File    : run_path.py
"""

import argparse
import os

import numpy as np
import csv
from utils import file_to_data, viterbi

if __name__ == '__main__':

    # add argument
    parser = argparse.ArgumentParser(description="calculate folding pathway.")
    parser.add_argument('-s', type=float, default=0.5, help='cluster TMscore.')
    parser.add_argument('-p', type=str, default='./', help='path of working directory.')
    parser.add_argument('-t', type=float, default=100, help='trajectory of stage2.')
    args = parser.parse_args()

    file_path = args.p
    trajectory = args.t

    # Get transition matrix and initial state probabilities
    state_tran = file_to_data(file_path, trajectory)
    init_state = [0 for i in range(len(state_tran))]
    init_state[0] = 1

    # Calculate the maximum path through the Viterbi algorithm and save it
    result = viterbi(state_tran, init_state, 10)
    save = file_path + "/pathway.txt"

    # Store the file, save is the incoming path
    with open(save, 'w', encoding='utf-8') as f:

        f.write("********************result****************************\n")
        f.write("*                                       *                                     *\n")
        f.write("*           Pathway               *          Probability         *\n")
        f.write("*                                       *                                     *\n")
        f.write("******************************************************\n")
        for i in range(len(result)):
            save_data = str(result[i])[1:-1]
            f.write(save_data + "\n")

    print(result)
