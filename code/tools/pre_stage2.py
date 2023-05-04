"""
!/usr/bin/env python
 -*- coding: utf-8 -*-
 @Time    : 2022/4/15 17:41
 @Author  : hzh
 @File    : pre_stage2.py
"""
import os
import argparse
from utils import cluster, rename_seed, multi_spicker

if __name__ == '__main__':
    # add argument
    parser = argparse.ArgumentParser(description="cluster combo.")
    parser.add_argument('-s', type=float, default=0.6, help='cluster TMscore.')
    parser.add_argument('-p', type=str, default='./', help='path of working directory.')
    args = parser.parse_args()

    print("==========================spicker=========================")
    multi_spicker(args.p)

    print("==========================cluster=========================")
    cluster(args.p, args.s)

    print("=======================update parameter========================")
