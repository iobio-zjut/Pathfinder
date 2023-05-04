"""
!/usr/bin/env python
 -*- coding: utf-8 -*-
 @Time    : 2022/3/22 19:22
 @Author  : hzh
 @File    : pre_stage2.py
"""
import copy

import csv
import multiprocessing
import os
import sys
import math

import numpy as np


def mv_file(c_name):
    """
    create file and get data
    :param c_name: path of cluster name defult = "./data"
    """
    os.system("mkdir " + c_name + "/cluster")
    index = os.listdir(c_name)
    index = [i for i in index if len(i) <= 2]
    c_nums = 1
    os.system("cp " + c_name + "/combo0.pdb " + c_name
              + "/cluster/combo0.pdb")
    for i in index:
        # get rst.dat
        os.system("cp " + c_name + "/" + i + "/rst.dat " + c_name
                  + "/cluster/rst.dat" + i)

        i_combo = os.listdir(c_name + "/" + i)
        for j in i_combo:
            # get combo information
            if j[:5] == 'combo' and j[-4:] == '.pdb':
                os.system("cp " + c_name + "/" + i +
                          "/" + j + " " + c_name +
                          "/cluster/combo" + str(c_nums) + ".pdb ")
                c_nums += 1


def TMscore(pdb_path, pdb1, pdb2):
    """
    compute pdb1 and pdb2 TMscore

    :param pdb_path: cluster pdb path
    :param pdb1: index of pdb1
    :param pdb2: index of pdb2
    :return: TMscore
    """

    # make sure tmscore file in current file
    tm_path = "./tools/TMscore"
    if not os.path.exists(tm_path):
        print("lack tmscore file")
        sys.exit(1)
    if os.path.exists(pdb_path + "/tmscore.log"):
        os.remove(pdb_path + "/tmscore.log")

    os.system("chmod 755 " + tm_path)
    # compute by c++ of TMscore
    os.system("nohup " + tm_path + " "
              + pdb_path + "/combo" + str(pdb1) + ".pdb" + " "
              + pdb_path + "/combo" + str(pdb2) + ".pdb"
              + " > " + pdb_path + "/tmscore.log 2>&1")

    # use nohup print info on log file
    with open(pdb_path + "/tmscore.log", 'r') as t_file:
        line = ''
        for i in t_file.readlines():
            if i and i != '\n' and i.split()[0] == 'TM-score':
                line = i.split()[2]
                break

        # TMscore file had remake to output only score value
        if line:
            os.remove(pdb_path + "/tmscore.log")
            return float(line)

    os.remove(pdb_path + "/tmscore.log")
    return -1


def get_combo(file_path):
    """
    get cluster information
    rst.dat contain the number of combo by spicker

    :param file_path: rst.dat path
    :return: [index,info]
    """
    # set the most value to native
    f_combo = [[0, 130000]]
    numbers = 1
    index = os.listdir(file_path)
    index = [i for i in index if len(i) <= 2]

    for n_rst in index:
        with open(file_path + "/cluster/rst.dat" + n_rst) as file:
            # get rst.dat about combo use number of spicker
            row = file.readlines()
            f_temp = os.listdir(file_path + "/" + n_rst)
            f_temp = [int(i[5:-4]) for i in f_temp if i[:5] == 'combo' and i[-4:] == '.pdb']
            eff = max(f_temp)
            lines = row[(21 + eff):(21 + 2 * eff)]

        # rename combo
        for data in lines:
            temp = [numbers, int(data.split()[1])]
            f_combo.append(temp)
            numbers += 1

    return f_combo


def max_cluster(c_data, c_score, c_path):
    """
    main cluster
    if TMscore > c_score, Use the larger number of clusters as the center

    :param c_path:
    :param c_data: combo information [index,number]
    :param c_score: cluster TM-score
    :return: index after cluster
    """
    pdb = []
    print("============================cluster=============================\n")
    print("all combo number: " + str(len(c_data)))
    # combo still be removed, use while avoid index problem
    while len(c_data):
        # use temp list compute TMscore to cluster
        temp_combo = c_data[:]

        for index in temp_combo:
            if index != c_data[0]:
                score = TMscore(c_path + "/cluster", c_data[0][0], index[0])

                if score > c_score:
                    c_data.remove(index)
                    # sum the seed number
                    c_data[0][1] += index[1]
                    # output the removed combo
                    with open(c_path + "/cluster/cluster.log", 'a')as t_f:
                        t_f.write(str(c_data[0][0]) + " and " + str(index[0]) + " " + str(score) + "\n")

        pdb.append(c_data[0])
        c_data.remove(c_data[0])

    # output the cluster combo
    pdb[0][1] -= 130000
    with open(c_path + "/cluster/cluster.log", 'a')as c_f:
        for i in range(len(pdb)):
            c_f.write(str(pdb[i][0]) + " numbers of cluster: " + str(pdb[i][1]) + "\n")
        c_f.write("=========================finish=========================\n")

    return pdb


def output_pdb(c_result, f_name, f_out):
    """
    output result
    :param c_result: result path
    :param f_name: result name
    :return: cluster numbers
    """
    os.system("mkdir " + f_out)
    p_nums = 0

    # cp the cluster combo to result
    for index in c_result:
        os.system("cp " + f_name + "combo" + str(index[0]) +
                  ".pdb " + f_out + "combo" + str(p_nums + 1) + ".pdb")
        p_nums += 1
    return p_nums


def cluster(p, s):
    """
    cluster main

    :param p: path of file
    :param s: TMscore of two combo, the parameter of cluster
    :return: no return
    """

    path = p + "data"
    # file initial and get data
    if not os.path.exists(p + 'cluster'):
        mv_file(path)

    # get combo data and rst data
    combo = get_combo(path)

    # sort combo by number of spicker
    combo.sort(key=(lambda x: x[1]), reverse=True)

    # Prevent multiple runs from overwriting results
    if os.path.exists(path + "/e_result"):
        os.system("rm -r " + path + "/e_result")

    # cluster and output
    nums = output_pdb(max_cluster(combo, s, path),
                      path + "/cluster/",
                      path + "/e_result/")

    rename_seed(p)

    # output stage2 info
    if nums > 15:
        update_par(p, 15)
    update_par(p, nums)


def update_par(p, nums):
    replace = "Seed_number" + "\t" + str(nums)
    data = []
    flag = False
    if not os.listdir(p + "data/e_result"):
        print("Seed number is zero!!!")
        sys.exit(1)

    with open(p + "tools/parameter", 'r+') as f1:
        data = f1.readlines()
        for i in range(len(data)):
            if data[i] and data[i] != '\n' and data[i].split()[0] == "Seed_number":
                data[i] = replace
                flag = True
    if flag:
        with open(p + "tools/parameter", 'w+') as f2:
            f2.writelines(data)
    else:
        with open(p + "tools/parameter", 'a') as f2:
            f2.write("\n" + replace)
    print("Seed numbers: " + str(nums))


def rename_seed(p):
    """
    mkdir stage2 file and rename seed

    :param p:path of file
    :return: no return
    """
    seed_data = os.listdir(p + "data/e_result")
    seed_data = [i for i in seed_data if i[:5] == 'combo' and i[-4:] == '.pdb']

    if os.path.exists(p + "HMM_data"):
        os.system("rm -r " + p + "HMM_data")
    os.system("mkdir " + p + "HMM_data")
    for i in seed_data:
        os.system("cp " + p + "data/e_result/"
                  + i + " " + p + "HMM_data/seed_"
                  + i[5:-4] + ".pdb")


def spicker(index):
    """
    spicker algorithm, Serial port.
    You can open comments if you want to parallelize.

    :param index: path of file
    :return: no return
    """
    os.chdir("./" + index)
    os.system("./spicker")
    os.chdir("../")

    print(index + " finished spicker")


def multi_spicker(p):
    """
    spicker algorithm, Serial port.
    You can open comments if you want to parallelize.

    :param p:path of file
    :return: no return
    """

    file_list = os.listdir(p + "data")
    file_list = len([i for i in file_list if len(i) <= 2])
    os.chdir(p + "data")
    index = 1
    last = file_list % 5

    while index < file_list - last + 1:
        p1 = multiprocessing.Process(target=spicker, args=(str(index),))
        p2 = multiprocessing.Process(target=spicker, args=(str(index + 1),))
        p3 = multiprocessing.Process(target=spicker, args=(str(index + 2),))
        p4 = multiprocessing.Process(target=spicker, args=(str(index + 3),))
        p5 = multiprocessing.Process(target=spicker, args=(str(index + 4),))

        p1.start()
        p2.start()
        p3.start()
        p4.start()
        p5.start()

        p1.join()
        p2.join()
        p3.join()
        p4.join()
        p5.join()

        index += 5

    while index < file_list + 1:
        spicker(str(index))
        index += 1

    os.chdir("../")


def file_to_data(data_path, data_num):
    """
    :param data_path: file path
    :keyword Convert documents into usable data
             and convert the input file to a correspondence matrix and save it.
    :param data_num: Final number of iterations
    :return: state_data: state transition matrix
             init_data: initial state probability
    """
    state_path = data_path + "/state"

    # 返回数据
    state_data = []
    flag = False
    with open(state_path) as f:
        read = csv.reader(f)
        for row in read:
            test_flag = row[0].split()
            if flag:
                # 对每一行数据暂存
                temp_data = []
                # 每一行数据求和
                sum_data = 0
                for i in range(len(test_flag) - 1):
                    temp_data.append(int(test_flag[i]))
                    sum_data += int(test_flag[i])
                if sum_data == 0:
                    state_data.append([0.001 for _ in range(len(test_flag) - 1)])
                else:
                    for i in range(len(test_flag) - 1):
                        # 对行求概率
                        temp_data[i] = float(test_flag[i]) / sum_data
                    state_data.append(temp_data)

            if test_flag[0] == "test_" + str(data_num):
                flag = True

    state_data = np.array(state_data)

    # 存储数据
    np.savetxt(data_path + "/state_tran.csv", state_data, delimiter=',', encoding="utf-8")

    return state_data


def viterbi(state, init_state, num):
    """
    :param state:  state transition matrix
    :param init_state: initial state probability
    :param num: length of path
    :return: state: folding pathway

    """
    # seed个数
    n_seed = len(init_state)
    # After transposition, the column data can be obtained directly through the index

    # initial probability
    prob = init_state

    # Since the file storage starts from 1, and the transfer matrix index starts from 0,
    # the following str(j+1) is the same
    path = [str(i) for i in range(1, n_seed + 1)]
    result = {}

    """
    维特比算法主体部分
        动态规划算法，减少时间复杂度,适合HMM模型预测最大路径
        初始状态应为native 从后往前推导
    """
    for i in range(num):

        temp_prob = [0] * n_seed
        temp_path = [''] * n_seed

        if np.all(np.array(state) == 0):
            break
        state = [i + 0.0001 for i in state]
        for j in range(n_seed):
            # Multiply the current probability by the probability of being transferred,
            # and only take the maximum value to get the probability of the current node
            temp = prob * state[j]
            temp_prob[j] = max(temp)

            # Get the index to get the transfer path, and output
            index = np.argmax(temp)
            if max(state[j]) < 0.01:
                temp_path[j] = str(path[j])
                continue
            # One-way
            state[index] *= 0.01
            # state[j][index] *= 0.1

            # Previous node data -> current node
            temp_path[j] = str(j + 1) + "->" + str(path[index])

        path = temp_path
        prob = temp_prob

    # Store the result in the dictionary (remove redundant data),
    # the key is the path, and the value is the maximum probability of the path
    for name, i in zip(path, range(n_seed)):

        # Multiple paths with the same name take the maximum value
        if name in result:
            result[name] = max(prob[i], result[name])
            result[name] = result[name]

        # Paths with probability 0 are not counted
        elif prob[i] == 0:
            continue
        elif len(path[i]) == 1:
            continue
        else:
            result[name] = prob[i]

    # Output Probability Normalization
    sum_result = sum(result.values())
    for i in result.keys():
        result[i] /= sum_result
        result[i] = round(result[i], 2)

    # Sort by probability
    result = sorted(result.items(), key=lambda x: x[1], reverse=True)

    return result[:5]


def read_pdb(f_path):
    """
        read .pdb file get CA,N,C to rebuild CB coordination
        CB is used to calculate contact

    :param f_path: path of pdb
    :return:
    """
    pdb_N = []
    pdb_CA = []
    pdb_C = []
    with open(f_path, 'r') as f1:
        reader = csv.reader(f1)
        for row in reader:
            row = row[0]
            if row.split()[0] == 'TER':
                break

            if row.split()[2] == 'N':
                temp = [float(row[31:38]),
                        float(row[39:46]),
                        float(row[47:54])]
                pdb_N.append(temp)

            if row.split()[2] == 'CA':
                temp = [float(row[31:38]),
                        float(row[39:46]),
                        float(row[47:54])]
                pdb_CA.append(temp)

            if row.split()[2] == 'C':
                temp = [float(row[31:38]),
                        float(row[39:46]),
                        float(row[47:54])]
                pdb_C.append(temp)

    return pdb_N, pdb_CA, pdb_C


def Contact_Order(pdb_path):
    """
        Contact order equation : (1 / N*L) * sum(sequence separation)
        sequence separation is the index of i and j distance in sequence

        In original paper contact use CB,and use CA insteal of CB,when
        there is no CB. So there is a little different from the data of
        the paper.

        'Plaxco K W, Simons K T, Baker D. Contact order,
        transition state placement and the refolding rates of single domain proteins[J].
        Journal of molecular biology, 1998, 277(4): 985-994.'

    :param pdb_N:
    :param pdb_CA:
    :param pdb_C:
    :return: CO: contact order
    """

    pdb_N, pdb_CA, pdb_C = read_pdb(pdb_path)

    ca = -0.58273431
    cb = 0.56802827
    cc = -0.54067466

    L = len(pdb_CA)
    deltaS = 0
    N = 0
    pdb_CB = []
    res_CO = []

    # rebuild CB
    for i in range(L):
        b = [pdb_CA[i][0] - pdb_N[i][0],
             pdb_CA[i][1] - pdb_N[i][1],
             pdb_CA[i][2] - pdb_N[i][2]]
        c = [pdb_C[i][0] - pdb_CA[i][0],
             pdb_C[i][1] - pdb_CA[i][1],
             pdb_C[i][2] - pdb_CA[i][2]]
        a = [b[1] * c[2] - b[2] * c[1],
             b[2] * c[0] - b[0] * c[2],
             b[0] * c[1] - b[1] * c[0]]
        temp = [ca * a[0] + cb * b[0] + cc * c[0] + pdb_CA[i][0],
                ca * a[1] + cb * b[1] + cc * c[1] + pdb_CA[i][1],
                ca * a[2] + cb * b[2] + cc * c[2] + pdb_CA[i][2]]
        pdb_CB.append(temp)

    for i in range(L):
        temp = 0
        for j in range(L):
            if i == j:
                continue
            x = pdb_CB[i][0] - pdb_CB[j][0]
            y = pdb_CB[i][1] - pdb_CB[j][1]
            z = pdb_CB[i][2] - pdb_CB[j][2]

            dist = math.sqrt(x * x + y * y + z * z)

            # if < 8A , contact
            if dist < 8:
                N += 1
                deltaS += abs(i - j)

            if dist < 20:
                temp += abs(i - j) / dist
        res_CO.append(temp)

    return deltaS / (L * N), res_CO
