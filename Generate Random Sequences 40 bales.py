# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 13:42:03 2021

@author: Dahui
"""

import random
import copy

sequenceListO = [[0,0],[0,0],[0,0],[1,0],[1,0],[1,0],[1,0],[1,0],[2,0],[2,0],[0,0],[0,0],[0,0],[1,0],[1,0],[1,0],[1,0],[1,0],[2,0],[2,0],\
                 [0,0],[0,0],[0,0],[1,0],[1,0],[1,0],[1,0],[1,0],[2,0],[2,0],[0,0],[0,0],[0,0],[1,0],[1,0],[1,0],[1,0],[1,0],[2,0],[2,0]]
randomSeed = 0
sequenceList = []
for i in range(5):
    sequenceListOC = copy.deepcopy(sequenceListO)
    sequenceListR = []
    indexRange = 39
    while indexRange >= 0:
        random.seed(randomSeed)
        locIndex = random.randint(0,indexRange)
        bale = sequenceListOC.pop(locIndex)
        sequenceListR.append(bale)
        indexRange -= 1
        randomSeed += 1
    sequenceList.append(sequenceListR)
for i in range(len(sequenceList)):
    print("Random Sequence " + str(i+1))
    print(sequenceList[i])
print("Random Sequence 1 Sort")
sequence1Sort = sorted(sequenceList[0], key=lambda x: x[1])
print(sequence1Sort)
print("Random Sequence 1 Sort1")
sequence1Sort1 = sorted(sequenceList[0], key=lambda x: x[0])
print(sequence1Sort1)
print(len(sequence1Sort))