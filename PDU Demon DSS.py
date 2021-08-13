#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 11:15:12 2020

@author: dahuiliu
"""

from gurobipy import *
import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#import math

def DataInput():
    #Parameter Define
    parameterDictionary = { #For blending process, equipment index binds to biomass type
                           'SLengthWidthHeight' : [], #l,w,h
                           'DDensity' : [], #d_sb
                           'TAverageDensity' : [], #d_isb
                           'DAverageDensityBlend' : [], #d_is for blending process
                           'TEquipmentCapacity' : [], #x_isb
                           'DEquipmentCapacityBlend' : [], #x_is for blending process
                           'DNumberOfBales' : [], #n_sb
                           'SDryMatterLoss' : [], #mu_i
                           'DBypassRatio' : [], #theta_sb
                           'SMassCapacityInventory' : [], #m_i
                           'SMassCapacityInventoryBlend' : [], #m_i for blending process
                           'SVolumnCapacityInventory' : [], #v_i
                           'SVolumnCapacityInventoryBlend' : [], #v_i for blending process
                           'DProcessingTimeOfBale' : [], #p_sb
                           'SAshContent' : [], #a_b
                           'SThermalContent' : [], #e_b
                           'SCarbohydrateContent' : [], #f_b
                           'STargetAshThermalCarbohydrate' : [], #a^*, e^*, f^*
                           'SReactorUpperLowerBigM' : []} #R^, R_, M

    #Read Parameter Data
    file = open(r'/Users/dahuiliu/Documents/BBP reactor feed rate/PDU Demon Data Basic.txt')
    listInput = file.readlines()
    tempList = []
    for i in range(len(listInput)):
        if listInput[i][0] == 'S' or listInput[i][0] == 'D' or listInput[i][0] == 'T':
            keyFirst = listInput[i][0] #Read dictionary keys
            keyName = listInput[i].strip()
        else:
            if keyFirst == 'S': #Read 1D list
                rowList = listInput[i].split()
                for j in range(len(rowList)):
                    parameterDictionary[keyName].append(rowList[j])
            elif keyFirst == 'D': #Read 2D list
                parameterDictionary[keyName].append(listInput[i].split())
            else: #Read 3D list
                if listInput[i][0] == 'N':
                    parameterDictionary[keyName].append(copy.deepcopy(tempList))
                    tempList[:] = []
                else:
                    tempList.append(listInput[i].split())
                #print(tempList)
    file.close()
    return parameterDictionary

def SequenceInitial(sequenceList, moistureIndex, biomassIndex, timeIndex, parameterDictionary):
    startAtTime = []
    for s in range(moistureIndex):
        startAtTimeRowCol = []
        for b in range(biomassIndex):
            startAtTimeRow = []
            for t in range(timeIndex):
                startAtTimeRow.append(0)
            startAtTimeRowCol.append(startAtTimeRow)
        startAtTime.append(startAtTimeRowCol)
    currentTime = 0
    for i in sequenceList:
        startAtTime[i[0]][i[1]][currentTime] = 1
        currentTime += int(parameterDictionary['DProcessingTimeOfBale'][i[0]][i[1]])
    return startAtTime

def MILPModel(equipmentIndex, moistureIndex, biomassIndex, equipmentIndexBlend, timeIndex, grinderSet, meteringbinSet, storeageSet, equipmentRemain, previousEquipments, previousEquipmentsBlend, previousEquipmentsReactor, parameterDictionary, volumnOfBale, startAtTime):
    #Define Model
    model = Model()
    #model.Params.MIPGap=0.0005

    #Define Variables
    outflowOfEquipment = model.addVars(equipmentIndex+1, moistureIndex, biomassIndex, timeIndex, lb=0.0, vtype=GRB.CONTINUOUS) #X_isbt
    outflowOfEquipmentBlend = model.addVars(equipmentIndexBlend, moistureIndex, timeIndex, lb=0.0, vtype=GRB.CONTINUOUS) #X_ist for blending process
    inventoryLevel = model.addVars(meteringbinSet, moistureIndex, biomassIndex, timeIndex, lb=0.0, vtype=GRB.CONTINUOUS) #M_isbt
    inventoryLevelBlend = model.addVars(storeageSet, moistureIndex, timeIndex, lb=0.0, vtype=GRB.CONTINUOUS) #M_ist for blending process
    velocityInFlow = model.addVars([0], moistureIndex, biomassIndex, timeIndex, lb=0.0, vtype=GRB.CONTINUOUS) #V_0sbt
    processAtTime = model.addVars(moistureIndex, biomassIndex, timeIndex, vtype=GRB.BINARY) #Z_sbt
    processAtTimeO = model.addVars(timeIndex, vtype=GRB.BINARY) #K_t
    upperReactor = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS)
    upperReactorL = model.addVars(timeIndex, lb=0.0, vtype=GRB.CONTINUOUS)
    #Non-negative constraints are also defined here

    #Set Objective
    objective = quicksum(processAtTime[s,b,t] for s in range(moistureIndex) for b in range(biomassIndex) for t in range(timeIndex)) + quicksum(processAtTimeO[t] for t in range(timeIndex))
    model.setObjective(objective, GRB.MINIMIZE)
    
    #Set Constraints
    #Capacity

    model.addConstrs(outflowOfEquipment[1,s,b,t] <= float(parameterDictionary['TEquipmentCapacity'][0][s][b]) * processAtTime[s,b,t] for s in range(moistureIndex) for b in range(biomassIndex) for t in range(timeIndex))
    model.addConstrs(outflowOfEquipment[i,s,b,t] <= float(parameterDictionary['TEquipmentCapacity'][i-1][s][b]) for i in grinderSet for s in range(moistureIndex) for b in range(biomassIndex) for t in range(timeIndex))
    model.addConstrs(quicksum(outflowOfEquipment[12,s,b,t] for s in range(moistureIndex) for b in range(biomassIndex)) <= quicksum(float(parameterDictionary['TEquipmentCapacity'][11][s][b]) * processAtTime[s,b,t] for s in range(moistureIndex) for b in range(biomassIndex)) for t in range(timeIndex))
    model.addConstrs(quicksum(outflowOfEquipmentBlend[i,s,t] for s in range(moistureIndex)) <= float(parameterDictionary['DEquipmentCapacityBlend'][i][0]) * processAtTimeO[t] for i in previousEquipmentsReactor for t in range(timeIndex)) #Equipment flow capacity for blending process

    model.addConstrs(quicksum(inventoryLevel[i,s,b,t] for s in range(moistureIndex) for b in range(biomassIndex)) <= float(parameterDictionary['SMassCapacityInventory'][0]) for i in meteringbinSet for t in range(timeIndex)) #Metering bin mass capacity. Index 0 is for metering bin
    model.addConstrs(quicksum(inventoryLevelBlend[i,s,t] for s in range(moistureIndex)) <= float(parameterDictionary['SMassCapacityInventoryBlend'][int((i-1)/3)]) for i in storeageSet for t in range(timeIndex)) #Storage bins mass capacity . Index 0 is for storage 1, 1 for 4, and etc.

    model.addConstrs(quicksum(inventoryLevel[i,s,b,t] / float(parameterDictionary['TAverageDensity'][0][s][b]) for s in range(moistureIndex) for b in range(biomassIndex)) <= float(parameterDictionary['SVolumnCapacityInventory'][0]) for i in meteringbinSet for t in range(timeIndex)) #Metering bin volumn capacity
    model.addConstrs(quicksum(inventoryLevelBlend[i,s,t] / float(parameterDictionary['DAverageDensityBlend'][int((i-1)/3)][s]) for s in range(moistureIndex)) <= float(parameterDictionary['SVolumnCapacityInventoryBlend'][int((i-1)/3)]) for i in storeageSet for t in range(timeIndex)) #Storage bins volumn capacity
    '''
    model.addConstrs(quicksum(outflowOfEquipmentBlend[i,s,t] for i in previousEquipmentsReactor for s in range(moistureIndex)) <= upperReactor for t in range(timeIndex))#Reactor rate needs to less than upper bound
    model.addConstrs(quicksum(outflowOfEquipmentBlend[i,s,t] for i in previousEquipmentsReactor for s in range(moistureIndex)) >= 0.9 * upperReactorL[t] for t in range(timeIndex)) #Reactor rate needs to great than lower bound
    model.addConstr(quicksum(outflowOfEquipmentBlend[i,s,t] for i in previousEquipmentsReactor for s in range(moistureIndex) for t in range(timeIndex)) >= 0.95 * quicksum(upperReactorL[t] for t in range(timeIndex)))
    '''
    #Operational
    model.addConstrs(quicksum(outflowOfEquipment[0,s,b,t] for t in range(timeIndex)) - volumnOfBale * float(parameterDictionary['DDensity'][s][b]) * float(parameterDictionary['DNumberOfBales'][s][b]) <= 0 for s in range(moistureIndex) for b in range(biomassIndex)) #Total inflow rate should less or equal to total amount privided for each moisture and biomass types
    '''
    model.addConstrs(quicksum((float(parameterDictionary['SAshContent'][int((i-2)/3)]) - float(parameterDictionary['STargetAshThermalCarbohydrate'][0])) * outflowOfEquipmentBlend[i,s,t] for i in previousEquipmentsReactor for s in range(moistureIndex)) <= 0 for t in range(timeIndex)) #Ash content needs to meet target level. Index 2 is for biomass type 0, 5 for 2, and etc.
    model.addConstrs(quicksum((float(parameterDictionary['SThermalContent'][int((i-2)/3)]) - float(parameterDictionary['STargetAshThermalCarbohydrate'][1])) * outflowOfEquipmentBlend[i,s,t] for i in previousEquipmentsReactor for s in range(moistureIndex)) >= 0 for t in range(timeIndex)) #Themal content needs to meet target level
    model.addConstrs(quicksum((float(parameterDictionary['SCarbohydrateContent'][int((i-2)/3)]) - float(parameterDictionary['STargetAshThermalCarbohydrate'][2])) * outflowOfEquipmentBlend[i,s,t] for i in previousEquipmentsReactor for s in range(moistureIndex)) >= 0 for t in range(timeIndex)) #Carbohydrate content needs to meet target level
    '''

    model.addConstr(quicksum(inventoryLevel[i,s,b,timeIndex-1] for i in meteringbinSet for s in range(moistureIndex) for b in range(biomassIndex)) + quicksum(inventoryLevelBlend[i,s,timeIndex-1] for i in storeageSet for s in range(moistureIndex)) <= 0)
    model.addConstrs(quicksum(velocityInFlow[0,s,b,tp] for tp in range(t,t+int(parameterDictionary['DProcessingTimeOfBale'][s][b]))) <=  float(parameterDictionary['SLengthWidthHeight'][0]) + (1-startAtTime[s][b][t]) * Msb for s in range(moistureIndex) for b in range(biomassIndex) for t in range(timeIndex-int(parameterDictionary['DProcessingTimeOfBale'][s][b])+1))
    model.addConstrs(quicksum(velocityInFlow[0,s,b,tp] for tp in range(t,t+int(parameterDictionary['DProcessingTimeOfBale'][s][b]))) >=  float(parameterDictionary['SLengthWidthHeight'][0]) - (1-startAtTime[s][b][t]) * Msb for s in range(moistureIndex) for b in range(biomassIndex) for t in range(timeIndex-int(parameterDictionary['DProcessingTimeOfBale'][s][b])+1))
    model.addConstrs(outflowOfEquipment[0,s,b,t] == (float(parameterDictionary['SLengthWidthHeight'][1]) * float(parameterDictionary['SLengthWidthHeight'][2]) * float(parameterDictionary['DDensity'][s][b])) * velocityInFlow[0,s,b,t] for s in range(moistureIndex) for b in range(biomassIndex) for t in range(timeIndex))

    model.addConstrs(quicksum(processAtTime[sp,bp,tp] for sp in range(moistureIndex) for bp in range(biomassIndex) for tp in range(t,t+int(parameterDictionary['DProcessingTimeOfBale'][s][b]))) - quicksum(processAtTime[s,b,tp] for tp in range(t,t+int(parameterDictionary['DProcessingTimeOfBale'][s][b]))) <= (1 - startAtTime[s][b][t]) * Msb for s in range(moistureIndex) for b in range(biomassIndex) for t in range(timeIndex-int(parameterDictionary['DProcessingTimeOfBale'][s][b])+1))
    model.addConstrs(processAtTimeO[t] <= processAtTimeO[t-1] for t in range(1,timeIndex))

    #Flow Balance
    model.addConstrs(outflowOfEquipment[i,s,b,t] == quicksum(outflowOfEquipment[p,s,b,t] for p in previousEquipments[i]) for i in equipmentRemain for s in range(moistureIndex) for b in range(biomassIndex) for t in range(timeIndex)) #Flow balance for "not specific" equipments
    model.addConstrs(outflowOfEquipment[i,s,b,t] == quicksum((1-float(parameterDictionary['SDryMatterLoss'][int((i-2)/4)])) * outflowOfEquipment[p,s,b,t] for p in previousEquipments[i]) for i in grinderSet for s in range(moistureIndex) for b in range(biomassIndex) for t in range(timeIndex)) #Flow balance for grinders. Index 2,6 is for information index 0,1

    model.addConstrs(outflowOfEquipmentBlend[i,s,t] == quicksum(outflowOfEquipmentBlend[p,s,t] for p in previousEquipmentsBlend[i]) for i in previousEquipmentsReactor for s in range(moistureIndex) for t in range(timeIndex)) #Flow balance for other equipments in blending process
    model.addConstrs(outflowOfEquipmentBlend[i,s,t] == quicksum(outflowOfEquipment[p,s,i/3,t] for p in previousEquipmentsBlend[i]) for i in range(0,equipmentIndexBlend,3) for s in range(moistureIndex) for t in range(timeIndex)) #Flow balance for the first b parallel equipments in blending process

    #Inventory Balance
    model.addConstrs(inventoryLevel[i,s,b,t] == inventoryLevel[i,s,b,t-1] + quicksum(outflowOfEquipment[p,s,b,t] for p in previousEquipments[i]) - outflowOfEquipment[i,s,b,t] for i in meteringbinSet for s in range(moistureIndex) for b in range(biomassIndex) for t in range(1,timeIndex)) #Inventory balance for time index 1 to the largest index for metering bin
    model.addConstrs(inventoryLevelBlend[i,s,t] == inventoryLevelBlend[i,s,t-1] + quicksum(outflowOfEquipmentBlend[p,s,t] for p in previousEquipmentsBlend[i]) - outflowOfEquipmentBlend[i,s,t] for i in storeageSet for s in range(moistureIndex) for t in range(1,timeIndex)) #Inventory balance for time index 1 to the largest index for storage bins
    model.addConstrs(inventoryLevel[i,s,b,0] == 0 + quicksum(outflowOfEquipment[p,s,b,0] for p in previousEquipments[i]) - outflowOfEquipment[i,s,b,0] for i in meteringbinSet for s in range(moistureIndex) for b in range(biomassIndex)) #Inventory balance for time index 0, where time index -1 should be the intial inventory (0) for metering bin
    model.addConstrs(inventoryLevelBlend[i,s,0] == 0 + quicksum(outflowOfEquipmentBlend[p,s,0] for p in previousEquipmentsBlend[i]) - outflowOfEquipmentBlend[i,s,0] for i in storeageSet for s in range(moistureIndex)) #Inventory balance for time index 0, where time index -1 should be the intial inventory (0) for storage bins

    #Linear
    model.addConstrs(upperReactorL[t] >= 0 for t in range(timeIndex))
    model.addConstrs(upperReactorL[t] <= 0.1 * processAtTimeO[t] for t in range(timeIndex))
    model.addConstrs(upperReactorL[t] <= upperReactor for t in range(timeIndex))
    model.addConstrs(upperReactorL[t] >= upperReactor + 0.1 * (processAtTimeO[t] - 1) for t in range(timeIndex))

    #Solve Optimization
    model.optimize()
    
    outflowOfEquipmentList, outflowOfEquipmentBlendList, inventoryLevelList, inventoryLevelBlendList, processAtTimeList, outflowOfEquipmentListP = [], [], [], [], [], []
    for s in range(moistureIndex):
        outflowOfEquipmentListM = []
        inventoryLevelListM = []
        processAtTimeListM = []
        outflowOfEquipmentListPM = []
        for b in range(biomassIndex):
            outflowOfEquipmentListMB = []
            inventoryLevelListMB = []
            processAtTimeListMB = []
            outflowOfEquipmentListPMB = []
            for t in range(timeIndex):
                outflowOfEquipmentListMB.append(outflowOfEquipment[0,s,b,t].x)
                inventoryLevelListMB.append(inventoryLevel[10,s,b,t].x)
                processAtTimeListMB.append(processAtTime[s,b,t].x)
                outflowOfEquipmentListPMB.append(outflowOfEquipment[12,s,b,t].x)
            outflowOfEquipmentListM.append(outflowOfEquipmentListMB)
            inventoryLevelListM.append(inventoryLevelListMB)
            processAtTimeListM.append(processAtTimeListMB)
            outflowOfEquipmentListPM.append(outflowOfEquipmentListPMB)
        outflowOfEquipmentList.append(outflowOfEquipmentListM)
        inventoryLevelList.append(inventoryLevelListM)
        processAtTimeList.append(processAtTimeListM)
        outflowOfEquipmentListP.append(outflowOfEquipmentListPM)
    for i in previousEquipmentsReactor:
        outflowOfEquipmentBlendListB = []
        for s in range(moistureIndex):
            outflowOfEquipmentBlendListBM = []
            for t in range(timeIndex):
                outflowOfEquipmentBlendListBM.append(outflowOfEquipmentBlend[i,s,t].x)
            outflowOfEquipmentBlendListB.append(outflowOfEquipmentBlendListBM)
        outflowOfEquipmentBlendList.append(outflowOfEquipmentBlendListB)
    for i in storeageSet:
        inventoryLevelBlendListB = []
        for s in range(moistureIndex):
            inventoryLevelBlendListBM = []
            for t in range(timeIndex):
                inventoryLevelBlendListBM.append(inventoryLevelBlend[i,s,t].x)
            inventoryLevelBlendListB.append(inventoryLevelBlendListBM)
        inventoryLevelBlendList.append(inventoryLevelBlendListB)
    timeNew = int(round(sum(processAtTimeO[t].x for t in range(timeIndex))))
    return outflowOfEquipmentList, outflowOfEquipmentBlendList, inventoryLevelList, inventoryLevelBlendList, processAtTimeList, timeNew, outflowOfEquipmentListP

def DataOutput(outflowOfEquipmentList, outflowOfEquipmentBlendList, inventoryLevelList, inventoryLevelBlendList, moistureIndex, biomassIndex, timeIndex, outflowOfEquipmentListP):
    fileOut = open(r'/Users/dahuiliu/Documents/BBP reactor feed rate/Output.txt','w')
    fileOut.write('RInputFlow')
    fileOut.write('\n')
    for s in range(moistureIndex):
        for b in range(biomassIndex):
            for t in range(timeIndex):
                fileOut.write(str(outflowOfEquipmentList[s][b][t]))
                fileOut.write(' ')
            fileOut.write('\n')
        fileOut.write('N')
        fileOut.write('\n')
    fileOut.write('ROutFlow')
    fileOut.write('\n')
    for i in range(biomassIndex):
        for s in range(moistureIndex):
            for t in range(timeIndex):
                fileOut.write(str(outflowOfEquipmentBlendList[i][s][t]))
                fileOut.write(' ')
            fileOut.write('\n')
        fileOut.write('N')
        fileOut.write('\n')
    fileOut.write('RMetering')
    fileOut.write('\n')
    for s in range(moistureIndex):
        for b in range(biomassIndex):
            for t in range(timeIndex):
                fileOut.write(str(inventoryLevelList[s][b][t]))
                fileOut.write(' ')
            fileOut.write('\n')
        fileOut.write('N')
        fileOut.write('\n')
    fileOut.write('RStorages')
    fileOut.write('\n')
    for i in range(biomassIndex):
        for s in range(moistureIndex):
            for t in range(timeIndex):
                fileOut.write(str(inventoryLevelBlendList[i][s][t]))
                fileOut.write(' ')
            fileOut.write('\n')
        fileOut.write('N')
        fileOut.write('\n')
    fileOut.write('RPelletFlow')
    fileOut.write('\n')
    for s in range(moistureIndex):
        for b in range(biomassIndex):
            for t in range(timeIndex):
                fileOut.write(str(outflowOfEquipmentListP[s][b][t]))
                fileOut.write(' ')
            fileOut.write('\n')
        fileOut.write('N')
        fileOut.write('\n')
    fileOut.close()
    
#Set and Index Setup
equipmentIndex, moistureIndex, biomassIndex = 12, 3, 4
equipmentIndexBlend = 3 * biomassIndex
timeIndex = 2400 #Max index for different sets
grinderSet = [2,6] #Set of grinders indices
meteringbinSet = [10] #Set of metering bin indices
storeageSet = [] #Set of storage bin (for different biomass) indices
for i in range(biomassIndex):
    storeageSet.append(1+i*3) #There are three equipments in each branch (equipments for one biomass in blending process) and storage bins are the second equipments.
equipmentRemain = [] #Set of equipments except some specific ones (grinders, metering bin, bypass equipments)
equipmentBeforeMetering = []
#equipmentAfterMetering = []
for i in range(1, equipmentIndex+1):
    if i not in grinderSet and i not in meteringbinSet:
        equipmentRemain.append(i)
for i in range(2, meteringbinSet[0]):
    if i not in grinderSet:
        equipmentBeforeMetering.append(i)
Msb = 1000

#Previous Equipments I_i
previousEquipments = [] #Set of predecessors of equipment
for i in range(equipmentIndex+1):
    previousEquipments.append([i-1])
previousEquipmentsBlend = [] #Set of predecessors of equipment in blending process
for i in range(equipmentIndexBlend):
    if i % 3 == 0: #First b parallel equipments in blending process and their predecessors
        previousEquipmentsBlend.append([equipmentIndex])
    else: #Other equipments in blending process and their predecessors
        previousEquipmentsBlend.append([i-1])
previousEquipmentsReactor = [] #Set of predecessors of reactor
for i in range(biomassIndex):
    previousEquipmentsReactor.append(2+i*3)

parameterDictionary = DataInput()
volumnOfBale = float(parameterDictionary['SLengthWidthHeight'][0]) * float(parameterDictionary['SLengthWidthHeight'][1]) * float(parameterDictionary['SLengthWidthHeight'][2]) #Bale volumn
countSum = 2

#1L-3M-1H 10C3-10C2-10S-10M
'''
sequenceList = [[0,2],[1,2],[1,1],[1,1],[2,2],[0,1],[1,2],[1,2],[1,1],[2,1],[0,2],[1,2],[1,1],[1,1],[2,2],[0,1],[1,2],[1,2],[1,1],[2,1],\
                [0,3],[1,3],[1,0],[1,0],[2,3],[0,0],[1,3],[1,3],[1,0],[2,0],[0,3],[1,3],[1,0],[1,0],[2,3],[0,0],[1,3],[1,3],[1,0],[2,0]]
'''
'''
sequenceList = [[0,2],[1,2],[1,2],[1,2],[2,2],[0,1],[1,1],[1,1],[1,1],[2,1],[0,3],[1,3],[1,3],[1,3],[2,3],[0,0],[1,0],[1,0],[1,0],[2,0],\
                [0,2],[1,2],[1,2],[1,2],[2,2],[0,1],[1,1],[1,1],[1,1],[2,1],[0,3],[1,3],[1,3],[1,3],[2,3],[0,0],[1,0],[1,0],[1,0],[2,0]]
'''
'''
sequenceList = [[0,2],[0,2],[1,2],[1,2],[1,2],[1,2],[1,2],[1,2],[2,2],[2,2],[0,1],[0,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[2,1],[2,1],\
                [0,3],[0,3],[1,3],[1,3],[1,3],[1,3],[1,3],[1,3],[2,3],[2,3],[0,0],[0,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[2,0],[2,0]]
'''
'''
sequenceList = [[2, 0], [2, 3], [1, 2], [1, 1], [0, 3], [2, 1], [0, 1], [1, 3], [0, 0], [1, 1], [0, 3], [1, 2], [1, 1], [1, 0], [1, 1], [2, 2], [1, 0], [1, 2], [1, 3], [1, 2],\
                [1, 3], [0, 0], [2, 2], [1, 3], [2, 0], [1, 2], [1, 2], [2, 1], [1, 0], [2, 3], [0, 2], [1, 0], [1, 1], [1, 1], [1, 0], [1, 3], [1, 0], [0, 1], [0, 2], [1, 3]]
'''
'''
sequenceList = [[1, 2], [1, 0], [2, 1], [2, 3], [0, 3], [1, 1], [1, 2], [1, 2], [1, 3], [1, 1], [1, 1], [2, 2], [1, 0], [0, 0], [1, 0], [2, 0], [2, 1], [1, 1], [1, 0], [1, 3],\
                [1, 2], [1, 0], [0, 2], [0, 0], [0, 2], [1, 3], [0, 3], [2, 2], [1, 2], [1, 1], [1, 3], [1, 2], [1, 3], [1, 1], [0, 1], [2, 3], [2, 0], [0, 1], [1, 3], [1, 0]]
'''

#2L-2M-1H 10C3-10C2-10S-10M
'''
sequenceList = [[0,3],[0,0],[1,3],[1,0],[2,3],[0,0],[0,3],[1,0],[1,3],[2,0],[0,3],[0,0],[1,3],[1,0],[2,3],[0,0],[0,3],[1,0],[1,3],[2,0],\
                [0,2],[0,1],[1,2],[1,1],[2,2],[0,1],[0,2],[1,1],[1,2],[2,1],[0,2],[0,1],[1,2],[1,1],[2,2],[0,1],[0,2],[1,1],[1,2],[2,1]]
'''
'''
sequenceList = [[0,3],[0,3],[1,3],[1,3],[2,3],[0,0],[0,0],[1,0],[1,0],[2,0],[0,2],[0,2],[1,2],[1,2],[2,2],[0,1],[0,1],[1,1],[1,1],[2,1],\
                [0,3],[0,3],[1,3],[1,3],[2,3],[0,0],[0,0],[1,0],[1,0],[2,0],[0,2],[0,2],[1,2],[1,2],[2,2],[0,1],[0,1],[1,1],[1,1],[2,1]]
'''
'''
sequenceList = [[0,3],[0,3],[0,3],[0,3],[1,3],[1,3],[1,3],[1,3],[2,3],[2,3],[0,0],[0,0],[0,0],[0,0],[1,0],[1,0],[1,0],[1,0],[2,0],[2,0],\
                [0,2],[0,2],[0,2],[0,2],[1,2],[1,2],[1,2],[1,2],[2,2],[2,2],[0,1],[0,1],[0,1],[0,1],[1,1],[1,1],[1,1],[1,1],[2,1],[2,1]]
'''
'''
sequenceList = [[2, 1], [2, 2], [1, 0], [1, 3], [0, 2], [2, 0], [0, 0], [1, 1], [0, 1], [1, 0], [0, 2], [0, 0], [1, 0], [1, 2], [1, 3], [2, 3], [1, 1], [0, 0], [0, 2], [0, 3],\
                [0, 1], [0, 1], [2, 3], [1, 1], [2, 1], [1, 0], [0, 3], [2, 0], [1, 2], [2, 2], [0, 3], [1, 1], [1, 3], [1, 3], [1, 2], [0, 2], [1, 2], [0, 0], [0, 3], [0, 1]]
'''
'''
sequenceList = [[1, 0], [1, 1], [2, 0], [2, 2], [0, 2], [1, 0], [0, 0], [0, 0], [0, 1], [1, 0], [1, 3], [2, 3], [1, 2], [0, 1], [1, 2], [2, 1], [2, 0], [1, 3], [1, 1], [0, 2],\
                [1, 0], [1, 2], [0, 3], [0, 1], [0, 3], [1, 1], [0, 2], [2, 3], [0, 3], [1, 3], [1, 1], [0, 3], [0, 2], [1, 3], [0, 0], [2, 2], [2, 1], [0, 0], [0, 1], [1, 2]]
'''

#1L-2M-2H 10C3-10C2-10S-10M
'''
sequenceList = [[0,3],[1,1],[1,3],[2,1],[2,3],[0,1],[1,3],[1,1],[2,3],[2,1],[0,3],[1,1],[1,3],[2,1],[2,3],[0,1],[1,3],[1,1],[2,3],[2,1],\
                [0,2],[1,0],[1,2],[2,0],[2,2],[0,0],[1,2],[1,0],[2,2],[2,0],[0,2],[1,0],[1,2],[2,0],[2,2],[0,0],[1,2],[1,0],[2,2],[2,0]]
'''
'''
sequenceList = [[0,3],[1,3],[1,3],[2,3],[2,3],[0,1],[1,1],[1,1],[2,1],[2,1],[0,2],[1,2],[1,2],[2,2],[2,2],[0,0],[1,0],[1,0],[2,0],[2,0],\
                [0,3],[1,3],[1,3],[2,3],[2,3],[0,1],[1,1],[1,1],[2,1],[2,1],[0,2],[1,2],[1,2],[2,2],[2,2],[0,0],[1,0],[1,0],[2,0],[2,0]]
'''
'''
sequenceList = [[0,3],[0,3],[1,3],[1,3],[1,3],[1,3],[2,3],[2,3],[2,3],[2,3],[0,0],[0,0],[1,0],[1,0],[1,0],[1,0],[2,0],[2,0],[2,0],[2,0],\
                [0,2],[0,2],[1,2],[1,2],[1,2],[1,2],[2,2],[2,2],[2,2],[2,2],[0,1],[0,1],[1,1],[1,1],[1,1],[1,1],[2,1],[2,1],[2,1],[2,1]]
'''
'''
sequenceList = [[2, 0], [2, 2], [1, 1], [1, 3], [0, 2], [2, 1], [0, 1], [1, 0], [0, 0], [2, 1], [0, 2], [1, 1], [2, 1], [1, 2], [2, 3], [2, 3], [2, 0], [1, 1], [1, 2], [1, 3],\
                [1, 0], [0, 0], [2, 3], [1, 0], [2, 0], [1, 1], [1, 3], [2, 1], [2, 2], [2, 2], [0, 3], [2, 0], [1, 3], [2, 3], [2, 2], [1, 2], [1, 2], [0, 1], [0, 3], [1, 0]]
'''

sequenceList = [[1, 1], [2, 0], [2, 1], [2, 2], [0, 2], [2, 1], [1, 1], [1, 1], [1, 0], [2, 1], [2, 3], [2, 3], [1, 2], [0, 0], [2, 2], [2, 0], [2, 1], [1, 3], [2, 0], [1, 2],\
                [1, 1], [2, 2], [0, 3], [0, 0], [0, 3], [1, 0], [0, 2], [2, 3], [1, 3], [2, 3], [1, 0], [1, 3], [1, 2], [1, 3], [0, 1], [2, 2], [2, 0], [0, 1], [1, 0], [1, 2]]

startAtTime = SequenceInitial(sequenceList, moistureIndex, biomassIndex, timeIndex, parameterDictionary)

while countSum !=0:
    outflowOfEquipmentList, outflowOfEquipmentBlendList, inventoryLevelList, inventoryLevelBlendList, processAtTimeList, timeIndex, outflowOfEquipmentListP = MILPModel(equipmentIndex, moistureIndex, biomassIndex, equipmentIndexBlend, timeIndex, grinderSet, meteringbinSet, storeageSet, equipmentRemain, previousEquipments, previousEquipmentsBlend, previousEquipmentsReactor, parameterDictionary, volumnOfBale, startAtTime)
    
    #Narrow Processing Time
    yList = []
    for s in range(moistureIndex):
        for b in range(biomassIndex):
            for t in range(timeIndex):
                if startAtTime[s][b][t] >= 1 - pow(10,-6):
                    yList.append([s,b,t])
    #print(yList)
    yListSort = sorted(yList, key=lambda x: x[2])
    #print(yListSort)
    count = []
    for i in range(len(yList)):
        count.append(0)
    for i in range(len(yListSort)):
        for t in range(yListSort[i][2],yListSort[i][2]+int(parameterDictionary['DProcessingTimeOfBale'][yListSort[i][0]][yListSort[i][1]])):
            #if processAtTime[yListSort[i][0],yListSort[i][1],t].x <= pow(10,-6):
            if processAtTimeList[yListSort[i][0]][yListSort[i][1]][t] <= pow(10,-6):
                count[i] += 1
    print(count)
    minReduce = []
    for s in range(moistureIndex):
        minReduceRow = []
        for b in range(biomassIndex):
            minReduceRow.append(1000)
        minReduce.append(minReduceRow)
    countSum = 0
    countTotal,countN,countLoc0,countLoc1 = 0,0,0,0
    for i in range(len(yListSort)):
        if i == 0:
            countN += 1
            countTotal += count[i]
            countLoc0 = yListSort[i][0]
            countLoc1 = yListSort[i][1]
        else:
            if countLoc0 == yListSort[i][0] and countLoc1 == yListSort[i][1]:
                countN += 1
                countTotal += count[i]
            else:
                averageR = int(countTotal/countN)
                if averageR <= minReduce[yListSort[i-1][0]][yListSort[i-1][1]]:
                    minReduce[yListSort[i-1][0]][yListSort[i-1][1]] = averageR
                countN = 1
                countTotal = count[i]
                countLoc0 = yListSort[i][0]
                countLoc1 = yListSort[i][1]
    averageR = int(countTotal/countN)
    if averageR <= minReduce[yListSort[len(yListSort)-1][0]][yListSort[len(yListSort)-1][1]]:
        minReduce[yListSort[len(yListSort)-1][0]][yListSort[len(yListSort)-1][1]] = averageR
    for s in range(moistureIndex):
        for b in range(biomassIndex):
            if int(parameterDictionary['DProcessingTimeOfBale'][s][b])-minReduce[s][b] >= 0:
                parameterDictionary['DProcessingTimeOfBale'][s][b] = str(int(parameterDictionary['DProcessingTimeOfBale'][s][b])-minReduce[s][b])
                countSum += minReduce[s][b]
    #startAtTime[:] = []
    startAtTime = SequenceInitial(sequenceList, moistureIndex, biomassIndex, timeIndex, parameterDictionary)

tempSum = 0
for s in range(moistureIndex):
    for b in range(biomassIndex):
        tempSum += int(parameterDictionary['DNumberOfBales'][s][b]) * int(parameterDictionary['DProcessingTimeOfBale'][s][b])
totalOutToReactor = 0
for i in range(biomassIndex):
    for s in range(moistureIndex):
        for t in range(timeIndex):
            totalOutToReactor += outflowOfEquipmentBlendList[i][s][t]
DataOutput(outflowOfEquipmentList, outflowOfEquipmentBlendList, inventoryLevelList, inventoryLevelBlendList, moistureIndex, biomassIndex, timeIndex, outflowOfEquipmentListP)
print(totalOutToReactor)
print(tempSum)
print(timeIndex)

