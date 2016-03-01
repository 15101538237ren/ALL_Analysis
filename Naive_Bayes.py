# -*- coding:utf-8 -*-
#  Example of Naive Bayes implemented from Scratch in Python
import csv
import random
import math
import os
class_hash={}
def loadCsv(filename,header=None,col_First=None):
    lines = csv.reader(open(filename, "rb"))
    dataset = list(lines)
    if header:
        del dataset[0]
    class_count=0
    for i in range(len(dataset)):
        for index,val in enumerate(dataset[i]):
            if col_First and index==0:
                continue
            #最后一列的话,将class_name转换为hash_value

            if index==len(dataset[i])-1:
                if val not in class_hash:
                    class_count=class_count+1
                    class_hash[val]=class_count
                val=class_hash[val]
            dataset[i][index]=float(val)
        if col_First:
            del dataset[i][0]
    return dataset
def splitDataset(dataset, splitRatio):
    trainSize = int(len(dataset) * splitRatio)
    trainSet = []
    copy = list(dataset)
    while len(trainSet) < trainSize:
        index = random.randrange(len(copy))
        trainSet.append(copy.pop(index))
    return [trainSet, copy]

def separateByClass(dataset):
    separated = {}
    for i in range(len(dataset)):
        vector = dataset[i]
        if (vector[-1] not in separated):
            separated[vector[-1]] = []
        separated[vector[-1]].append(vector)
    return separated

def mean(numbers):
    return sum(numbers)/float(len(numbers))

def stdev(numbers):
    avg = mean(numbers)
    variance = sum([pow(x-avg,2) for x in numbers])/float(len(numbers)-1)
    return math.sqrt(variance)

def summarize(dataset):
    summaries = [(mean(attribute), stdev(attribute)) for attribute in zip(*dataset)]
    del summaries[-1]
    return summaries

def summarizeByClass(dataset):
    separated = separateByClass(dataset)
    summaries = {}
    for classValue, instances in separated.iteritems():
        summaries[classValue] = summarize(instances)
    return summaries

def calculateProbability(x, mean, stdev):
    if stdev==0.0:
        stdev=math.pow(10,-5)
    exponent = math.exp(-(math.pow(x-mean,2)/(2*math.pow(stdev,2))))
    return (1 / (math.sqrt(2*math.pi) * stdev)) * exponent

def calculateClassProbabilities(summaries, inputVector):
    probabilities = {}
    for classValue, classSummaries in summaries.iteritems():
        probabilities[classValue] = 1
        for i in range(len(classSummaries)):
            mean, stdev = classSummaries[i]
            x = inputVector[i]
            probabilities[classValue] *= calculateProbability(x, mean, stdev)
    return probabilities

def predict(summaries, inputVector):
    probabilities = calculateClassProbabilities(summaries, inputVector)
    bestLabel, bestProb = None, -1
    for classValue, probability in probabilities.iteritems():
        if bestLabel is None or probability > bestProb:
            bestProb = probability
            bestLabel = classValue
    return bestLabel

def getPredictions(summaries, testSet):
    predictions = []
    for i in range(len(testSet)):
        result = predict(summaries, testSet[i])
        predictions.append(result)
    return predictions

def getAccuracy(testSet, predictions):
    correct = 0
    for i in range(len(testSet)):
        if testSet[i][-1] == predictions[i]:
            correct += 1
    return (correct/float(len(testSet))) * 100.0

def main():

    n_fold=10
    splitRatio = 1.0-1.0/n_fold

    input_dir="data"+os.sep+"methylation_data"

    for filename in os.listdir(input_dir):
        if os.path.isfile(os.path.join(input_dir,filename)):
            dataset = loadCsv(os.path.join(input_dir,filename),True,True)
            i=1
            while i<=n_fold:
                trainingSet, testSet = splitDataset(dataset, splitRatio)
                print('Split {0} rows into train={1} and test={2} rows').format(len(dataset), len(trainingSet), len(testSet))
                # prepare model
                summaries = summarizeByClass(trainingSet)
                # test model
                predictions = getPredictions(summaries, testSet)
                accuracy = getAccuracy(testSet, predictions)
                print('{0}:Accuracy for {1}-th run is {2}%').format(filename,i,accuracy)
                i=i+1
main()