'''
@author: Garvin
'''
from numpy import *
import matplotlib.pyplot as plt
import os,csv
class_hash={}
def loadDataSet(fileName, delim=','):
  fr = open(fileName)
  stringArr = [line.strip().split(delim) for line in fr.readlines()]
  datArr = [map(float,line) for line in stringArr]
  return mat(datArr)

def percentage2n(eigVals,percentage):
    sortArray=sort(eigVals)
    sortArray=sortArray[-1::-1]
    arraySum=sum(sortArray)
    tmpSum=0
    num=0
    for i in sortArray:
        tmpSum+=i
        num+=1
        if tmpSum>=arraySum*percentage:
            return num
def pca(dataMat,percentage=0.99):
  meanVals = mean(dataMat, axis=0)
  meanRemoved = dataMat - meanVals #remove mean
  covMat = cov(meanRemoved, rowvar=0)
  eigVals,eigVects = linalg.eig(mat(covMat))
  topNfeat=percentage2n(eigVals,percentage)
  eigValInd = argsort(eigVals)			#sort, sort goes smallest to largest
  eigValInd = eigValInd[-1:-(topNfeat+1):-1]
  redEigVects = eigVects[:,eigValInd]	   #reorganize eig vects largest to smallest
  lowDDataMat = meanRemoved * redEigVects#transform data into new dimensions
  reconMat = (lowDDataMat * redEigVects.T) + meanVals
  return lowDDataMat, reconMat

def plotBestFit(dataSet1,dataSet2):
  dataArr1 = array(dataSet1)
  dataArr2 = array(dataSet2)
  n = shape(dataArr1)[0]
  n1=shape(dataArr2)[0]
  xcord1 = []; ycord1 = []
  xcord2 = []; ycord2 = []
  xcord3=[];ycord3=[]
  j=0
  for i in range(n):

      xcord1.append(dataArr1[i,0]); ycord1.append(dataArr1[i,1])
      xcord2.append(dataArr2[i,0]); ycord2.append(dataArr2[i,1])
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.scatter(xcord1, ycord1, s=30, c='red', marker='s')
  ax.scatter(xcord2, ycord2, s=30, c='green')

  plt.xlabel('X1'); plt.ylabel('X2');
  plt.show()

if __name__=='__main__':
    input_dir="data"+os.sep+"methylation_data"
    filename="methy_gene_list_ALL_and_Normal2.csv"
    #dataset = loadCsv(os.path.join(input_dir,filename),True,True)
    mata=loadDataSet(os.path.join(input_dir,filename))
    lowDDataMat,reconMat= pca(mata, 0.99)
    print "haha"