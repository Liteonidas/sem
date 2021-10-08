from __future__ import division, print_function
import os
import re
import random
import time
import binascii
from bisect import bisect_right
from heapq import heappop, heappush


numBands = 20
numRows = 10
numHashes = numRows * numBands

numDocs = 1000
# numDocs = 2500
# numDocs = 9950
File = "./data/artykul_" + str(numDocs) + ".txt"


#    Convert Documents To Sets of Shingles


print("Shingling ...")


curShingleID = 0
docsAsShingleSets = {}


f = open(File, "r")

docNames = []

t0 = time.time()

totalShingles = 0

for i in range(0, numDocs):
    words = f.readline().split(" ")
    docID = words[0]
    # print(docID)

    docNames.append(docID)
    del words[0]
    shinglesInDoc = set()
    for index in range(0, len(words) - 2):
        
        shingle = words[index] + " " + \
            words[index + 1] + " " + words[index + 2]

    # for index in range(0, len(words) - 5):
    #     shingle = words[index] + " " + \
    #         words[index + 1] + " " + words[index + 2] + " " + \
    #         words[index + 3] + " " + words[index + 4] + " " + \
    #         words[index + 5]

    # for index in range(0, len(words) - 8):
    #     shingle = words[index] + " " + \
    #         words[index + 1] + " " + words[index + 2] + " " + \
    #         words[index + 3] + " " + words[index + 4] + " " + \
    #         words[index + 5] + " " + words[index + 6] + " " + \
    #         words[index + 7] + " " + words[index + 8]

        crc = binascii.crc32(str.encode(shingle)) & 0xffffffff
        shinglesInDoc.add(crc)
    docsAsShingleSets[docID] = shinglesInDoc
    totalShingles = totalShingles + (len(words) - 2)
f.close()

# with open('tablica1.txt', 'w') as file:
#     file.write(str(docsAsShingleSets['t120']))
t1 = time.time() - t0
print('\nShingling ' + str(numDocs) +
      ' artyk trwal %.2f sek.' % t1)

print('\nSrednia k-shingle: %.2f' % (totalShingles / numDocs))


#----------------Define Triangle Matrices-----------



numElems = int(numDocs * (numDocs - 1) / 2)


JSim = [0 for x in range(numElems)]
estJSim = [0 for x in range(numElems)]


# Define a function to map a 2D matrix coordinate into a 1D index.
def getTriangleIndex(i, j):
    # If i == j that's an error.
    if i == j:
        sys.stderr.write("Can't access triangle matrix with i == j")
        sys.exit(1)
    # If j < i just swap the values.
    if j < i:
        temp = i
        i = j
        j = temp

    k = int(i * (numDocs - (i + 1) / 2.0) + j - i) - 1

    return k


#------------------Calculate Jaccard Similarities--------------------


if numDocs < 2500:
    print("\nJaccard ")
    t0 = time.time()
    for i in range(0, numDocs):
      s1 = docsAsShingleSets[docNames[i]]
      for j in range(i + 1, numDocs):
        s2 = docsAsShingleSets[docNames[j]]
        JSim[getTriangleIndex(i, j)] = (len(s1.intersection(s2)) / len(s1.union(s2)))
    elapsed = (time.time() - t0)

    print("\nJaccard trwal %.2fsec" % elapsed)

del JSim


#--------------Generate MinHash Signatures---------------


t0 = time.time()

maxShingleID = 2 ** 32 - 1
nextPrime = 4294967311

def genHashFunc(k):
    randList = []
    while k > 0:
        randIndex = random.randint(0, maxShingleID)
        while randIndex in randList:
            randIndex = random.randint(0, maxShingleID)
        randList.append(randIndex)
        k = k - 1

    return randList


coeffA = genHashFunc(numHashes)
coeffB = genHashFunc(numHashes)

print('\nGenerowanie MinHash')

signatures = []

for docID in docNames:

    shingleIDSet = docsAsShingleSets[docID]
    signature = []

    for i in range(0, numHashes):

        minHashCode = nextPrime + 1
        for shingleID in shingleIDSet:
            hashCode = (coeffA[i] * shingleID + coeffB[i]) % nextPrime
            if hashCode < minHashCode:
                minHashCode = hashCode
        signature.append(minHashCode)
    signatures.append(signature)
elapsed = (time.time() - t0)

print("\nGenerowanie MinHash trwalo %.2fsek" % elapsed)


#-----------LSH Buckets------------------------------


print('\nGenerowanie bucketow')

t0 = time.time()
buckets = []

for i in range(0, numBands):
    bucket = {}
    for signature in signatures:
        hash = 0
        index = signatures.index(signature)
        for j in range(0, numRows):
            hash = (hash * 2654435761 + signature[numRows * i + j])
            hash = hash % (1<<32)
        if hash not in bucket.keys():
            bucket[hash] = set()
        bucket[hash].add(docNames[index])
    buckets.append(bucket)
    del bucket


elapsed = (time.time() - t0)
            
print("\nGenerowanie bucketow zajelo %.2fsec" % elapsed)


#------------LSH select bucket-------------------


candidateDocs = {}
dicNames = {}


for bucket in buckets:
    for hash in bucket:
        if (len(bucket[hash]) > 1):
            
            for docID in bucket[hash]:
                if docID not in candidateDocs.keys():
                    candidateDocs[docID] = []
                candidateDocs[docID] = list(set(bucket[hash]) | set(candidateDocs[docID]))
                candidateDocs[docID].remove(docID)
# Remove duplicated document
for docID in list(candidateDocs):
    for id in list(candidateDocs[docID]):
        if id in candidateDocs:
            if docID in candidateDocs[id]:
                candidateDocs[id].remove(docID)
                # candidateDocs.pop(docID)
                # del candidateDocs[docID]
    
# Delete empty list
for docID in list(candidateDocs):
    if (len(candidateDocs.get(docID)) < 1):
        # candidateDocs.pop(docID)
        del candidateDocs[docID]

c = 1
for bucket in buckets:
    print("Bucket" + str(c))
    c += 1
    for hash in bucket:
        if (len(bucket[hash]) > 1):
            print("hash: " + str(hash) + str(bucket[hash]))

# for bucket in buckets:
#     for hashSet in bucket:
#         nBuck = len(bucket[hashSet]
#         if (nBuck > 1):
#             for doc in bucket(hashSet):
#                 # bucketHash jest slownikiem zbiorow
#                 index = bucket(hashSet).index(doc)
#                 if (index + 1) < nBuck
#                 ind2 = index + 1
#                 for doc + 1 in bucket(hashSet):

#             print("hash: " + str(hash) + str(bucket[hash]))
# docsAsShingleSets[docID]
# JacSimCandidateCheck
count = 0
for docID in candidateDocs:
    count += 1
    # similarDocs.add(candidateDocs[docID].values())
    print(docID + "=> " + str(candidateDocs[docID]))
    # print(candidateDocs[docID])
print(count)
# for i in similarDocs:
#     print(i)


