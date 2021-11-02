from __future__ import division, print_function
import random
import time
import binascii
import sys
import argparse

#    Convert Documents To Sets of Shingles


def shingling(file, numDocs):
    print("\nShingling ...")

    # curShingleID = 0
    docsAsShingleSets = {}
    global docNames

    f = open(file, "r")

    t0 = time.time()

    totalShingles = 0

    for i in range(0, numDocs):
        words = f.readline().split(" ")
        docID = words[0]
        # print(docID)

        docNames.append(docID)
        del words[0]
        shinglesInDoc = set()
        if numK == 3:
            for index in range(0, len(words) - 2):
                shingle = words[index] + " " + \
                    words[index + 1] + " " + words[index + 2]
                crc = binascii.crc32(str.encode(shingle)) & 0xffffffff
                shinglesInDoc.add(crc)
        elif numK == 6:
            for index in range(0, len(words) - 5):
                shingle = words[index] + " " + \
                    words[index + 1] + " " + words[index + 2] + " " + \
                    words[index + 3] + " " + words[index + 4] + " " + \
                    words[index + 5]
                crc = binascii.crc32(str.encode(shingle)) & 0xffffffff
                shinglesInDoc.add(crc)
        elif numK == 9:
            for index in range(0, len(words) - 8):
                shingle = words[index] + " " + \
                    words[index + 1] + " " + words[index + 2] + " " + \
                    words[index + 3] + " " + words[index + 4] + " " + \
                    words[index + 5] + " " + words[index + 6] + " " + \
                    words[index + 7] + " " + words[index + 8]
                crc = binascii.crc32(str.encode(shingle)) & 0xffffffff
                shinglesInDoc.add(crc)
        else:
            sys.stderr.write("Wrong variable k!")
            sys.exit(1)

        docsAsShingleSets[docID] = shinglesInDoc
        if numK == 3:
            totalShingles = totalShingles + (len(words) - 2)
        elif numK == 6:
            totalShingles = totalShingles + (len(words) - 5)
        elif numK == 9:
            totalShingles = totalShingles + (len(words) - 8)
    f.close()
    t1 = time.time() - t0
    print('\nShingling ' + str(numDocs) +
          ' article took %.2f sec.' % t1)

    print('\nAverage k-shingle: %.2f' % (totalShingles / numDocs))
    return docsAsShingleSets


# ----------------Define Triangle Matrices-----------

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


# ------------------Calculate Jaccard Similarities--------------------

# This function check the time taken by created Jaccard similarity
# Matrix is deleted after comparision and not used

def jacSimilarities(shingleSets):
    JSim = [0 for x in range(numElems)]
    print("\nBuilding index Jaccard ")
    t0 = time.time()
    for i in range(0, numDocs):
        s1 = shingleSets[docNames[i]]
        for j in range(i + 1, numDocs):
            s2 = shingleSets[docNames[j]]
            JSim[getTriangleIndex(i, j)] = (
                len(s1.intersection(s2)) / len(s1.union(s2)))
    elapsed = (time.time() - t0)
    print("\nCompare Jaccarda took %.2fsec" % elapsed)
    # return JSim
    # del JSim
    # return JSim[getTriangleIndex(i, j)]


# --------------Generate MinHash Signatures---------------


def genMinHashSig(shingleSets, numHashes):
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

    print('\nGenerating MinHash')

    signatures = []

    for docID in docNames:

        shingleIDSet = shingleSets[docID]
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

    print("\nGenarating MinHash took %.2fsec" % elapsed)
    return signatures


# -----------Compare MinHash signatures------------------------------


def compareMinHash(signatures, numHashes):
    estJSim = [0 for x in range(numElems)]
    t0 = time.time()
    # For each of the test documents...
    for i in range(0, numDocs):
        signature1 = signatures[i]
        for j in range(i + 1, numDocs):
            signature2 = signatures[j]
            count = 0
            for k in range(0, numHashes):
                count = count + (signature1[k] == signature2[k])
            estJSim[getTriangleIndex(i, j)] = (count / numHashes)
    elapsed = (time.time() - t0)
    print("\nComparing MinHash signatures took %.2fsec" % elapsed)
    return estJSim


def displayComparisionMinHash(estJSim, threshold):
    print("                   Est. J   Act. J")
    for i in range(0, numDocs):
        for j in range(i + 1, numDocs):
            estJ = estJSim[getTriangleIndex(i, j)]
            if estJ > threshold:
                s1 = docsAsShingleSets[docNames[i]]
                s2 = docsAsShingleSets[docNames[j]]
                J = (len(s1.intersection(s2)) / len(s1.union(s2)))
                print("  %5s --> %5s   %.2f     %.2f" %
                      (docNames[i], docNames[j], estJ, J))


# -----------LSH Buckets------------------------------

def genBuckets(signatures, numBands, numRows):
    print('\nGenerating buckets')

    t0 = time.time()
    buckets = []

    for i in range(0, numBands):
        bucket = {}
        for signature in signatures:
            hash = 0
            index = signatures.index(signature)
            for j in range(0, numRows):
                hash = (hash * 2654435761 + signature[numRows * i + j])
                hash = hash % (1 << 32)
            if hash not in bucket.keys():
                bucket[hash] = set()
            bucket[hash].add(docNames[index])
        buckets.append(bucket)
        del bucket

    elapsed = (time.time() - t0)

    print("\nGenerating buckets took %.2fsec" % elapsed)
    return buckets

# ------------LSH select bucket-------------------


def selectCandidatePair(buckets):
    candidateDocs = {}
    for bucket in buckets:
        for hash in bucket:
            if (len(bucket[hash]) > 1):
                for docID in bucket[hash]:
                    if docID not in candidateDocs.keys():
                        candidateDocs[docID] = []
                    candidateDocs[docID] = list(
                        set(bucket[hash]) | set(candidateDocs[docID]))
                    candidateDocs[docID].remove(docID)
    # Remove duplicated document
    for docID in list(candidateDocs):
        for id in list(candidateDocs[docID]):
            if id in candidateDocs:
                if docID in candidateDocs[id]:
                    candidateDocs[id].remove(docID)

    # Delete empty list
    for docID in list(candidateDocs):
        if (len(candidateDocs.get(docID)) < 1):
            del candidateDocs[docID]

    return candidateDocs


def displayBuckets():
    c = 1
    for bucket in buckets:
        print("Bucket" + str(c))
        c += 1
        for hash in bucket:
            if (len(bucket[hash]) > 1):
                print("hash: " + str(hash) + str(bucket[hash]))


def countCandidatePair(candidateDocs):
    count = 0
    for docID in candidateDocs:
        count += 1
        print(docID + "=> " + str(candidateDocs[docID]))
    print(count)


if __name__ == '__main__':
    numBands = 20
    numRows = 10
    numHashes = numRows * numBands
    threshold = 0.7
    numK = 3
    numDocs = 1000
    File = "./data/artykul_" + str(numDocs) + ".txt"
    docNames = []
    # Number of elements triangle matrix
    numElems = int(numDocs * (numDocs - 1) / 2)
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("-jaccard", action="store_true",
                        help="calculate Jaccard similarity")
    parser.add_argument("-minhash", action="store_true",
                        help="calculate MinHash similarity")
    parser.add_argument("-t", type=float, action="store",
                        help="set threshold minhash; use with -minhash(default \
                             is 0.7)")
    group.add_argument("-s", "--small", action="store_true",
                       help="number of documents=1000(default is small=1000)")
    group.add_argument("-m", "--medium", action="store_true",
                       help="number of documents=2500(default is small=1000)")
    group.add_argument("-l", "--large", action="store_true",
                       help="number of documents=9950(default is small=1000)")
    parser.add_argument(
        "-k", type=int, choices=[3, 6, 9], default=3,
        help="pick how long must be k-shingle(default is 3)")
    parser.add_argument("-bands", type=int, action='store',
                        default=20, help="number of bands(default is 20)")
    parser.add_argument("-rows", type=int, action='store', default=10,
                        help="number of rows in bands(default is 10)")
    parser.add_argument("-dspHash", action="store_true",
                        help="display hash in buckets")
    args = parser.parse_args()

    if args.small:
        numDocs = 1000
    elif args.medium:
        numDocs = 2500
    elif args.large:
        numDocs = 9950

    if args.t is not None:
        if (args.t > 0 and args.t < 1):
            if args.minhash:
                print(f"Run with parameter threshold = {args.t}")
            else:
                sys.stderr.write("Option -t run only with parameter -minhash!")
                sys.exit(1)
        elif (args.t <= 0 or args.t >= 1):
            sys.stderr.write(
                "Threshold must be within the range (0, 1).\n\
Use -h or --help for help.")
            sys.exit(1)

    print(f"\nProgram running on {numDocs} docs.")

    if args.bands:
        numBands = args.bands
        print(f"\nSet number of bands = {numBands}")

    if args.rows:
        numRows = args.rows
        print(f"set number of rows = {numRows}")

    if args.k == 3:
        numK = 3
    elif args.k == 6:
        numK = 6
    elif args.k == 9:
        numK = 9

    print(f"\nRun with parameter k = {numK}")

    docsAsShingleSets = shingling(File, numDocs)

    if args.jaccard:
        jacSimilarities(docsAsShingleSets)

    signatures = genMinHashSig(docsAsShingleSets, numHashes)
    if args.minhash:
        estJSim = compareMinHash(signatures, numHashes)
        displayComparisionMinHash(estJSim, threshold)

    buckets = genBuckets(signatures, numBands, numRows)
    candidateDocs = selectCandidatePair(buckets)
    if args.dspHash:
        displayBuckets()

    countCandidatePair(candidateDocs)
