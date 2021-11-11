from __future__ import division, print_function
import random
import time
import binascii
import sys
import argparse

#    Convert Documents To Sets of Shingles


def shingling(file, numDocs):
    print("\nShingling ...")

    # Create a dictionary of the articles, mapping the article identifier
    # to the list of shingle IDs that appear in the document.
    docsAsShingleSets = {}
    global docNames

    f = open(file, "r")

    t0 = time.time()

    totalShingles = 0

    for i in range(0, numDocs):
        # Read all of the words (they are all on one line)
        # and split them by white space.
        words = f.readline().split(" ")
        # Get article ID.
        docID = words[0]

        # Maintain a list of all document IDs.
        docNames.append(docID)
        del words[0]

        # 'shinglesInDoc' will hold all of the unique shingle IDs
        # present in the current document. If a shingle ID occurs multiple
        # times in the document, it will only appear once in the set.
        shinglesInDoc = set()
        if numK == 3:
            # For each word in the document...
            for index in range(0, len(words) - 2):
                # Construct the shingle text by combining three words together.
                shingle = words[index] + " " + \
                    words[index + 1] + " " + words[index + 2]
                # Hash the shingle to a 32-bit integer.
                crc = binascii.crc32(str.encode(shingle)) & 0xffffffff

                # Add the hash value to the list of shingles
                # for the current document. Note that set objects will only
                # add the value to the set if the set
                # doesn't already contain it.
                shinglesInDoc.add(crc)
        elif numK == 6:
            for index in range(0, len(words) - 5):
                # Construct the shingle text by combining six words together.
                shingle = words[index] + " " + \
                    words[index + 1] + " " + words[index + 2] + " " + \
                    words[index + 3] + " " + words[index + 4] + " " + \
                    words[index + 5]
                crc = binascii.crc32(str.encode(shingle)) & 0xffffffff
                shinglesInDoc.add(crc)
        elif numK == 9:
            for index in range(0, len(words) - 8):
                # Construct the shingle text by combining nine words together.
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

        # Store the completed list of shingles
        # for this document in the dictionary.
        docsAsShingleSets[docID] = shinglesInDoc

        # Count the number of shingles across all documents.
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

    # Calculate the index within the triangular array.
    # This indexing scheme is taken from pg. 223 of:
    # http://infolab.stanford.edu/~ullman/mmds/ch6.pdf
    # But I adapted it for a 0-based index.
    k = int(i * (numDocs - (i + 1) / 2.0) + j - i) - 1

    return k


# ------------------Calculate Jaccard Similarities--------------------

# This function check the time taken by created Jaccard similarity
# Matrix is deleted after comparision and not used

def jacSimilarities(shingleSets):
    # Initialize empty list to store the similarity values.
    # 'JSim' will be for the actual Jaccard Similarity values.
    JSim = [0 for x in range(numElems)]
    print("\nBuilding index Jaccard ... ")
    t0 = time.time()
    # For every document pair...
    for i in range(0, numDocs):
        # Retrieve the set of shingles for document i.
        s1 = shingleSets[docNames[i]]
        for j in range(i + 1, numDocs):
            # Retrieve the set of shingles for document j.
            s2 = shingleSets[docNames[j]]
            # Calculate and store the actual Jaccard similarity.
            JSim[getTriangleIndex(i, j)] = (
                len(s1.intersection(s2)) / len(s1.union(s2)))
    # Calculate the elapsed time (in seconds)
    elapsed = (time.time() - t0)
    print("\nCompare Jaccarda took %.2fsec" % elapsed)
    # Delete the Jaccard Similarities, because it's a big matrix.
    del JSim


# --------------Generate MinHash Signatures---------------


def genMinHashSig(shingleSets, numHashes):
    t0 = time.time()

    # Record the maximum shingle ID that we assigned.
    maxShingleID = 2 ** 32 - 1
    # Largest prime number above 'maxShingleID'.
    nextPrime = 4294967311

    # Our random hash function will take the form of:
    # h(x) = (a*x + b) % c
    # Where 'x' is the input value, 'a' and 'b' are random coefficients,
    # and 'c' is a prime number just greater than maxShingleID.

    # Generate a list of 'k' random coefficients for the random hash functions,
    # while ensuring that the same value does not appear multiple times in the
    # list.
    def genHashFunc(k):
        # Create a list of 'k' random values.
        randList = []
        while k > 0:
            # Get a random shingle ID.
            randIndex = random.randint(0, maxShingleID)
            # Ensure that each random number is unique.
            while randIndex in randList:
                randIndex = random.randint(0, maxShingleID)
            # Add the random number to the list.
            randList.append(randIndex)
            k = k - 1

        return randList
    # For each of the 'numHashes' hash functions,
    # generate a different coefficient 'a' and 'b'.
    coeffA = genHashFunc(numHashes)
    coeffB = genHashFunc(numHashes)

    print('\nGenerating MinHash ...')

    # List of documents represented as signature vectors
    signatures = []

    # Rather than generating a random permutation of all possible shingles,
    # we'll just hash the IDs of the shingles that are
    # actually in the document, then take the lowest resulting
    # hash code value. This corresponds to the index of the first shingle
    # that you would have encountered in the random order.

    # For each document...
    for docID in docNames:

        # Get the shingle set for this document.
        shingleIDSet = shingleSets[docID]
        # The resulting minhash signature for this document.
        signature = []

        # For each of the random hash functions...
        for i in range(0, numHashes):
            # For each of the shingles actually in the document,
            # calculate its hash code using hash function 'i'.

            # Track the lowest hash ID seen. Initialize 'minHashCode'
            # to be greater than the maximum
            # possible value output by the hash.
            minHashCode = nextPrime + 1

            # For each shingle in the document...
            for shingleID in shingleIDSet:
                # Evaluate the hash function.
                hashCode = (coeffA[i] * shingleID + coeffB[i]) % nextPrime
                # Track the lowest hash code seen.
                if hashCode < minHashCode:
                    minHashCode = hashCode
            # Add the smallest hash code value as
            # component number 'i' of the signature.
            signature.append(minHashCode)
        # Store the MinHash signature for this document.
        signatures.append(signature)
    elapsed = (time.time() - t0)

    print("\nGenarating MinHash took %.2fsec" % elapsed)
    return signatures


# -----------Compare MinHash signatures------------------------------


def compareMinHash(signatures, numHashes):
    # Initialize empty list to store the similarity values.
    # 'estJSim' will be for the estimated Jaccard Similarities
    # found by comparing the MinHash signatures.
    estJSim = [0 for x in range(numElems)]
    t0 = time.time()
    # For each of the test documents...
    for i in range(0, numDocs):
        # Get the MinHash signature for document i.
        signature1 = signatures[i]
        # For each of the other test documents...
        for j in range(i + 1, numDocs):
            # Get the MinHash signature for document j.
            signature2 = signatures[j]
            count = 0
            # Count the number of positions in the minhash signature
            # which are equal.
            for k in range(0, numHashes):
                count = count + (signature1[k] == signature2[k])
            estJSim[getTriangleIndex(i, j)] = (count / numHashes)
    elapsed = (time.time() - t0)
    print("\nComparing MinHash signatures took %.2fsec" % elapsed)
    return estJSim


def displayComparisionMinHash(estJSim, threshold):
    print("                   Est. J   Act. J")
    # For each of the document pairs...
    for i in range(0, numDocs):
        for j in range(i + 1, numDocs):
            # Retrieve the estimated similarity value for this pair.
            estJ = estJSim[getTriangleIndex(i, j)]
            # If the similarity is above the threshold...
            if estJ > threshold:
                # Calculate the actual Jaccard similarity for validation.
                s1 = docsAsShingleSets[docNames[i]]
                s2 = docsAsShingleSets[docNames[j]]
                J = (len(s1.intersection(s2)) / len(s1.union(s2)))

                print("  %5s --> %5s   %.2f     %.2f" %
                      (docNames[i], docNames[j], estJ, J))


# -----------LSH Buckets------------------------------

def genBuckets(signatures, numBands, numRows):
    print('\nGenerating buckets ...')

    t0 = time.time()
    # Create a list of buckets.
    buckets = []
    # For each band...
    for i in range(0, numBands):
        # Create bucket which is dictionary.
        bucket = {}
        # For every document(list signatures) in list of documents...
        for signature in signatures:
            # Set the hash = 0
            hash = 0
            # Take the index of document [0..., numDocs - 1]
            index = signatures.index(signature)
            # For every hash in band...
            for j in range(0, numRows):
                # Knuth's Multiplicative Hashing
                hash = (hash * 2654435761 + signature[numRows * i + j])
                # For 32 bit word size
                hash = hash % (1 << 32)
            # If hash bucket doesn't exist then create it.
            if hash not in bucket.keys():
                bucket[hash] = set()
            # Add LSH hash as a docID to bucket.
            bucket[hash].add(docNames[index])
            # Store the LSH bucket signatures
        buckets.append(bucket)
        del bucket

    elapsed = (time.time() - t0)

    print("\nGenerating buckets took %.2fsec" % elapsed)
    return buckets

# ------------LSH select bucket-------------------


def selectCandidatePair(buckets):
    # Create a dictionary of the candidate pairs.
    candidateDocs = {}
    # For every band...
    for bucket in buckets:
        # For every bucket...
        for hash in bucket:
            # If bucket have more than 1 item.
            if (len(bucket[hash]) > 1):
                # For every docID in bucket.
                for docID in bucket[hash]:
                    # If docID doesn't exist in bucket name(key),
                    # create bucket.
                    if docID not in candidateDocs.keys():
                        candidateDocs[docID] = []
                    # Union sets with other bands
                    candidateDocs[docID] = list(
                        set(bucket[hash]) | set(candidateDocs[docID]))
                    # Remove docID which is name of bucket
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
