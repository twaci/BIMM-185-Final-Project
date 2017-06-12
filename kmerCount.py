def KmerCount(file, kv):

    infile = open(file,'r')
    genome = infile.read()

    kmerList = {}
    k = int(kv)

    for i in range(0, len(genome) - k +1):
        kmer = str(genome[i:i+k])
        if kmer[len(kmer)-2:len(kmer)] == "GG":
            if kmer not in kmerList:
                kmerList[kmer] = 1
            else:
                kmerList[kmer] = kmerList[kmer] + 1

    outfile = open('ecoli_22mers.txt', 'w')

    for index, i in enumerate(kmerList):
        if kmerList[i] == 1:
            outfile.write(i)
            if index != len(kmerList) - 1:
                outfile.write("\n")

    outfile.close()

KmerCount("only_sequence.fasta", 25)
