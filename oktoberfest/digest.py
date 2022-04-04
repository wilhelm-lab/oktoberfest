from __future__ import print_function

import sys
import csv
import itertools
import collections

import numpy as np

def main(args):
    args = parseArgs()    
    # python digest.py --enzyme trypsinp --cleavages 0 --fasta /home/matthewt/data/MaxLFQ_benchmark/ups.fasta --prosit_input ~/data/Prosit_test/ups_prosit_input.csv
    if args.prosit_input:
        writer = getTsvWriter(args.prosit_input, delimiter = ',')
        writer.writerow("modified_sequence,collision_energy,precursor_charge".split(","))

        print(args.prosit_input)
        print(args.fasta)
        prositInputFileWithProteins = args.prosit_input.replace(".csv", "_with_proteins.csv")
        writerWithProteins = getTsvWriter(prositInputFileWithProteins, delimiter = ',')
        writerWithProteins.writerow("modified_sequence,collision_energy,precursor_charge,protein".split(","))
        
        pre, not_post = getCleavageSites(args.enzyme)    
        for peptide, proteins in getPeptideToProteinMap(
                    args.fasta, 
                    db = 'concat', 
                    digestion = args.digestion, 
                    min_len = args.min_length, 
                    max_len = args.max_length, 
                    pre = pre, 
                    not_post = not_post, 
                    miscleavages = args.cleavages, 
                    methionineCleavage = True, 
                    specialAAs = list(args.special_aas), 
                    useHashKey = False).items():
            if validPrositPeptide(peptide):
                for charge in [2,3,4]:
                    writer.writerow([peptide, 30, charge])
                    writerWithProteins.writerow([peptide, 30, charge, proteins[0]])

    # python digest.py --fasta /media/kusterlab/internal_projects/active/Mouse_proteome/stuff/10090_UP000000589_UniProtKB_Mouse_CanIso_2018_03_27.fasta --peptide_protein_map ../data/fasta/10090_UP000000589_UniProtKB_Mouse_CanIso_2018_03_27.peptide_to_protein_map.txt            
    if args.peptide_protein_map:
        with open(args.peptide_protein_map + '.params.txt', 'w') as f:
            f.write(" ".join(sys.argv))
        writer = getTsvWriter(args.peptide_protein_map, delimiter = '\t')
        
        pre, not_post = getCleavageSites(args.enzyme)    
        for peptide, proteins in getPeptideToProteinMap(args.fasta, db = 'concat', digestion = args.digestion, min_len = args.min_length, max_len = args.max_length, pre = pre, not_post = not_post, miscleavages = args.cleavages, methionineCleavage = True, specialAAs = list(args.special_aas), useHashKey = False).items():
            writer.writerow([peptide, ";".join(proteins)])
    
    if args.ibaq_map:
        writer = getTsvWriter(args.ibaq_map, delimiter = '\t')
        
        numPeptidesPerProtein = getNumIbaqPeptidesPerProtein(args)
        
        for protein, numPeptides in numPeptidesPerProtein.items():
            writer.writerow([protein, numPeptides])
        
    #writeProteinToGeneMap(fastaFile, outputFile)


def validPrositPeptide(peptide):
    return len(peptide) <= 30 and "U" not in peptide and "X" not in peptide


def parseArgs():
    import argparse
    apars = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('--fasta', default="/media/kusterlab/internal_projects/active/ProteomeTools/TMT/datasets/HumanYeast_Mix/Homo_Sapiens.fasta", metavar = "F", required = False,
                                         help='''Fasta file used as input
                                                    ''')
    
    apars.add_argument('--prosit_input', default="prosit_input.csv", metavar = "M", required = False,
                                         help='''Path to file where to write the prosit input file.
                                                    ''')
                                                    
    apars.add_argument('--peptide_protein_map', default=None, metavar = "M", required = False,
                                         help='''Write mapping from peptides to all its proteins to 
                                                         the specified file.
                                                    ''')
    
    apars.add_argument('--ibaq_map', default=None, metavar = "M", required = False,
                                         help='''Write number of peptides per protein to the specified 
                                                         file that meet the iBAQ criteria
                                                         (6 <= pepLen <= 30, no miscleavages).
                                                    ''')
                                                                                            
    addArguments(apars)
                                                                                                    
    # ------------------------------------------------
    args = apars.parse_args()
    
    return args

    
def addArguments(apars):    
    apars.add_argument('-e', '--enzyme', default = "trypsin", metavar='E',
                                         help='''Type of enzyme "no_enzyme","elastase","pepsin",
                                                         "proteinasek","thermolysin","chymotrypsin",
                                                         "lys-n","lys-c","arg-c","asp-n","glu-c","trypsin",
                                                         "trypsinp".
                                                    ''')        
                                                        
    apars.add_argument('-c', '--cleavages', default = 2, metavar='C', type=int,
                                         help='''Number of allowed miss cleavages used in the search 
                                                         engine (Only valid when using option -F).
                                                    ''')
    
    apars.add_argument('-l', '--min-length', default = 7, metavar='L', type=int,
                                         help='''Minimum peptide length allowed used in the search 
                                                         engine (Only valid when using option -F).
                                                    ''')
    
    apars.add_argument('-t', '--max-length', default = 60, metavar='L', type=int,
                                         help='''Maximum peptide length allowed used in the search 
                                                         engine (Only valid when using option -F).
                                                    ''')                                     
    
    apars.add_argument('--special-aas', default = 'KR', metavar='S', 
                                         help='''Special AAs that MaxQuant uses for decoy generation.
                                                    ''')
    
    apars.add_argument('--digestion', default = 'full', metavar='D', 
                                         help='''Digestion mode ('full', 'semi' or 'none').
                                                    ''')


def writeProteinToGeneMap(fastaFile, outputFile):
    writer = csv.writer(open(outputFile, 'w'), delimiter = '\t')
    for proteinName, _ in readFastaTide(fastaFile, db = "target"):
        proteinId = proteinName.split("|")[1]
        geneId = proteinName.split("|")[2].split(" ")[0]
        writer.writerow([proteinId, geneId])



parseUntilFirstSpace = lambda x : x.split(" ")[0]
parseProteinNameFunc = lambda x : " ".join(x.split(" OS=")[0].split(" ")[1:])
parseGeneNameFunc = lambda x : x.split(" GN=")[1].split(" ")[0] if "GN=" in x else ""
#parseId = lambda x : x.replace(" ", "") # default MaxQuant 
#parseId = lambda x : x.split("|")[1] # extract Uniprot ID, for PrDB runs (mouse proteome, Dongxue original)

def parseUniProtId(fastaId):
    proteinId = parseUntilFirstSpace(fastaId)
    if "|" in proteinId:
        return proteinId.split("|")[1]
    else:
        return proteinId

def readFastaProteins(filePath, db = "concat", parseId = parseUntilFirstSpace, parseProteinName = parseProteinNameFunc, parseGeneName = parseGeneNameFunc):
    name, seq = None, []
    with open(filePath, 'r') as fp:
        for line in itertools.chain(fp, [">"]):
            line = line.rstrip()
            if line.startswith(">"):
                if db in ["target", "concat"]:
                    yield (parseId(line[1:]), parseProteinName(line[1:]), parseGeneName(line[1:]), line[1:])
                
                if db in ["decoy", "concat"]:
                    yield ("REV__" + parseId(line[1:]), "", "REV__" + parseGeneName(line[1:]), "")


def getProteinAnnotations(fastaFile, parseId):
    proteinAnnotations = dict()
    
    if not fastaFile:
        return proteinAnnotations
    
    for protein, proteinName, geneName, fastaHeader in readFastaProteins(fastaFile, parseId = parseId):
        if protein not in proteinAnnotations:
            proteinAnnotations[protein] = (proteinName, geneName, fastaHeader)
    
    return proteinAnnotations


def hasGeneNames(proteinAnnotations):
    counts = sum(1 for _, geneName, _ in proteinAnnotations.values() if len(geneName) > 0)
    return counts / len(proteinAnnotations) > 0.5


def readFastaTide(filePath, db = "target", parseId = parseUntilFirstSpace):
    readFastaMaxQuant(filePath, db, parseId, specialAAs = [], decoyPrefix = 'decoy_')
 

def readFastaMaxQuant(filePath, db = "target", parseId = parseUntilFirstSpace, specialAAs = ['K', 'R'], decoyPrefix = 'REV__'):
    if db not in ["target", "decoy", "concat"]:
        sys.exit("unknown db mode: %s" % db)
    
    hasSpecialAAs = len(specialAAs) > 0
    name, seq = None, []
    with open(filePath, 'r') as fp:
        for line in itertools.chain(fp, [">"]):
            line = line.rstrip()
            if line.startswith(">"):
                if name: 
                    seq = "".join(seq)
                    if db in ["target", "concat"]:
                        yield (name, seq)
                    
                    if db in ["decoy", "concat"]:
                        revSeq = seq[::-1]
                        if hasSpecialAAs:
                            revSeq = swapSpecialAAs(revSeq, specialAAs)
                        yield (decoyPrefix + name, revSeq)
                    
                if len(line) > 1:
                    name, seq = parseId(line[1:]), []
            else: seq.append(line)


#from . import digestfast
#readFasta = digestfast.readFastaMaxQuant
readFasta = readFastaMaxQuant


# swaps the specialAAs with its preceding amino acid, as is done in MaxQuant
# e.g. specialAAs = ['R', 'K'] transforms ABCKDEFRK into ABKCDERKF
def swapSpecialAAs(seq, specialAAs):
    seq = list(seq)
    for i in range(1, len(seq)):
        if seq[i] in specialAAs:
            swapPositions(seq, i, i-1)
    seq = "".join(seq)
    return seq


def swapPositions(seq, pos1, pos2):         
    seq[pos1], seq[pos2] = seq[pos2], seq[pos1] 


def getProteinIds(filePath):
    proteinIds = list()
    for proteinId, _ in readFasta(filePath):
        proteinIds.append(proteinId)
    return set(proteinIds)


def getProteinSequences(filePath, parseId):
    proteinSequences = dict()
    for proteinId, proteinSeq in readFasta(filePath, db = 'concat', parseId = parseId):
        if proteinId not in proteinSequences: # keep only first sequence per identifier
            proteinSequences[proteinId] = proteinSeq
    return proteinSequences


def filterFastaFile(fastaFile, filteredFastaFile, proteins):
    with open(filteredFastaFile, 'w') as f:
        for prot, seq in readFasta(fastaFile):
            if prot in proteins:
                f.write('>' + prot + '\n' + seq + '\n')
                #f.write('>decoy_' + prot + '\n' + seq[::-1] + '\n')


def getPeptides(fastaFile, db = "concat", min_len = 6, max_len = 50, pre = ['K', 'R'], not_post = ['P'], digestion = 'full', miscleavages = 0, methionineCleavage = True):
    for protein, seq in readFasta(fastaFile, db):
        for peptide in getDigestedPeptides(seq, min_len, max_len, pre, not_post, digestion, miscleavages, methionineCleavage):
            yield peptide


#@profile
def getDigestedPeptides(seq, min_len = 6, max_len = 50, pre = ['K', 'R'], not_post = ['P'], digestion = 'full', miscleavages = 0, methionineCleavage = True):
    if digestion == 'none':
        yield from nonSpecificDigest(seq, min_len, max_len)
    elif digestion == 'semi':
        yield from semiSpecificDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)
    else:
        yield from fullDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)


def nonSpecificDigest(seq, min_len, max_len):
    lenS = len(seq)
    for i in range(lenS + 1):
        for j in range(i + min_len, min(lenS+1, i + max_len + 1)):
            if j <= lenS:
                yield seq[i:j]


def semiSpecificDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage):
    lenS, starts = len(seq), [0]
    methionineCleavage = methionineCleavage and seq[0] == "M"
    length_accepted = lambda x : x >= min_len and x <= max_len
    
    for i in range(lenS + 1):
        isCleavageSite = (seq[min([lenS-1,i])] in pre and seq[min([lenS-1,i+1])] not in not_post)
        isMethionineCleavageSite = (i == 0 and methionineCleavage)
        if i == lenS or isCleavageSite or isMethionineCleavageSite:
            # peptides with enzymatic C-terminal (both enzymatic and non-enzymatic N-terminal)
            start = starts[0]
            for j in range(start, min([i+1, lenS])):
                lenP = min([i, lenS - 1]) - j + 1
                if length_accepted(lenP):
                    yield (seq[j : i + 1])
            starts.append(i + 1)
            methionineCleaved = int(starts[0] == 0 and methionineCleavage)
            if len(starts) > miscleavages + 1 + methionineCleaved or i == lenS:        
                starts = starts[1 + methionineCleaved:]
        else: # peptides with non enzymatic C-terminal
            for start in starts:
                lenP = i - start + 1
                if length_accepted(lenP) and i + 1 not in starts:
                    yield (seq[start : i + 1])


def fullDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage):
    lenS, starts = len(seq), [0]
    methionineCleavage = methionineCleavage and seq[0] == "M"
    length_accepted = lambda x : x >= min_len and x <= max_len
    
    cleavageSites = [0] if methionineCleavage else []    
    cleavageSites.extend([i for i in range(lenS) if seq[i] in pre and seq[min([lenS-1,i+1])] not in not_post])
    cleavageSites.append(lenS)
    for i in cleavageSites:
        for start in starts:
            lenP = i - start + 1
            if length_accepted(lenP):
                yield (seq[start : i + 1])
        starts.append(i + 1)
        methionineCleaved = int(starts[0] == 0 and methionineCleavage)
        if len(starts) > miscleavages + 1 + methionineCleaved:        
            starts = starts[1 + methionineCleaved:]


def get_peptide_to_protein_map(args, parseId):
    if args.fasta:
        print("In silico protein digest for peptide-protein mapping")
        pre, not_post = getCleavageSites(args.enzyme)
        if args.enzyme == "no_enzyme":
            args.digestion = "none"
        
        peptideToProteinMap = getPeptideToProteinMap(
                args.fasta, db = 'concat', digestion = args.digestion, 
                min_len = args.min_length, max_len = args.max_length, 
                pre = pre, not_post = not_post, 
                miscleavages = args.cleavages, methionineCleavage = True, 
                specialAAs = list(args.special_aas), useHashKey = (args.digestion == "none"), parseId = parseId)
    elif args.peptide_protein_map:
        print("Loading peptide to protein map")
        peptideToProteinMap = getPeptideToProteinMapFromFile(args.peptide_protein_map, useHashKey = False)
    else:
        sys.exit("No fasta or peptide to protein mapping file detected, please specify either the --fasta or --peptide_protein_map flags")
    
    return peptideToProteinMap


#function_to_be_profiled = profile(digestfast.getDigestedPeptides) 
#function_to_be_profiled = profile(digestfast.readFastaMaxQuant)
#@profile
def getPeptideToProteinMap(fastaFile, db = "concat", min_len = 6, max_len = 52, pre = ['K', 'R'], not_post = ['P'], digestion = 'full', miscleavages = 2, methionineCleavage = True, useHashKey = False, specialAAs = ['K', 'R'], parseId = parseUntilFirstSpace):
    peptideToProteinMap = collections.defaultdict(list)
    proteinToSeqMap = dict()
    for proteinIdx, (protein, seq) in enumerate(readFasta(fastaFile, db, parseId, specialAAs = specialAAs)):
        if (proteinIdx+1) % 10000 == 0:
            print("Digesting protein", proteinIdx+1)
        seenPeptides = set()
        proteinToSeqMap[protein] = seq
        #for peptide in digestfast.getDigestedPeptides(seq, min_len, max_len, pre, not_post, digestion, miscleavages, methionineCleavage):
        for peptide in getDigestedPeptides(seq, min_len, max_len, pre, not_post, digestion, miscleavages, methionineCleavage):
            peptide = peptide
            if useHashKey:
                hashKey = peptide[:6]
            else:
                hashKey = peptide
            if hashKey not in seenPeptides:
                seenPeptides.add(hashKey)
                peptideToProteinMap[hashKey].append(protein)
    
    if useHashKey:
        return (peptideToProteinMap, proteinToSeqMap)
    else:
        return peptideToProteinMap


def getPeptideToProteinMapFromFile(peptideToProteinMapFile, useHashKey = False):
    if useHashKey:
        print("Hash key not supported yet, continuing without hash key...")
        useHashKey = False
    peptideToProteinMap = collections.defaultdict(list)
    reader = getTsvReader(peptideToProteinMapFile)
    for i, row in enumerate(reader):
        if (i+1) % 1000000 == 0:
            print("Processing peptide", i+1)
        
        peptide, proteins = row[0], row[1].split(";")
        if useHashKey:
            sys.exit("Hash key not supported yet...")
            hashKey = peptide[:6]
        else:
            hashKey = peptide
        for protein in proteins:
            peptideToProteinMap[hashKey].append(protein)
    return peptideToProteinMap

    
def getProteins(peptideToProteinMap, peptide):
    peptide = peptide#.replace("I", "L")
    if len(peptideToProteinMap) == 2:
        hashKey = peptide[:6]
        proteins = list()
        if hashKey in peptideToProteinMap[0]:
            for protein in peptideToProteinMap[0][hashKey]:
                #TODO: This does not work correctly for full or partial digestion, since we might find the peptide with the wrong number of enzymatic terminals
                if peptide in peptideToProteinMap[1][protein]:
                    proteins.append(protein)
            proteins = sorted(proteins)
        #else:
        #    print("WARNING: Could not find peptide " + peptide + " in fasta database")
        return proteins
    else:
        return peptideToProteinMap.get(peptide, [])


def getAllProteins(peptideToProteinMap):
    seenProteins = set()
    if len(peptideToProteinMap) == 2:
        for _, proteins in peptideToProteinMap[0].items():
            for protein in proteins:
                if protein not in seenProteins:
                    seenProteins.append(protein)
    else:
        for _, proteins in peptideToProteinMap.items():
            for protein in proteins:
                if protein not in seenProteins:
                    seenProteins.append(protein)
    return list(seenProteins)


def getIbaqPeptideToProteinMap(args):
    pre, not_post = getCleavageSites(args.enzyme)
    return getPeptideToProteinMap(args.fasta, db = 'concat', digestion = 'full', min_len = max([6, args.min_length]), max_len = min([30, args.max_length]), pre = pre, not_post = not_post, miscleavages = 0, methionineCleavage = False, specialAAs = list(args.special_aas))


def getNumIbaqPeptidesPerProtein(args):
    peptideToProteinMapIbaq = getIbaqPeptideToProteinMap(args)
    return getNumPeptidesPerProtein(peptideToProteinMapIbaq)


def getNumPeptidesPerProtein(peptideToProteinMap):
    numPeptidesPerProtein = collections.defaultdict(int)
    for peptide, proteins in peptideToProteinMap.items():
        for protein in proteins:
            numPeptidesPerProtein[protein] += 1
    
    return numPeptidesPerProtein

    
def getCleavageSites(enzyme):
    if enzyme == "trypsinp":
        pre = ['K', 'R']
        not_post = []
    elif enzyme == "trypsin":
        pre = ['K', 'R']
        not_post = ['P']
    elif enzyme == "no_enzyme":
        pre = []
        not_post = []
    elif enzyme == "chymotrypsin":
        pre = ['F', 'W', 'Y', 'L']
        not_post = ['P']
    elif enzyme == "proteinasek":
        pre = ['A', 'E', 'F', 'I', 'L', 'T', 'V', 'W', 'Y']
        not_post = []
    elif enzyme == "elastase":
        pre = ['L', 'V', 'A', 'G']
        not_post = ['P']
    elif enzyme == "lys-c":
        pre = ['K']
        not_post = ['P']
    elif enzyme == "arg-c":
        pre = ['R']
        not_post = ['P']
    elif enzyme == "glu-c":
        pre = ['E']
        not_post = ['P']
    elif enzyme == 'v8-de':
        pre = ['N', 'D', 'E', 'Q']
        not_post = ['P']
    else:
        sys.exit("Enzyme", enzyme, "not implemented yet")
    
    return pre, not_post

    
def isEnzymatic(aa1, aa2, pre = ['K', 'R'], not_post = ['P'], methionineCleavage = True):
    return aa1 == "-" or aa2 == "-" or (aa1 in pre and aa2 not in not_post) or (methionineCleavage and aa1 == "M")


def hasMiscleavage(seq, pre = ['K', 'R'], not_post = ['P']):
    for i in range(len(seq) - 1):
        if isEnzymatic(seq[i], seq[i+1], pre, not_post):
            return True
    return False


def getTsvReader(filename, delimiter = '\t'):
    # Python 3
    if sys.version_info[0] >= 3:
        return csv.reader(open(filename, 'r', newline = ''), delimiter = delimiter)
    # Python 2
    else:
        return csv.reader(open(filename, 'rb'), delimiter = delimiter)


def getTsvWriter(filename, delimiter = '\t'):
    # Python 3
    if sys.version_info[0] >= 3:
        return csv.writer(open(filename, 'w', newline = ''), delimiter = delimiter)
    # Python 2
    else:
        return csv.writer(open(filename, 'wb'), delimiter = delimiter)
        
if __name__ == "__main__":
    main(sys.argv[1:])
