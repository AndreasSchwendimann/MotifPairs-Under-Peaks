import sys
class Range:
	def __init__(self,start,end):
		self.start = start
		self.end = end

	def __str__(self):
		return 'Start={0}, End={1}'

class Pair:
	def __init__(self,start1,end1,start2,end2,firstM,motSeq):
		self.start1 = start1
		self.end1 = end1
		self.start2 = start2
		self.end2 = end2
		self.firstM = firstM
		self.motSeq = motSeq

	def __str__(self):
		return 'Start1={0}, End1={1}, Start2={2}, End2={3}, First Motif={4}, Motif Sequence={5}'

class Hit:
	def __init__(self,Mstart,Mend,Pstart,Pend):
		self.Mstart = Mstart
		self.Mend = Mend
		self.Pstart = Pstart
		self.Pend = Pend
	def __str__(self):
		return 'Motif range:{0},{1}; Peak range:{2},{3}'

def readBED(File):
	validChrom = ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']
	BEDInfo = {}
	with open(File) as f:
		line = f.readlines()
		for l in line:
			elements = l.split("\t")
			if elements[0] != "chr":
				chrom = elements[0]
				start = int(elements[1])
				end = int(elements[2])
				if chrom in validChrom:
					if chrom not in BEDInfo.keys():
						BEDInfo[chrom] = [Range(start,end)]
					else:
						BEDInfo[chrom].append(Range(start,end))

	return(BEDInfo)

def readPairs(File):
	Chrom = range(1,25)
	validChrom = ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']
	PairInfo = {}
	with open(File) as f:
		line = f.readlines()
		for l in line:
			elements = l.split(",")
			chrom = validChrom[int(elements[0])-1]
			start1 = int(elements[2])
			end1 = int(elements[3])
			start2 = int(elements[4])
			end2 = int(elements[5])
			firstM = int(elements[1])
			motSeq = elements[7]
			if chrom in validChrom:
				if chrom not in PairInfo.keys():
					PairInfo[chrom] = [Pair(start1,end1,start2,end2,firstM,motSeq)]
				else:
					PairInfo[chrom].append(Pair(start1,end1,start2,end2,firstM,motSeq))
	return(PairInfo)

def main():
	## load BED files for peaks and Motifs
	# 1: peaks as bed file
	# 2: REAREs in format: chr,firstMotif,start1,end1,start2,end2
	# 3: savefile label
	print("Reading of BED files ...")
	print("------- Settings -------")
	peakFile = sys.argv[1]
	peaksChr = readBED(peakFile)
	print("> Peaks: "+peakFile)

	MotifFile = sys.argv[2]
	MotifChr = readPairs(MotifFile)
	print("> Motifs: "+MotifFile)

	## Output file to save hits in
	HitFile = sys.argv[3]+".csv"
	Hits = open(HitFile,"w")

	HitFile1 = sys.argv[3]+"_onlyARE"+".csv"
	HitFile2 = sys.argv[3]+"_onlyRE"+".csv"
	Hits1 = open(HitFile1,"w")
	Hits2 = open(HitFile2,"w")

	print("> Output-Files:")
	print(">> both motifs under peak: "+HitFile)
	print(">> only motif1 under peak: "+HitFile1)
	print(">> only motif2 under peak: "+HitFile2)

	## go through each Chromosome
	print()
	print("------- Analysis -------")
	validChrom = ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']
	for c in validChrom:
		print("Analysis of "+c)
		peakSet = peaksChr[c]
		MotifSet = MotifChr[c]
		counter = 0
		print("# of peaks:\t"+str(len(peakSet)))
		print("# of motifs:\t"+str(len(MotifSet)))
		for r1 in MotifSet:
			for r2 in peakSet:
				# both motifs fully below peak
				if r2.start <= r1.start1 and r2.start <= r1.start2 and r1.end1 <= r2.end and r1.end2 <= r2.end:
					counter += 1
					Hits.write(c+","+str(r2.start)+","+str(r2.end)+","+str(r1.start1)+","+str(r1.end2)+","+r1.motSeq)
				# only the first motif of pair fully below peak
				elif r2.start <= r1.start1 and r2.end >= r1.end1 and (r1.start2 > r2.end or r1.end2 < r2.start):
					if r1.firstM == 1:
						Hits1.write(c+","+str(r2.start)+","+str(r2.end)+","+str(r1.start1)+","+str(r1.end2)+","+r1.motSeq)
					elif r1.firstM == 2:
						Hits2.write(c+","+str(r2.start)+","+str(r2.end)+","+str(r1.start1)+","+str(r1.end2)+","+r1.motSeq)
				# only the second motif of pair fully below peak
				elif r2.start <= r1.start2 and r2.end >= r1.end2 and (r1.end1 < r2.start or r1.start1 > r2.end):
					if r1.firstM == 1:
						Hits1.write(c+","+str(r2.start)+","+str(r2.end)+","+str(r1.start1)+","+str(r1.end2)+","+r1.motSeq)
					elif r1.firstM == 2:
						Hits2.write(c+","+str(r2.start)+","+str(r2.end)+","+str(r1.start1)+","+str(r1.end2)+","+r1.motSeq)
		print("# of hits:\t"+str(counter))

	Hits.close()
	Hits1.close()
	Hits2.close()

main()
