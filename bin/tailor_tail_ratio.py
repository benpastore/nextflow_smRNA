
import sys
import os

def calc_tail_ratio(tail_bed, perfect_bed, outfile) : 

	"""
	A function to calculate the ratio of tailed to perfectly matched reads, rpm MUST be in the 5th column 

	Input : tail_bed, perfect_bed and outfile name

	Output : tailed rpm, perfect rpm, ratio

	"""

	tailed = 0
	with open(tail_bed, 'r') as f : 
		for line in f : 
			info = line.strip().split("\t")
			rpm = float(info[4])
			tailed += rpm
	f.close() 

	perfect = 0
	with open(perfect_bed, 'r') as f : 
		for line in f : 
			info = line.strip().split("\t")
			rpm = float(info[4])
			perfect += rpm
	f.close() 

	ratio = 100 * ( tailed / (tailed + perfect) )

	if os.path.exists(outfile) : 
		os.remove(outfile)
	
	out = open(outfile, 'a')

	out.write(f"{tailed}\t{perfect}\t{ratio}\n")
	out.close()

def main() : 

	tail_bed = sys.argv[1]
	perfect_bed = sys.argv[2]
	outfile = sys.argv[3]
	calc_tail_ratio(tail_bed, perfect_bed, outfile)

if __name__ == "__main__" : 

	main()



