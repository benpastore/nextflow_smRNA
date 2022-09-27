
import sys
import os

def calc_ratio_by_tail(tailed_bed, perfect_bed, outfile) : 

	"""
	A function to calculate tailing ratio by tail

	Input : tailed_bed, perfect_bed, outfielname 

	Output : tail ratio by tail
	
	"""

	perfect = 0
	with open(perfect_bed, 'r') as f : 
		for line in f : 
			info = line.strip().split("\t")
			rpm = float(info[4])
			perfect += rpm
	f.close() 

	tailed = 0
	with open(tailed_bed, 'r') as f : 
		for line in f : 
			info = line.strip().split("\t")
			rpm = float(info[4])
			tailed += rpm
	f.close() 

	total = perfect + tailed

	tail_dict = {}
	with open(tailed_bed, 'r') as f : 
		for line in f : 
			info = line.strip().split("\t")
			rpm = float(info[4])
			tail = str(info[6])

			if len( set(tail) ) == 1 : 
				key = list(set(tail))[0]

				if key not in tail_dict.keys() : 
					tail_dict[key] = 100*(rpm / total)
				else : 
					tail_dict[key] += 100*(rpm / total)
			
			else : 
				key = 'other'
				
				if key not in tail_dict.keys() : 
					tail_dict[key] = 100*(rpm / total)
				else : 
					tail_dict[key] += 100*(rpm / total)

	f.close() 

	if os.path.exists(outfile) : 
		os.remove(outfile)
	
	out = open(outfile, 'a')

	for key, value in tail_dict.items() : 

		out.write(f"{key}\t{value}\n")

	out.close() 


def main() : 

	tail_bed = sys.argv[1]
	perfect_bed = sys.argv[2]
	outfile = sys.argv[3]

	calc_ratio_by_tail(tail_bed, perfect_bed, outfile)

if __name__ == "__main__" : 

	main()







	