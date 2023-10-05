# -----------------------------------------------------------------------------------------------
# 
#
#
# Jonathan Holmes
#
# Permission to use, copy, modify, and/or distribute this software or any part thereof for any
# purpose with or without fee is hereby granted provided that:
#     (1) the original author is credited appropriately in the source code
#         and any accompanying documentation
# and (2) that this requirement is included in any redistribution.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL
# THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
# NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
# OF THIS SOFTWARE.
#
# E-mail: jh795@leicester.ac.uk
# ------------------------------------------------------------------------------------------------
# 
#
# 		Generates ON state data from PSanalysis results files
#
# 		python3 GeneStateExtraction.py 
#
# ------------------------------------------------------------------------------------------------


## NOTE: This code currently only works for polyN repeat tracts no functionality for di/tri/tetra/etc has been implemented here.

# Input file
file = "Example_File.txt"

# generates output file
output = open(file.replace(".txt","_output.txt"),"w")

# Functions

# This holds the tract length data for ON states for each gene tested: change this to suit input genes
def state(gene,tract):
	track_data = {"capa":"11","cj0031":"9","cj0046":"11","cj1296":"10","cj1318":"11","cj1321":"10","cj1421":"9","cj1422":"9","cj1426":"10","cj1429":"10"}

	status = 0

	if tract == track_data[gene]:
		status = status + 1 

	return str(status)

# The slippage resolved from single colony data is held here, to calculate slippage refer to: https://www.sciencedirect.com/science/article/pii/S2215016123003886
def resolve_slippage(new_line):

	prev_slip = {"cj1429":1.396491035,"cj1296":1.067412549,"cj1321":5.949057501,"cj1426":9.548569328,"cj0031":2.524190745,"cj1421":3.651462177,"cj1318":4.870101745,"cj0046":2.955134747,"cj1422":0.353581791,"capa":0.201128939}
	post_slip = {"cj1429":5.170100643,"cj1296":0.227146598,"cj1321":16.69435917,"cj1426":6.983350521,"cj0031":14.13813188,"cj1421":9.812889814,"cj1318":6.188596157,"cj0046":24.0502116,"cj1422":4.811977806,"capa":0.539651872}
	prev_value =  str(float(new_line[5])*((100 - prev_slip[new_line[1]])/100))
	post_value =  str(float(new_line[7])*((100 - post_slip[new_line[1]])/100))

	new_line[5] = prev_value 
	new_line[7] = post_value


	return new_line




## removes header
header = open(file).read().split("\n")[0].split("\t")[:-1]

# Holder output data
output_matrix = []

# Holds known header variables
output_header = ["Well","Gene","Major_Tract","Major_Area","Minor_Tract_Prev","Minor_Area_Prev","Minor_Tract_Post","Minor_Area_Post"]
output.write("\t".join(output_header) + "\n")

# Runs through PSA input file
for line in open(file).read().split("\n")[1:]:
		array = line.split("\t")
		if len(array) != 1:
		
			
			new_line = [array[0],array[1],state(array[1],array[11]),array[4]]
			
			## Checks for major and minor peaks

			if array[6] != "N":
				minor_prev = [state(array[1],str(int(array[11])-1)),array[6]]
			if array[9] != "N":
				minor_post = [state(array[1],str(int(array[11])+1)),array[9]]
			if array[6] == "N":
				minor_prev = [state(array[1],str(int(array[11])-1)),"0"]
			if array[9] == "N":
				minor_post = [state(array[1],str(int(array[11])+1)),"0"]			

			
			new_line = new_line + minor_prev + minor_post

			# Resolves slippage to modulate peak heights for post a prior sizes

			corrected_line = resolve_slippage(new_line)
	
			
			# Converts tract lengths to ON/OFF states		

			ON = 0
			OFF = 0

			for i in range(0,len(corrected_line)):
				if corrected_line[i] == "0":
					OFF = OFF + float(corrected_line[i + 1])
				if corrected_line[i] == "1":
					ON = ON + float(corrected_line[i + 1])
			ON_per = (ON/(ON + OFF))*100
			
			corrected_line.append(str(ON_per))

			# appends results to output
	
			output_matrix.append("\t".join(corrected_line))

			


# Writes results to output file

output.write("\n".join(output_matrix))
output.close()











