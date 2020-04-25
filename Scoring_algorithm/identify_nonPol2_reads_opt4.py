from __future__ import division, print_function
import numpy as np
import sys

def scoring_PV(input_table="/home/shriyas/cifs-lab/Shriya/Samples/Masa_RIP_terminator/read_Pol_assignment_scoring/scoring/top_200000_input.bed", coverage_table="input_chr_loci_counts_nonzero.bed",output_name="scores.bed"):
    '''This function takes the following input parameters:
        1) table name of the mapped dedup sorted bed file with 7columns:1-chr#,2-start,3-stop,4-read name,5-random,6-strand,7-gene overlap binary
      	2) table name of the file with the genome-wide coverage of each dataset (Col0,Nrpe1): 1-chr#, 2-position, 3-strand, 4-Col0 count, 5-nrpe1 count
        '''
    reads_list = [] ##reading the reads bed file as a list
    with open(input_table)as f:
        for line in f:
            reads_list.append(line.strip().split())

    cov_list = [] ##reading the genome coverage file as a list
    with open(coverage_table)as f:
        for line in f:
            cov_list.append(line.strip().split())
    j=0
    j2=0
    scores=np.zeros(len(reads_list)) ## declaring output scores table with each row representing the score for each read
    for i in range(len(reads_list)): ## loop to iterate through each line of read file
        col0=0 ## variable for calculating the col0 coverage over the read
        nrpe1=0 ## variable for calculating the nrpe1 coverage over the read
        opp_col0=0 ## variable for calculating the col0 coverage on opposite strand of the read
        opp_nrpe1=0 ## variable for calculating the nrpe1 coverage on opposite strand of the read
        count1=int(reads_list[i][2])-int(reads_list[i][1]) ## number of bases the read spans (same orientation)
	while j<len(cov_list):
	    if int(reads_list[i][0])==int(cov_list[j][0]) and int(reads_list[i][1])<=int(cov_list[j][1]) and int(reads_list[i][2])>=int(cov_list[j][1]):  ## check for chr, coverage position between read start and stop
                if reads_list[i][5]==cov_list[j][2]: ## check if coverage locus and read have same strandedness
                    col0+=int(cov_list[j][3]) ## count coverage on same strand of read in Col0
                    nrpe1+=int(cov_list[j][4]) ## count coverage on same strand of read in nrpe1
                else: ## for coverage on opposite strand compared to read
                    opp_col0+=int(cov_list[j][3]) ## count coverage on opposite strand of read in Col0
                    opp_nrpe1+=int(cov_list[j][4]) ## count coverage on opposite strand of read in nrpe1
		if (j-1)>=0 and int(reads_list[i][1])>int(cov_list[j-1][1]): ##
		    j2=j
		elif (j-1)>=0 and int(reads_list[i][0])<int(cov_list[j-1][0]):
		    j2=j
		elif j==0:
		    j2=0
	    elif int(reads_list[i][0])==int(cov_list[j][0]) and int(reads_list[i][1])<=int(cov_list[j][1]) and int(reads_list[i][2])<int(cov_list[j][1]):
	        j=j2
                break
            j+=1
        if count1>0:
            avg_col0=col0/count1
            avg_nrpe=nrpe1/count1
	    avg_col0_opp=opp_col0/count1
            avg_nrpe_opp=opp_nrpe1/count1
            if avg_nrpe<=1:
                scores[i]+=2
            if avg_col0>(4*avg_nrpe):
                scores[i]+=2
            if avg_col0_opp>(2*avg_nrpe_opp):
                scores[i]+=2
            if int(reads_list[i][6])>0:
                scores[i]-=5
            if (avg_col0_opp)<=(0.1*avg_col0):
		scores[i]-=3
    return scores

def main():
    input_table=sys.argv[1] #1st argument: input tsv file name
    coverage_table=sys.argv[2]
    output_name=sys.argv[3]
    np.savetxt(output_name, scoring_PV(input_table, coverage_table,output_name), fmt='%i')

if __name__ == '__main__':
    main()
