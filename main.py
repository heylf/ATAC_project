import sys

from scipy import stats

import statsmodels.stats.multitest as multi
import numpy as np
import time
import multiprocessing
import argparse
import os

B_VALUE_LIST = dict()
B_ID_LIST = dict()
B_PEAK_COORDS = dict()
NUM_B_PEAKS = 0

# Function to obtain the number of lines in a file.
def get_line_count(filepath):
    file = open(filepath, "r")
    count = 0
    for line in file:
        count += 1
    file.close()
    return count

# Function calculates the regulatory potential for a region.
def reg_potential(a_coords, b_coords, max_dist):
    # Get the coordinates for region A.
    start_a = int(a_coords[1])
    end_a = int(a_coords[2])

    # Get the coordniates for region B.
    start_b = int(b_coords[1])
    end_b = int(b_coords[2])

    # Intitialize potential.
    potential = 0.0

    # Caluclate the shortest distance between A and B.
    dist = np.min([np.abs(start_a - start_b), np.abs(end_a - start_b), np.abs(start_a - end_b), np.abs(end_a - end_b)])

    # Calculate the relative distance based on the user defined maximal distance.
    relative_dist = dist / max_dist

    # Region has to be closer than the maximal distance.
    if ( relative_dist <= 1.0 ):
        potential = np.exp(-(0.5 + 4 * relative_dist))

    return potential

# Function to do a parallel investigation of a region A to all other regions B.
def init_peak_analysis(i, a_data, potential_dict, significance_dict, t, max_dist, strand):

    print(i)

    # Initialize potential
    potential = 0.0

    # Get coordinates for A.
    coords_a = a_data[0]
    chr_a = coords_a[0]
    strand_a = coords_a[3]

    a_values = a_data[1]

    if ( a_values.count(a_values[0]) == len(a_values) ):
        print("[ERROR] Something is wrong with you list A. An entry in list A has constant values for all samples.")

    # Go over all regions of set B.
    for j in range(0, NUM_B_PEAKS):

        # Get coordinates of B
        coords_b = B_PEAK_COORDS[B_ID_LIST[j]]
        chr_b = coords_b[0]
        strand_b = coords_b[3]

        # If strandness does not matter.
        if ( not strand ):
            strand_a = strand_b

        # Check if both region are on the same chr and strand.
        if ( chr_a == chr_b and strand_a == strand_b ):

            b_values = B_VALUE_LIST[j]

            corr = 0.0
            pval = 1.0

            # Check if the valuesi B changes or if it is just empty or constant over all samples.
            if ( b_values.count(b_values[0]) != len(b_values) ):
                # Calculate the correlation and pvalue between region A and B.
                corr, pval = stats.spearmanr(a_values, B_VALUE_LIST[j])

            # Check if the correlation and pval are not defined.
            if ( np.isnan(corr) ):
                corr = 0.0
            if ( np.isnan(pval) ):
                pval = 1.0

            # Only consider regions B if they pass the pval threshold.
            if ( pval < t ):
                potential += reg_potential(coords_a, coords_b, max_dist)

    # Insert potential and significance in the dictionary.
    potential_dict[i] = potential
    significance_dict[i] = (1 - stats.norm.cdf(potential))


if __name__ == "__main__":

    ####################
    ##   ARGS INPUT   ##
    ####################

    tool_description = """
        The tool can calculate the regulation potential for a reagion A (e.g., ATAC-Seq peaks) to all other regions of a set B
        (e.g., RNA-Seq peaks). 
        """

    # parse command line arguments
    parser = argparse.ArgumentParser(description=tool_description,
                                     usage='%(prog)s [-h] [options] -a *.tabular -b *.tabular -c *.bed -d *.bed -e *.tabular',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # version
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 0.0")

    # positional arguments
    parser.add_argument(
        "-a", "--regions_a_file_path",
        metavar='*.tabular',
        required=True,
        help="Path to the counts for the regions of set A in tabular format (e.g., ATAC-Seq peaks).")
    parser.add_argument(
        "-b", "--regions_b_file_path",
        metavar='*.tabular',
        required=True,
        help="Path to the counts for the regions of set B in tabular format (e.g., RNA-Seq peaks).")
    parser.add_argument(
        "-c", "--regions_a_coords_path",
        metavar='*.bed',
        required=True,
        help="Path to cooridnates of regions of the set A in bed format (e.g., ATAC-Seq peaks).")
    parser.add_argument(
        "-d", "--regions_b_coords_path",
        metavar='*.bed',
        required=True,
        help="Path to cooridnates of regions of the set B in bed format (e.g., RNA-Seq peaks).")
    parser.add_argument(
        "-e", "--blastn_short_path",
        metavar='*.tabular',
        required=True,
        help="Path to BLAST file that contains the TF motifs for your regions in set A.")

    # optional arguments
    parser.add_argument(
        "-o", "--output_folder",
        metavar='path/',
        default=os.getcwd(),
        help="Write results to this path. [Default: Operating Path]")
    parser.add_argument(
        "-t", "--threads",
        metavar='int',
        default=1,
        help="Number of threads for parallelization. [Default: 1]")
    parser.add_argument(
        "--pval_thresh",
        metavar='float',
        default=0.05,
        help="P-value threshold to consider region of set B to calculate the regulatory potential of a region A [Default: 0.05]")
    parser.add_argument(
        "--max_dist",
        metavar='int',
        default=100000,
        help="Maximal distance for the calculation of the regulatory potential of a region A [Default: 100kb]")
    parser.add_argument(
        "--stranded",
        action='store_true',
        help="Set this option if your data is stranded.")

    print("[START]")

    args = parser.parse_args()

    #atac_peaks_path = "/home/florian/Documents/ATAC_Project/atac_peaks.tabular"
    #regions_a_file_path = "/home/florian/Documents/ATAC_Project/atac_peaks_subset.tabular"
    #regions_b_file_path = "/home/florian/Documents/ATAC_Project/rna_seq.tabular"

    output_path = "/home/florian/Documents/ATAC_Project/"

    # Get number of peaks of region set A and B.
    num_a_peaks = get_line_count(args.regions_a_file_path) - 1
    NUM_B_PEAKS = get_line_count(args.regions_b_file_path) - 1

    print("[NOTE] Analysing {} of regions A".format(num_a_peaks))
    print("[NOTE] Analysing {} of regions B".format(NUM_B_PEAKS))

    a_value_list = [ [0.0] ] * num_a_peaks
    a_id_list = [""] * num_a_peaks

    # Get ids and values of regions set A.
    regions_a_file = open(args.regions_a_file_path, "r")
    headline_a = regions_a_file.readline()
    for i in range(0, num_a_peaks):
        line_atac = regions_a_file.readline()
        x = line_atac.strip("\n").split("\t")
        a_id_list[i] = x[0]
        a_value_list[i] = [float(val) for val in x[1:]]
    regions_a_file.close()

    # Get ids and values of regions set B.
    regions_b_file = open(args.regions_b_file_path, "r")
    headline_b = regions_b_file.readline()
    for i in range(0, NUM_B_PEAKS):
        line_rna = regions_b_file.readline()
        y = line_rna.strip("\n").split("\t")
        B_ID_LIST[i] = y[0]
        B_VALUE_LIST[i] = [float(val) for val in y[1:]]
    regions_b_file.close()

    id_delimiter = "|"

    # Get coordinates for region set A.
    #regions_a_coords_path = "/home/florian/Documents/ATAC_Project/atac_peaks.bed"
    region_a_coords_file = open(args.regions_a_coords_path, "r")

    regions_a_coords = dict()
    regions_a_coords_to_id = dict()

    for line in region_a_coords_file:
        values = line.strip("\n").split("\t")
        if ( values[3] not in regions_a_coords ):
            coords = [values[0], values[1], values[2], values[5]]
            regions_a_coords[values[3]] = coords
            regions_a_coords_to_id[id_delimiter.join(coords)] = values[3]
    region_a_coords_file.close()

    # Get coordinates for region set B.
    #regions_b_coords_path = "/home/florian/Documents/ATAC_Project/mm10_annotation.bed"
    regions_b_coords_file = open(args.regions_b_coords_path, "r")

    for line in regions_b_coords_file:
        values = line.strip("\n").split("\t")
        if ( values[3] not in B_PEAK_COORDS ):
            B_PEAK_COORDS[values[3]] = [values[0], values[1], values[2], values[5]]
    regions_b_coords_file.close()

    # Read blastn file if you want to check also if ATAC-Peaks cover specific transcription factors.
    #blastn_short_path = "/home/florian/Documents/ATAC_Project/blastn_short.tabular"
    blastn_file = open(args.blastn_short_path, "r")

    num_motif_dict = dict()

    # Go over every line in the blast file.
    for line in blastn_file:
        # Get the coordinates of the regions with a TF motif. I will use the coordinates as ids. Becuase blastn
        # Does not uses saves the ids.
        values = line.strip("\n").split("\t")
        coords = values[0].split(id_delimiter)
        coords = coords[1:]

        # Check if strandness is relevant.
        if ( not args.stranded ):
            coords[-1] = "."

        # Add up the number of motifs for a region in the regions set A.
        if ( id_delimiter.join(coords) in regions_a_coords_to_id ):
            id = regions_a_coords_to_id[id_delimiter.join(coords)]
            if( id not in num_motif_dict ):
                num_motif_dict[id] = 0
            else:
                num_motif_dict[id] += 1

    blastn_file.close()

    # Create manager for parallel processing.
    manager = multiprocessing.Manager()
    potential_dict = manager.dict()
    significance_dict = manager.dict()

    print("[NOTE] Using {} threads.".format(args.threads))
    pool = multiprocessing.Pool(int(args.threads))

    print("[NOTE] Filter for a correlation with a p-value of {}.".format(args.pval_thresh))
    print("[NOTE] Setting a maximal distance of {}.".format(args.max_dist))

    init_peak_analysis(39, [regions_a_coords[a_id_list[39]], a_value_list[39]], potential_dict, significance_dict, args.pval_thresh, args.max_dist, args.stranded)
    sys.exit(1)

    start = time.time()
    # Analyse any region in the region set A.
    for i in range(0, num_a_peaks):
        atac_data = [regions_a_coords[a_id_list[i]], a_value_list[i]]
        try:
            pool.apply_async(init_peak_analysis, args=(i, atac_data, potential_dict, significance_dict,
                                                       args.pval_thresh, args.max_dist, args.stranded))
        except:
            print("[ERROR 1] Failed Analysis")

    pool.close()
    pool.join()

    end = time.time()
    print(end - start)

    # Do a p-vlaue correction for multiple hypothesis testing.
    corrected_pvals = multi.multipletests(significance_dict.values(), alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

    # List of improtant ATAC_peaks:
    # All regions A with high correlation --> calculate from them the regulatory potential for the atac peak
    # pvalue from the cdf of a normal distribution
    potential_file = open(output_path + "/potential_file.tsv", "w")
    potential_file.write("ATAC_Peak\tRegulatory_Potential\tp-val\tcorrected_p-val\tNum_Motifs\n")

    for i in range(0, len(a_id_list)):
        peak_id = a_id_list[i]
        if ( peak_id in num_motif_dict ):
            potential_file.write(("{}\t{}\t{}\t{}\t{}\n").format(peak_id, potential_dict[i], significance_dict[i],
                                                             corrected_pvals[1][i], num_motif_dict[peak_id]))
        else:
            potential_file.write(("{}\t{}\t{}\t{}\t{}\n").format(peak_id, potential_dict[i], significance_dict[i],
                                                             corrected_pvals[1][i], 0))

    potential_file.close()

    print("[END]")