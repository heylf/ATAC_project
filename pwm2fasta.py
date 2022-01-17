
import numpy as np


print("[START]")

pcm_path = "/home/florian/Documents/ATAC_Project/HOCOMOCOv11_core_pcms_MOUSE_mono.txt"
output_path = "/home/florian/Documents/ATAC_Project/"

pcm_file = open(pcm_path, "r")
seq_file = open(output_path + "sequences.fa", "w")

nucleotides = ["A", "C" , "G", "T"]
motif_seq = ""
first_seq = 0

for line in pcm_file:
    if ( line.startswith(">") ):
        if ( first_seq != 0 ):
            seq_file.write(motif_seq + "\n")
            motif_seq = ""
        else:
            first_seq = 1
        seq_file.write(line)
    else:
        values = line.strip("\n").split("\t")
        values = [float(x) for x in values]
        sum_count = np.sum(values)
        probs = [ x/sum_count for x in values ]
        motif_seq += nucleotides[int(np.argmax(probs))]


pcm_file.close()
seq_file.close()