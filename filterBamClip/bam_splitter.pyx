import pysam
from pysam.calignmentfile cimport AlignmentFile, AlignedSegment
import re
import numpy as np
from numpy cimport ndarray
from cpython cimport bool


cpdef bool qualify_aln(AlignedSegment aln):
    '''
    Check if alignment is properly mapped
    '''
    cdef:
        bool qualify
        int flag

    flag = aln.flag
    qualify = (flag == 99) or (flag == 147) or (flag == 163) or (flag == 83)
    return qualify


cpdef int split_bam(str in_bam, str outputprefix):

    cdef:
        AlignmentFile inbam, out_bam_1, out_bam_2
        AlignedSegment aln
        int count_R1 = 0, count_R2 = 0

    r1_out_bam_file = outputprefix + '_R1.bam'
    r2_out_bam_file = outputprefix + '_R2.bam'
    with pysam.AlignmentFile(in_bam,'rb') as inbam:
        with pysam.Samfile(r1_out_bam_file,'wb',template = inbam) as out_bam_1,\
                pysam.Samfile(r2_out_bam_file,'wb',template = inbam) as out_bam_2:
            for aln in inbam:
                if not aln.is_unmapped:
                    if aln.is_read1:
                        out_bam_1.write(aln)
                        count_R1 += 1
                    if aln.is_read2:
                        out_bam_2.write(aln)
                        count_R2 += 1
    print 'Done splitting'
    print 'Read 1: %i' %(count_R1)
    print 'Read 2: %i' %(count_R2)
    return 0

cpdef ndarray split_cigar(cigar_string):
    '''
    split cigar string to numpy array
    return:
        [ [list of numbers],
          [list of cigar operators correspongs to the numbers] ]
    '''
    cdef:
        ndarray cigar_array

    cigar_numbers = re.findall(r'[0-9]+', cigar_string)
    cigar_operator = re.findall(r'[A-Z]', cigar_string)
    cigar_array = np.array([cigar_numbers, cigar_operator])
    return cigar_array

cdef int check_aln(int output_count, AlignmentFile outbam, AlignedSegment aln,
                float single_end_thresh, float both_end_thresh):
    '''
    Compare soft clip length and threshold and write
    '''
    cdef:
        int total_clipped, seq_len
        ndarray cigar_array, all_soft_clipped
        int max_single_clipped
        float bet, set

    seq_len = len(aln.query_sequence)
    set = (single_end_thresh * seq_len)
    bet = (both_end_thresh * seq_len)
    cigar_array = split_cigar(str(aln.cigarstring))
    all_soft_clipped = np.array(cigar_array[0][cigar_array[1]=='S'],dtype=np.int16)
    total_clipped = all_soft_clipped.sum()
    max_single_clipped =  all_soft_clipped.max()
    single_pass = abs(max_single_clipped) <  set
    both_pass =  abs(total_clipped) < bet
    if single_pass and both_pass and not aln.is_unmapped:
        outbam.write(aln)
        output_count += 1
    return output_count

cpdef int filter_bam(str in_bam, str out_bam, float single_end_thresh, float both_end_thresh):
    '''
    This function filter softclipped sequence by tresholds regarding to the sequence length
    '''
    cdef:
        AlignmentFile inbam, outbam
        AlignedSegment aln
        ndarray cigar_array, all_soft_clipped
        int output_count

    with pysam.AlignmentFile(in_bam,'rb') as inbam:
        with pysam.Samfile(out_bam,'wb',template = inbam) as outbam:
            for count, aln in enumerate(inbam):
                if qualify_aln(aln):
                    if 'S' not in aln.cigarstring:
                        outbam.write(aln)
                        output_count += 1
                    else:
                        check_aln(output_count, outbam, aln, single_end_thresh, both_end_thresh)
                if count % 1000000 == 0 and count != 0:
                    print 'Parsed %i alignments' %(count)
    return output_count
