#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "librk.h"

int main (int argc, char **argv)
{
    
    if ((argc!=8)&(argc!=9)) {
        fprintf(stderr, "\nRabin-Karp implementation for nucleotides query search with substitutions. Requires pre-computed hashes as input. Outputs a tsv file with hits of query in the library of entries.\n\n");
        fprintf(stderr, "Usage:\t%s [query.bin] [library.bin] [entry_seq_len] [query_seq_len] [min_shift_left] [min_shift_right] [allowed_mismatches] (right_to_left, optional) (report_single_hit, optional) > output.tsv\n\n", argv[0]);
        fprintf(stderr, "\t query.bin : string filename \n\t\t input binary hash for query/queries, compiled by fasta2hash/fastq2hash from librk. All queries should be the same length, corresponds to encoded_query_seq_len!\n");
        fprintf(stderr, "\t library.bin : string filename \n\t\t input binary hash for library of sequence entries, compiled by fasta2hash/fastq2hash from librk. All entries should be the same length, same as entry_seq_len!\n");
        fprintf(stderr, "\t query_seq_len : int \n\t\t the length of the queries.\n");
        fprintf(stderr, "\t entry_seq_len : int \n\t\t read sequence length of target library.\n");
        fprintf(stderr, "\t min_shift_left : int \n\t\t minimal allowed distance from the left side of the entry where the search will be started.\n");
        fprintf(stderr, "\t min_shift_right : int \n\t\t minimal allowed distance from the right side of the entry after which the search will be finished.\n");
        fprintf(stderr, "\t allowed_mismatches : int \n\t\t allowed mismatches per hit.\n");
        fprintf(stderr, "\t right_to_left : int (optional) \n\t\t if 1 (true), then search from right to left (from the end of read to its start), 0 is default.\n");
        fprintf(stderr, "\t report_single_hit : int (optional) \n\t\t if 1 (true), then a single hit per entry will be reported for each sequence with dump output if the hit was not is_found (reports 99999 as the position and number of mismatches).\n\t\t if 0 (false), then report all but valid hits\n\n");
        return 0;
        }

    printf("#entry_index\tquery_index\tmismatches\tstart_hit\tend_hit\n");

    // Reading input parameters
    unsigned int entry_seq_len, query_seq_len, allowed_mis, add_mis, report_single_hit, right_to_left;
    int shift_max_right, shift_max_left;

    sscanf(argv[3],"%u",&query_seq_len);
    sscanf(argv[4],"%u",&entry_seq_len);
    sscanf(argv[5],"%i",&shift_max_left);
    sscanf(argv[6],"%i",&shift_max_right);
    sscanf(argv[7],"%u",&allowed_mis);
    if (argc>=9) {
        sscanf(argv[8],"%u",&right_to_left);
    } else {
        right_to_left = 0;
    }
    if (argc==10) {
        sscanf(argv[9],"%u",&report_single_hit);
    } else {
        report_single_hit = 1;
    }

    char * input_query  = argv[1], * input_library = argv[2];
//
    // Open and read file with total reads and with query library (simultaneously, NB)
    unsigned int query_length = 3+(query_seq_len-1)/8; // ((rna_complementary_length+3).intdiv(8))+1
    unsigned int entry_length = 3+(entry_seq_len-1)/8;
    FILE *library_bin, *query_bin;
    library_bin = fopen(input_library, "rb");
    query_bin = fopen(input_query, "rb");

    unsigned int buf[entry_length+1], buf_query[query_length+1], * array_shift_l[10], * array_shift_r[10];
    unsigned int c, i, j, oligo_len;
    int is_hit=0, shift, is_found = 0;
    unsigned int shift_l, shift_s, target_region_len;

    for (i=0; fread(buf, 4, entry_length-1, library_bin)>0; i++)
    {
        fread(buf_query, 4, query_length-1, query_bin);
        is_found = 0;

        oligo_len = query_length-2; // in bins
        for (shift_l=1; shift_l<=4; shift_l++)
        {
            shift_right(shift_l, buf, &(array_shift_r[shift_l]), entry_length-2);
            shift_left( shift_l, buf, &(array_shift_l[shift_l]), entry_length-2);
        }

//        if (right_to_left>0)
//        {
//            for (shift=shift_max_right; shift<=shift_max_left; shift--)
//            {
//                if (abs(shift)%8 == 0)
//                {
//                    shift_l = abs(shift)/8;
//                    if (shift >= 0)
//                    {
//                        target_region_len = ( (entry_seq_len-shift)>8*oligo_len ) ? oligo_len-1 : (entry_seq_len-shift)/8+1;
//                        is_hit = find_(buf,   buf_query, target_region_len, &c, shift_l, 0, 2, allowed_mis);
//                        if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
//                    }
//                    else
//                    {
//                        is_hit = find_(buf,   buf_query, oligo_len, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
//                        if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
//                    }
//
//                    if (is_hit) is_found = 1;
//                }
//                else if (abs(shift)%8 < 4)
//                {
//                    shift_l = abs(shift)/8;
//                    shift_s = abs(shift)%8;
//                    if (shift >= 0)
//                    {
//                        target_region_len = ( (entry_seq_len-shift)>8*oligo_len ) ? oligo_len-1 : (entry_seq_len-shift)/8;
//                        is_hit = find_(array_shift_l[shift_s], buf_query, target_region_len+1, &c, shift_l+1, 0, 2, allowed_mis);
//                        if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
//                    }
//                    else
//                    {
//                        is_hit = find_(array_shift_r[shift_s], buf_query, oligo_len+2, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
//                        if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
//                    }
//
//                    if (is_hit) is_found = 1;
//                }
//                else
//                {
//                    shift_l = (abs(shift)/8);
//                    shift_s = 8-abs(shift)%8;
//
//                    if (shift >= 0)
//                    {
//                        target_region_len = ( (entry_seq_len-shift)>8*oligo_len ) ? oligo_len : (entry_seq_len-shift)/8+1;
//                        is_hit = find_(array_shift_r[shift_s], buf_query, target_region_len, &c, shift_l+1, 0, 2, allowed_mis);
//                        if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
//                    }
//                    else
//                    {
//                        is_hit = find_(array_shift_l[shift_s], buf_query, oligo_len, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
//                        if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
//                    }
//
//                    if (is_hit) is_found = 1;
//                }
//                if (is_found&report_single_hit) break;
//            }
//        }
//        else // left to right search
//        {
                for (shift=shift_max_left; shift<=shift_max_right; shift++)
                {
                    if (abs(shift)%8 == 0)
                    {
                        shift_l = abs(shift)/8;
                        if (shift >= 0)
                        {
                            target_region_len = ( (entry_seq_len-shift)>8*oligo_len ) ? oligo_len : (entry_seq_len-shift)/8+1;
                            is_hit = find_(buf,   buf_query, target_region_len, &c, shift_l, 0, 2, allowed_mis);
                            if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
                        }
                        else
                        {
                            is_hit = find_(buf,   buf_query, oligo_len, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
                            if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
                        }
                        if (is_hit) is_found = 1;
                    }
                    else if (abs(shift)%8 < 4)
                    {
                        shift_l = abs(shift)/8;
                        shift_s = abs(shift)%8;
                        if (shift >= 0)
                        {
                            target_region_len = ( (entry_seq_len-shift)>8*oligo_len ) ? oligo_len : (entry_seq_len-shift)/8+1;
                            is_hit = find_(array_shift_l[shift_s], buf_query, target_region_len, &c, shift_l+1, 0, 2, allowed_mis);
                            if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
                        }
                        else
                        {
                            is_hit = find_(array_shift_r[shift_s], buf_query, oligo_len+2, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
                            if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
                        }
                        if (is_hit) is_found = 1;
                    }
                    else
                    {
                        shift_l = (abs(shift)/8);
                        shift_s = 8-abs(shift)%8;
                        if (shift >= 0)
                        {
                            target_region_len = ( (entry_seq_len-shift)>8*oligo_len ) ? oligo_len : (entry_seq_len-shift)/8+1;
                            is_hit = find_(array_shift_r[shift_s], buf_query, target_region_len, &c, shift_l+1, 0, 2, allowed_mis);
                            if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
                        }
                        else
                        {
                            is_hit = find_(array_shift_l[shift_s], buf_query, oligo_len, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
                            if (is_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+query_seq_len<=entry_seq_len ? shift+query_seq_len : entry_seq_len );
                        }
                        if (is_hit) is_found = 1;
                    }
                    if (is_found&report_single_hit) break;
                }
//        }
        // Free the pointer of array buffer:
        for (shift_l=1; shift_l<=4; shift_l++)
        {
            free(array_shift_r[shift_l]);
            free(array_shift_l[shift_l]);
        }
        if ((!is_found)&report_single_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, 99999, 99999, 99999);
    }
    fclose( library_bin );
    fclose( query_bin );
    return 0;
}
