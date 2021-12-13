#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "librk.h"

int main (int argc, char **argv)
{
    
    if ((argc!=10)&(argc!=11)) {
        fprintf(stderr, "\nRabin-Karp implementation for nucleotides query search with substitutions. Requires pre-computed hashes as input. Outputs a tsv file with hits of query in the library of entries.\n\n");
        fprintf(stderr, "Usage:\t%s query.bin library.bin encoded_query_seq_len library_seq_len encoded_library_seq_len min_shift_left min_shift_right allowed_mismatches reported_len (report_single_hit) > output.tsv\n\n", argv[0]);
        fprintf(stderr, "\t query.bin : string filename \n\t\t input binary hash for query/queries, compiled by fasta2hash/fastq2hash from librk. All queries should be the same length, corresponds to encoded_query_seq_len!\n");
        fprintf(stderr, "\t library.bin : string filename \n\t\t input binary hash for library of sequence entries, compiled by fasta2hash/fastq2hash from librk. All entries should be the same length, same as library_seq_len!\n");
        fprintf(stderr, "\t encoded_query_seq_len : int \n\t\t sequence length of queries converted to length of binary hash.\n");
        fprintf(stderr, "\t library_seq_len : int \n\t\t read sequence length of target library.\n");
        fprintf(stderr, "\t encoded_library_seq_len : int \n\t\t sequence length of target library converted to length of binary hash. The numbers for common Illumina lengths are: [80:12, 101:15, 125:18, 151:21, 133:19, 251:34]\n");
        fprintf(stderr, "\t min_shift_left : int \n\t\t minimal allowed distance from the left side of the entry where the search will be started.\n");
        fprintf(stderr, "\t min_shift_right : int \n\t\t minimal allowed distance from the right side of the entry after which the search will be finished.\n");
        fprintf(stderr, "\t allowed_mismatches : int \n\t\t allowed mismatches per hit.\n");
        fprintf(stderr, "\t reported_len : int \n\t\t the length of the hit for reporting (might be larger).\n");
        fprintf(stderr, "\t report_single_hit : int (optional) \n\t\t if 1 (true), then a single hit per entry will be reported for each sequence with dump output if the hit was not found (reports 99999 as the position and number of mismatches).\n\t\t if 0 (false), then report all but valid hits\n\n");
        return 0;
        }

    printf("#entry_index\tquery_index\tmismatches\tstart_hit\tend_hit\n");

    // Reading input parameters
    unsigned int library_seq_len, library_seq_len_norm, query_seq_len_norm, len_queries, allowed_mis, add_mis;
    int shift_max_right, shift_max_left, report_len;
    int report_single_hit;

    sscanf(argv[3],"%u",&query_seq_len_norm);
    sscanf(argv[4],"%u",&library_seq_len  );
    sscanf(argv[5],"%u",&library_seq_len_norm   );
    sscanf(argv[6],"%i",&shift_max_left );
    sscanf(argv[7],"%i",&shift_max_right);
    sscanf(argv[8],"%u",&allowed_mis    );
    sscanf(argv[9],"%u",&report_len    );
    if (argc!=11) {
        sscanf(argv[10],"%d",&report_single_hit);
    } else {
        report_single_hit = 1;
    }

    char * input_query  = argv[1],
    * input_library = argv[2];

    // Opening and reading file with total reads and with query library (simultaneously, NB)
    FILE * bin_lib, *bin_query;
    bin_lib = fopen(input_library, "rb");
    bin_query = fopen(input_query, "rb");
    
    unsigned int ccc[library_seq_len_norm+1], ccc_query[library_seq_len_norm+1], * array_shift_l[10], * array_shift_r[10];
    unsigned int c;
    int res;
    unsigned int i, j;
    
    unsigned int oligo_len;
    int shift;
    unsigned int shift_l, shift_s, comparison_length;
    
    int found = 0;
    
    for (i=0; fread(ccc, 4, library_seq_len_norm-1, bin_lib)>0; i++)
    {
        
        fread(ccc_query, 4, query_seq_len_norm, bin_query); //
        
        found = 0;

        
            oligo_len = query_seq_len_norm; // in bins
            
            for (shift_l=1; shift_l<=4; shift_l++)
            {
                shift_right(shift_l, ccc, &(array_shift_r[shift_l]), library_seq_len_norm-2);
                shift_left( shift_l, ccc, &(array_shift_l[shift_l]), library_seq_len_norm-2);
            }
            
            for (shift=shift_max_left; shift<=shift_max_right; shift++)
            {
                if (abs(shift)%8 == 0)
                {
                    
                    shift_l = abs(shift)/8;
                    
                    if (shift >= 0)
                    {
                        comparison_length = ( (library_seq_len-shift)>8*oligo_len ) ? oligo_len-1 : (library_seq_len-shift)/8+1;
                        res = find_(ccc,   ccc_query, comparison_length, &c, shift_l, 0, 2, allowed_mis);
                        if (res) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=library_seq_len ? shift+report_len : library_seq_len );
                    }
                    else
                    {
                        res = find_(ccc,   ccc_query, oligo_len, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
                        if (res) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=library_seq_len ? shift+report_len : library_seq_len );
                    }
                    
                    if (res) found = 1;
                }
                else if (abs(shift)%8 < 4)
                {
                    shift_l = abs(shift)/8;
                    shift_s = abs(shift)%8;
                    
                    if (shift >= 0)
                    {
                        comparison_length = ( (library_seq_len-shift)>8*oligo_len ) ? oligo_len-1 : (library_seq_len-shift)/8;
                        res = find_(array_shift_l[shift_s], ccc_query, comparison_length+1, &c, shift_l+1, 0, 2, allowed_mis);
                        if (res) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=library_seq_len ? shift+report_len : library_seq_len );
                    }
                    else
                    {
                        res = find_(array_shift_r[shift_s], ccc_query, oligo_len+2, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
                        if (res) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=library_seq_len ? shift+report_len : library_seq_len );
                    }
                    
                    if (res) found = 1;
                }
                else
                {
                    shift_l = (abs(shift)/8);
                    shift_s = 8-abs(shift)%8;
                    
                    if (shift >= 0)
                    {
                        comparison_length = ( (library_seq_len-shift)>8*oligo_len ) ? oligo_len : (library_seq_len-shift)/8+1;
                        res = find_(array_shift_r[shift_s], ccc_query, comparison_length, &c, shift_l+1, 0, 2, allowed_mis);
                        if (res) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=library_seq_len ? shift+report_len : library_seq_len );
                    }
                    else
                    {
                        res = find_(array_shift_l[shift_s], ccc_query, oligo_len, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
                        if (res) printf("%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=library_seq_len ? shift+report_len : library_seq_len );
                    }
                    
                    if (res) found = 1;
                }
                
                if (found&report_single_hit) break;
            }
            for (shift_l=1; shift_l<=4; shift_l++)
            {
                free(array_shift_r[shift_l]);
                free(array_shift_l[shift_l]);
            }
        
    
    
        if ((!found)&report_single_hit) printf("%u\t%u\t%u\t%i\t%i\n", i, j, 99999, 99999, 99999);
    }
    
    return 0;
    
}
