#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compare.h"

int main (int argc, char **argv)
{
    
    if (argc!=12) {printf("Usage:\n%s bridge_library.bin fasta.bin (mode:1,0) (fasta_seq_len) (fasta_seq_len_norm) (n_queries) (shift_r) (shift_l) (allowed_mis) (add_mis) (report_len) > output.txt\n",argv[0]);return 0;}
    
    printf("This is output of program '%s' processing file '%s'\n",argv[0],argv[1]);
    
    // Reading input parameters
    unsigned int fasta_seq_len, fasta_seq_len_norm, query_seq_len_norm, len_queries, allowed_mis, add_mis;
    int shift_max_right, shift_max_left, report_len;
    int mode;

    sscanf(argv[3],"%d",&mode);
    sscanf(argv[4],"%u",&fasta_seq_len  );
    sscanf(argv[5],"%u",&fasta_seq_len_norm   );
    sscanf(argv[6],"%u", &query_seq_len_norm);
    sscanf(argv[7],"%i",&shift_max_left );
    sscanf(argv[8],"%i",&shift_max_right);
    sscanf(argv[9],"%u",&allowed_mis    );
    sscanf(argv[10],"%u",&add_mis       );
    sscanf(argv[11],"%u",&report_len    );
    
    char * input_query  = argv[1],
    * input_library = argv[2];

    // Opening and reading file with total reads and with query library (simultaneously, NB)
    FILE * bin_lib, *bin_query;
    bin_lib = fopen(input_library, "rb");
    bin_query = fopen(input_query, "rb");
    
    unsigned int ccc[fasta_seq_len_norm+1], ccc_query[fasta_seq_len_norm+1], * array_shift_l[10], * array_shift_r[10];
    unsigned int c;
    int res;
    unsigned int i, j;
    
    unsigned int oligo_len;
    int shift;
    unsigned int shift_l, shift_s, comparison_length;
    
    int found = 0;
    
    for (i=0; fread(ccc, 4, fasta_seq_len_norm-1, bin_lib)>0; i++)
    {
        
        fread(ccc_query, 4, query_seq_len_norm, bin_query); //
        
        found = 0;

        
            oligo_len = query_seq_len_norm; // in bins
            
            for (shift_l=1; shift_l<=4; shift_l++)
            {
                shift_right(shift_l, ccc, &(array_shift_r[shift_l]), fasta_seq_len_norm-2);
                shift_left( shift_l, ccc, &(array_shift_l[shift_l]), fasta_seq_len_norm-2);
            }
            
            for (shift=shift_max_left; shift<=shift_max_right; shift++)
            {
                if (abs(shift)%8 == 0)
                {
                    
                    shift_l = abs(shift)/8;
                    
                    if (shift >= 0)
                    {
                        comparison_length = ( (fasta_seq_len-shift)>8*oligo_len ) ? oligo_len-1 : (fasta_seq_len-shift)/8+1;
                        res = find_(ccc,   ccc_query, comparison_length, &c, shift_l, 0, 2, allowed_mis);
                        if (res) printf("\t%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=fasta_seq_len ? shift+report_len : fasta_seq_len );
                    }
                    else
                    {
                        res = find_(ccc,   ccc_query, oligo_len, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
                        if (res) printf("\t%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=fasta_seq_len ? shift+report_len : fasta_seq_len );
                    }
                    
                    if (res) found = 1;
                }
                else if (abs(shift)%8 < 4)
                {
                    shift_l = abs(shift)/8;
                    shift_s = abs(shift)%8;
                    
                    if (shift >= 0)
                    {
                        comparison_length = ( (fasta_seq_len-shift)>8*oligo_len ) ? oligo_len-1 : (fasta_seq_len-shift)/8;
                        res = find_(array_shift_l[shift_s], ccc_query, comparison_length+1, &c, shift_l+1, 0, 2, allowed_mis);
                        if (res) printf("\t%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=fasta_seq_len ? shift+report_len : fasta_seq_len );
                    }
                    else
                    {
                        res = find_(array_shift_r[shift_s], ccc_query, oligo_len+2, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
                        if (res) printf("\t%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=fasta_seq_len ? shift+report_len : fasta_seq_len );
                    }
                    
                    if (res) found = 1;
                }
                else
                {
                    shift_l = (abs(shift)/8);
                    shift_s = 8-abs(shift)%8;
                    
                    if (shift >= 0)
                    {
                        comparison_length = ( (fasta_seq_len-shift)>8*oligo_len ) ? oligo_len : (fasta_seq_len-shift)/8+1;
                        res = find_(array_shift_r[shift_s], ccc_query, comparison_length, &c, shift_l+1, 0, 2, allowed_mis);
                        if (res) printf("\t%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=fasta_seq_len ? shift+report_len : fasta_seq_len );
                    }
                    else
                    {
                        res = find_(array_shift_l[shift_s], ccc_query, oligo_len, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
                        if (res) printf("\t%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=fasta_seq_len ? shift+report_len : fasta_seq_len );
                    }
                    
                    if (res) found = 1;
                }
                
                if (found&mode) break;
            }
            for (shift_l=1; shift_l<=4; shift_l++)
            {
                free(array_shift_r[shift_l]);
                free(array_shift_l[shift_l]);
            }
        
    
    
        if ((!found)&mode) printf("\t%u\t%u\t%u\t%i\t%i\n", i, j, 500, 500, 500);
    }
    
    return 0;
    
}
