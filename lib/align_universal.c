#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compare.h"

int read_queries(char *input_query, unsigned int n_queries, unsigned int len_query, unsigned int ***ff, unsigned int **oligo_lens_ff, unsigned int *len_ff)
{
    FILE * p_ff;
    p_ff = fopen(input_query,"rb");
    
    *oligo_lens_ff = (unsigned int *) malloc (n_queries * sizeof(unsigned int));
    *ff = (unsigned int **) malloc (n_queries * sizeof(unsigned int *));
    
    unsigned int i, oligo_buf[2], oligo_flag, ind_buf;
    
    for (i = 0; i < n_queries; i++)
    {
        (*ff)[i] = (unsigned int *) malloc (len_query * sizeof(unsigned int));
    }
    
    ind_buf = 0;
    oligo_flag = fread(oligo_buf,4,1,p_ff);
    
    for (i=0; oligo_flag>0; i++)
    {
        
        if (ind_buf==0)
        {
            for (ind_buf=0; ind_buf<len_query; ind_buf++)
            {
                (*ff)[i][ind_buf] = 0x00000000;
            }
            ind_buf = 0;
        }
        while (oligo_flag>0)
        {
            if (oligo_buf[0]==0x20000000)
            {
                oligo_flag = fread(oligo_buf, 4, 1, p_ff);
                (*oligo_lens_ff)[i] = ind_buf;
                ind_buf = 0;
                break;
            }
            else
            {
                (*ff)[i][ind_buf] = oligo_buf[0];
                oligo_flag = fread(oligo_buf, 4, 1, p_ff);
                ind_buf++;
            }
        }
    }
    *len_ff = i;
    return 0;
}

int main (int argc, char **argv)
{

    if (argc!=12) {printf("Red-C Rabin-Karp implementation for adapters and . Usage:\n%s bridge_library.bin fasta.bin (mode:2,1,0) (fasta_seq_len) (seq_len_norm) (n_queries) (shift_r) (shift_l) (allowed_mis) (add_mis) (report_len) > output.txt\n",argv[0]);return 0;}
    
    printf("This is output of program '%s' processing file '%s'\n",argv[0],argv[1]);
    
    // Reading input parameters
    unsigned int fasta_seq_len, seq_len_norm, n_queries, allowed_mis, add_mis;
		int shift_max_right, shift_max_left, report_len;
    int mode;

    sscanf(argv[3],"%d",&mode);
    sscanf(argv[4],"%u",&fasta_seq_len  );
    sscanf(argv[5],"%u",&seq_len_norm   );
    sscanf(argv[6],"%u",&n_queries      );
    sscanf(argv[7],"%i",&shift_max_left );
    sscanf(argv[8],"%i",&shift_max_right);
    sscanf(argv[9],"%u",&allowed_mis    );
    sscanf(argv[10],"%u",&add_mis       );
    sscanf(argv[11],"%u",&report_len    );
    
    char * input_query  = argv[1],
    * input_library = argv[2];

    // Opening and reading sequences from ff
    
    unsigned int **ff, *oligo_lens_ff, len_ff, len_query=seq_len_norm;
    
    read_queries(input_query, n_queries+1, len_query, &ff, &oligo_lens_ff, &len_ff);
    
    // Opening and reading file with total reads
    FILE * bin;
    bin = fopen(input_library, "rb");
    
    unsigned int ccc[len_query+1], * array_shift_l[10], * array_shift_r[10];
    unsigned int c;
    int res;
    unsigned int i, j;
    
    unsigned int oligo_len;
    int shift;
    unsigned int shift_l, shift_s, comparison_length;
    
    int found = 0;
    
    for (i=0; fread(ccc, 4, len_query-1, bin)>0; i++)
    {
        
        found = 0;
        
        for (j=0; j<len_ff-1; j++)
        {
            oligo_len = oligo_lens_ff[j];
            
            for (shift_l=1; shift_l<=4; shift_l++)
            {
                shift_right(shift_l, ccc, &(array_shift_r[shift_l]), len_query-2);
                shift_left( shift_l, ccc, &(array_shift_l[shift_l]), len_query-2);
            }
            
            for (shift=shift_max_right; shift>=shift_max_left; shift--)
            {
                if (abs(shift)%8 == 0)
                {
                    
                    shift_l = abs(shift)/8;
                    
                    if (shift >= 0)
                    {
                        comparison_length = ( (fasta_seq_len-shift)>8*oligo_len ) ? oligo_len : (fasta_seq_len-shift)/8+1;
                        res = find_(ccc,   ff[j], comparison_length, &c, shift_l, 0, 2, allowed_mis);
                        if (res) printf("\t%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=fasta_seq_len ? shift+report_len : fasta_seq_len );
                    }
                    else
                    {
                        res = find_(ccc,   ff[j], oligo_len, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
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
                        comparison_length = ( (fasta_seq_len-shift)>8*oligo_len ) ? oligo_len : (fasta_seq_len-shift)/8;
                        res = find_(array_shift_l[shift_s], ff[j], comparison_length+1, &c, shift_l+1, 0, 2, allowed_mis);
                        if (res) printf("\t%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=fasta_seq_len ? shift+report_len : fasta_seq_len );
                    }
                    else
                    {
                        res = find_(array_shift_r[shift_s], ff[j], oligo_len+2, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
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
                        res = find_(array_shift_r[shift_s], ff[j], comparison_length, &c, shift_l+1, 0, 2, allowed_mis);
                        if (res) printf("\t%u\t%u\t%u\t%i\t%i\n", i, j, c, shift, shift+report_len<=fasta_seq_len ? shift+report_len : fasta_seq_len );
                    }
                    else
                    {
                        res = find_(array_shift_l[shift_s], ff[j], oligo_len, &c, oligo_len-shift_l-1, oligo_len-1, 1, allowed_mis);
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
            if (found&mode) break;
        }
       
        if ((!found)&mode) printf("\t%u\t%u\t%u\t%i\t%i\n", i, j, 500, 500, 500);
    }
    
    return 0;
    
}
