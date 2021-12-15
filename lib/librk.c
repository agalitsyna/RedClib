// Library for Rabin-Karp C implementation of string search, optimized for nucleotide sequences
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Defining binary constants for shifting and mismatch counting
#define B1 0x11111111
#define U  0xffffffff
#define M1 0x55555555 // binary: 0101 0101 0101 ... (unity each second bit)
#define M2 0x33333333 // binary: 0011 0011 0011 ... (two zeros, two ones)
#define M4 0x0f0f0f0f // binary: 0000 1111 0000 ... (4 zeros, 4 ones)
#define H01 0x01010101 // binary: 0000 0001 0000 0001 ... 

unsigned int find_mismatches( unsigned int * seq1, unsigned int * seq2, unsigned int n , unsigned int * mismatches)
{
  unsigned int c,i,j;

  *mismatches = 0;
  for (i=0; i<n; i++)
  {
    c = seq1[i]^seq2[i];
    c = (c&B1)|((c>>1)&B1)|((c>>2)&B1);
    c -= (c >> 1) & M1;             //put count of each 2 bits into those 2 bits
    c = (c & M2) + ((c >> 2) & M2); //put count of each 4 bits into those 4 bits 
    c = (c + (c >> 4)) & M4;        //put count of each 8 bits into those 8 bits 
    c = (c * H01)>>24;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
    *mismatches = *mismatches + c;
  }
  return 1;
}

unsigned int find_substitutions_in_region( unsigned int * seq1, unsigned int * seq2, unsigned int n, unsigned int * v, unsigned int seq1_from, unsigned int seq2_from)
{
  unsigned int c,i,j;
  *v = 0;
  for (i=seq1_from, j=seq2_from; i<seq1_from+n; i++, j++)
  {
      c = seq1[i]^seq2[j];
      c = (c&B1)|((c>>1)&B1)|((c>>2)&B1);
      c -= (c >> 1) & M1;             //put count of each 2 bits into those 2 bits
      c = (c & M2) + ((c >> 2) & M2); //put count of each 4 bits into those 4 bits 
      c = (c + (c >> 4)) & M4;        //put count of each 8 bits into those 8 bits 
      c = (c * H01)>>24;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
      *v = *v + c;
  }
  return 1;
}

unsigned int find_substitutions_in_region_inversed( unsigned int * seq1, unsigned int * seq2, unsigned int n, unsigned int * v, unsigned int seq1_to, unsigned int seq2_to)
{
  unsigned int c,i,j;
  for (i=seq1_to, j=seq2_to; i>seq1_to-n; --i, --j) {printf("1: %x 2: %x\n",seq1[i],seq2[j]);}
  *v = 0;
  for (i=seq1_to, j=seq2_to; i>seq1_to-n; --i, --j)
  {
      c = seq1[i]^seq2[j];
      c = (c&B1)|((c>>1)&B1)|((c>>2)&B1);
      c -= (c >> 1) & M1;             //put count of each 2 bits into those 2 bits
      c = (c & M2) + ((c >> 2) & M2); //put count of each 4 bits into those 4 bits 
      c = (c + (c >> 4)) & M4;        //put count of each 8 bits into those 8 bits 
      c = (c * H01)>>24;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
      *v = *v + c;
  }
  return 1;
}

unsigned int find_mask( unsigned int * seq1, unsigned int * seq2, unsigned int n, unsigned int * v, unsigned int seq1_to, unsigned int seq2_to, unsigned int mode, unsigned int max_subst,unsigned int mask_start, unsigned int mask_end)
{
/*
seq1 - 1st sequence to compare
seq2 - 2nd
n - number of bits to compare
v - pointer of unsigned int to save count
seq1_to and seq2_to - start position for comparison
mode - 1: from end to start; 2: from start to end
max_subst - when to break search
mask_start - how many chars should not be compared from the beginning of sequences
mask_end - how many chars should not be compared from the end of sequences
*/
  unsigned int c,i,j;
  *v = 0;
  if (mode == 1)
  {
    for (i=seq1_to, j=seq2_to; i>seq1_to-n; --i, --j)
    {
        c = seq1[i]^seq2[j];
        c = (c&B1)|((c>>1)&B1)|((c>>2)&B1);
        c -= (c >> 1) & M1;             //put count of each 2 bits into those 2 bits
        c = (c & M2) + ((c >> 2) & M2); //put count of each 4 bits into those 4 bits 
        c = (c + (c >> 4)) & M4;        //put count of each 8 bits into those 8 bits 
        c = (c * H01)>>24;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
        *v = *v + c;
        if (*v > max_subst) {//printf("%i, i:%i, j:%i\n",*v,i,j);
          *v = 10000; return 0;}
    }
  }
  else if (mode == 2)
  {
    for (i=seq1_to, j=seq2_to; i<seq1_to+n; i++, j++)
    {
        c = seq1[i]^seq2[j];
        c = (c&B1)|((c>>1)&B1)|((c>>2)&B1);
        c -= (c >> 1) & M1;             //put count of each 2 bits into those 2 bits
        c = (c & M2) + ((c >> 2) & M2); //put count of each 4 bits into those 4 bits 
        c = (c + (c >> 4)) & M4;        //put count of each 8 bits into those 8 bits 
        c = (c * H01)>>24;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
        *v = *v + c;
        if (*v > max_subst) {//printf("%i, i:%i, j:%i\n",*v,i,j);
          *v = 10000; return 0;}
    }
  }
  return 1;
}

unsigned int align_seq( unsigned int * seq1, unsigned int * seq2_for, unsigned int * seq2_rev,
                        unsigned int n, unsigned int seq1_to, unsigned int seq2_to,
                        unsigned int mode, unsigned int max_subst,
                        unsigned int * start, unsigned int * end)
{
/*
seq1 - 1st sequence to compare
seq2 - 2nd
n - number of bits to compare
v - pointer of unsigned int to save count
seq1_to and seq2_to - start position for comparison
mode - 1: from end to start; 2: from start to end
max_subst - when to break search
mask_start - how many chars should not be compared from the beginning of sequences
mask_end - how many chars should not be compared from the end of sequences
*/

  return 1;
}


const char *byte_to_binary(int x)
{
    static char b[40];
    b[0] = '\0';

    unsigned int z;
    for (z = 0x80000000; z > 0; z >>= 1)
    {
        strcat(b, ((x & z) == z) ? "1" : "0");
    }

    return b;
}

unsigned int interpret_buf(int buf)
{
    unsigned int i;
    char c;

    char s[10];
    s[0] = '\0';

    int l[8] = {0x00000007,0x00000070,0x00000700,0x00007000,0x00070000,0x00700000,0x07000000,0x70000000};
    int hex_A = 0x00000001, hex_T = 0x00000004, hex_C = 0x00000003, hex_G = 0x00000007;


    for (i=0; i<8; i++)
    {
        c = (buf & l[(7-i)]) >> 4*(7-i); // Generating the right order
        strcat(s, (c==hex_A) ? "A" : (c==hex_T) ? "T" : (c==hex_G) ? "G" : (c==hex_C) ? "C" : "N");
    }
    printf("Start interpreting: %s \t -> \t %s\n", byte_to_binary(buf), s);

    return 0;
}


unsigned int find_(unsigned int * seq1,
                   unsigned int * seq2,
                   unsigned int n,
                   unsigned int * v,
                   unsigned int seq1_to,
                   unsigned int seq2_to,
                   unsigned int mode,
                   unsigned int max_subst)
{
/*
seq1 - 1st sequence to compare
seq2 - 2nd
n - number of bits to compare
v - pointer of unsigned int to save count
seq1_to and seq2_to - start position for comparison
mode - 1: from end to start; 2: from start to end
max_subst - when to break search
*/
  unsigned int c,  u1, u2;
  int i,j, limit;
  *v = 0;

  if (mode == 1)
  {
    limit = ((seq1_to>=n) ? (seq1_to-n):0);
    for (i=seq1_to, j=seq2_to; i>=limit; --i, --j)
    {
        c = seq1[i]^seq2[j];

        u1 = U^seq1[i];
        u2 = U^seq2[j];
        u1 = (u1&(u1>>1)&(u1>>2)&B1);
        u2 = (u2&(u2>>1)&(u2>>2)&B1);

        c = ((c&B1)|((c>>1)&B1)|((c>>2)&B1))&(~(u1|u2));
        c -= (c >> 1) & M1;             //put count of each 2 bits into those 2 bits
        c = (c & M2) + ((c >> 2) & M2); //put count of each 4 bits into those 4 bits
        c = (c + (c >> 4)) & M4;        //put count of each 8 bits into those 8 bits
        c = (c * H01)>>24;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
        *v = *v + c;

        if (*v > max_subst) {//printf("%i, i:%i, j:%i\n",*v,i,j);
          return 0;}

    }
  }
  else if (mode == 2) // Add negative seq1_to and seq2_to? No, just add iteration over 2nd shift
  {
    if (n==0)
    {
        n=1;
    }
    for (i=seq1_to, j=seq2_to; i<seq1_to+n; i++, j++)
    {
        c = seq1[i]^seq2[j];

        u1 = U^seq1[i];
        u2 = U^seq2[j];
        u1 = (u1&(u1>>1)&(u1>>2)&B1);
        u2 = (u2&(u2>>1)&(u2>>2)&B1);

        c = ((c&B1)|((c>>1)&B1)|((c>>2)&B1))&(~(u1|u2));
        c -= (c >> 1) & M1;             //put count of each 2 bits into those 2 bits
        c = (c & M2) + ((c >> 2) & M2); //put count of each 4 bits into those 4 bits
        c = (c + (c >> 4)) & M4;        //put count of each 8 bits into those 8 bits
        c = (c * H01)>>24;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
        *v = *v + c;

        if (*v > max_subst) {return 0;}
    }
  }
  return 1;
}

unsigned int shift_right(unsigned int x,
                         unsigned int * seq,
                         unsigned int ** seq_mod,
                         unsigned int n )
{

  unsigned int i;
  unsigned int b[8] = {0x0000000f,0x000000ff,0x00000fff,0x0000ffff,0x000fffff,0x00ffffff,0x0fffffff,0xffffffff};
  (*seq_mod) = (unsigned int *)malloc((n+2)*sizeof(int));

  (*seq_mod)[0] = ((~b[x-1])&seq[0])>>(x*4);
  for (i=1;i<n;i++)
  {
    (*seq_mod)[i] = (((~b[x-1])&seq[i])>>(x*4))|((seq[i-1]&b[x-1])<<((8-x)*4));
  }
  (*seq_mod)[n] = ((seq[n-1]&b[x-1])<<((8-x)*4));
  return 1;
}

unsigned int shift_left(unsigned int x,
                        unsigned int * seq,
                        unsigned int ** seq_mod,
                        unsigned int n ) // modifying only last/fist byte
{
  unsigned int i;
  unsigned int b[8] = {0xf0000000,0xff000000,0xfff00000,0xffff0000,0xfffff000,0xffffff00,0xfffffff0,0xffffffff};
  (*seq_mod) = (unsigned int *)malloc((n+3)*sizeof(int));

  (*seq_mod)[0] = ((b[x-1])&seq[0])>>((8-x)*4);
  for (i=0;i<n-1;i++)
  {
    (*seq_mod)[i+1] = (((~b[x-1])&seq[i])<<(x*4))|((seq[i+1]&b[x-1])>>((8-x)*4));
  }

    (*seq_mod)[n] = ((~b[x-1])&seq[n-1])<<(x*4);
  return 1;
}
