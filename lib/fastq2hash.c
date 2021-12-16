#include <stdio.h>
#include "librk.h"

#define B1 0x11111111
#define U6 0x66666666

unsigned int main (int argc, char **argv)
{
  // Arguments parser
  if (argc!=3) {printf("Convert FASTQ to binary hash  (each character A/T/G/C/N is encoded as 4 bits). Case insensitive. \nUsage:\n\n%s filename.fasta filename.bin\n",argv[0]);return 0;}

  FILE * fastq;
  fastq = fopen(argv[1],"r");

  FILE * bin;
  bin = fopen(argv[2],"wb");

  int chunk_size=5000; // Buffer size
  int name_length=5000; // Maximal length of read name
  char s[chunk_size+1];  // Buffer string for the current fastq line
  int i, c, buf=0x00000000,n, u1;
  int l[8] = {0x00000007,0x00000070,0x00000700,0x00007000,0x00070000,0x00700000,0x07000000,0x70000000};

  while(1)
  {
    fgets(s,name_length,fastq); // Getting line with name
    fgets(s,chunk_size,fastq);  // Getting line with sequence
    if (feof(fastq)){break;}    // Break out if reach the end

    n=0;                        // Iterating over the characters in line
    for (i=0; (s[i]!='\n')&&(s[i]!='\0'); i++)
    {
      c = s[i];
      c = (c & l[0])<<4*(7-i%8); // Converting character to binary mask
      buf = buf | c;             // Summing up binary mask with previous result
      if (i%8==7) 
      {
        u1 = ~(U6^buf);
        u1 = (u1&(u1>>1)&(u1>>2)&B1);
        u1 = u1|(u1<<1)|(u1<<2)|(u1<<3);
        buf = buf&(~u1);

        fwrite(&buf,4,1,bin); 
        buf=0x00000000;
        n++;
      }   // The limit of the buffer is reached
    }
    i--;
    if ((i%8!=7)) 
    {
      fwrite(&buf,4,1,bin); 
      n++; 
    } // Wringing the last forgotten characters

    buf=0x20000000;
    fwrite(&buf,4,1,bin);
    buf=0x00000000;
    n++;
    fgets(s,chunk_size,fastq);
    fgets(s,chunk_size,fastq);
  }
  buf=0x50000000;
  //fwrite(&buf,4,1,bin);         // Intentional empty sequence

  return 0;
}