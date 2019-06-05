unsigned int find_mismatches( unsigned int * cc1, unsigned int * cc2, unsigned int n , unsigned int * v);
unsigned int shift_right( unsigned int x, unsigned int * cc, unsigned int ** cc_mod, unsigned int n );
unsigned int shift_left( unsigned int x, unsigned int * cc, unsigned int ** cc_mod, unsigned int n );
unsigned int find_mismatches_in_region( unsigned int * cc1, unsigned int * cc2, unsigned int n, unsigned int * v, unsigned int cc1_from, unsigned int cc2_from);
unsigned int find_mismatches_in_region_inversed( unsigned int * cc1, unsigned int * cc2, unsigned int n, unsigned int * v, unsigned int cc1_to, unsigned int cc2_to);

unsigned int find_( unsigned int * cc1,  // seq1
                    unsigned int * cc2,  //seq2
                    unsigned int n,      //number of bits to compare
                    unsigned int * v,    // pointer to save count
                    unsigned int cc1_to, unsigned int cc2_to, // start positions to compare
                    unsigned int mode,   // mode - 1: from end to start; 2: from start to end
                    unsigned int max_subst // when to break search
                    );

unsigned int find_mask( unsigned int * cc1, unsigned int * cc2,
                        unsigned int n, unsigned int * v, unsigned int cc1_to,
                        unsigned int cc2_to, unsigned int mode, unsigned int max_subst,
                        unsigned int mask_start, unsigned int mask_end);

unsigned int align_seq( unsigned int * seq1,
                        unsigned int * seq2_for,
                        unsigned int * seq2_rev,
                        unsigned int n,             // number of bits to compare
                        unsigned int seq1_to, unsigned int seq2_to, // start position for comparison
                        unsigned int mode,
                        unsigned int max_subst,
                        unsigned int,
                        unsigned int * start, unsigned int * end
                        );

unsigned int interpret_buf(int buf);

const char *byte_to_binary(int x);