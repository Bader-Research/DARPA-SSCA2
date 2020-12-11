#include "simple.h"
#include <unistd.h>
#include "symbreak.h"

#define DO_SERIAL 0

#define DO_R 1
#define DO_S 1
#define DO_C 1

#define MIN_SIZE (1 << 16)
#define MAX_SIZE (1 << 25)
#define NUM_ITER        5

#define DEBUG 0

#define _INP_RANDOM   0
#define _INP_SEQ      1
#define _INP_CYCLIC02 2
#define _INP_CYCLIC04 3
#define _INP_CYCLIC08 4
#define _INP_CYCLIC16 5
#define _INP_CYCLIC32 6
#define _INP_CYCLIC64 7

char *inputLabel(const int input) {
  char *a;
  switch(input) {
  case _INP_RANDOM:   a = "  R"; break;
  case _INP_SEQ:      a = "  S"; break;
  case _INP_CYCLIC02: a = "C02"; break;
  case _INP_CYCLIC04: a = "C04"; break;
  case _INP_CYCLIC08: a = "C08"; break;
  case _INP_CYCLIC16: a = "C16"; break;
  case _INP_CYCLIC32: a = "C32"; break;
  case _INP_CYCLIC64: a = "C64"; break;
  default: a = "E";
  }
  return(a);
}

/* ---------------------------------------------------------------------------- 
   graph_color.cpp  
   
   Beginning of Ja Ja 3-coloring algorithm on shared-memory SMP. 
   pfh 
   
   Note: Tabs at 4, indent style BSD 
   
   To print out nicely, enscript -2rG --pretty-print -T4 
    However, enscript still does not get the tabs correct, sorry. 
   
   Change log: 
   7/99 Re-write after NT ate the previous version in a blue screen crash. 
   MS sucks. 
   Now back to coding under linux. 
   
   9/99 Ported to DEC unix. Only had to add the RAND_MAX hack.  
   See notes in header file for details. 
   3 cheers for portable code! 
   
   As of 8/99, serial version working and verified. Now to parallelize 
   this bad boy. 
   
   As of 10/10/99, correct on 256 or fewer points, all routines parallel. 
   Removing some of the now-unused debug code, etc. 
   Still need a bit of code to time the actual algorithm running. 
   
   10-12-99 or so, all code working, added loops to run par/serial code 
   N times and time the results. 
    
   10-17-99 Rewrote serial version to near optimal, added std. dev. calc. 
   from numerical recipes, timings on Alpha 8400 show parallel code never 
   faster...oops. 

   10-19-99 Incorporating DAB improvements to see if we can speed it up. 
   Also fixed N==1 sigma bug. 
   Machine is pretty busy, load > 6, but around 4M the parallel wins now! 
*/

/* ---------------------------------------------------------------------------- */


/* Defining this disables status printouts and assert macros */
#if DEBUG
#undef NDEBUG
#else
#define NDEBUG
#endif


#if DEBUG
/* ---------------------------------------------------------------------------- */
/* Util fn for debugging, flushes stdout after each. */
/* Disabled if NDEBUG defined, which also disables assert macros. */
void debug_print(const char *str, THREADED)
{

#if !defined(NDEBUG)
    on_one_thread
    {
	fprintf(stdout,"%s", (char *) str);
	fflush(stdout);
    }
#endif

    return;
}
#endif

/* ---------------------------------------------------------------------------- */
/* DA Bader's SMP counting sort, modified to sort structs instead of */
/* just ints.  */
/* 10-9-99 Modifying this to return the final histogram - need it */
/* for the final pardo in the Ja Ja algorithm. */
/* Obscure commentary: You Are Not Expected to Understand This. */

/****************************************************/
/* R (range)      must be a multiple of SMPS */
/* q (elems/proc) must be a multiple of SMPS */
/****************************************************/
void all_countsort_smp(int q,
		       s_struct *lKey,
		       s_struct *lSorted,
		       int *pdf,
		       int R,
		       int bitOff, 
		       int m,
		       THREADED)
{
    register int
	i,
	j,
	k,
        last, temp,
	offset;
    
    int *myHisto,
        *mhp,
        *mps,
        *psHisto,
        *allHisto,
	*ptr;

    myHisto  = (int *)node_malloc(THREADS*R*sizeof(int), TH);
    psHisto  = (int *)node_malloc(THREADS*R*sizeof(int), TH);

    mhp = myHisto + MYTHREAD*R;

    for (k=0 ; k<R ; k++)
	mhp[k] = 0;
    
    pardo(k, 0, q, 1)
	mhp[bits(lKey[k].c,bitOff,m)]++;

    node_Barrier();

    pardo(k, 0, R, 1) 
	{
	    last = psHisto[k] = myHisto[k];
	    for (j=1 ; j<THREADS ; j++) 
	    {
		temp = psHisto[j*R + k] = last + myHisto[j*R +  k];
		last = temp;
	    }
	}

    allHisto = psHisto+(THREADS-1)*R;
    
    node_Barrier();

    offset = 0;

    mps = psHisto + (MYTHREAD*R);
    for (k=0 ; k<R ; k++) 
    {
	mhp[k]  = (mps[k] - mhp[k]) + offset;
	offset += allHisto[k];
    }
    
    /* pfh - as per discussions with DAB, at this point we can */
    /* score the final PDF. Use pardo to copy it to a smaller */
    /* array that we will pass back to the caller. */

    /* This j term is constant; have all threads compute it before they */
    /* run the loop. Moved the barrier just above to after this  */
    /* calculation. */

    /* As per 10-19, replaced with DAB's tweaked version */
    ptr = myHisto + (THREADS - 1)*R;
    node_Barrier();

    pardo(i, 0, R, 1)
	pdf[i] = *(ptr + i);

    node_Barrier();

    /* Back to stock code */
    pardo(k, 0, q, 1) 
	{
	    j = bits(lKey[k].c,bitOff,m);

	    /* 
	       pfh mod - Replace shallow with deep copy 

	       Original code: 
	       lSorted[mhp[j]] = lKey[k];

	    */
	    lSorted[mhp[j]].c = lKey[k].c;
	    lSorted[mhp[j]].v = lKey[k].v;
	    mhp[j]++;
	}

    node_Barrier();

    node_free(psHisto, TH);
    node_free(myHisto, TH);
}

/* ---------------------------------------------------------------------------- */
/* Simple debug routine to print a graph, inorder or sequential or both. */
void print_graph(bool inorder, bool sequential, const int *S, const int num_pts)
{
    int i, current;

    if(inorder)
    {
      /* Start at first node in array */
	current = 0;

	fprintf(stdout,"\nInorder traversal: \n");
	for(i = 0; i < num_pts ; i++)
	{
	    fprintf(stdout,"%d->", S[current]);
	    current = S[current];
	}
    }

    if(sequential)
    {
	fprintf(stdout,"\nSequential listing:\n");

	for(i = 0; i < num_pts; i++)
	    fprintf(stdout,"%d ", S[i]);
    }

    fprintf(stdout,"\nEnd of graph listing.\n");
    return;
}

/* ---------------------------------------------------------------------------- */
/* Expanded version, prints colors & predecessors as well */
void print_full_graph(const int *S, const int *P, const int *C, const int num_pts)
{
    int cur_vertex;
    int i;


    cur_vertex = 0;

    cur_vertex = P[S[cur_vertex]];

    fprintf(stdout,"Complete graph:\n(Pred)<-[ curr ]->(succ) color\n\n");

    for(i = 0; i < num_pts; i++)
    {
	fprintf(stdout,"(%2d)<-[%2d]->(%2d), color %2d\n",
	       P[cur_vertex], cur_vertex,S[cur_vertex], C[cur_vertex]);

	cur_vertex = S[cur_vertex];
    }

    fprintf(stdout,"\nEnd of graph listing.\n");
    return;
}

/* ---------------------------------------------------------------------------- */
/* Test function - count the unique colors in a given colormap */
/* as a simple check of algorithm efficacy. */
/* To do: Check predecessor/successor colors also. */
/* Note that this routine starts at array location zero, since the */
/* colormaps use zero also. */
int count_colors(const int *colors, const int num_pts)
{
    int total_used;
    int *ref_cnt;
    int i;


    ref_cnt = (int *) malloc(sizeof(int) * num_pts);
    assert(ref_cnt != NULL);

    for(i = 0; i < num_pts; i++)
	ref_cnt[i] = 0;

    /* Walk the color list, carefully incrementing ref counts */
    for(i = 0; i < num_pts; i++)
    {
	assert(colors[i] >= 0);
	assert(colors[i] <= num_pts);

	ref_cnt[colors[i]]++;
    }

    /* Compute total used */
    total_used = 0;
    for(i = 0; i < num_pts; i++)
    {
	if(ref_cnt[i] > 0)
	    total_used++;
    }

    free(ref_cnt);

    return(total_used);
}

/* ---------------------------------------------------------------------------- */
/* Function from D Bader to permute a given array.  */
/* This does an in-place modification of the argument, and expects an array of n+1  */
/* points from 0 to n. */
void permute_array(int *S, const int num_pts)
{
    int		i, rnum, itmp;
    double  dtmp;
	
    /* 
       Seed the random-number generator with current time so that
       the numbers will be different every time we run.
       Code snipped from MSDN.

       Per DAB, use srandom and random instead of rand/srand. Not in 
       C standard, but better PRNG.
    */

    /* srandom should return zero if init'd ok */
#if 0
    assert((srandom(time(NULL)) == 0));
#else
    srandom(time(NULL));
#endif

    for(i = num_pts - 1; i >= 1; i--)
    {
      /* Get random number from 0 to 1; map to [0, i] */
	dtmp = (double) random() / RAND_MAX;
	dtmp = dtmp * i;
	rnum = (int) floor(dtmp);

	assert(rnum >= 0);
	assert(rnum <  num_pts);

	/* Swap s[i], s[j] */
	itmp = S[i];
	S[i] = S[rnum];
	S[rnum] = itmp;
    }

    return;
}

/* ---------------------------------------------------------------------------- */
/* Generate a graph, in the format of Ja Ja's successor array. */
/* Caller responsible for freeing the returned vector. */
int *generate_array(const int num_pts, const int input, THREADED)
{
    int *result;
    int i;
    
    result = (int *) node_malloc(sizeof(int) * num_pts , TH);

    switch(input) {
    case _INP_RANDOM : {
      int	*list;
      list = (int *) node_malloc(sizeof(int) * num_pts , TH);

      pardo(i, 0, num_pts, 1) {
	list[i] = i;
      }

      node_Barrier();
    
      on_one_thread {
	/* Permute the array */
	permute_array(list, num_pts);
	
	/* Pad last entry with first - wraparound padding */
	/* list[num_pts] = list[0]; */
	    
	/* Algorithm, courtesy of Chris Hurlburt, to convert  */
	/* the permuted vertex list into a valid successor array. */
	for(i = 0; i < num_pts-1 ; i++)
	  result[list[i]] = list[i + 1];
	result[list[num_pts-1]] = list[0];
      }
    
      /* Free permuted vector listing */
      node_free(list, TH);
    }
    break;
    case _INP_SEQ : {
      on_one_thread {
	/* Algorithm, courtesy of Chris Hurlburt, to convert  */
	/* the permuted vertex list into a valid successor array. */
	for(i = 0; i < num_pts-1 ; i++)
	  result[i] = i + 1;
	result[num_pts-1] = 0;
      }
    }
    break;
    case _INP_CYCLIC02 : {
      on_one_thread {
	for (i=0; i < num_pts - 2 ; i++)
	  result[i] = i + 2;
	for (i=0 ; i<2 ; i++)
	  result[num_pts-2+i] = (i+1) % 2;
      }
    }
    break;
    case _INP_CYCLIC04 : {
      on_one_thread {
	for (i=0; i < num_pts - 4 ; i++)
	  result[i] = i + 4;
	for (i=0 ; i<4 ; i++)
	  result[num_pts-4+i] = (i+1) % 4;
      }
    }
    break;
    case _INP_CYCLIC08 : {
      on_one_thread {
	for (i=0; i < num_pts - 8 ; i++)
	  result[i] = i + 8;
	for (i=0 ; i<8 ; i++)
	  result[num_pts-8+i] = (i+1) % 8;
      }
    }
    break;
    case _INP_CYCLIC16 : {
      on_one_thread {
	for (i=0; i < num_pts - 16 ; i++)
	  result[i] = i + 16;
	for (i=0 ; i<16 ; i++)
	  result[num_pts-16+i] = (i+1) % 16;
      }
    }
    break;
    case _INP_CYCLIC32 : {
      on_one_thread {
	for (i=0; i < num_pts - 32 ; i++)
	  result[i] = i + 32;
	for (i=0 ; i<32 ; i++)
	  result[num_pts-32+i] = (i+1) % 32;
      }
    }
    break;
    case _INP_CYCLIC64 : {
      on_one_thread {
	for (i=0; i < num_pts - 64 ; i++)
	  result[i] = i + 64;
	for (i=0 ; i<64 ; i++)
	  result[num_pts-64+i] = (i+1) % 64;
      }
    }
    break;
    default: fprintf(stderr,"ERROR: No input selected\n"); exit(-1);
    }
    
    /* Return S vector */
    return(result);
}


/* ---------------------------------------------------------------------------- */
/* Function to generate a predecessor list from an S-list. */
/* As allways, caller frees returned list. */
/* Serial algorithm is link dragging, order N. */
/* Parallel version uses pardo for N/P, but the memory references */
/* are so randomized as to parallelize poorly...or so I suspect.  */
/* This could stand some investigation. */
void S_to_P(const int *S, int *P, const int num_pts, THREADED)
{
    int i;
    
    /* node_Barrier(); */

    pardo(i, 0, num_pts, 1)
	P[S[i]] = i;

    /* node_Barrier(); */

    return;
}

/* ---------------------------------------------------------------------------- */
/* Function to return the LSB where two integers differ. */
/* Assumes 8-bit bytes, does not assume ints are any particular size,  */
/* can probably be made faster - this is linear in the */
/* number of bits. Returns -1 if numbers are equal. */
int lsb_diff(const int in_1, const int in_2)
{
    const int num_bits = sizeof(int) * 8;

    int tmp1, tmp2;
    int	xor;
    int i;

    /* Suggestion from DAB - do XOR outside loop */
    xor = (in_1 ^ in_2);

    for(i = 0; i < num_bits; i++)
    {
      /* Shift right so we can see LSB of interest */
	tmp1 = (xor >> i);

	/* Mask off all except lowest bit */
	tmp2 = (tmp1 & 0x1);

	if(tmp2 == 1)
	    return(i);
    }

    /* If got here, all bits equal. Gack. */
    return(-1);
}

/* ---------------------------------------------------------------------------- */
/* Util function, picks the color for a vertex from the rules of Ja Ja,  */
/* page 80 */
/* Only needs to know the colors of its neighbors, pass as const */
/* 10-19-99, new version optimized by DAB for speed, a little less readable */
/* but much less work required. */
INLINE int color_vertex(const int pred, const int succ)
{
    register int retval;

    if((pred + succ) == 1)
	retval = 2;
    else
	retval = ((!pred) || (!succ));

#if DEBUG
    assert(retval >= 0);
    assert(retval <= 2);
    assert(retval != pred);
    assert(retval != succ);
#endif

    return(retval);
}

/* ---------------------------------------------------------------------------- */
/* Verify a colormap by inorder traversal */
/* False for bad graph, true for good. Does not need to check both P/S arrays,  */
/* but does so anyway to catch moron errors. */
bool verify_colormap(const int *S, const int *P, const int *C, const int num_pts)
{
    int current = 0;
    int i;

    for(i = 0; i < num_pts; i++)
    {
	if((C[P[current]] == C[current]) || (C[S[current]] == C[current]))
	    return(false);
	
	current = S[current];
    }

    return(true);
}
    
/* ---------------------------------------------------------------------------- */
/* Serial version */
/* As above, caller must free the returned array. */
void serial_three_color(const int *S, const int num_pts, int *init_colors)
{
#if 0
  int *init_colors = NULL;
#endif
    int idx;
    int current_vertex = 0;
    
#if 0
    /* Allocate memory */
    init_colors = (int *) malloc(sizeof(int) * num_pts);
    assert(init_colors != NULL);
#endif
    
    /* ------- */
    /* Coloring algorithm from JaJa pg75, my mod follows. */
    for(idx = 0; idx < num_pts-1 ; idx++)
    {
	init_colors[current_vertex] = (idx & 0x1);
	current_vertex = S[current_vertex];
    }

    /* My mod - color the last point 2, requires no decisions */
    init_colors[current_vertex] = 2;
    
#if 0
    return(init_colors);
#else
    return;
#endif
}

/* ---------------------------------------------------------------------------- */
/* Algorithm 2.10 from ja ja for the 3-coloring of a simple cycle. */
/* As above, caller must free the returned array. */
/* Assumes eight-bit bytes. */
void three_color(const int *S, int *second_colors, const int num_pts, THREADED)
{
  /* Node-malloc'd memory */
  int *init_colors = NULL;
    int *pdf = NULL;

    int *P = NULL;
    s_struct *bag_in = NULL;
    s_struct *bag_out = NULL;

    /* Auto (stack) variables */
    int cur_color, num_elements, start_idx, end_idx;
    int i, idx, k;
    int up_limit;
    /* As of 10-19, these are local/auto to each thread. */
    int cur_vertex,
	pred, 
	succ;
#if 0
    double t0;
#endif

    /* Allocate memory for colormaps */
    init_colors = (int *) node_malloc(sizeof(int) * num_pts, TH);
    second_colors = (int *) node_malloc(sizeof(int) * num_pts, TH);
    
    /* Predecessor array */
    P = (int *) node_malloc(sizeof(int) * num_pts, TH);

    /* Sort structs - 2-ints, vertex & color */
    bag_in = (s_struct *) node_malloc(sizeof(s_struct) * num_pts, TH);
    bag_out = (s_struct *) node_malloc(sizeof(s_struct) * num_pts, TH);

    /* 256 ints, holding PDF of sorted color array */
    pdf = (int *) node_malloc(sizeof(int) * NUM_HBINS, TH);

#if 0
    t0 = get_seconds();
#endif
    
    /* --------------------------------------------------- */
    /* Step one: Init each color to its index */
#if DEBUG
    debug_print("\nStep 1: Init the colormap.", TH);
#endif
    
    node_Barrier();

    pardo(i, 0, num_pts, 1)
	init_colors[i] = i;

    /* --------------------------------------------------- */
    /* Step two: Compute the reduced colormap */
#if DEBUG
    debug_print("\nStep 2: Applying algorithm 2.9...", TH);
#endif
    
    /* Run the Ja Ja algorithm 2.9 to reduce the colormap */
    node_Barrier();

    pardo(i, 0, num_pts, 1)
    {
	k = lsb_diff(init_colors[i], init_colors[S[i]]);

	second_colors[i] = (k << 1) + kth_lsb(init_colors[i], k);
    }

    node_Barrier();

    /* Done with the initial colormap */
    node_free(init_colors, TH);

#if DEBUG
    /* If NDEBUG defined, skip checks and printouts */
#if !defined(NDEBUG)

    /* DDT run diags on one CPU */
    on_one_thread
    {
	fprintf(stdout,"done. Now %d colors.",
	       count_colors(second_colors, num_pts));

	if(num_pts <= SIZE_TOOBIG)
	{
	    fprintf(stdout,"\nColormap after 2.9: ");
	    for(i = 0; i < num_pts; i++)
		fprintf(stdout,"%d ", second_colors[i]);
	    
	    fprintf(stdout,"\n");
	}
	
    }
#endif
#endif
    
    /* --------------------------------------------------- */
    /* Ancillary step - generate predecessor array */
#if DEBUG
    debug_print("\nStep 3a: Generating P array...", TH);
#endif
    
    node_Barrier();
    S_to_P(S, P, num_pts, TH);

    /* Now, step 3 of Ja Ja 2.10 -- sort graph by colors */
    /* First, copy data into vector of new structs. */
#if DEBUG
    debug_print("\nStep 3b: Copying data into struct via pardo...", TH);
#endif
    
    node_Barrier();
    pardo(i, 0, num_pts, 1)
    {
	bag_in[i].c = second_colors[i];
	bag_in[i].v = i;
    }

    /* ---------------------------------------------------     */
#if DEBUG
    debug_print("\nStep 3c: Sorting graph by color... ", TH);
#endif
    
    /*
      Use DAB's sort - since we have up_limit as the max color, we
      only need to sort the LS byte - factor of eight speedup!
      Once this is debugged, we can lower the bound even more - 
      we only need ~40 colors, so could lower the bound to 64.

      Note the pointer arithmetic: sorting code uses 0 to n-1,
      so send it base+1 as starting point.
    */
    node_Barrier();

    all_countsort_smp(num_pts, bag_in, bag_out, pdf,
		      NUM_HBINS,  0, NUM_HBITS, TH);

    /* Max color that is possible from step 2.9 */
    up_limit = 2 * (int) ceil((log( (double) num_pts) / log(2.0)));

    /* Done with input list, free soonest */
    node_free(bag_in, TH);

    /* --------------------------------------------------- */
    /* Step 4 begins */
#if DEBUG
    debug_print("\nStep 4: 3-coloring (Alg 2.10)...", TH);
#endif
    
    /* Set current color to starting color */
    cur_color = 3;

    /* Main loop - pick a color, recolor all nodes present */
    while(cur_color <= up_limit)
    {
      /* Extract start/stop indices from PDF */
	start_idx = pdf[cur_color - 1] ;
	end_idx = pdf[cur_color] - 1 ;
	num_elements = (end_idx - start_idx) + 1;

	/* Make sure all nodes sync'd w/correct range */
	node_Barrier();

#if DEBUG
#if !defined(NDEBUG)
	/* debug printout if graph small enough */
	on_one_thread
	{
	    if(num_pts <= SIZE_TOOBIG)
		fprintf(stdout,"\nColor: %2d start_idx: %d end_idx: %d elements: %d",
		       cur_color, start_idx, end_idx, num_elements);

	    if((num_elements <= 0) && (num_pts <= SIZE_TOOBIG))
		fprintf(stdout," - skipping.");

	    fflush(stdout);
	}
#endif
#endif
	/* Test if this color isn't present */
	if(num_elements <= 0)
	{
	    cur_color++;
	    continue;
	}

	/* Ending condition - out of vertices to process */
	if(start_idx > num_pts)
	    break;

	/* Loop over all vertices of this color, in parallel,  */
	/* and recolor them. */
	pardo(idx, start_idx, (end_idx+1), 1)
	{
	  /* Look up each threads' vertex */
	    cur_vertex = bag_out[idx].v;

	    /* From there, get the neighbors */
	    /* This is cache hell, but it is read-only */
	    pred = second_colors[P[cur_vertex]];
	    succ = second_colors[S[cur_vertex]];

	    /* Use these to color the current vertex */
	    /* This is worse, cache-wise, and all writing */
	    second_colors[cur_vertex] = color_vertex(pred, succ);
	}

	/* Go to next color */
	cur_color++;
    }

    /* Threads sync here after breaking out of the above while loop */
    node_Barrier();

#if 0
    t0 = get_seconds() - t0;
    on_one_thread ffprintf(stdout,stdout,"Time: %9.6f\n",t0);
#endif
    
#if DEBUG
    /* If NDEBUG defined, skip all checks */
#if !defined(NDEBUG)

    /* Do graph diagnostics on cpu zero only */
    on_one_thread
    {
	debug_print("Done!\n\n", TH);
	
	i = count_colors(second_colors, num_pts);
	
	fprintf(stdout,"%d colors after algorithm 2.10 - ", 
	       count_colors(second_colors, num_pts), i);
	
	if(i > 3)
	{
	    fprintf(stdout,"Error: should have 3 colors!\n", i);

	    fprintf(stdout,"Last vertex is colored %d\n\n", second_colors[num_pts]);
	}
	else
	    fprintf(stdout,"Correct.\n");
	
	
	if(num_pts <= SIZE_TOOBIG)
	    print_full_graph(S, P, second_colors, num_pts);

	debug_print("\nVerifying colormap...", TH);

	if(verify_colormap(S, P, second_colors, num_pts) == false)
	    debug_print("Error in colormap!!\n", TH);
	else
	    debug_print("OK\n", TH);
    }
#endif
#endif
    /* Housecleaning */
    node_free(pdf, TH);
    node_free(P, TH);
    node_free(bag_out, TH);

    /* Make sure we're all done */
    node_Barrier();
    
    return;
}

/* ---------------------------------------------------------------------------- */
/* Script function - run the algorithm at a given size for a given
   number of iterations, and return results & stats for same.
*/
void run_the_code(const int num_pts, const int num_iterations, 
		  double *min_time,
		  double *max_time,
		  double *avg_time,
		  double *sigma,
		  bool do_serial,
		  const int input,
		  THREADED)
{
    int idx;
    double start_time, run_time, ttl_time;
    int *data_set = NULL;
    int *color_map = NULL;
    double *run_times = NULL;
    double dtmp, sum;

    ttl_time = 0.0;

    /* Set min, max to values such that will be replaced on the first run */
    *min_time = (double) num_pts * 10.0;
    *max_time = -0.01;

    /* Allocate memory to hold runtimes - don't know how much */
    /* we will need until runtime. */
    on_one_thread
    {
	run_times = (double *) malloc(sizeof(double) * num_iterations);
	assert(run_times != NULL);
    }

    for(idx = 0; idx < num_iterations; idx++)
    {
      /* Generating the data set is not counted in the runtime. */
	data_set = generate_array(num_pts, input, TH);

	/* Do we want serial or parallel version? */
	if(do_serial == true)
	{
	    on_one_thread
	    {
	      color_map = (int *) malloc(sizeof(int) * (num_pts + 1));
	      assert(color_map != NULL);
	        
		start_time = get_seconds();
		
		serial_three_color(data_set, num_pts, color_map);

		run_time = get_seconds() - start_time;

		free(color_map);
	    }
	}
	else /* Parallel code */
	{
#if 0
	  all_Barrier(TH);
#else
	  node_Barrier();
#endif

	    start_time = get_seconds();
	
	    three_color(data_set, color_map, num_pts, TH);
	}

	run_time = get_seconds() - start_time;

	/* Save time for stats later */
	on_one_thread
	    run_times[idx] = run_time;

        /* Update min, max, total times - do these on one thread. */
	on_one_thread
	{
	    if(run_time < *min_time)
		*min_time = run_time;
	    
	    if(run_time > *max_time)
		*max_time = run_time;
	    
	    ttl_time += run_time;
	}

	/* Housecleaning */
	node_free(data_set, TH);

	/* Serial code uses malloc, so different cleanup */
	if(do_serial == true) {
#if 0
	    on_one_thread
		free(color_map);
#endif
	}
	else
	    node_free(color_map, TH);
    }

    /* Compute stats */
    on_one_thread
    {
      /* Mean */
	*avg_time = (double) (ttl_time / num_iterations);

	/* Compute std deviation */
	sum = 0.0;

	for(idx = 0; idx < num_iterations; idx++)
	{
	    dtmp = (run_times[idx] - *avg_time);
	    sum += (dtmp * dtmp);
	}

	/* Sigma undefined for N==1 */
	if(num_iterations > 1)
	    *sigma = sqrt(sum) * (1.0 / (num_iterations - 1.0));
	else
	    *sigma = *avg_time;

	free(run_times);
    }

#if 0
    all_Barrier(TH);
#else
    node_Barrier();
#endif
    
    return;
}


/* ---------------------------------------------------------------------------- */
/* Main - create data, run the algorithm, save the stats. */
void *symbreak(const int input, int cur_size, int num_iterations, int testrun, THREADED)
{
    double ser_avg, par_avg, ser_min, par_min, ser_max, par_max, 
      ser_sig, par_sig;

#if DO_SERIAL
    /* Serial first */
    run_the_code(cur_size,
		 num_iterations,
		 &ser_min, &ser_max, &ser_avg, &ser_sig, true, input, TH);

    if (!testrun) {
      on_one_thread {
	fprintf(stdout,"%3s\t%4d\tSer\t%12d\t%9.6f\t%9.6f\t%9.6f\t%9.6f\n",
		inputLabel(input),
		THREADS,
		cur_size,
		ser_min,
		ser_max,
		ser_avg,
		ser_sig);
	fflush(stdout);
      }
    }
#endif

    /* Parallel timings */
    run_the_code(cur_size,
		 num_iterations,
		 &par_min, &par_max, &par_avg, &par_sig, false, input, TH);

    if (!testrun) {
      on_one_thread {
	fprintf(stdout,"%3s\t%4d",
		inputLabel(input),
		THREADS);
#if DO_SERIAL
	fprintf(stdout,"\tPar");
#endif
	fprintf(stdout,"\t%12d\t%9.6f\t%9.6f\t%9.6f\t%9.6f",
		cur_size,
		par_min,
		par_max,
		par_avg,
		par_sig);

#if DO_SERIAL
	/* Div by zero if size too small - oops! */
	if(par_min < 0.0001)
	  par_min = 0.0001;
	  fprintf(stdout,"\t%9.6f",ser_min/par_min);
#endif

	  fprintf(stdout,"\n");
	  fflush(stdout);
      }
    }
		     
    return;
}

void *SIMPLE_main(THREADED) {
  const int min_size = MIN_SIZE;
  const int max_size = MAX_SIZE;
  const int num_iterations = NUM_ITER;

  int cur_size;

  on_one_thread {
    fprintf(stdout,"Parallel and serial three-coloring.\n");
    fprintf(stdout,"Min size: %d Max size: %d Iterations per size: %d\n\n",
	    min_size, max_size, num_iterations);

    fprintf(stdout,"INP\t   T");
#if DO_SERIAL
    fprintf(stdout,"\tP/S");
#endif
    fprintf(stdout,"\tSize\tMinimum\tMaximum\tAverage\tSigma");
#if DO_SERIAL
    fprintf(stdout,"\tSpeedup");
#endif
    fprintf(stdout,"\n");
    fflush(stdout);
  }

  for (cur_size = min_size; cur_size <= max_size; cur_size *= 2) {
  
#if DO_R
    symbreak(_INP_RANDOM,   cur_size, num_iterations, 0, TH);
    symbreak(_INP_RANDOM,   cur_size, num_iterations, 1, TH);
#endif
#if DO_S
    symbreak(_INP_SEQ,      cur_size, num_iterations, 0, TH);
    symbreak(_INP_SEQ,      cur_size, num_iterations, 1, TH);
#endif
#if DO_C
    symbreak(_INP_CYCLIC02, cur_size, num_iterations, 0, TH);
    symbreak(_INP_CYCLIC02, cur_size, num_iterations, 1, TH);
#endif
#if DO_C
    symbreak(_INP_CYCLIC04, cur_size, num_iterations, 0, TH);
    symbreak(_INP_CYCLIC04, cur_size, num_iterations, 1, TH);
#endif
#if DO_C
    symbreak(_INP_CYCLIC08, cur_size, num_iterations, 0, TH);
    symbreak(_INP_CYCLIC08, cur_size, num_iterations, 1, TH);
#endif
#if DO_C
    symbreak(_INP_CYCLIC16, cur_size, num_iterations, 0, TH);
    symbreak(_INP_CYCLIC16, cur_size, num_iterations, 1, TH);
#endif
#if DO_C
    symbreak(_INP_CYCLIC32, cur_size, num_iterations, 0, TH);
    symbreak(_INP_CYCLIC32, cur_size, num_iterations, 1, TH);
#endif
#if DO_C
    symbreak(_INP_CYCLIC64, cur_size, num_iterations, 0, TH);
    symbreak(_INP_CYCLIC64, cur_size, num_iterations, 1, TH);
#endif
  }
  
  return;
}
