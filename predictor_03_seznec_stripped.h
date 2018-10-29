#ifdef INCLUDEPRED
/*** Stripped all ISL components ***/

/*
Code has been succesively derived from the tagged PPM predictor simulator from Pierre Michaud, the OGEHL predictor simulator from by André Seznec, the TAGE predictor simulator from  André Seznec and Pierre Michaud

*/

#include <inttypes.h>
#include <math.h>
#define SHARINGTABLES		// let us share the physical among several logic predictor tables
#define INITHISTLENGTH		// uses the "best" history length we found

#define NHIST 15		// 15 tagged tables + 1 bimodal table
#define LOGB 15			// log of number of entries in bimodal predictor
#define HYSTSHIFT 2		// sharing an hysteris bit between 4 bimodal predictor entries
#define LOGG (LOGB-4)		// initial definition was with 2K entries per tagged tables

//The Best  set of history lengths
int m[NHIST + 1] = {0, 3, 8, 12, 17, 33, 35, 67, 97, 138, 195, 330, 517, 1193, 1741, 1930};


#ifndef INITHISTLENGTH
#define MINHIST 8		// shortest history length
#define MAXHIST 2000		// longest history length
#endif

#define PHISTWIDTH 16		// width of the path history


#ifndef SHARINGTABLES
#define TBITS 6
#define MAXTBITS 15
#endif

#define CWIDTH 3		// predictor counter width on the tagged tables
#define HISTBUFFERLENGTH 4096	// we use a 4K entries history buffer to store the branch history


// utility class for index computation
// this is the cyclic shift register for folding
// a long global history into a smaller number of bits; see P. Michaud's PPM-like predictor at CBP-1
class folded_history
{
public:
  unsigned comp;
  int CLENGTH;
  int OLENGTH;
  int OUTPOINT;

  folded_history ()
  {
  }

  void init (int original_length, int compressed_length)
  {
    comp = 0;
    OLENGTH = original_length;
    CLENGTH = compressed_length;
    OUTPOINT = OLENGTH % CLENGTH;
  }

  void update (uint8_t * h, int PT)
  {
    comp = (comp << 1) | h[PT & (HISTBUFFERLENGTH - 1)];
    comp ^= h[(PT + OLENGTH) & (HISTBUFFERLENGTH - 1)] << OUTPOINT;
    comp ^= (comp >> CLENGTH);
    comp &= (1 << CLENGTH) - 1;
  }
};


class bentry			// TAGE bimodal table entry
{
public:
  int8_t hyst;
  int8_t pred;
  bentry ()
  {
    pred = 0;
    hyst = 1;
  }
};

class gentry			// TAGE global table entry
{
public:
  int8_t ctr;
  uint16_t tag;
  int8_t u;
  gentry ()
  {
    ctr = 0;
    tag = 0;
    u = 0;
  }
};

int8_t USE_ALT_ON_NA;		// "Use alternate prediction on newly allocated":  a 4-bit counter  to determine whether the newly allocated entries should be considered as  valid or not for delivering  the prediction
int TICK, LOGTICK;		//control counter for the smooth resetting of useful counters
int phist;			// use a path history as on  the OGEHL predictor
uint8_t ghist[HISTBUFFERLENGTH];
int Fetch_ptghist;
int Fetch_ptghistos;
int Fetch_phist;		//path history
int Fetch_phistos;		//path os history
folded_history Fetch_ch_i[NHIST + 1];	//utility for computing TAGE indices
folded_history Fetch_ch_t[2][NHIST + 1];	//utility for computing TAGE tags
int Retire_ptghist;
int Retire_phist;		//path history
folded_history Retire_ch_i[NHIST + 1];	//utility for computing TAGE indices
folded_history Retire_ch_t[2][NHIST + 1];	//utility for computing TAGE tags


//For the TAGE predictor
bentry *btable;			//bimodal TAGE table
gentry *gtable[NHIST + 1];	// tagged TAGE tables
int TB[NHIST + 1];		// tag width for the different tagged tables
int logg[NHIST + 1];		// log of number entries of the different tagged tables

int GI[NHIST + 1];		// indexes to the different tables are computed only once
int GTAG[NHIST + 1];		// tags for the different tables are computed only once
int BI;				// index of the bimodal table

bool pred_taken;		// prediction
bool alttaken;			// alternate  TAGEprediction
bool tage_pred;			// TAGE prediction
bool LongestMatchPred;
int HitBank;			// longest matching bank
int AltBank;			// alternate matching bank


class my_predictor
{
public:
  my_predictor (void)
  {
    USE_ALT_ON_NA = 0;
    LOGTICK = 8;

    TICK = 0;
    Fetch_phist = 0;
    Retire_phist = 0;
    for (int i = 0; i < HISTBUFFERLENGTH; i++)
      ghist[0] = 0;
    Fetch_ptghist = 0;
    Fetch_ptghistos = 0;
#ifndef INITHISTLENGTH
    m[1] = MINHIST;
    m[NHIST] = MAXHIST;
    for (int i = 2; i <= NHIST; i++)
    {
      m[i] = (int) (((double) MINHIST *
                     pow ((double) (MAXHIST) / (double) MINHIST,
                          (double) (i - 1) / (double) ((NHIST - 1)))) + 0.5);
    }
    for (int i = 2; i <= NHIST; i++)
      if (m[i] <= m[i - 1] + 2)
        m[i] = m[i - 1] + 2;
#endif

#ifndef SHARINGTABLES
    for (int i = 1; i <= NHIST; i++)
      TB[i] = TBITS + (i - 1);
    TB[1]++;
    for (int i = 1; i <= NHIST; i++)

      if (TB[i] > MAXTBITS)
        TB[i] = MAXTBITS;
    // log2 of number entries in the tagged components
    for (int i = 1; i <= 3; i++)
      logg[i] = LOGG;
    for (int i = 4; i <= 6; i++)
      logg[i] = LOGG + 1;
    for (int i = 7; i <= 10; i++)
      logg[i] = LOGG;
    for (int i = 10; i <= NHIST; i++)
      logg[i] = LOGG - 1;
    for (int i = 1; i <= NHIST; i++)
    {
      gtable[i] = new gentry[1 << (logg[i])];
    }
#endif


    //initialisation of the functions for index and tag computations

    for (int i = 1; i <= NHIST; i++)
    {
      Fetch_ch_i[i].init (m[i], (logg[i]));
      Fetch_ch_t[0][i].init (Fetch_ch_i[i].OLENGTH, TB[i]);
      Fetch_ch_t[1][i].init (Fetch_ch_i[i].OLENGTH, TB[i] - 1);
    }

    for (int i = 1; i <= NHIST; i++)
    {
      Retire_ch_i[i].init (m[i], (logg[i]));
      Retire_ch_t[0][i].init (Retire_ch_i[i].OLENGTH, TB[i]);
      Retire_ch_t[1][i].init (Retire_ch_i[i].OLENGTH, TB[i] - 1);
    }

    //allocation of the all predictor tables
    btable = new bentry[1 << LOGB];
  }


  // index function for the bimodal table

  int bindex (uint32_t pc)
  {
    return ((pc) & ((1 << (LOGB)) - 1));
  }


  // the index functions for the tagged tables uses path history as in the OGEHL predictor
  //F serves to mix path history
  int F (int A, int size, int bank)
  {
    int A1, A2;
    A = A & ((1 << size) - 1);
    A1 = (A & ((1 << logg[bank]) - 1));
    A2 = (A >> logg[bank]);
    A2 =
      ((A2 << bank) & ((1 << logg[bank]) - 1)) + (A2 >> (logg[bank] - bank));
    A = A1 ^ A2;
    A = ((A << bank) & ((1 << logg[bank]) - 1)) + (A >> (logg[bank] - bank));
    return (A);
  }

  // gindex computes a full hash of pc, ghist and phist
  int gindex (unsigned int pc, int bank, int hist, folded_history * ch_i)
  {
    int index;
    int M = (m[bank] > PHISTWIDTH) ? PHISTWIDTH : m[bank];
    index =
      pc ^ (pc >> (abs (logg[bank] - bank) + 1)) ^
      ch_i[bank].comp ^ F (hist, M, bank);
    return (index & ((1 << (logg[bank])) - 1));
  }

  //  tag computation
  uint16_t gtag (unsigned int pc, int bank, folded_history * ch0,
                 folded_history * ch1)
  {
    int tag = pc ^ ch0[bank].comp ^ (ch1[bank].comp << 1);
    return (tag & ((1 << TB[bank]) - 1));
  }

  // up-down saturating counter
  void ctrupdate (int8_t & ctr, bool taken, int nbits)
  {
    if (taken)
    {
      if (ctr < ((1 << (nbits - 1)) - 1))
        ctr++;
    }
    else
    {
      if (ctr > -(1 << (nbits - 1)))
        ctr--;
    }
  }

  bool getbim ()
  {
    return (btable[BI].pred > 0);
  }

  // update  the bimodal predictor: a hysteresis bit is shared among 4 prediction bits
  void baseupdate (bool Taken)
  {
    int inter = (btable[BI].pred << 1) + btable[BI >> HYSTSHIFT].hyst;
    if (Taken)
    {
      if (inter < 3)
        inter += 1;
    }
    else if (inter > 0)
      inter--;
    btable[BI].pred = inter >> 1;
    btable[BI >> HYSTSHIFT].hyst = (inter & 1);
  };


  //  TAGE PREDICTION: same code at fetch or retire time but the index and tags must recomputed

  void Tagepred ()
  {
    HitBank = 0;
    AltBank = 0;
    //Look for the bank with longest matching history
    for (int i = NHIST; i > 0; i--)
    {
      if (gtable[i][GI[i]].tag == GTAG[i])
	  {
	    HitBank = i;
	    break;
	  }
    }
    //Look for the alternate bank
    for (int i = HitBank - 1; i > 0; i--)
    {
      if (gtable[i][GI[i]].tag == GTAG[i])
	  {

	    AltBank = i;
	    break;
	  }
    }
    //computes the prediction and the alternate prediction
    if (HitBank > 0)
    {
      if (AltBank > 0)
        alttaken = (gtable[AltBank][GI[AltBank]].ctr >= 0);
      else
        alttaken = getbim ();
      LongestMatchPred = (gtable[HitBank][GI[HitBank]].ctr >= 0);
      //if the entry is recognized as a newly allocated entry and
      //USE_ALT_ON_NA is positive  use the alternate prediction
      if ((USE_ALT_ON_NA < 0)
          || (abs (2 * gtable[HitBank][GI[HitBank]].ctr + 1) > 1))
        tage_pred = LongestMatchPred;
      else
        tage_pred = alttaken;
    }
    else
    {
      alttaken = getbim ();
      tage_pred = alttaken;
      LongestMatchPred = alttaken;

    }
  }

  //compute the prediction

  bool predict_brcond (unsigned int pc, uint16_t brtype)
  {
    if (brtype & IS_BR_CONDITIONAL)
    {
      // computes the TAGE table addresses and the partial tags
      for (int i = 1; i <= NHIST; i++)
	  {
	    GI[i] = gindex (pc, i, Fetch_phist, Fetch_ch_i);
	    GTAG[i] = gtag (pc, i, Fetch_ch_t[0], Fetch_ch_t[1]);
	  }
      BI = pc & ((1 << LOGB) - 1);
      Tagepred ();

      pred_taken = tage_pred;
    }
    return pred_taken;
  }


  //  UPDATE  FETCH HISTORIES   + spec update of the loop predictor + Update of the IUM
  void FetchHistoryUpdate (uint32_t pc, uint16_t brtype, bool taken,
                           uint32_t target)
  {
    HistoryUpdate (pc, brtype, taken, target, Fetch_phist, Fetch_ptghist,
                   Fetch_ch_i, Fetch_ch_t[0], Fetch_ch_t[1]);
  }

  void HistoryUpdate (uint32_t pc, uint16_t brtype, bool taken,
                      uint32_t target, int &X, int &Y, folded_history * H,
                      folded_history * G, folded_history * J)
  {
    //special treatment for indirects and returnd: inhereited from the indirect branch predictor submission
    int maxt = (brtype & IS_BR_INDIRECT) ? 4 : 1;
    if (brtype & IS_BR_CALL)
      maxt = 5;


    int T = ((target ^ (target >> 3) ^ pc) << 1) + taken;
    int PATH = pc;
    for (int t = 0; t < maxt; t++)
    {
      bool TAKEN = (T & 1);
      T >>= 1;
      bool PATHBIT = (PATH & 1);
      PATH >>= 1;
      //update  history
      Y--;
      ghist[Y & (HISTBUFFERLENGTH - 1)] = TAKEN;
      X = (X << 1) + PATHBIT;
      X = (X & ((1 << PHISTWIDTH) - 1));
      //prepare next index and tag computations for user branchs
      for (int i = 1; i <= NHIST; i++)
	  {

	    H[i].update (ghist, Y);
	    G[i].update (ghist, Y);
	    J[i].update (ghist, Y);
	  }
    }

    //END UPDATE  HISTORIES
  }

  // PREDICTOR UPDATE

  void update_brcond (uint32_t pc, uint16_t brtype, bool taken,
                      uint32_t target)
  {

    if (brtype & IS_BR_CONDITIONAL)
    {
      // Recompute the prediction

      //Note that on a real hardware processor, one would avoid this recomputation to save an extra read port on the branch predictor. One  would  have to propagate informations needed for update with the instruction: HitBank, updated through the STATCOR or not, etc.

      //Recompute the prediction (the indices and tags) with the Retire history.
      for (int i = 1; i <= NHIST; i++)
	  {
	    GI[i] = gindex (pc, i, Retire_phist, Retire_ch_i);

	    GTAG[i] = gtag (pc, i, Retire_ch_t[0], Retire_ch_t[1]);
	  }
      BI = pc & ((1 << LOGB) - 1);
      Tagepred ();

      {
        bool ALLOC = ((tage_pred != taken) & (HitBank < NHIST));

        // try to allocate a  new entries only if TAGE prediction was wrong

        if (HitBank > 0)
	    {
          // Manage the selection between longest matching and alternate matching
          // for "pseudo"-newly allocated longest matching entry

	      bool PseudoNewAlloc =
            (abs (2 * gtable[HitBank][GI[HitBank]].ctr + 1) <= 1);
          // an entry is considered as newly allocated if its prediction counter is weak
	      if (PseudoNewAlloc)
          {
            if (LongestMatchPred == taken)
              ALLOC = false;
            // if it was delivering the correct prediction, no need to allocate a new entry
            //even if the overall prediction was false
            if (LongestMatchPred != alttaken)
              ctrupdate (USE_ALT_ON_NA, (alttaken == taken), 4);
          }
	    }

        //Allocate entries on mispredictions
        if (ALLOC)
	    {

          /* for such a huge predictor allocating  several entries is better*/
	      int T = 3;
	      for (int i = HitBank + 1; i <= NHIST; i += 1)
          {
            if (gtable[i][GI[i]].u == 0)
		    {
		      gtable[i][GI[i]].tag = GTAG[i];
		      gtable[i][GI[i]].ctr = (taken) ? 0 : -1;
		      gtable[i][GI[i]].u = 0;
		      TICK--;
		      if (TICK < 0)
                TICK = 0;
		      if (T == 0)
                break;
		      i += 1;
		      T--;
		    }
            else
              TICK++;
          }
	    }
        //manage the u  bit
        if ((TICK >= (1 << LOGTICK)))
	    {
	      TICK = 0;
          // reset the u bit
	      for (int i = 1; i <= NHIST; i++)
            for (int j = 0; j < (1 << logg[i]); j++)
              gtable[i][j].u >>= 1;

	    }

        //update the prediction

        if (HitBank > 0)
	    {
	      ctrupdate (gtable[HitBank][GI[HitBank]].ctr, taken, CWIDTH);
          // acts as a protection
	      if ((gtable[HitBank][GI[HitBank]].u == 0))
          {
            if (AltBank > 0)
              ctrupdate (gtable[AltBank][GI[AltBank]].ctr, taken,
                         CWIDTH);
            if (AltBank == 0)
              baseupdate (taken);
          }
	    }
        else
          baseupdate (taken);
        // update the u counter
        if (HitBank > 0)
          if (LongestMatchPred != alttaken)
	      {
            if (LongestMatchPred == taken)
            {
              if (gtable[HitBank][GI[HitBank]].u < 1)
                gtable[HitBank][GI[HitBank]].u++;

            }
	      }
        //END PREDICTOR UPDATE

      }

    }

    //  UPDATE RETIRE HISTORY
    HistoryUpdate (pc, brtype, taken, target, Retire_phist, Retire_ptghist,
                   Retire_ch_i, Retire_ch_t[0], Retire_ch_t[1]);
  }

};
#endif
