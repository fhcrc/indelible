// Please uncomment this line if compiling source code on Windows machine
// #define WINDOWS

/*
 * INDELible V1.03
 *  "A comprehensive and flexible simulator of molecular sequence evolution"
 *  Copyright (C) 2009 William Fletcher
 *
 *  If using this software please cite the following reference:
 *
 *  Fletcher, W. and Yang, Z. 2009.
 *      "INDELible: a flexible simulator of biological sequence evolution."
 *      Mol. Biol. and Evol. (submitted on Friday 13th March 2009).
 *
 *  If you need to contact me with any queries, questions, problems or (god-forbid) bugs
 *  then please go to http://abacus.gene.ucl.ac.uk/software/indelible/bug.php
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  Please visit http://www.gnu.org/licenses/ for a full copy of the license.
 */

///////////////////////////////////////
// debugging options

//#define printoverrideok  // this allows the use of the same filename for different blocks

//if(globalthisrep==6877+1 || globalthisrep==5775+1)

//#define inheritancycheck	// used to check if inheritancy of insertions is correct, needs
// a guide tree with a max. of 5 branches depth root to tip

int fixedsize = 3;                      // fixed size of indels if inheritancy check is used.

//#define rippartitions		// used to rip partitioned datasets in to separate datasets for each partition for checking purposes.

//#define splittree			// used to split 24 taxa dataset into 3 x 8 taxa datasets t o see if subtrees evolve properly with different models

#define printmaxindelsize	// just prints out the new maximum indel length as it changes
// for picking up bugs that steal lots of memory etc here

// then if there is an error - it can be repeated and tracked.
int globalthisrep;

#ifdef printmaxindelsize
int maxindelsize = 0;
#endif

//#define INDELLOG			// lists all indel events, where they occurred, changes etc.  for debugging indel routines.

// N.B. the pinv test crashes when used with partitions, and will be incorrect with partitions or non-stationary models.
int inspinvcount, corepinvcount;

//#define timewatch			// times each step of simulation, setup, evolving, then printing. keep disabled as prints
// lots of junk when simulations with a large number of replicates are performed.
//#define checksitesclass
/////////////////////////////////////////

//#include <cmath>



// needed for time monitoring

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * // example use of time monitoring
 *
 * clock_t start, finish;
 * double  duration;
 *
 * // Measure the duration of an event.
 * printf( "Time to do %ld empty loops is ", i );
 * start = clock();
 * while( i-- )
 *    ;
 * finish = clock();
 * duration = (double)(finish - start) / CLOCKS_PER_SEC;
 * printf( "%2.1f seconds\n", duration );
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <string>
#include <vector>
#include <math.h>
#include <sstream>
#include <algorithm>
#include <map>
#include <assert.h>

#ifdef WINDOWS
#include "windows.h"
#else
#include "limits.h"
#endif

// #include "control.cpp"
#include "models.h"
#include "control.h"
#include "randoms.h"
#include "paml.h"

using namespace std;

enum class Event { error=-1, core_sub=0, ins_sub=1, core_ins=2, ins_ins=3, core_del=4, ins_del=5 };

// used for timing
clock_t start, finish, startofprinting, startofevolution, endofprinting, endofblock, startofblock;
double  duration;


// this annoying warning always pops up in Visual C++ V6.0, it is a confirmed bug in the compiler.
#pragma warning(disable:4786)

// this one too - it warns of exceptions for using <string> class etc
#pragma warning(disable:4530)

ofstream *results  = new ofstream;      // pointers to output streams for file output.
ofstream *results2 = new ofstream;      // pointers to output streams for file output.
ofstream *results3 = new ofstream;      // pointers to output streams for file output.

int oldindellength;

#ifdef INDELLOG
ofstream indellog("indellog.txt");
#endif
int         lasttreeseed = 0;
vector<int> mymods;
//////////////////////////////////////
vector<vector<int> > mycodes;

// pointer to a function used to produce P(t) from Q for method 1.
// pointer needed to choose between functions for reversible and non-reversible substitution models
// i.e. UNREST for NUCLEOTIDE uses a different function than all other reversible models.

vector<vector<double> > (*Pt)(vector<vector<double> > Qvec, vector<double>& basefreqs, double t);


class site
{
    // every inserted site is represented as an instance of "site"
public:
    int    base;                // state at that site
    double rate;                // relative substitution rate used for e.g. gamma models
    int    siteclass;           // site class for site, used for codon site-models or discrete gamma
    double subrate;             // substitution rate for site that depends on state held in base (for method 2)
    double timeleft;            // if an insertion occurs at distance x into branch of length y then timeleft = y - x
    // used by method 1 to distinguish from sites that were inserted on different branches

    site(int b, double ra, int si, double r, double t)
    {
	base      = b;
	rate      = ra;
	siteclass = si;
	subrate   = r;
	timeleft  = t;
    }
};

class insert
{
    // every site k in the core sequence has an instance of "insert" i_k say
    // that mirrors it to hold information about inserted sequences.

public:
    vector<site> insertvec;     // all sites that are inserted after core position k are stored here
    int          length;        // number of sites in insertvec
    double       subrate;       // total substitution rate for all sites in insertvec

    insert(vector<site> mytempvec, int mylength, double mysubrate)
    {
	subrate   = mysubrate;
	length    = mylength;
	insertvec = mytempvec;
    }
};


class SUMS
{
    // an instance of SUMS is used in every simulation.
    // for method 1 the vector<double> are never used.
    // for method 2 they are used to choose site where substitutions occur
    // both methods use the vector<int> for choosing where indels occur.
    /*
     * maintain maps of relative probabilities of insertion at each
     * position in a sequence.  
     *
     * In order to select positions according to a multinomial
     * probability distribution, we maintain a cumulative sum of
     * probabilities associated with each site.  
     * Given the cumulative sum of probabilities associated with each site (cumsum) and a 
     * random number (r), identify the highest position at which cumsum < r.
     * That will be the randomly chosen position at which in the insertion should occur.
     *
     * To make this process fast when a sequencs has millions of sites,  
     * maintain a hierarchical summation.  Each level in the hierarchy
     * covers 10 times the number of sites at the next lowest level.
     *
     */

public:
    vector< vector<int> > IIlevels;	//insertion rates for each particular site in inserted sites.
    vector< vector<int> > CIlevels;	//insertion rates for each particular site in core sequence.
    vector< vector<double> > CSlevels;  //substitution rates for each particular site in core sequence.
    vector< vector<double> > ISlevels;  //substitution rates for each particular site in inserted sites.

    SUMS() : IIlevels(11), CIlevels(11), CSlevels(11), ISlevels(11) {
	
    }

    void clear()
    {
	// simply clears all vectors above in an instnce of SUMS.
	// used when the model changes from branch to branch
	for (auto& level : IIlevels)	level.clear();
	for (auto& level : CIlevels)	level.clear();
	for (auto& level : CSlevels)	level.clear();
	for (auto& level : ISlevels)	level.clear();
    }
};



class RATES
{
    // an instance of RATES is used in each simulation to track and store changes in
    // substitution, insertion and deletion rates for inserted sites and core sites.
    // the total rate is used to generate exponential waiting times
    // the other rates are used to choose between different types of event, and
    // to choose whether the event takes place in inserted or "core" sites

public:
    double totalrate;           // total event rate for all sites and all events
    double corerate;            // total event rate for all core sites and all events
    double coresubrate;         // total substitution rate for all core sites
    double coreinsertrate;      // total insertion rate for all core sites
    double coredeleterate;      // total deletion rate for all core sites
    double insrate;             // total event rate for all inserted sites and all events
    double inssubrate;          // total substitution rate for all inserted sites
    double insinsertrate;       // total insertion rate for all inserted sites
    double insdeleterate;       // total deletion rate for all inserted sites

    // at the root inslength = 0, and the others all equal the root length specified.
    int rootlength;             //
    int corelength;             // this starts as rootlength but decreases with every core site deleted
    int inslength;              //
    int totallength;            //
    int partitionlength;        // ... fill in later.

    RATES()
    {
	totalrate      = 0;
	corerate       = 0;
	coresubrate    = 0;
	coreinsertrate = 0;
	coredeleterate = 0;

	insrate       = 0;
	inssubrate    = 0;
	insinsertrate = 0;
	insdeleterate = 0;

	rootlength      = 0;
	corelength      = 0;
	inslength       = 0;
	totallength     = 0;
	partitionlength = 0;
    }
};

////////////////////////////////////////
vector<model>        mymodels;
vector<int>          siteclasses;
vector<vector<int> > siteclassmodels;


vector<vector<double> > totalratevecs;

bool           noinsertionsyet;
vector<double> ratevec;
vector<int>    siteclassvec;

int  mypos;
bool weareon;

#ifdef checksitesclass
ofstream gout("g.txt");
#endif

model          *m;
siteclass      *s;
branchclass    *b;
partitionclass *p;
evolve         *e;
Tree           *treeC;

void changezipfrandoms()
{
    Zq1    = (*m).q1;
    Zq2    = (*m).q2;
    ZHx0   = (*m).Hx0;
    Zs     = (*m).s;
    ZHimax = (*m).Himax;
}


vector<vector<int> >    insPOS;                 // for use inside TinsPOS
vector<vector<int> >    sequencesINT;           // for use inside TsequencesINT
vector<vector<insert> > insINT;                 // for use inside TinsINT

vector<vector<vector<int> > > TinsPOS;          // holds *all* inserted site position information in a block.
vector<vector<vector<int> > > TsequencesINT;    // holds *all* core sequence information in a block.
// outside vector is for each partition
// middle vector is positions for each node on the tree for the partition
// internal vector is actual sequence information.

vector<vector<vector<insert> > > TinsINT;               // holds *all* core insertion information in a block.
// outside vector is for each partition
// middle vector is positions for each node on the tree for the partition
// internal vector represents insertion information positions of insertions relative to core.
// internal vector is an actual insertion


//Jump Vector for use when choosing new bases
vector<vector<double> > JVec;

// to store codon site class frequencies
vector<double> Csitefreqs;

// to store codon model information from control file
vector<vector<double> > codoninfo;

// vector of jump matrices for each site class
vector<vector<vector<double> > > matrixJVsites;

// vector of different matrixJVsites for each branch category
vector<vector<vector<vector<double> > > > matrixJVbranches;

vector<double> GDfreq;          // used in gamma rate variation
vector<double> GDrates;         // used in gamma rate variation

double substitutioncount;       // counts the number of substitution events in a block
double insertioncount;          // counts the number of insertion events in a block
double insertiontotlength;      // tracks cumulative insertion length in a block
double deletioncount;           // counts the number of deletion events in a block
double deletiontotlength;       // tracks cumulative deletion length in a block

string originaltree;
string VersionNumber = "1.03";


int linecount      = 0;         // used to identify the line in the control file where an error occurs
int deletionlength = 0;         // tracks number of bases needing to be deleted when switching between core sequence and insertions

//int breakonerror	= 1;	// should be left as 1 unless debugging

vector<double> defaultAAbasefreqs, AAbasefreqs;   // Protein base frequencies vectors
vector<double> defaultDNAbasefreqs, DNAbasefreqs; // DNA base frequencies vectors
vector<double> defaultCbasefreqs, Cbasefreqs;     // DNA base frequencies vectors

// Valid characters for input and output, used for translation to and from integer representation
string DU[4] = { "T", "C", "A", "G" };
string DL[4] = { "t", "c", "a", "g" };
string AU[20] = { "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V" };
string AL[20] = { "a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v" };


string CDU[64] =
{
    "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
    "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
    "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
    "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"
};

string CDL[64] =
{
    "ttt", "ttc", "tta", "ttg", "tct", "tcc", "tca", "tcg", "tat", "tac", "taa", "tag", "tgt", "tgc", "tga", "tgg",
    "ctt", "ctc", "cta", "ctg", "cct", "ccc", "cca", "ccg", "cat", "cac", "caa", "cag", "cgt", "cgc", "cga", "cgg",
    "att", "atc", "ata", "atg", "act", "acc", "aca", "acg", "aat", "aac", "aaa", "aag", "agt", "agc", "aga", "agg",
    "gtt", "gtc", "gta", "gtg", "gct", "gcc", "gca", "gcg", "gat", "gac", "gaa", "gag", "ggt", "ggc", "gga", "ggg"
};

vector<string> myDUletters(DU, DU + 4);
vector<string> myDLletters(DL, DL + 4);
vector<string> myAUletters(AU, AU + 20);
vector<string> myALletters(AL, AL + 20);
vector<string> myCDUletters(CDU, CDU + 64);
vector<string> myCDLletters(CDL, CDL + 64);


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A lot of the following variables are now defunct and not used - but some are still in use - I just haven't had the time to check which.  //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<string> commandblockvec;                         // used to store a command block to append output files (e.g. MrBayes or PAUP block)
vector<string> block;                                   // this will be used to store information about one particular block
vector<string> taxanames;                               // this stores the real taxon names and internal nodel labels
vector<string> taxaspacenames;                          // formatted taxanames
vector<string> partonetaxanames;                        // this stores the real taxon names and internal nodel labels for the first partition
// - used for comparison to make sure ordering of taxon names is same across partitions
// - and to make sure that the total number of taxa/taxa names is same across partitions
vector<vector<string> > commandblocks;  // this will store the command blocks in total, as read from the control file

int NgamCAT;                            // Number of categories for the discrete gamma rate variation
int medianORmean;                       // 1 = use medians, 0 = use means, to represent categories in discrete gamma rate variation
int rootlength;                         // length of root sequence to be created or when given.
int reps;                               // number of repetitions
int ntaxa;                              // number of taxa
int mymodel;                            // substitution model number
int startnodelabel;                     // used in manipulation of guide tree

int dashlength;                         // used to calculate length of true alignment;
int totallength;

int partitionsdone;             // offset for block titles etc
int partitionnumber;            // number of partitions for each repetition
int partitionstogo;             // number of partitions left to do.
int geneticCode;                // used for translation between amino acids and codons.

// which bool values reset, which ones don;t?

bool isitsites;
bool isitbranches;


bool insertionBFfrommodel;              // If true insertion base frequencies is same as model base frequencies, otherwise they come from ((*m).insertfreqs).
bool ECMbasefreqOVERRIDE;               // Used for +F version of ECM codon models.  set to true, model base frequencies come from AAbasefreqs
bool AAbasefreqOVERRIDE;                // Used for +F protein model.  set to true, model base frequencies come from AAbasefreqs
bool modelbasefreqs;                    // set to true means that root base frequencies come from the model, if root sequence is created.
bool rootTOmodelOVERRIDE;               // if true, then model stationary base frequencies are determined from the given root sequence.
bool rootseqgiven;                      // true if root sequence given, if false then rootTOmodelOVERRIDE above has no effect.
bool ratechange;                        // whether changes to the site rate variation commands occurs in a block
bool codonrates;                        // whether to used codon rates
bool paupblock;                         // whether to print PAUP block or MrBayes block
bool warningprint;                      // whether two taxa names are the same in the guide tree & a warning has been printed
bool printonlyroot;                     // whether to print ancestral sequences.
bool sites;                             // for defining sites and sites-branches models for codons
bool branches;                          // for defining branches and branches-sites models for codons

double pinv;                            // proportion of invariant sites, 0 for none, 1 for all.
double alpha;                           // alpha value for gamma rate distribution
//	double (*m).instodel;			// ratio of insertions to deletions e.g. 0.48 will give 48 insertions to every 52 deletions on average.
//	double (*m).indeltosub;			// ratio of indels to substitutions, e.g. 0.01 gives 1 insertion/deletion event for every 99 substitutions on average
//	double georate;				// rate parameter for the geometric distribution controlling indel lengths
double currenttreelength;               // amount of treelength in a given repetition that has been completed
double treelength;                      // total treelength to be evolved in a given repetition, used with currenttreelength to calculate percentage complete.

string filenamestub;                    // this is a unique file name for every block appended with the repetition number for the output files
string paupblockname;                   // filename for paupblock.
string nonamestree;                     // all used in tree manipulation.
string origtreewithnodes;
string nonamestreewithnodes;
string workingtree;
string guidetree;

//	vector<double> ratevec;				// Used to store information on continuous/discrete gamme rate variation, & proportion of invariant sites
vector<double> usermodel;                               // Used for user-defined Amino Acid substitution model
//	vector<double> rootbasefreqs;		// base frequencies used to create root sequence if not using model base frequencies
//	vector<double> INSERTIONbasefreqs;	// base frequencies used to choose inserted bases if not using model base frequencies
//	vector<double> myrates;				// these are the rates used when estimating exponential waiting times

//	vector<string> rootseqCHAR;			// alphabetical representation of root sequence
//	vector<int> rootseqINT;				// numerical representation of root sequence
vector<int> partitionposition;                  // partitionposition is used to correctly print sequence information when using partitions with different guide trees

double codonrate[3] = { 1, 1, 1 };


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double smalldiff = 0.0000000001;  //for comparing floating point values



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void printtitle(ofstream& LOG, string& title)
{
    // prints out block specific title for print settings function
    int mysize = title.size();

    LOG << endl;
    LOG << "+";
    for (int tg = -2; tg < mysize; tg++) {
	LOG << "-";
    }
    LOG << "+" << endl;
    LOG << "| " << title << " |" << endl;
    LOG << "+";
    for (int tg2 = -2; tg2 < mysize; tg2++) {
	LOG << "-";
    }
    LOG << "+" << endl;
    LOG << endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void format(string& tss, string& css, string& ass, string& gss, double t, double c, double a, double g)
{
/*
 * string format(double Value, int nPrecision)
 * {
 *      stringstream s;
 *      s.width(nPrecision);
 *      s <<  Value;
 *       return s.str();
 * }
 */


    stringstream ts, cs, as, gs;

//	int maxsize;

    if (t == 1.0) {
	tss = "1.0";
    } else if (t == 0.0) {
	tss = "0.0";
    } else {
	ts << t;
	ts.width(6);
	tss = ts.str();
    }                                                                                             //cout<<tss<<endl;}
    if (c == 1.0) {
	css = "1.0";
    } else if (c == 0.0) {
	css = "0.0";
    } else {
	cs << c;
	cs.width(6);
	css = cs.str();
    }                                                                                             //cout<<css<<endl;}
    if (a == 1.0) {
	ass = "1.0";
    } else if (a == 0.0) {
	ass = "0.0";
    } else {
	as << a;
	as.width(6);
	ass = as.str();
    }                                                                                             //cout<<ass<<endl;}
    if (g == 1.0) {
	gss = "1.0";
    } else if (g == 0.0) {
	gss = "0.0";
    } else {
	gs << g;
	gs.width(6);
	gss = gs.str();
    }                                                                                             //cout<<gss<<endl;}


/*
 *      stringstream ts,cs,as,gs;
 *      int maxsize;
 *
 *      if(t==1.0) tss="1.0"; else if(t==0.0) tss="0.0"; else {ts<<t; tss=ts.str(); maxsize=tss.size();                         }//cout<<tss<<endl;}
 *      if(c==1.0) css="1.0"; else if(c==0.0) css="0.0"; else {cs<<c; css=cs.str(); if(maxsize<css.size()) maxsize=css.size();  }//cout<<css<<endl;}
 *      if(a==1.0) ass="1.0"; else if(a==0.0) ass="0.0"; else {as<<a; ass=as.str(); if(maxsize<ass.size()) maxsize=ass.size();  }//cout<<ass<<endl;}
 *      if(g==1.0) gss="1.0"; else if(g==0.0) gss="0.0"; else {gs<<g; gss=gs.str(); if(maxsize<gss.size()) maxsize=gss.size();  }//cout<<gss<<endl;}
 *
 *      for(int jh1=0; jh1<maxsize+1; jh1++)
 *      {
 *              if(jh1>tss.size()) tss+="0";
 *              if(jh1>css.size()) css+="0";
 *              if(jh1>ass.size()) ass+="0";
 *              if(jh1>gss.size()) gss+="0";
 *      }
 */
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void chooseNEWbase(int& newbase, vector<double> JVecRow)
{
    // choose new base according to jump matrix
    myrand = mtrand1();

    for (int i = 0; i < JVecRow.size(); i++) {
	if (myrand < JVecRow.at(i)) {
	    newbase = i;
	    break;
	}
    }
}


//////////////////////////////////////////////////////////////////////////


int returnNewSiteClass(/*vector<int> &siteclasscount, */ vector<double> siteprops, int po)
{
    if (siteprops.empty()) {
	cout << endl << endl << "ERROR 1 in returnNEWsiteClass" << endl << endl;
	return 0;
    }
    double rand      = mtrand1(); //cout<<rate
    int    siteclass = -1;
    for (int i1 = 0; i1 < siteprops.size(); i1++) {
	if (rand < siteprops.at(i1)) {
	    siteclass = i1;
	    break;
	}
    }
    //cout<<" "<<rate<<endl;
    if (siteclass == -1) {
	cout << endl << endl << "ERROR 2 in returnNEWsiteClass" << endl << endl;
    }

    //(siteclasscount.at(siteclass))++;

    //if(po<50) cout<<po<<"  "<<siteclass<<endl;
    return siteclass;
}


///////////////
double returnNewSiteRate(int pos)
{
    double rate;

    // this function sets the site rates for evolution of the core sequence.
    // inserted bases will inherit the rate of the base  where the insertion occurred.

    if (model_type == model::Type::codon) {
	rate = mtrand1();       //cout<<rate;
	for (int i1 = 0; i1 < ((*m).cumfreqs).size(); i1++) {
	    if (rate < ((*m).cumfreqs).at(i1)) {
		rate = i1;
		break;
	    }
	}
	//cout<<" "<<rate<<endl;

	return rate;
    }

    double multiplier = 1;
    if (model_type == model::Type::nucleotide) {
	int codonpos = pos % 3;
	multiplier = ((*m).codonrates)[codonpos];
    }

    if (codonrates) {
    } else if (((*m).alpha == 0) && ((*m).pinv == 0)) {  // no gamma, no codon rates, no proportion invariant
	// constant rates across whole sequence
	rate = multiplier;
    } else {     // proportion invariant with either no gamma, discrete gamma, or continuous gamma
	if ((*m).alpha == 0) {
	    //this is proportion invariant
	    if (mtrand1() < (*m).pinv) {
		rate = 0;
	    } else {
		rate = multiplier / (1 - ((*m).pinv));
	    }
	} else if (((*m).alpha > 0) && ((*m).ngamcat == 0)) {
	    //Continuous Gamma + Proportion invariant
	    if (mtrand1() < (*m).pinv) {
		rate = 0;
	    } else {
		rate = rndgamma(((*m).alpha) / (((*m).alpha) * (1 - ((*m).pinv))));
	    }
	} else if (((*m).alpha > 0) && ((*m).ngamcat > 0)) {
	    //Discrete Gamma + Proportion invariant

	    if (mtrand1() > (*m).pinv) {
		rate = GDrates.at((int)(mtrand1() * ((*m).ngamcat))) / (1 - ((*m).pinv));
	    } else {
		rate = 0;
	    }
	} else {
	    (*LOG) << "SetSite Rates error 1 - alpha or ngamcat is less than zero" << endl;
	}
    }

    return rate;
}


double returnsitesubrate(double rate, int siteclass, int currentbase)
{
    if (model_type == model::Type::codon) {
	return (((*m).myratesvec).at(siteclass)).at(currentbase);
    } else {
	return rate * (((*m).myrates).at(currentbase));
    }
}


int returnsiteclass()
{
    int    i1;
    double rate = mtrand1();

    for (i1 = 0; i1 < ((*m).cumfreqs).size(); i1++) {
	if (rate <= ((*m).cumfreqs).at(i1)) {
	    break; //{ ratevec.push_back(i1);   break;}
	}
    }

    return i1;
}


double returnsiterate(int& rateclass, bool core)
{
    if ((*m).continuousgamma) {
	double rand = mtrand1();

	rateclass = 0;

	if (rand < (*m).pinv) {
	    if (core) {
		corepinvcount++;
	    } else {
		inspinvcount++;
	    }
	    return 0;
	} else {
	    return rndgamma(((*m).alpha)) / (((*m).alpha) * (1 - ((*m).pinv)));
	}
    } else {
	rateclass = returnsiteclass();
	double myrate = ((*m).Rrates).at(rateclass);

	if (myrate == 0) {
	    if (core) {
		corepinvcount++;
	    } else {
		inspinvcount++;
	    }
	}
	return myrate;
    }
}


void SetSiteRates2(RATES& rates, int repnumber, int mysize)
{
    siteclassvec.clear();
    ratevec.clear();
    siteclassvec.reserve(rates.rootlength);
    ratevec.reserve(rates.rootlength);

    if (model_type == model::Type::codon) {
	if (fixtrueproportions) {
	    // this option was just for me to use in one of my papers.
	    siteclassvec.push_back(0);
	    for (int x = 0; x < 10; x++) {
		for (int y = 0; y < 20; y++) {
		    siteclassvec.push_back(x);
		}
	    }

	    for (int t1 = 0; t1 < rates.rootlength; t1++) {
		//	siteclassvec.push_back(returnsiteclass2());
		ratevec.push_back(0);
	    }
	    cout << rates.rootlength << "  " << ratevec.size() << "  " << siteclassvec.size() << endl;
	} else {
	    for (int t1 = 0; t1 < rates.rootlength; t1++) {
		siteclassvec.push_back(returnsiteclass());
		ratevec.push_back(0);
	    }
	}
    } else {
	for (int t1 = 0; t1 < rates.rootlength; t1++) {
	    int myclass;
	    ratevec.push_back(returnsiterate(myclass, true));
	    siteclassvec.push_back(myclass);
	}
    }
}


//////////////////////////////////////////

vector<double> makecumfreqs(vector<double> basefreqs)
{
    vector<double> boundaries;
    double         s = 0;
    for (int yh = 0; yh < basefreqs.size(); yh++) {
	s += basefreqs.at(yh); /*cout<<"  W "<<basefreqs.at(yh)<<"  A  "<<s<<endl;*/ boundaries.push_back(s);
    }

    if (boundaries.back() - 1 > 0.00000001) {
	cout << "make seq error in sum to 1 on boundaries " << boundaries.back() << " " << boundaries.back() + 1 << " " << boundaries.back() - 1 << endl;
    }

    return boundaries;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void makeseq2(int length, vector<site>& myseq, double& sdiff, double timeleft)
{
    // this assigns to "myseq" a sequence of length "length" subject to the base frequencies in "basefreqs" in integer form

    vector<double> boundaries2;

    boundaries2 = makecumfreqs((*m).insertfreqs);

    for (int yh2 = 0; yh2 < length; yh2++) {
	myrand = mtrand1();
	bool wearenotdone = true;

	for (int yh3 = 0; yh3 < boundaries2.size(); yh3++) {
	    if (myrand <= boundaries2.at(yh3)) {
		wearenotdone = false;

		double rate      = -1;
		int    siteclass = -1;

		if (model_type == model::Type::codon) {
		    siteclass = returnsiteclass();
		} else {
		    rate = returnsiterate(siteclass, false);
		}

		double sitesubrate;
		if (oldmethod) {
		    sitesubrate = returnsitesubrate(rate, siteclass, yh3);
		    sdiff      += sitesubrate;
		}

		myseq.push_back(site(yh3, rate, siteclass, sitesubrate, timeleft));

		break;
	    }
	}
	if (wearenotdone) {
	    cout << "NO BASE PICKED in makeseq2, last boundaries was "
		 << boundaries2.back()
		 << " myrand was " << myrand
		 << " sdiff is " << sdiff << endl;
	}
    }
}


// Create a random sequence according to the base frequency distribution.
// Operates on both model base frequencies or branch base frequencies, whatever that means... -csw

void makeseq(int length2, vector<int>& myseq)
{
    // this assigns to "myseq" a sequence of length "length" subject to the base frequencies in "basefreqs" in integer form

    int length = length2 - 1;

    vector<double> boundaries2;
    if (isitbranches) {
	boundaries2 = makecumfreqs((*b).rootbasefreqs);
    } else {
	boundaries2 = makecumfreqs((*m).rootbasefreqs);
    }

    myseq.push_back(-1);
    for (int yh2 = 0; yh2 < length; yh2++) {
	myrand = mtrand1();

	for (int yh3 = 0; yh3 < boundaries2.size(); yh3++) {
	    if (myrand <= boundaries2.at(yh3)) {
		myseq.push_back(yh3);
		break;
	    }
	}
    }
}


void fromprintseq(vector<string>& seq, vector<int>& myseq)
{
    // translates alphabetical protein/dna sequence into integer representation
    myseq.clear();
    myseq.reserve(seq.size());

    vector<string> myletters, myletters2;
    int            mysize;

    if (model_type == model::Type::nucleotide) {
	mysize     = 4;
	myletters  = myDUletters;
	myletters2 = myDLletters;
    } else if (model_type == model::Type::aminoacid) {
	mysize     = 20;
	myletters  = myAUletters;
	myletters2 = myALletters;
    } else {
	cout << "error in fromprintseq" << endl;
    }


    for (int gv1 = 0; gv1 < seq.size(); gv1++) {
	string bb        = seq.at(gv1);
	int    minimatch = -1;

	for (int gv3 = 0; gv3 < myletters.size(); gv3++) {
	    if (bb == myletters.at(gv3)) {
		minimatch = gv3;
		break;
	    }
	}

	if (minimatch == -1) {
	    for (int gv4 = 0; gv4 < myletters2.size(); gv4++) {
		if (bb == myletters2.at(gv4)) {
		    minimatch = gv4;
		    break;
		}
	    }
	}

	if (minimatch == -1) {
	    cout << "ERROR in fromprintseq" << endl << " ";
	} else {
	    myseq.push_back(minimatch);
	}
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
string nowhitespace(string guidetree)
{
    //takes white space out of tree and checks number of right and left parentheses
    char   c;
    int    bracketright = 0, bracketleft = 0;
    string guidetree2S;

    for (int ghy11 = 0; ghy11 < guidetree.size(); ghy11++) {
	// get rid of white space in tree
	c = guidetree[ghy11];
	if (c != ' ') {
	    guidetree2S += c;
	}

	if (c == ')') {
	    bracketright++;
	}
	if (c == '(') {
	    bracketleft++;
	}
    }

    if (bracketleft != bracketright) {
	cout << "ERROR ON GUIDE TREE - mismatch on number of brackets" << endl << " ";
    }

    return guidetree2S;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

int gettaxanames(int partitionnumber, string guidetree, vector<string>& taxanames)
{
    // called from the settreeup function

    // gets taxon names from guide tree and puts into vector
    // also truncates filenames to 10 characters if output in Phylip truncated is chosen
    // final vector entry is tree with the label integer n replacing taxon name placed in nth position of taxanames vector

    char   c3 = 'q', c2 = 'q', c1 = 'q';
    string taxaname, tempstring, guidetree3;
    int    taxanameon = 0, maxnamesize = 0;

    taxanames.clear();
    taxanames.push_back("ROOT");

    for (int jh1 = 0; jh1 < guidetree.size(); jh1++) { // get taxanames and replace them in tree;
	c3 = c2;
	c2 = c1;
	c1 = guidetree[jh1];
	if (taxanameon == 0) {
	    if ((c2 == '(') || (c2 == ',')) {
		if (c1 != '(') {
		    taxanameon = 1;
		    taxaname   = "";
		}
	    }
	}
	if ((taxanameon == 1) && (c1 == ':')) {
	    taxanameon = 2;
	}
	if (taxanameon == 1) {
	    taxaname += c1;
	}
	if (taxanameon == 2) {
	    taxanameon = 0;
	    stringstream fd;

	    int minimatch = -1;

	    if (partitionnumber == 0) {
		fd << taxanames.size();
		string fddd = fd.str();
		guidetree3 += fddd;
	    } else {
		for (int gb2 = 0; gb2 < partonetaxanames.size(); gb2++) {
		    if (taxaname == partonetaxanames.at(gb2)) {
			minimatch = gb2;
			break;
		    }
		}

		if (minimatch == -1) {
		    stringstream f;
		    f << partitionnumber + 1;
		    string r = f.str();
		    controlerrorprint2("[PARTITIONS]", (*p).name, "", "Error with taxon " + taxaname + "\nIt was found in tree for partition " + r + " but not found in tree for partition 1", "");
		    (*LOG) << guidetree << endl;

		    //			(*LOG)<<endl<<endl<<" ERROR: Taxon "<<taxaname<<" in guide tree for partition "<<partitionnumber+1<<" not found in guide tree for partition 1"<<endl<<endl;
		    //			cout<<endl<<endl<<" ERROR: Taxon "<<taxaname<<" in guide tree for partition "<<partitionnumber+1<<" not found in guide tree for partition 1"<<endl<<endl;
		    return -1;
		} else {
		    fd << minimatch;
		    string fddd = fd.str();
		    guidetree3 += fddd;
		}
	    }

	    // truncate filename if necessary
	    if (phylipnametruncate && (taxaname.size() > 10)) {
		string tempstring2 = taxaname;
		taxaname = "";
		for (int fv = 0; fv < 10; fv++) {
		    taxaname += tempstring2[fv];
		}
	    }
	    if (phylipnametruncate && (taxaname.size() < 10)) {
		for (int fv = taxaname.size(); fv < 10; fv++) {
		    taxaname += " ";
		}
	    }

	    taxanames.push_back(taxaname);
	    if (taxaname.size() > maxnamesize) {
		maxnamesize = taxaname.size();
	    }
	}
	if (taxanameon == 0) {
	    guidetree3 += c1;
	}
    }

    taxanames.push_back(guidetree3);

    int partnowsize = taxanames.size() - 1;

    // check whether two taxa have the same name
    for (int gn2 = 0; gn2 < partnowsize; gn2++) {
	for (int gn3 = 0; gn3 < partnowsize; gn3++) {
	    if (gn3 != gn2) {
		if (taxanames.at(gn2) == taxanames.at(gn3)) {
		    stringstream f;
		    f << partitionnumber + 1;
		    string r = f.str();

		    if (phylipnametruncate) {
			controlerrorprint2("[TREE]", (*treeC).name, "", "Two taxa in partition " + r + " tree have the same name after truncation for PHYLIP format: " + (taxanames.at(gn2)), "");

			(*LOG) << guidetree << endl;
			//(*LOG)<<endl<<endl<<" ERROR: two taxa in the partition "<<partitionnumber+1<<" guidetree have the same name after truncation for PHYLIP format: "<<taxanames.at(gn2)<<endl<<endl;

			return -1;
		    } else {
			controlerrorprint2("[TREE]", (*treeC).name, "", "Two taxa in partition " + r + " tree have the same name: " + (taxanames.at(gn2)), "");
			(*LOG) << guidetree << endl;
			//(*LOG)<<endl<<endl<<" ERROR: two taxa in the partition "<<partitionnumber+1<<" guidetree have the same name: "<<taxanames.at(gn2)<<endl<<endl;
			return -1;
		    }
		}
	    }
	}
    }

    if (partitionnumber == 0) {
	partonetaxanames = taxanames;
	taxaspacenames   = taxanames;

	// add white space to taxa names so that output lines up
	for (int gn1 = 0; gn1 < partnowsize; gn1++) {
	    string tempname = taxanames.at(gn1);

	    if (!phylipnametruncate) {
		for (int fg = tempname.size(); fg < maxnamesize + 5; fg++) {
		    tempname += " ";
		}
	    }

	    taxaspacenames.at(gn1) = tempname;
	}
    } else {
	// if performing a partitioned analysis check for consistency of taxa names and numbers
	int partonesize = partonetaxanames.size() - 1;

	if (partonesize != partnowsize) {
	    stringstream f;
	    f << partitionnumber + 1;
	    string       r = f.str();
	    stringstream g;
	    g << partonesize - 1;
	    string       s = g.str();
	    stringstream h;
	    h << partnowsize - 1;
	    string t = h.str();

	    controlerrorprint2("[TREE]", (*treeC).name, "", "Guide tree for partition 1 has " + s + " taxa\nGuide tree for partition " + r + " has " + t + " taxa.", "");

//			(*LOG)<<endl<<endl<<" ERROR: Guide tree for partition 1 has "<<partonesize-1<<" taxa\n";
//			(*LOG)<<            "        Guide tree for partition "<<partitionnumber+1<<" has "<<partnowsize-1<<" taxa"<<endl<<endl;
	    return -1;
	}
    }

    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
string addnodestostring(string mytree, int& nodelabel)
{
    // places internal node labels on tree
    char c1 = 'q', c = 'q';

    string workingtree;

    for (int gh11 = 0; gh11 < mytree.size(); gh11++) {
	c1 = c;
	c  = mytree[gh11];
	if ((c1 == ')') && (c == ':')) {
	    nodelabel++;
	    workingtree += 'N';
	    stringstream g;
	    g << nodelabel;
	    string h = g.str();
	    workingtree += h;
	}

	if ((c != ' ') && (c != ';')) {
	    workingtree += c;
	}
    }

    return workingtree;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getlabelstring(int& label, double& length, string& mystring)
{
    //gets tip label name and final length from terminal branch - used in evolvebranch
    string labelstring, lengthstring;
    int    num = 0;
    char   c   = mystring[num];

    while (c != ':') {
	labelstring += c;
	num++;
	c = mystring[num];
    }

    c = mystring[num];
    num++;

    while (num < mystring.size()) {
	c = mystring[num];
	num++;
	lengthstring += c;
    }

    length = atof(lengthstring.c_str());
    label  = atoi(labelstring.c_str());
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void getmynewbits(int& label, double& length, string& mystring, vector<string>& mynewbits)
{
    //returns length to next node and strings of each daughter branch from that node - used in evolvebranch
    string labelstring, lengthstring;
    int    bracketlevel = 1;
    string mynewstring;

    for (int pi1 = 1; pi1 < mystring.size(); pi1++) {
	char c = mystring[pi1];

	if (bracketlevel == -3) {
	    lengthstring += c;
	}
	if ((bracketlevel == -2) && (c == ':')) {
	    bracketlevel = -3;
	}
	if (bracketlevel == -2) {
	    labelstring += c;
	}
	if ((bracketlevel == -1) && (c == 'N')) {
	    bracketlevel = -2;
	}

	if ((bracketlevel == 1) && (c == ',')) {
	    mynewbits.push_back(mynewstring);
	    mynewstring = "";
	} else {
	    if (bracketlevel > 0) {
		mynewstring += c;
	    }
	}
	if (c == '(') {
	    bracketlevel++;
	}
	if (c == ')') {
	    bracketlevel--;
	}

	if (bracketlevel == 0) {
	    string mynewstring2;
	    for (int yhb = 0; yhb < mynewstring.size() - 1; yhb++) {
		mynewstring2 += mynewstring[yhb];
	    }

	    mynewbits.push_back(mynewstring2);
	    mynewstring  = "";
	    bracketlevel = -1;
	}
    }

    length = atof(lengthstring.c_str());
    label  = atoi(labelstring.c_str());
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<string> MAKEmyCAUletters(int partition, int taxon)
{
    //int gencode=((*p).geneticcodes).at(partition);
    int gencode = (mycodes.at(partition)).at(taxon);

    vector<string> output;
    for (int ig = 0; ig < 65; ig++) {
	string s;
	s += GeneticCodeTable[gencode][ig];
	output.push_back(s);
    }                                                                                                  //cout<<s<<endl;}

    return output;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<string> MAKEmyCALletters(int partition, int taxon)
{
    int gencode = (mycodes.at(partition)).at(taxon);

    vector<string> output;
    for (int ig = 0; ig < 65; ig++) {
	string s; 

	s += ::tolower(GeneticCodeTable[gencode][ig]);
	output.push_back(s);
    }   

    return output;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void makeprintseqINT(int partition, int expectedcount, vector<int>& inspos, vector<int>& tempintseq, vector<insert>& insertstuff, int whichseq, ofstream& results3, int whichresults)
{
    int currentcount = 0, seqnow;

    // leaf and ancestral printing is kept seperate for speed

    // translates core sequence to alphabetical form, includes insertion and prints to file.
    string now, blank = "-", insblank = "-";

    if (markdeletedinsertions) {
	insblank = "*";
    }

    vector<string> myletters, myletters2;

    if (model_type == model::Type::nucleotide) {
	myletters = myDUletters;
	if (insertaslowercase) {
	    myletters2 = myDLletters;
	} else {
	    myletters2 = myletters;
	}
    } else {
	if (model_type == model::Type::aminoacid) {
	    myletters = myAUletters;
	    if (insertaslowercase) {
		myletters2 = myALletters;
	    } else {
		myletters2 = myletters;
	    }
	} else {
	    if (model_type == model::Type::codon) {
		if (printcodonsasDNA) {
		    myletters = myCDUletters;
		    if (insertaslowercase) {
			myletters2 = myCDLletters;
		    } else {
			myletters2 = myletters;
		    }
		} else {
		    myletters = MAKEmyCAUletters(partition, whichseq);
		    if (insertaslowercase) {
			myletters2 = MAKEmyCALletters(partition, whichseq);
		    } else {
			myletters2 = myletters;
		    }
		}
		blank    = "---";
		insblank = "---";
		if (markdeletedinsertions) {
		    insblank = "***";
		}
	    } else {
		cout << "/nERROR in makeprintseq" << endl << " ";
	    }
	}
    }

    // ancestral sequences


    //eternal link
    if (inspos.at(0) != -1) {
	insert *ii = &(insertstuff.at(inspos.at(0)));

	//if((*ii).length==0) continue;

	vector<site> *s = &((*ii).insertvec);

	for (int j = 0; j < (*s).size(); j++) {
	    seqnow = ((*s).at(j)).base;

	    if (seqnow == -1) {
		results3 << insblank;
		currentcount++;
	    }
	    else {
		results3 << myletters2.at(seqnow);
		currentcount++;
	    }
	}
    }


    for (int i = 1; i < tempintseq.size(); i++) {
	seqnow = tempintseq.at(i);

	if (seqnow == -1) {
	    results3 << blank;
	    currentcount++;
	}
	else {
	    results3 << myletters.at(seqnow);
	    currentcount++;
	}


	if (inspos.at(i) != -1) {
	    insert *ii = &(insertstuff.at(inspos.at(i)));

	    //if((*ii).length==0) continue;

	    vector<site> *s = &((*ii).insertvec);

	    for (int j = 0; j < (*s).size(); j++) {
		seqnow = ((*s).at(j)).base;

		if (seqnow == -1) {
		    results3 << insblank;
		    currentcount++;
		}
		else {
		    results3 << myletters2.at(seqnow);
		    currentcount++;
		}
	    }
	}
    }
//if(currentcount!=expectedcount) cout<<endl<<"ERROR in length of sequence in this file"<<endl<<"currentcount was "<<currentcount<<" and expectedcount is "<<expectedcount<<endl;
}


void makeprintseqLEAF(int partition, int expectedcount, vector<int>& inspos, vector<int>& tempintseq,
		      vector<insert>& insertstuff, int whichseq, ofstream& results, ofstream& results2, int whichresults)
{
    int currentcount = 0, seqnow;
    // leaf and ancestral printing is kept seperate for speed

    // translates core sequence to alphabetical form, includes insertion and prints to file.
    string now, blank = "-", insblank = "-";

    if (markdeletedinsertions) {
	insblank = "*";
    }

    vector<string> myletters, myletters2;

    if (model_type == model::Type::nucleotide) {
	myletters = myDUletters;
	if (insertaslowercase) {
	    myletters2 = myDLletters;
	} else {
	    myletters2 = myletters;
	}
    } else {
	if (model_type == model::Type::aminoacid) {
	    myletters = myAUletters;
	    if (insertaslowercase) {
		myletters2 = myALletters;
	    } else {
		myletters2 = myletters;
	    }
	} else {
	    if (model_type == model::Type::codon) {
		if (printcodonsasDNA) {
		    myletters = myCDUletters;
		    if (insertaslowercase) {
			myletters2 = myCDLletters;
		    } else {
			myletters2 = myletters;
		    }
		} else {
		    myletters = MAKEmyCAUletters(partition, whichseq);
		    if (insertaslowercase) {
			myletters2 = MAKEmyCALletters(partition, whichseq);
		    } else {
			myletters2 = myletters;
		    }
		}
		blank    = "---";
		insblank = "---";
		if (markdeletedinsertions) {
		    insblank = "***";
		}
	    } else {
		cout << "/nERROR in makeprintseq" << endl << " ";
	    }
	}
    }
    // Leaf sequences

    //eternal link
    if (inspos.at(0) != -1) {
	insert *ii = &(insertstuff.at(inspos.at(0)));

	//	if((*ii).length==0) continue;

	vector<site> *s = &((*ii).insertvec);

	for (int j = 0; j < (*s).size(); j++) {
	    seqnow = ((*s).at(j)).base;

	    if (seqnow == -1) {
		results << insblank;
		currentcount++;
	    }
	    else {
		now = myletters2.at(seqnow);
		results << now;
		results2 << now;
		currentcount++;
	    }
	}
    }



    for (int i = 1; i < tempintseq.size(); i++) {
	seqnow = tempintseq.at(i);

	if (seqnow == -1) {
	    results << blank;
	    currentcount++;
	}
	else {
	    now = myletters.at(seqnow);
	    results << now;
	    results2 << now;
	    currentcount++;
	}

	if (inspos.at(i) != -1) {
	    insert *ii = &(insertstuff.at(inspos.at(i)));

	    //	if((*ii).length==0) continue;

	    vector<site> *s = &((*ii).insertvec);

	    for (int j = 0; j < (*s).size(); j++) {
		seqnow = ((*s).at(j)).base;

		if (seqnow == -1) {
		    results << insblank;
		    currentcount++;
		}
		else {
		    now = myletters2.at(seqnow);
		    results << now;
		    results2 << now;
		    currentcount++;
		}
	    }
	}
    }
//if(currentcount!=expectedcount) cout<<endl<<"ERROR in length of sequence in this file"<<endl<<"currentcount was "<<currentcount<<" and expectedcount is "<<expectedcount<<endl;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PrintProgress(int blocknumber, int totalblock, int repnumber, int totalrep, int percent, string extra)
{
    // screen output of progress when running
    string extra2 = "                          ";

    if (extra != "") {
	extra2 = extra;
    }
    //cout<<"\r";
    cout << endl;
    cout << "   Block " << blocknumber << "/" << totalblock << "  rep ";
    if ((repnumber < 10) && (totalrep > 9)) {
	cout << " ";
    }
    if ((repnumber < 100) && (totalrep > 99)) {
	cout << " ";
    }
    if ((repnumber < 1000) && (totalrep > 999)) {
	cout << " ";
    }
    cout << repnumber << "/" << totalrep << "  ";
    if (percent < 100) {
	cout << " ";
    }
    if (percent < 10) {
	cout << " ";
    }
    if (extra == "") {
	cout << percent << " % of treelength done.  ";
    } else {
	cout << extra2;
    }
}


void PrintProgress3(int blocknumber, int totalblock, int repnumber, int totalrep, int percent, string extra)
{
    // screen output of progress when running
    string extra2 = "                          ";

    if (extra != "") {
	extra2 = extra;
    }
    //cout<<"\r";
    cout << endl;
    cout << "   Block " << blocknumber << "/" << totalblock << "  rep ";
    if ((repnumber < 10) && (totalrep > 9)) {
	cout << " ";
    }
    if ((repnumber < 100) && (totalrep > 99)) {
	cout << " ";
    }
    if ((repnumber < 1000) && (totalrep > 999)) {
	cout << " ";
    }
    cout << repnumber << "/" << totalrep << "  ";
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PrintProgress2(int blocknumber, int numberofevolveblocks, int thisrep, int reps, int printcount)
{
    // screen output of progress when running
//	stringstream fd2;  fd2<<total-1; string fd2bit=fd2.str();

    stringstream fd1;

    fd1 << printcount;
    string fd1bit = fd1.str();

    string mystring = "Printing sequence " + fd1bit;



    /*
     * string mystring="Printing sequence ";
     * if(total+1>9 && printcount<10) mystring+=" ";
     * if(total+1>99 && printcount<100) mystring+=" ";
     * if(total+1>999 && printcount<1000) mystring+=" ";
     *
     * string endbit=fd1bit+" of "+fd2bit+"    ";
     * mystring+=endbit;
     */

//	cout<<"\t"<<mystring<<endl;

    PrintProgress(blocknumber, numberofevolveblocks, thisrep, reps, 100, mystring);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void sumrates(RATES& rates)
{
/*
 *      rates.totallength = rates.corelength + rates.inslength;
 *
 *
 * //specified this way when raw
 *      rates.totalrate = rates.coresubrate	+ rates.inssubrate	+(  (rates.totallength -1)*((*m).insertrate + (*m).deleterate)  )
 *
 + (  ((*m).deleterate)*((*m).dellength)  )  + (*m).insertrate - (*m).deleterate;
 +
 */

/*
 * //but these are specified in buildsums functions
 *      x1=rates.coreinsertrate; x2=((*m).insertrate)*(rates.corelength); diff=x2-x1; if(diff<0) diff=-diff; if(diff>0.00001) {cout<<"ERROR in core insertionrate in buildsums"<<endl; cout<<rates.coreinsertrate<<"  "<<(*m).insertrate<<"  "<<rates.corelength<<"  "<<((*m).insertrate)*(rates.corelength)<<endl;}
 *      x1=rates.coredeleterate; x2=((*m).deleterate)*(rates.corelength-1); diff=x2-x1; if(diff<0) diff=-diff; if(diff>0.00001) {cout<<"ERROR in core deletionrate in buildsums"<<endl;  cout<<rates.coredeleterate<<"  "<<(*m).deleterate<<"  "<<rates.corelength<<"  "<<((*m).deleterate)*(rates.corelength)<<endl;}
 *      x1=rates.insinsertrate;  x2=((*m).insertrate)*(rates.inslength);  diff=x2-x1; if(diff<0) diff=-diff; if(diff>0.00001) {cout<<"ERROR in ins insertionrate in buildsums"<<endl;  cout<<rates.insinsertrate <<"  "<<(*m).insertrate<<"  "<<rates.inslength<<"  "<<((*m).insertrate)*(rates.inslength)<<endl;}
 *      x1=rates.insdeleterate;  x2=((*m).deleterate)*(rates.inslength);  diff=x2-x1; if(diff<0) diff=-diff; if(diff>0.00001) {cout<<"ERROR in ins deletionrate in buildsums"<<endl;   cout<<rates.insdeleterate <<"  "<<(*m).deleterate<<"  "<<rates.inslength<<"  "<<((*m).deleterate)*(rates.inslength)<<endl;}
 *
 * (rates.coredeleterate)+=(  ((*m).deleterate)*(((*m).dellength)-1)  );
 *
 */

    //so we can just say:

    rates.totalrate = rates.coresubrate + rates.inssubrate + rates.coredeleterate + rates.insdeleterate + rates.coreinsertrate + rates.insinsertrate;

}


////////////////////


int buildsumsold(RATES& rates, SUMS& sums, vector<int>& fromseq, vector<insert>& fromins, vector<int>& inspos)
{
    int i, j, y = 0, size = rates.rootlength - 1, lastinslength = rates.inslength, lastcorelength = rates.corelength;
    double x = 0, z = 0;

    rates.insinsertrate  = 0;
    rates.insdeleterate  = 0;
    rates.inssubrate     = 0;
    rates.inslength      = 0; // not used in this method
    rates.corelength     = 0; // not used in this method
    rates.coreinsertrate = 0;
    rates.coredeleterate = 0;
    rates.coresubrate    = 0;

    sums.clear();

    (sums.IIlevels.at(0)).reserve(rates.rootlength);
    (sums.CIlevels.at(0)).reserve(rates.rootlength);
    (sums.CSlevels.at(0)).reserve(rates.rootlength);
    (sums.ISlevels.at(0)).reserve(rates.rootlength);


    site *s;
    insert *ii;
    vector<site> *ss;

    for (i = 0; i < rates.rootlength; i++) {
	if (inspos.at(i) == -1) {
	    (sums.ISlevels.at(0)).push_back(0);
	    (sums.IIlevels.at(0)).push_back(0);
	} else {
	    //deal with insertion information.

	    ii = &(fromins.at(inspos.at(i)));

	    ss = &((*ii).insertvec);

	    int templength = (*ss).size();

	    (*ii).subrate = 0;
	    (*ii).length  = 0;

	    for (j = 0; j < templength; j++) {
		s = &((*ss).at(j));

		if (((*s).base) == -1) {
		    (*s).subrate = 0;
		} else {
		    ((*ii).length)++;
		    (*s).subrate = returnsitesubrate((*s).rate, (*s).siteclass, (*s).base);

		    ((*ii).subrate) += (*s).subrate;
		}
	    }

	    rates.insinsertrate += (((*m).insertrate) * ((*ii).length));
	    rates.insdeleterate += (((*m).deleterate) * ((*ii).length));

	    (sums.IIlevels.at(0)).push_back((*ii).length);
	    (sums.ISlevels.at(0)).push_back((*ii).subrate);

	    //end of dealing with insertion information
	}

	if (i == 0) {
	    rates.coreinsertrate += (*m).insertrate;
	    (sums.CIlevels.at(0)).push_back(1);
	    (sums.CSlevels.at(0)).push_back(0);
	} else {
	    int test = fromseq.at(i);

	    if (test == -1) {
		(sums.CIlevels.at(0)).push_back(0);
		(sums.CSlevels.at(0)).push_back(0);
	    } else {
		(sums.CSlevels.at(0)).push_back(returnsitesubrate(ratevec.at(i), siteclassvec.at(i), test));

		rates.coreinsertrate += (*m).insertrate;
		rates.coredeleterate += (*m).deleterate;

		(sums.CIlevels.at(0)).push_back(1);
	    }
	}
    }

    // Accumulate sums according to log-level
    for (int loglevel = 0; loglevel < sums.IIlevels.size()-1; loglevel++) {
	auto& level = sums.IIlevels.at(loglevel);
	int y = 0;
	int size = level.size();
	for (int i = 0; i < size; i++) {
	    y += level.at(i);
	    if ((i % 10 == 9) || (i == size-1)) {
		sums.IIlevels.at(loglevel+1).push_back(y);
		y = 0;
	    }
	}
    }

    // Accumulate sums according to log-level
    for (int loglevel = 0; loglevel < sums.ISlevels.size()-1; loglevel++) {
	auto& level = sums.ISlevels.at(loglevel);
	double y = 0;
	int size = level.size();
	for (int i = 0; i < size; i++) {
	    y += level.at(i);
	    if ((i % 10 == 9) || (i == size-1)) {
		sums.ISlevels.at(loglevel+1).push_back(y);
		y = 0;
	    }
	}
    }

    // sum the last level to get to total core length
    for (auto y : *(sums.ISlevels.rbegin())) 
	rates.inssubrate += y;



    // Accumulate CIsums according to log-level
    for (int loglevel = 0; loglevel < sums.CIlevels.size()-1; loglevel++) {
	auto const& level = sums.CIlevels.at(loglevel);
	int y = 0;
	int size = level.size();
	for (int i = 0; i < size; i++) {
	    y += level.at(i);
	    if ((i % 10 == 9) || (i == size-1)) {
		sums.CIlevels.at(loglevel+1).push_back(y);
		y = 0;
	    }
	}
    }

    // sum the last level to get to total core length
    for (auto y : *(sums.CIlevels.rbegin())) 
	rates.corelength += y;

    // Accumulate CSsums according to log-level
    for (int loglevel = 0; loglevel < sums.CSlevels.size()-1; loglevel++) {
	auto const& level = sums.CSlevels.at(loglevel);
	double y = 0;
	int size = level.size();
	for (int i = 0; i < size; i++) {
	    y += level.at(i);
	    if ((i % 10 == 9) || (i == size-1)) {
		sums.CSlevels.at(loglevel+1).push_back(y);
		y = 0;
	    }
	}
    }

    // sum the last level to get to total core subrate length
    rates.coresubrate = 0;  // temporary during transition
    for (auto y : *(sums.CSlevels.rbegin())) 
	rates.coresubrate += y;

    rates.corerate = rates.coreinsertrate + rates.coredeleterate;

    rates.insrate = rates.insinsertrate + rates.insdeleterate;

    double x1, x2, diff;

    x1   = rates.coreinsertrate;
    x2   = ((*m).insertrate) * double((rates.corelength));
    diff = x2 - x1;
    if (diff < 0) {
	diff = -diff;
    }
    if (diff > 0.001) {
	cout << diff << " " << x1 << " " << x2 << " ERROR in core insertionrate in buildsums" << endl;
	cout << rates.coreinsertrate << "  " << (*m).insertrate << "  " << rates.corelength << "  " << ((*m).insertrate) * double((rates.corelength)) << endl;
    }
    x1   = rates.coredeleterate;
    x2   = ((*m).deleterate) * double((rates.corelength - 1));
    diff = x2 - x1;
    if (diff < 0) {
	diff = -diff;
    }
    if (diff > 0.001) {
	cout << diff << " " << x1 << " " << x2 << " ERROR in core deletionrate in buildsums" << endl;
	cout << rates.coredeleterate << "  " << (*m).deleterate << "  " << rates.corelength << "  " << ((*m).deleterate) * double((rates.corelength)) << endl;
    }
    x1   = rates.insinsertrate;
    x2   = ((*m).insertrate) * double((rates.inslength));
    diff = x2 - x1;
    if (diff < 0) {
	diff = -diff;
    }
    if (diff > 0.001) {
	cout << diff << " " << x1 << " " << x2 << " ERROR in ins insertionrate in buildsums" << endl;
	cout << rates.insinsertrate << "  " << (*m).insertrate << "  " << rates.inslength << "  " << ((*m).insertrate) * double((rates.inslength)) << endl;
    }
    x1   = rates.insdeleterate;
    x2   = ((*m).deleterate) * double((rates.inslength));
    diff = x2 - x1;
    if (diff < 0) {
	diff = -diff;
    }
    if (diff > 0.001) {
	cout << diff << " " << x1 << " " << x2 << " ERROR in ins deletionrate in buildsums" << endl;
	cout << rates.insdeleterate << "  " << (*m).deleterate << "  " << rates.inslength << "  " << ((*m).deleterate) * double((rates.inslength)) << endl;
    }

    if (lastinslength != rates.inslength) {
	cout << "ERROR in total inserted sites length in buildsums " << endl << "lastinslength was " << lastinslength << "  and rates.inslength was " << rates.inslength << endl;
    }

    (rates.coredeleterate) += (((*m).deleterate) * (((*m).delmeansize) - 1));
    sumrates(rates);

    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int buildsumsnew(RATES& rates, SUMS& sums, vector<int>& fromseq, vector<insert>& fromins, vector<int>& inspos)
{
    int size = rates.rootlength - 1, lastinslength = rates.inslength, lastcorelength = rates.corelength;
    double x = 0;

    rates.insinsertrate  = 0;
    rates.insdeleterate  = 0;
    rates.inssubrate     = 0;   // not used in this method
    rates.inslength      = 0;
    rates.corelength     = 0;
    rates.coreinsertrate = 0;
    rates.coredeleterate = 0;
    rates.coresubrate    = 0;   // not used in this method

    sums.clear();

    (sums.IIlevels.at(0)).reserve(rates.rootlength);
    (sums.CIlevels.at(0)).reserve(rates.rootlength);

    site *s;
    insert *ii;
    vector<site> *ss;

    for (int i = 0; i < rates.rootlength; i++) {
	if (inspos.at(i) == -1) {
	    (sums.IIlevels.at(0)).push_back(0);
	} else {
	    //deal with insertion information.

	    ii = &(fromins.at(inspos.at(i)));

	    ss = &((*ii).insertvec);

	    int templength = (*ss).size();

	    (*ii).length = 0;

	    for (int j = 0; j < templength; j++) {
		s = &((*ss).at(j));

		if (((*s).base) != -1) {
		    ((*ii).length)++;
		}
	    }

	    //(rates.inslength)+=((*ii).length);

	    rates.insinsertrate += (((*m).insertrate) * ((*ii).length));
	    rates.insdeleterate += (((*m).deleterate) * ((*ii).length));

	    (sums.IIlevels.at(0)).push_back((*ii).length);

	    //end of dealing with insertion information
	}

	//deal with core information

	if (i == 0) {
	    rates.coreinsertrate += (*m).insertrate;
	    (sums.CIlevels.at(0)).push_back(1);
	} else {

	    if (fromseq.at(i) == -1) {
		(sums.CIlevels.at(0)).push_back(0);
	    } else {
		rates.coreinsertrate += (*m).insertrate;
		rates.coredeleterate += (*m).deleterate;
		(sums.CIlevels.at(0)).push_back(1);
	    }
	}
    }

    // Accumulate sums for insertion at insert site levels (IIlevels)
    for (int loglevel = 0; loglevel < sums.IIlevels.size()-1; loglevel++) {
	vector<int>& level = sums.IIlevels.at(loglevel);
	int y = 0;
	int size = level.size();
	for (int i = 0; i < size; i++) {
	    y += level.at(i);
	    if ((i % 10 == 9) || (i == size-1)) {
		sums.IIlevels.at(loglevel+1).push_back(y);
		y = 0;
	    }
	}
    }

    // sum the last level to get to total insertion length
    for (auto y : *(sums.IIlevels.rbegin())) 
	rates.inslength += y;


    // Accumulate sums for insertion at core site levels (CIlevels)
    for (int loglevel = 0; loglevel < sums.CIlevels.size()-1; loglevel++) {
	vector<int>& level = sums.CIlevels.at(loglevel);
	int y = 0;
	int size = level.size();
	for (int i = 0; i < size; i++) {
	    y += level.at(i);
	    if ((i % 10 == 9) || (i == size-1)) {
		sums.CIlevels.at(loglevel+1).push_back(y);
		y = 0;
	    }
	}
    }
    
    // sum the last level to get to total core insertion length
    for (auto y : *(sums.CIlevels.rbegin())) 
	rates.corelength += y;


//		cout<<"           QWE 4"<<endl;
    rates.corerate = rates.coreinsertrate + rates.coredeleterate;

    rates.insrate = rates.insinsertrate + rates.insdeleterate;

    double x1, x2, diff;

    x1   = rates.coreinsertrate;
    x2   = ((*m).insertrate) * double((rates.corelength));
    diff = x2 - x1;
    if (diff < 0) {
	diff = -diff;
    }
    if (diff > 0.001) {
	cout << diff << " " << x1 << " " << x2 << " ERROR in core insertionrate in buildsums" << endl;
	cout << rates.coreinsertrate << "  " << (*m).insertrate << "  " << rates.corelength << "  " << ((*m).insertrate) * double((rates.corelength)) << endl;
    }
    x1   = rates.coredeleterate;
    x2   = ((*m).deleterate) * double((rates.corelength - 1));
    diff = x2 - x1;
    if (diff < 0) {
	diff = -diff;
    }
    if (diff > 0.001) {
	cout << diff << " " << x1 << " " << x2 << " ERROR in core deletionrate in buildsums" << endl;
	cout << rates.coredeleterate << "  " << (*m).deleterate << "  " << rates.corelength << "  " << ((*m).deleterate) * double((rates.corelength)) << endl;
    }
    x1   = rates.insinsertrate;
    x2   = ((*m).insertrate) * double((rates.inslength));
    diff = x2 - x1;
    if (diff < 0) {
	diff = -diff;
    }
    if (diff > 0.001) {
	cout << diff << " " << x1 << " " << x2 << " ERROR in ins insertionrate in buildsums" << endl;
	cout << rates.insinsertrate << "  " << (*m).insertrate << "  " << rates.inslength << "  " << ((*m).insertrate) * double((rates.inslength)) << endl;
    }
    x1   = rates.insdeleterate;
    x2   = ((*m).deleterate) * double((rates.inslength));
    diff = x2 - x1;
    if (diff < 0) {
	diff = -diff;
    }
    if (diff > 0.001) {
	cout << diff << " " << x1 << " " << x2 << " ERROR in ins deletionrate in buildsums" << endl;
	cout << rates.insdeleterate << "  " << (*m).deleterate << "  " << rates.inslength << "  " << ((*m).deleterate) * double((rates.inslength)) << endl;
    }

    if (lastinslength != rates.inslength) {
	cout << "ERROR in total inserted sites length in buildsums " << endl << "lastinslength was " << lastinslength << "  and rates.inslength was " << rates.inslength << endl;
    }

//	double dd; dd=lastcorelength - rates.corelength; if(dd<0) dd=-dd;
//	if(dd>0.0001) cout<<"ERROR in total core sequence length in buildsums "<<endl<<"lastcorelength was "<<lastcorelength<<"  and rates.corelength was "<<rates.corelength<<endl;

    //if(lastcorelength != rates.corelength-1) cout<<"ERROR in total core sequence length in buildsums "<<endl<<"lastcorelength was "<<lastcorelength<<"  and rates.corelength was "<<rates.corelength<<endl;

    (rates.coredeleterate) += (((*m).deleterate) * (((*m).delmeansize) - 1));

    sumrates(rates);

//	cout<<rates.totalrate<<"  is the TOTAL rate"<<endl;

    return 0;
}

// convert an enumeration value to integer (for purposes of printing, etc.)
// http://stackoverflow.com/a/11421471/1135316
//
template <typename Enumeration>
auto as_integer(Enumeration const value)
    -> typename std::underlying_type<Enumeration>::type
{
    return static_cast<typename std::underlying_type<Enumeration>::type>(value);
}

// select a position for subsitution events in inserted sites
//
int findpos_ins_sub(Event event, vector<int>& updatepositions, double unirand, const SUMS& sum, double& S)
{
    // substitution in inserted sites
    updatepositions.clear();
    S = 0;
    int pos = -1;
    double s = 0;
    int j = 0;
    for (int loglevel = sum.ISlevels.size()-1; loglevel >= 0; loglevel--) {
	auto const & level = sum.ISlevels.at(loglevel);
	pos = -1;
	for (int i = j; i < level.size(); i++) {
	    double s = level.at(i);
	    if (unirand <= s + S) {
		j = 10 * i;
		updatepositions.push_back(i);
		pos = i;
		break;
	    } else {
		S += s;
	    }
	}
	if (pos == -1) {
	    cout << "ERROR in findpos_ins_sub level " << loglevel << " at event " << as_integer(event) << endl;
	    return -1;
	}
    }

    {
	double s = sum.ISlevels.begin()->at(pos);
	assert(!(S - unirand > 0));
	assert(!(unirand - s - S > 0));
    }

    return pos;
}


// select a position for substitution events in core sites
//
int findpos_core_sub(Event event, vector<int>& updatepositions, double unirand, const SUMS& sum, double& S)
{
    //event 0	substitution in core
    // in core sequence sites: 0 for substitution, 2 for insertion, 4 for deletion
    // in inserted sites:      1 for substitution, 3 for insertion, 5 for deletion

    int j = 0;
    int pos = -1;

    updatepositions.clear();
    S = 0;
    for (int loglevel = sum.CSlevels.size()-1; loglevel >= 0; loglevel--) {
	auto const & level = sum.CSlevels.at(loglevel);
	pos = -1;
	for (int i = j; i < level.size(); i++) {
	    double s = level.at(i);
	    if (unirand <= s + S) {
		j = 10 * i;
		updatepositions.push_back(i);
		pos = i;
		break;
	    } else {
		S += s;
	    }
	}
	if (pos == -1) {
	    cout << "ERROR in findpos_core_sub level " << loglevel << " at event " << as_integer(event) << endl;
	    // return -1;
	}
    }

    {
	double s = sum.CSlevels.begin()->at(pos);
	assert(!(S - unirand > 0));
	assert(!(unirand - s - S > 0));
    }
    // not needed now as simply prevent random numbers that are "EXACTLY ZERO" to be used here for that reason.
    //if(pos==0) pos=1; // guards against pseudo-random values of unirand that are exactly zero to machine precision
    //  - in this case it is possible to choose the "imaginary" eternal link position
    //  which we do not want to be possible!  this will not happend if unirand is non-zero.


    return pos;
}

// select a position for indel events in core sites
//
int findpos_core_indel(Event event, vector<int>& updatepositions, int mypos, const SUMS& sum, int& SI)
{
    // this is for event 2 or 4
    // in core sequence sites: 0 for substitution, 2 for insertion, 4 for deletion
    // in inserted sites:      1 for substitution, 3 for insertion, 5 for deletion

    int j = 0, si = 0;
    int pos = -1;

    updatepositions.clear();

    SI = 0;
    for (auto rit = sum.CIlevels.rbegin(); rit != sum.CIlevels.rend(); rit++) {
	auto const & level = *rit;
	pos = -1;
	for (int i = j; i < level.size(); i++) {
	    int si = level.at(i);
	    if (mypos <= si + SI) {
		j = 10 * i;
		updatepositions.push_back(i);
		pos = i;
		break;
	    } else {
		SI += si;
	    }
	}
	if (pos == -1) {
	    cout << "ERROR in findpos_core_indel level " << sum.CIlevels.rend() - rit << " at event " << as_integer(event) << endl;
	    return -1;
	}
    }
    return pos;
}


// select a position for indel events in inserted sites
//
int findpos_ins_indel(Event event, vector<int>& updatepositions, int mypos, SUMS& sum, int& SI)
{
    // this is for event 3 or 5
    //event numbers	// in core sequence sites: 0 for substitution, 2 for insertion, 4 for deletion
    // in inserted sites:                          1 for substitution, 3 for insertion, 5 for deletion

    updatepositions.clear();

    SI = 0;

    int i, j, si = 0;

    j = 0;
    SI = 0;
    int pos = -1;
    for (int loglevel = sum.IIlevels.size()-1; loglevel >= 0; loglevel--) {
	auto const & level = sum.IIlevels.at(loglevel);
	pos = -1;
	for (int i = j; i < level.size(); i++) {
	    int si = level.at(i);
	    if ((mypos <= si + SI) && (si + SI != 0)) {
		j = 10 * i;
		updatepositions.push_back(i);
		pos = i;
		break;
	    } else {
		SI += si;
	    }
	}
	if (pos == -1) {
	    cout << "ERROR in findpos_ins_indel level " << loglevel << " at event " << as_integer(event) << endl;
	    // return -1;
	}
    }

    return pos;
}


// Randomly select and event according to probabilities in `rates`

Event chooseevent(RATES& rates)
{
    int event = -1;     // in core sequence sites: 0 for substitution, 2 for insertion, 4 for deletion
    			// in inserted sites:      1 for substitution, 3 for insertion, 5 for deletion
    double rand = 1;

    while (rand == 1 || rand == 0) {
	rand = mtrand1();
    }  // prevents against pseudo-random numbers that are exactly equal to 1 or 0 to machine precision

    rand *= (rates.totalrate);

    // in core sequence sites: 0 for substitution, 2 for insertion, 4 for deletion
    // in inserted sites:      1 for substitution, 3 for insertion, 5 for deletion
    double temprate = rates.coreinsertrate;
    if (rand < temprate) {
	return Event::core_ins;  	// core sites insertion
    } 

    temprate += (rates.insinsertrate);
    if (rand < temprate) {
	return Event::ins_ins;       // inserted sites insertion
    } 

    temprate += (rates.coredeleterate);
    if (rand < temprate) {
	return Event::core_del;       // core sites deletion
    } 

    temprate += (rates.insdeleterate);
    if (rand < temprate) {
	return Event::ins_del;       // inserted sites deletion
    } 

    temprate += (rates.coresubrate);
    if (rand < temprate) {
	return Event::core_sub;       // core sites substitution
    }

    temprate += (rates.inssubrate);
    if (rand < temprate) {
	return Event::ins_sub;       // inserted sites substitution
    } 

    cerr << "ERROR IN CHOOSEEVENT - random number was " << rand << "  and rates.totalrate was " << rates.totalrate << " and temprate was " << temprate << endl;
    return Event::error;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

int updatesubsums0(vector<int> updatepositions, double sdiff, SUMS& sums)
{
    // substitution in core sequence
    assert(updatepositions.size() == sums.CSlevels.size());

    int i = updatepositions.size()-1;
    for (auto &level : sums.CSlevels) {
	level.at(updatepositions.at(i)) += sdiff;
	i--;
    }

    return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

int updatesubsums1(vector<int> updatepositions, double sdiff, SUMS& sums)
{
    // substitution in inserted sites
    assert(updatepositions.size() == sums.ISlevels.size());

    int i = updatepositions.size()-1;
    for (auto &level : sums.ISlevels) {
	level.at(updatepositions.at(i)) += sdiff;
	i--;
    }

    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int updateindelsums24(vector<int> updatepositions, int idiff, SUMS& sums)
{
    // indel in core sequence
    assert(updatepositions.size() == sums.CIlevels.size());

    int i = updatepositions.size()-1;
    for (auto &level : sums.CIlevels) {
	level.at(updatepositions.at(i)) += idiff;
	i--;
    }

    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int updateindelsums35(vector<int> updatepositions, int idiff, SUMS& sums)
{
    // indel in inserted sites

    assert(updatepositions.size() == sums.IIlevels.size());

    int i = updatepositions.size()-1;
    for (auto &level : sums.IIlevels) {
	level.at(updatepositions.at(i)) += idiff;
	i--;
    }

    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void dodeletion(int inspos, SUMS& sums, RATES& rates, vector<int>& newseqINT,
		vector<int>& updatepositions, vector<insert>& insINT2, vector<int>& insPOS2, int endlength, int indellength)
{
    //inspos = -1 is deletion beginning in core.
    //          otherwise deletion starts at position inspos in the "insert" at that core position currentev.

    // event 4 is deletion in core, event 5 is deletion in inserted sites

    int coreindellength = 0, insindellength = 0;
    //12345567
    site *s;
    insert *ii;

    vector<site> *ss;

    deletioncount++;

    int idiff = 0;

    double csdiff = 0, isdiff = 0, diffs = 0;


    int currentev = updatepositions.at(10);


#ifdef INDELLOG
    int tttt;
    if (inspos == -1) {
	tttt = 1;
	indellog << "deletion in core\tlength\t" << indellength << "\tcore pos\t" << currentev << "\t";
    } else {
	tttt = 2;
	indellog << "deletion in ins\tlength\t" << indellength << "\tcore pos\t" << currentev << "\tins pos\t" << inspos << "\t";
    }
#endif

    while (indellength > 0 && currentev != endlength) {
	if (inspos == -1) {
	    if (newseqINT.at(currentev) > -1) {
#ifdef INDELLOG
		if (tttt == 2) {
		    tttt = -1;
		    indellog << "\tthat went out to core" << endl;
		}
		indellog << " " << newseqINT.at(currentev) << " ";
#endif

		indellength--;

		newseqINT.at(currentev) = -1;

		if (oldmethod) {
		    diffs = -((sums.CSlevels.at(0)).at(currentev));

		    csdiff += diffs;

		    updatesubsums0(updatepositions, diffs, sums);
		}

		coreindellength--;

		updateindelsums24(updatepositions, -1, sums);
	    }

	    inspos = 0;
	}

	if (insPOS2.at(currentev) != -1) {
	    //ii=&(insINT2.at(currentev));
	    ii = &(insINT2.at(insPOS2.at(currentev)));

	    ss = &((*ii).insertvec);

	    //`ddiff=0;
	    idiff = 0;
	    diffs = 0;

	    while (indellength > 0 && inspos < (*ss).size()) {       //(*ii).length)
		s = &((*ss).at(inspos));
		inspos++;

		if ((*s).base == -1) {
		    continue;
		}


#ifdef INDELLOG
		if (tttt == 1) {
		    tttt = -1;
		    indellog << "\tthat went into an insertion" << endl;
		}
		indellog << " " << (*s).base << " ";
#endif

		(*s).base = -1;

		if (oldmethod) {
		    diffs -= ((*s).subrate);

		    ((*s).subrate) = 0;
		}

		idiff--;
		indellength--;
	    }


	    if (oldmethod) {
		((*ii).subrate) += diffs;

		isdiff += diffs;

		updatesubsums1(updatepositions, diffs, sums);
	    }

	    //indellength+=idiff;

	    ((*ii).length) += idiff;

	    insindellength += idiff;

	    updateindelsums35(updatepositions, idiff, sums);
	}

	inspos = -1;
	currentev++;
	(updatepositions.at(10))++;

	if (currentev % 10 == 0) {
	    (updatepositions.at(9))++;
	} else {
	    continue;
	}
	if (currentev % 100 == 0) {
	    (updatepositions.at(8))++;
	} else {
	    continue;
	}
	if (currentev % 1000 == 0) {
	    (updatepositions.at(7))++;
	} else {
	    continue;
	}
	if (currentev % 10000 == 0) {
	    (updatepositions.at(6))++;
	} else {
	    continue;
	}
	if (currentev % 100000 == 0) {
	    (updatepositions.at(5))++;
	} else {
	    continue;
	}
	if (currentev % 1000000 == 0) {
	    (updatepositions.at(4))++;
	} else {
	    continue;
	}
	if (currentev % 10000000 == 0) {
	    (updatepositions.at(3))++;
	} else {
	    continue;
	}
	if (currentev % 100000000 == 0) {
	    (updatepositions.at(2))++;
	} else {
	    continue;
	}
	if (currentev % 1000000000 == 0) {
	    (updatepositions.at(1))++;
	} else {
	    continue;
	}
	//if(currentev%10000000000==0)(updatepositions.at(0))++; else continue;
    }

    (rates.coreinsertrate) += (coreindellength * ((*m).insertrate));
    (rates.coredeleterate) += (coreindellength * ((*m).deleterate));

    (rates.insinsertrate) += (insindellength * ((*m).insertrate));
    (rates.insdeleterate) += (insindellength * ((*m).deleterate));

    (rates.corelength) += coreindellength;
    (rates.inslength)  += insindellength;

    (rates.coresubrate) += csdiff;
    (rates.inssubrate)  += isdiff;

    sumrates(rates);

#ifdef INDELLOG
    indellog << endl;
#endif
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void doinsertion(double timeleft, int inspos, SUMS& sums, RATES& rates, vector<int>& newseqINT, vector<int>& updatepositions, vector<insert>& insINT2,
		 vector<vector<insert> >& insINT, int label, vector<vector<int> >& insPOS, vector<int>& insPOS2, int indellength)
{
    //inspos = -1 is insertion beginning in core.
    //          otherwise insertion starts at position inspos in the "insert" at that core position currentev.

    // event 2 is insertion in core, event 3 is insertion in inserted sites
    // 2 is never really used.
    // if "2" i.e. an insertion in core then insertion comes at beginning of vector at <insert> for that core positions.

    site *s;
    insert *ii;

    vector<site> *ss;
    vector<insert> *ivec;

    double ddiff = 0, sdiff = 0;


    (rates.insinsertrate) += (indellength * ((*m).insertrate));
    (rates.insdeleterate) += (indellength * ((*m).deleterate));

    (rates.inslength) += indellength;

    dashlength += indellength;
    insertioncount++;

    vector<site> insertseq;

    makeseq2(indellength, insertseq, sdiff, timeleft);

    int currentev = updatepositions.at(10);

    if (oldmethod) {
	updatesubsums1(updatepositions, sdiff, sums);
    }


    updateindelsums35(updatepositions, indellength, sums);

//	vector<site> blank; for(int j=0; j<indellength; j++) blank.push_back(site(-1,-1,-1,-1,-1));

    vector<site> blank = insertseq;
    for (int j = 0; j < indellength; j++) {
	s             = &(blank.at(j));
	(*s).base     = -1;
	(*s).subrate  = 0;
	(*s).timeleft = -1;
    }

    if (insPOS2.at(currentev) == -1) {
	insPOS2.at(currentev) = insINT2.size();
	insINT2.push_back(insert(insertseq, indellength, sdiff));

#ifdef INDELLOG
	indellog << "new insertion\tlength\t" << indellength << "\tcore pos\t" << currentev << "\t";
	for (int j = 0; j < indellength; j++) {
	    s = &(insertseq.at(j));
	    indellog << (*s).base;
	}
	indellog << endl;
#endif
    } else {
	ii = &(insINT2.at(insPOS2.at(currentev)));

	ss = &((*ii).insertvec);

	((*ii).length) += indellength;

	((*ii).subrate) += sdiff;

	if (inspos == -1) {
	    (*ss).insert((*ss).end(), insertseq.begin(), insertseq.end());
	} else {
	    (*ss).insert((*ss).begin() + inspos, insertseq.begin(), insertseq.end());
	}


#ifdef INDELLOG
	indellog << "insertion in insertion\tlength\t" << indellength << "\tcore pos\t" << currentev << "\tinsert pos\t" << inspos << endl;
	for (int j = 0; j < indellength; j++) {
	    s = &(insertseq.at(j));
	    indellog << (*s).base;
	}
	indellog << endl;
#endif
    }

    for (int i = 0; i < insINT.size(); i++) {
	if (i == label) {
	    continue;
	}

	ivec = &(insINT.at(i));
	if ((*ivec).empty()) {
	    continue;
	}

	vector<int> *myinsPOS = &(insPOS.at(i));

	int *temp = &((*myinsPOS).at(currentev));

	if ((*temp) == -1) {
	    (*temp) = (*ivec).size();

	    (*ivec).push_back(insert(blank, 0, 0));
	} else {
	    ii = &((*ivec).at(*temp));

	    ss = &((*ii).insertvec);

	    if (inspos == -1) {
		(*ss).insert((*ss).end(), blank.begin(), blank.end());
	    } else {
		(*ss).insert((*ss).begin() + inspos, blank.begin(), blank.end());
	    }
	}
    }

    (rates.inssubrate) += sdiff;

    sumrates(rates);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void func(int branchID, double branchlength, RATES& rates, vector<int>& newseqINT, vector<insert>& insINT2,
	  SUMS& sums, vector<vector<insert> >& insINT, int label, vector<int>& insPOS2, vector<vector<int> >& insPOS)
{
    //cout<<"Q "<<rates.totalrate<<endl;

    /* EVENT */                 // in core sequence sites: 0 for substitution, 2 for insertion, 4 for deletion
    // in inserted sites:      1 for substitution, 3 for insertion, 5 for deletion


//	cout<<"NEW BRANCH"<<endl;
    vector<int> updatepositions;

    double currenttime = branchlength, lastlength;

    double exprand = expmyrand(rates.totalrate);

    int currentbase;
    Event event;
    double testlength = currenttime - 0.01;

    int goodcount = 0, badcount = 0;

    site *s;
    insert *ii;
    vector<site> *ss;

    if (oldmethod || ((*m).insertrate != 0) || ((*m).deleterate != 0)) {
	while (currenttime > exprand) {
	    //cout<<currenttime<<"\t"<<exprand<<"\t"<<currenttime-exprand<<endl;

	    currenttime -= exprand;   //(((*m).Rrates).at(  ratevec.at(currentev)  ));


	    double unirand = 0, sdiff = 0, dsumreached = 0, ratepos;
	    int intrand    = 0, sumreached = 0, indellength, currentev = -99, j, pos = -1, siteclass;

	    // in core sequence sites: 0 for substitution, 2 for insertion, 4 for deletion
	    // in inserted sites:      1 for substitution, 3 for insertion, 5 for deletion
	    event = chooseevent(rates);

	    if (event != Event::core_sub && event != Event::ins_sub) {   // if NOT a substition....
		// choose indel length

		if (event == Event::core_ins || event == Event::ins_ins) { // insertion
		    indellength         = (*m).insrandomsize((*m).insI, (*m).insD, (*m).insV);
		    insertiontotlength += indellength;	// length of insertion
		} else {          // deletion
		    indellength        = (*m).delrandomsize((*m).delI, (*m).delD, (*m).delV);
		    deletiontotlength += indellength;	//length of deletion
		}

#ifdef printmaxindelsize
		if (maxindelsize < indellength) {
		    cout << "************************************" << endl;
		    cout << "max was " << maxindelsize << " now it is " << indellength << endl;
		    cout << "************************************" << endl;
		    maxindelsize = indellength;
		}
#endif


		int takeoff = 0;
		int inspos  = -1;

		if (event == Event::core_ins) {  // insertion into core sites
		    // find position of insertion.
		    double rand = 1;
		    while (rand == 1 || rand == 0) {
			rand = mtrand1();
			assert(rand != 1 && rand != 0);
		    } // prevents against pseudo-random numbers that are exactly equal to 1 or 0 to machine precision

		    intrand = int((rates.corelength) * rand );  //position of insertion in core sequence
		    currentev = findpos_core_indel((Event) event, updatepositions, intrand, sums, sumreached);

		} else if (event == Event::core_del) { // deletion from core sites

		    // if rates.corelength is 1 and rates.inslength=0 then every position in the core sequence has been deleted.
		    // "core deletions" still have a non-zero positive probability due to the possibility of
		    // deletions occuring before the core sequence and therefore deleting some of the beginning
		    // of the core sequence.  If there is nothing to delete we just ignore this.

		    if ((rates.corelength == 1) && (rates.inslength == 0)) {
			continue;
		    }

		    // e.g. Real sequence has length 100, at positions 1,2,...,99,100 but rates.corelength is 101.
		    //		Position 0 is a "fictitious" position in core sequence to allow "eternal link" for insertions.
		    //	    Let indellength be 5,
		    //		This means rates.corelength-1 +indellength-1 = 104
		    //		So if x is the possible random numbers of int( (rates.corelength-1 +indellength-1)*mtrand1() )
		    //		 then x can be any of 0,1,....,103
		    //		So if intrand is the possible random numbers of int( (rates.corelength-1 +indellength-1)*mtrand1() )+1-indellength
		    //		 then intrand can be any of -4,-3,-2,-1,0,1,2,....,99
		    //		IF intrand<0 then indellength+=intrand, and intrand=0, and inspos=0;
		    //		  i.e. intrand=-4, means *indellength is 1
		    //		  i.e. intrand=-3, means *indellength is 2
		    //		  i.e. intrand=-2, means *indellength is 3
		    //		  i.e. intrand=-1, means *indellength is 4
		    //		  so deletion begins at the beginning of the insertion before "real position" 1,
		    //          and has new length *indellength.
		    //		BUT IF intrand>=0, then we have 0<=intrand<=99
		    //		  and we know that position 0 is not allowed and that real sequence is 1,2,...,100
		    //        so intrand must be increased by 1 to get 1<=intrand<=100

		    double rand = 1;
		    while (rand == 1 || rand == 0) {
			rand = mtrand1();
			assert(rand != 1 && rand != 0);
		    } // prevents against pseudo-random numbers that are exactly equal to 1 or 0 to machine precision
		    oldindellength = indellength;
		    intrand        = int((rates.corelength - 1 + indellength - 1) * rand ) + 1 - indellength; //position of deletion in core sequence

		    if (intrand < 0) {
			indellength += intrand;
			currentev    = 0;
			inspos       = 0;
			updatepositions.assign(11, 0);
		    } else {
			intrand++;
			currentev = findpos_core_indel((Event) event, updatepositions, intrand, sums, sumreached);
		    }
		} else {	// indel at insertion sites
		    assert(event == Event::ins_ins || event == Event::ins_del);

		    double rand = 1;
		    intrand = 0;
		    while (rand == 1) {
			rand    = mtrand1();
			intrand = int((rates.inslength) * rand );
		    }

		    currentev = findpos_ins_indel(event, updatepositions, intrand, sums, sumreached);
		}

		if ((event == Event::ins_ins) || (event == Event::ins_del)) { // indel in inserted sites
		    // find exact position
		    int test = 0;

		    ii = &(insINT2.at(insPOS2.at(currentev)));

		    ss = &((*ii).insertvec);

		    for (j = 0; j < (*ss).size(); j++) {
			s = &((*ss).at(j));

			if ((*s).base != -1) {
			    if (intrand == 0) {
				inspos = j;
				break;
			    } else {
				test++;
			    }
			}
			if (intrand == test + sumreached) {
			    inspos = j;
			    break;
			}
		    }
		    if (inspos == -1) {
			cout << "ERROR 2 IN EVENT " << as_integer(event) << " IN funcnew()" << endl;
			cout << rates.inslength << "  2 " << intrand << "  " << currentev << "  " << inspos << "  " << test << "  " << sumreached << "  " << test + sumreached << endl;

			(*LOG) << "ERROR 2 IN EVENT " << as_integer(event) << " IN funcnew()" << endl;
			(*LOG) << rates.inslength << "  2 " << intrand << "  " << currentev << "  " << inspos << "  " << test << "  " << sumreached << "  " << test + sumreached << endl;
		    }
		}

		// in core sequence sites: 0 for substitution, 2 for insertion, 4 for deletion
		// in inserted sites:      1 for substitution, 3 for insertion, 5 for deletion
		if ((event == Event::core_ins) || (event == Event::ins_ins)) {
		    doinsertion(currenttime, inspos, sums, rates, newseqINT, updatepositions, insINT2, insINT, label, insPOS, insPOS2, indellength);
		} else if ((event == Event::core_del) || (event == Event::ins_del)) {
		    dodeletion(inspos, sums, rates, newseqINT, updatepositions, insINT2, insPOS2, newseqINT.size() /*rates.rootlength*/, indellength);
		}

		sumrates(rates);

	    } else if (event == Event::core_sub) {	// substitution in core sequence

		substitutioncount++;

		double blahrand = 1;
		while (blahrand == 1 || blahrand == 0) {
		    blahrand = mtrand1();
		}                                                                                 // prevents against pseudo-random numbers that are exactly equal to 1 or 0 to machine precision

		unirand = (rates.coresubrate) * blahrand;

		currentev = findpos_core_sub(event, updatepositions, unirand, sums, dsumreached);

		if (currentev == -1) {
		    double yh = 0;
		    for (auto fg : *(sums.CSlevels.rbegin()))
			yh += fg;
		    // I'm not sure we ever visit this line....-csw
		    cout << endl << endl << " ERROR ERROR 0  yh total is " << yh << " as compared to " << rates.coresubrate << " rates and " << unirand << " unirand." << endl << "yh - rates " << yh - rates.coresubrate << " yh-unirand " << yh - unirand << " unirand - rates " << unirand - rates.coresubrate << endl;
		}

		currentbase = newseqINT.at(currentev);

		ratepos = ratevec.at(currentev);

		siteclass = siteclassvec.at(currentev);


		sdiff -= returnsitesubrate(ratepos, siteclass, currentbase);


		chooseNEWbase(currentbase, (((*m).Jvecs).at(siteclass)).at(currentbase));                               // choose new base

		newseqINT.at(currentev) = currentbase;


		sdiff += returnsitesubrate(ratepos, siteclass, currentbase);


		updatesubsums0(updatepositions, sdiff, sums);

		(rates.coresubrate) += sdiff;
		(rates.totalrate) += sdiff;

	    } else if (event == Event::ins_sub) {	//substitution in inserted sites

		substitutioncount++;

		double blahrand = 1;
		while (blahrand == 1 || blahrand == 0) {
		    blahrand = mtrand1();
		}                                                                                 // prevents against pseudo-random numbers that are exactly equal to 1 or 0 to machine precision

		unirand = (rates.inssubrate) * blahrand;

		currentev = findpos_ins_sub(event, updatepositions, unirand, sums, dsumreached);

		if (currentev == -1) {
		    // I have never seen execution reach this point...-csw
		    
		    double yh = 0;
		    for (auto fg : *(sums.ISlevels.rbegin())) 
			yh += fg;
		    cout << endl << endl << " ERROR ERROR 0  yh total is " << yh << " as compared to " << rates.coresubrate << " rates and " << unirand << " unirand." << endl << "yh - rates " << yh - rates.inssubrate << " yh-unirand " << yh - unirand << " unirand - rates " << unirand - rates.inssubrate << endl;
		}

		ii = &(insINT2.at(insPOS2.at(currentev)));

		ss = &((*ii).insertvec);

		int templength = (*ss).size();

		double odsumreached = dsumreached;

		for (j = 0; j < templength; j++) {
		    s = &((*ss).at(j));

		    dsumreached += ((*s).subrate);

		    if (unirand < dsumreached) {
			pos = j;

			currentbase = (*s).base;

			ratepos   = (*s).rate;
			siteclass = (*s).siteclass;

			sdiff -= returnsitesubrate(ratepos, siteclass, currentbase);

			chooseNEWbase(currentbase, (((*m).Jvecs).at(siteclass)).at(currentbase));                               // choose new base

			(*s).subrate = returnsitesubrate(ratepos, siteclass, currentbase);

			sdiff += ((*s).subrate);

			(*s).base = currentbase;

			updatesubsums1(updatepositions, sdiff, sums);

			((*ii).subrate) += sdiff;
			(rates.inssubrate) += sdiff;
			(rates.totalrate) += sdiff;                   //+idiff+ddiff;

			break;
		    }
		}

		if (pos == -1) {
		    badcount++;
		} else {
		    goodcount++;
		}
	    } else {
		cout << "ERROR in choosing event for func" << endl;
	    }

	    exprand = expmyrand(rates.totalrate); 
	}
    }

    //do substitutions



    if (badcount != 0) {
	cout << "ERROR on branch with label " << label << " the goodcount is " << goodcount << " and the badcount is " << badcount << endl; //else cout<<"IT IS OK on branch with label "<<label<<endl;
    }
    if ((model_type == model::Type::nucleotide) && ((*m).modelnumber == 16)) {
	Pt = &matexp;
    } else {
	Pt = &PMatQRev;
    }

    if (!oldmethod) {
	///////////////////////////////////////////////////////////////////////////////////////////////
	// new method substitutions
	//////////////////////////////////////////////////////////////////////////////////////////////////

	if (model_type != model::Type::codon) {
	    if ((*m).continuousgamma) {
		//type is 1 or 2 and it is continuous gamma

		vector<vector<double> > PMat;

		for (int core = 1; core < rates.rootlength; core++) {
		    currentbase = newseqINT.at(core);

		    if ((ratevec.at(core) == 0) || (currentbase == -1)) {
			continue;
		    }


		    PMat = Pt((*m).Qvec, (*m).basefreqs, branchlength * ratevec.at(core));

		    chooseNEWbase(currentbase, PMat.at(currentbase));

		    newseqINT.at(core) = currentbase;
		}

		for (int ins = 0; ins < insINT2.size(); ins++) {
		    ss = &((insINT2.at(ins)).insertvec);

		    for (int j = 0; j < (*ss).size(); j++) {
			s = &((*ss).at(j));

			if (((*s).base != -1) && ((*s).rate != 0)) {
			    double lastleft = (*s).timeleft;

			    if (lastleft == -1) {
				PMat = Pt((*m).Qvec, (*m).basefreqs, branchlength * ((*s).rate));

				chooseNEWbase((*s).base, PMat.at((*s).base));
			    } else {
				PMat = Pt((*m).Qvec, (*m).basefreqs, lastleft * ((*s).rate));

				chooseNEWbase((*s).base, PMat.at((*s).base));

				(*s).timeleft = -1;
			    }
			}
		    }
		}
	    }            //end of: type is 1 or 2 and it is continuous gamma
	    else {
		//type is 1 or 2 and it is not continuous gamma

		vector<vector<vector<double> > > PMat2, PMat, blank;
		vector<vector<double> > blank2;
		blank.assign((*m).numberofsiteclasses, blank2);

		for (int b = 0; b < (*m).numberofsiteclasses; b++) {
		    PMat.push_back(Pt((*m).Qvec, (*m).basefreqs, (((*m).Rrates).at(b)) * branchlength));
		}

		lastlength = -1;

		for (int core = 1; core < rates.rootlength; core++) {
		    currentbase = newseqINT.at(core);

		    if ((ratevec.at(core) == 0) || (currentbase == -1)) {
			continue;
		    }

		    chooseNEWbase(currentbase, (PMat.at(siteclassvec.at(core))).at(currentbase));
		    newseqINT.at(core) = currentbase;
		}

		for (int ins = 0; ins < insINT2.size(); ins++) {
		    ss = &((insINT2.at(ins)).insertvec);

		    for (int j = 0; j < (*ss).size(); j++) {
			s = &((*ss).at(j));

			if (((*s).base != -1) && ((*s).rate != 0)) {
			    if ((*s).timeleft == -1) {
				chooseNEWbase((*s).base, (PMat.at((*s).siteclass)).at((*s).base));
			    } else if ((*s).timeleft == lastlength) {
				if ((PMat2.at((*s).siteclass)).empty()) {
				    PMat2.at((*s).siteclass) = Pt((*m).Qvec, (*m).basefreqs, ((*s).rate) * lastlength);
				}

				chooseNEWbase((*s).base, (PMat2.at((*s).siteclass)).at((*s).base));

				(*s).timeleft = -1;
			    } else {
				lastlength = (*s).timeleft;

				PMat2 = blank;

				PMat2.at((*s).siteclass) = Pt((*m).Qvec, (*m).basefreqs, ((*s).rate) * lastlength);

				chooseNEWbase((*s).base, (PMat2.at((*s).siteclass)).at((*s).base));

				(*s).timeleft = -1;
			    }
			}
		    }
		}
	    }
	} else {
	    vector<vector<vector<double> > > PMat2, PMat, blank;
	    vector<vector<double> > blank2;
	    blank.assign((*m).numberofsiteclasses, blank2);

	    for (int nc = 0; nc < (*m).numberofsiteclasses; nc++) {
		PMat.push_back(Pt(((*m).Qvecs).at(nc), (*m).basefreqs, branchlength));
	    }

	    lastlength = -1;

	    //vector<int> fgh; fgh.assign(5,0);
	    for (int core = 1; core < rates.rootlength; core++) {
		currentbase = newseqINT.at(core);

		if (currentbase == -1) {
		    continue;
		}

		//int ghp=siteclassvec.at(core);  (fgh.at(ghp))++;
		//int ghp=fgh.at(siteclassvec.at(core));  fgh.at(siteclassvec.at(core))=ghp+1;

		chooseNEWbase(currentbase, (PMat.at(siteclassvec.at(core))).at(currentbase));

		newseqINT.at(core) = currentbase;
	    }


	    for (int ins = 0; ins < insINT2.size(); ins++) {
		ss = &((insINT2.at(ins)).insertvec);

//				cout<<"(*ss).size() "<<(*ss).size()<<endl;

		for (int j = 0; j < (*ss).size(); j++) {
		    s = &((*ss).at(j));

		    if ((*s).base != -1) {
			//(fgh.at((*s).siteclass))++;
			if ((*s).timeleft == -1) {
			    chooseNEWbase((*s).base, (PMat.at((*s).siteclass)).at((*s).base));
			} else if ((*s).timeleft == lastlength) {
			    if ((PMat2.at((*s).siteclass)).empty()) {
				PMat2.at((*s).siteclass) = Pt(((*m).Qvecs).at((*s).siteclass), (*m).basefreqs, lastlength);
			    }

			    chooseNEWbase((*s).base, (PMat2.at((*s).siteclass)).at((*s).base));

			    (*s).timeleft = -1;
			} else {
			    lastlength = (*s).timeleft;

			    PMat2 = blank;

			    PMat2.at((*s).siteclass) = Pt(((*m).Qvecs).at((*s).siteclass), (*m).basefreqs, lastlength);

			    chooseNEWbase((*s).base, (PMat2.at((*s).siteclass)).at((*s).base));

			    (*s).timeleft = -1;
			}
		    }
		}
	    }
	}
    }     // end of new method substitutions if

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int findnext(vector<double> currenttime)
{
    double max = 0;
    int maxpos = -1;

    for (int yg = 0; yg < currenttime.size(); yg++) {
	double x = currenttime.at(yg);
	if (x > max) {
	    max    = x;
	    maxpos = yg;
	}
    }
    return maxpos;
}


////////////////////
// took "&" off of insPOS and insINT.  there is problem. on e.g. ((A,B)AB,(C,D)CD);  root to AB is ok, and AB to A is ok, but AB to B is wrong as insPOS and insINT changed

void evolvebranch(vector<vector<int> >& insPOS, /*vector<int> &oldinsPOSin,*/ RATES rates, int parentID, int& branchID, int totalpartitions, int partitionnumber, vector<vector<insert> >& insINT,
                  /* vector<insert> &oldinsIN,*/ int parentlabel, vector<vector<int> >& sequencesINT, string mystring,                 /* vector<int>  &oldseqIN, */
		  int blocknumber, int repnumber, int totalblock, int totalrep, SUMS sums)
{
    // this function performs the main structure of evolving a particular branch.
    // if the branch is a tip it simply evolves the terminal branch length
    // if the branch is internal it evolves to the next node and recursively calls itself on the daughter branches

    int mymode = 0, ev = 0, label, newID = 0;



    int skip     = 0;
    int myadjust = 0;

    double length;

    vector<string> mynewbits;

    if (mystring[0] != '(') {
	mymode = -1;                                                                            // this means it is not a terminal branch e.g. (1:0.1,2:0.2):0.3 instead of A:0.1
    }
    if (mymode == -1) {
	getlabelstring(label, length, mystring);                        // in this case the length would be 0.1 and the label would be "1"
    } else {
	getmynewbits(label, length, mystring, mynewbits);               // in this case the length would be 0.3 and mynewbits would contain "1:0.1" and "2:0.2"
    }
    int siteclass = -1;


    // calculating tree length percent remaining to be evolved
    currenttreelength -= length;

    int percent = int(100 * (treelength - currenttreelength) / treelength);


#ifdef INDELLOG
    indellog << "*********************************************************************************" << endl;
    indellog << " branch label\t" << label << "\tparent:\t" << parentlabel << endl;
    indellog << "*********************************************************************************" << endl;
#endif


    sequencesINT.at(label) = sequencesINT.at(parentlabel);

    insINT.at(label) = insINT.at(parentlabel);

    insPOS.at(label) = insPOS.at(parentlabel);

    /*
     * sequencesINT.at(label)=oldseqIN;
     *
     * insINT.at(label)=oldinsIN;
     *
     * insPOS.at(label)=oldinsPOSin;
     */

    //if(isitbranches) b=&(totalbranches.at( ((*p).mbsposvec).at(partitionnumber) ));

    //cout<<"LABEL "<<label<<"    PARENTLABEL "<<parentlabel<<endl;



    if (isitbranches) {
	//branches model OR  sites-branches model

	int oldoldpos = ((*b).modelpositions).at(branchID);


	branchID++;
	newID = branchID;


	double oldalpha         = (*m).alpha;
	bool oldcontinuousgamma = (*m).continuousgamma;

	//cout<<"old alpha is "<<oldalpha<<endl;


	int oldpos = ((*b).modelpositions).at(parentID);


	int mypos = ((*b).modelpositions).at(branchID);

	//cout<<totalmodels.size()<<"  "<<mypos<<endl;
	m = &(totalmodels.at(mypos));


	if ((oldpos != mypos) || (oldoldpos != mypos)) {
	    changezipfrandoms();

	    if (((*m).continuousgamma || oldcontinuousgamma) && (oldalpha != (*m).alpha)) {
		bool nothing = false;
		if ((*m).alpha == 0) {
		    nothing = true;                                                         //  if(((*m).Rrates).size()>0) nothing=true;
		}
		for (int t1 = 0; t1 < rates.rootlength; t1++) {
		    double rand = mtrand1();

		    if (rand < (*m).pinv) {
			ratevec.at(t1) = 0;
		    } else {
			if (nothing) {
			    ratevec.at(t1) = 1 / (1 - ((*m).pinv));
			} else {
			    ratevec.at(t1) = rndgamma(((*m).alpha)) / (((*m).alpha) * (1 - ((*m).pinv)));
			}
		    }
		}



		vector<insert> *ins = &(insINT.at(label));

		for (int t2 = 0; t2 < (*ins).size(); t2++) {
		    insert *ii = &((*ins).at(t2));

		    vector<site> *s = &((*ii).insertvec);

		    for (int t3 = 0; t3 < (*s).size(); t3++) {
			site *ss = &((*s).at(t3));

			if ((*ss).base >= 0) {
			    double rand = mtrand1();

			    if (rand < (*m).pinv) {
				(*ss).rate = 0;
			    } else {
				if (nothing) {
				    (*ss).rate = 1 / (1 - ((*m).pinv));
				} else {
				    (*ss).rate = rndgamma(((*m).alpha)) / (((*m).alpha) * (1 - ((*m).pinv)));
				}
			    }
			}
		    }
		}
	    } else if ((model_type != model::Type::codon) && (oldalpha != (*m).alpha)) {
		mymods.push_back(mypos);

		for (int t1 = 0; t1 < rates.rootlength; t1++) {
		    ratevec.at(t1) = ((*m).Rrates).at(siteclassvec.at(t1));
		}


		vector<insert> *ins = &(insINT.at(label));


		for (int t2 = 0; t2 < (*ins).size(); t2++) {
		    insert *ii = &((*ins).at(t2));

		    vector<site> *s = &((*ii).insertvec);

		    for (int t3 = 0; t3 < (*s).size(); t3++) {
			site *ss = &((*s).at(t3));

			if ((*ss).base >= 0) {
			    (*ss).rate = ((*m).Rrates).at((*ss).siteclass);
			}
		    }
		}
	    }
	}         //end of oldpos mypos if bracket
    }             //end of isitbranches if bracket
    else {
	branchID = 0;
    }

    if (oldmethod) {
	buildsumsold(rates, sums, sequencesINT.at(label), insINT.at(label), insPOS.at(label));
    } else {
	buildsumsnew(rates, sums, sequencesINT.at(label), insINT.at(label), insPOS.at(label));
    }



    func(branchID, length, rates, sequencesINT.at(label), insINT.at(label), sums, insINT, label, insPOS.at(label), insPOS);



//////////////////////////////////////////////////////////////////////////////////////////////////


    deletionlength = 0;
    // deletionlength is global int.  simply resetting values so that it does not carry over from
    // the end of one sequence to the beginning of the first daughter sequence.



//	string fccc; if(totalpartitions>1) {stringstream fc; fc<<partitionnumber+1; string fcc=fc.str();  fccc="Partition "+fcc+"                         ";}
//	PrintProgress( blocknumber,  totalblock,  repnumber,  totalrep,  percent,fccc);	// print progress, percentage of treelength complete

    if (mymode != -1) {
	// mymode = -1 when terminal branch
	// so this if bracket recursively calls evolvebranch on the daughter branches of the current branch
	for (int asdf3 = 0; asdf3 < mynewbits.size(); asdf3++) {
	    string mynextbit = mynewbits.at(asdf3);
	    evolvebranch(insPOS, /*insPOS.at(label),*/ rates, newID, branchID, totalpartitions, partitionnumber, insINT, /*insINT.at(label),*/ label, sequencesINT, mynextbit,
	                 /*	sequencesINT.at(label),*/blocknumber, repnumber, totalblock, totalrep, sums);
	}
    }

    (mycodes.at(partitionnumber)).at(label) = (*m).geneticcode;
}  // end of evolvebranch function


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int settreeup(int partitionnumber, string& originaltree, string& nonamestree, string& origtreewithnodes, string& nonamestreewithnodes, int& startnodelabel, vector<string>& taxanames, int& ntaxa)
{
    // called at the end of the call control function

    // this function takes the guide tree and creates a version with numbers instead of taxanames
    // and adds interior node labels to a copy of the original tree and to a copy of the numbered tree
    // a taxon name gets replace by "label" and an interior node is labelled N+"label" where "label"
    // is a number corresponding to the position of that nodes' sequence in sequencesINT

    origtreewithnodes = nonamestree = nonamestreewithnodes = "";
    int error   = 0;
    int printon = 0;

    int bracketright = 0, bracketleft = 0;

    originaltree = nowhitespace(originaltree);


    if (printon == 1) {
	cout << originaltree << endl << "dgfweg " << endl << " wegweg";
    }

    error = gettaxanames(partitionnumber, originaltree, taxanames);

    ntaxa = taxanames.size() - 2;
    if (printon == 1) {
	for (int gn1 = 0; gn1 < taxanames.size(); gn1++) {
	    cout << gn1 << " " << taxanames.at(gn1) << endl << " ";
	}
    }
    nonamestree = taxanames.at(taxanames.size() - 1);
    taxanames.pop_back();

    startnodelabel    = taxanames.size() - 1;
    origtreewithnodes = addnodestostring(originaltree, startnodelabel);

    startnodelabel++;

    origtreewithnodes += "ROOT;";

    startnodelabel       = taxanames.size() - 1;
    nonamestreewithnodes = addnodestostring(nonamestree, startnodelabel);

    startnodelabel++;
    nonamestreewithnodes += 'N';
    stringstream sd1;
    sd1 << startnodelabel;
    string fv1 = sd1.str();
    nonamestreewithnodes += fv1;
    nonamestreewithnodes += ';';
    //nonamestreewithnodes+="ROOT;";


    // add interior node labels to the taxa names list
    for (int gn2 = taxanames.size(); gn2 < startnodelabel + 1; gn2++) {
	stringstream gbh;
	gbh << gn2;
	string gbhs   = gbh.str();
	string metaxa = "N" + gbhs;
	taxanames.push_back(gbhs);
    }



    return error;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void sorttreeout(vector<vector<int> >& insPOS, RATES& rates, vector<vector<insert> >& insINT, int mynodelabel, string workingtree, int currentrep, int& ntaxa, vector<string>& taxanames,
		 vector<vector<int> >& sequencesINT, int blocknumber, int repnumber, int totalblock, int totalrep, int partitionnumber, int totalpartitions,
		 SUMS& sums)

{
    // this function simply splits the tree into branches from the root and sends them off to be evolved
    // one by one
    char c             = 'Q';
    int mybracketlevel = 0;
    string abranch;


    int mylim = workingtree.size() - 1;
    int currentbranchnumber = 0;

    for (int yu1 = 1; yu1 < mylim + 1; yu1++) {
	// goes through guide tree and seperates into branches from root
	c = workingtree[yu1];
	if (((c == ',') && (mybracketlevel == 0)) || (yu1 == mylim)) {
	    // evolves a branch from root once parsing of that branch is complete
	    evolvebranch(insPOS, rates, currentbranchnumber, currentbranchnumber, totalpartitions, partitionnumber, insINT,
			 0, sequencesINT, abranch, blocknumber, repnumber, totalblock, totalrep, sums);

	    abranch = "";
	} else {
	    abranch += c;
	}

	if (c == ')') {
	    mybracketlevel--;
	} else if (c == '(') {
	    mybracketlevel++;
	}
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void controlerrorprint(int blocknumber, string instring, int linecount, string myline)
{
    // this functions provides a framework for a standard output of an error to the screen and log file
    // mostly it is just formatting the border of the error box and white space etc

    cout << "\r ERROR occurred in block " << blocknumber << ".  See below for details or consult LOG.txt                       ";

    vector<string> toprint;

    string startline = "ERROR in command line number ";
    stringstream fd1;
    fd1 << linecount;
    string fd1s = fd1.str();
    startline += fd1s;
    startline += " of control file:";

    toprint.push_back(startline);

    string tempstring;
    char c;
    int themaxsize = startline.size();

    for (int j0 = 0; j0 < instring.size(); j0++) {
	c = instring[j0];
	if (c == '\n') {
	    toprint.push_back(tempstring);
	    if (themaxsize < tempstring.size()) {
		themaxsize = tempstring.size();
	    }
	    tempstring = "";
	} else {
	    tempstring += c;
	}
    }
    toprint.push_back(tempstring);
    if (themaxsize < tempstring.size()) {
	themaxsize = tempstring.size();
    }

    string endline = "Last Input was: ";
    endline += myline;
    if (themaxsize < endline.size()) {
	themaxsize = endline.size();
    }
    toprint.push_back(endline);

    cout << endl << endl;
    (*LOG) << endl;

    for (int i0 = 0; i0 < toprint.size(); i0++) {
	string tempstring2 = toprint.at(i0);
	for (int h1 = tempstring2.size(); h1 < themaxsize; h1++) {
	    tempstring2 += " ";
	}
	toprint.at(i0) = tempstring2;
    }

    cout << endl << " +";
    (*LOG) << endl << "+";
    for (int i1 = 0; i1 < themaxsize + 2; i1++) {
	cout << "-";
	(*LOG) << "-";
    }
    cout << "+" << endl;
    (*LOG) << "+" << endl;

    for (int i2 = 0; i2 < toprint.size(); i2++) {
	cout << " | " << toprint.at(i2) << " |" << endl;
	(*LOG) << "| " << toprint.at(i2) << " |" << endl;
    }

    cout << " +";
    (*LOG) << "+";
    for (int i3 = 0; i3 < themaxsize + 2; i3++) {
	cout << "-";
	(*LOG) << "-";
    }
    cout << "+" << endl << endl;
    (*LOG) << "+" << endl << endl;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void printresultsout(int currentrep, string filenamestub, int ntaxa, int mylength, vector<string> taxanames, int blocknumber, int numberofevolveblocks, int reps, int numberofpartitions)
{
//	cout<<endl<<endl<<" Printing leaf sequences"<<endl;


    // this function provides the skeleton structure of printing results to file.
    // the actual translation of integer sequences to alphabetical ones, and the
    // printing to file of these sequences is performed by makeprintseq


    int Mlength         = mylength - numberofpartitions + dashlength;
    int lastlengthcount = Mlength;

    if (printcodonsasDNA && (model_type == model::Type::codon)) {
	Mlength *= 3;
    }

    //L<<Mlength; length=L.str();

    string filename1 = filenamestub + "_TRUE_";         // true alignment
    string filename2 = filenamestub + "_";              // true unaligned sequences
    string filename3 = filenamestub + "_ANCESTRAL_";    // Ancestral sequences
    string endbit    = ".";

    if (fileperrep) {
	stringstream asdd;
	asdd << currentrep;
	string currentrepS = asdd.str();

	// choose file extension depending on desired output type
	if (outputtype == Output::fasta) {
	    endbit += fastaextension;
	} else if (outputtype == Output::phylip) {
	    endbit += phylipextension;
	} else if (outputtype == Output::nexus) {
	    endbit += nexusextension;
	}
	filename1 += currentrepS;
	filename1 += endbit;
	filename2 += currentrepS;
	filename2 += ".";
	filename2 += fastaextension;       //filename2+=endbit;   // unaligned sequences always output in FASTA format
	filename3 += currentrepS;
	filename3 += endbit;

	(*results).close();
	(*results2).close();
	(*results).clear();
	(*results2).clear();

	(results)->open(filename1.c_str());
	(results2)->open(filename2.c_str());

	(*results) << paupstart;
    }


#ifdef INDELLOG
    (*results) << "*********************************************************************************" << endl;
    (*results) << origtreewithnodes << endl;
    (*results) << "*********************************************************************************" << endl;
#endif

    if (outputtype == Output::nexus) {
	// header for nexus file format
	(*results) << "#NEXUS" << endl << endl << "BEGIN DATA;" << endl << "   DIMENSIONS NTAX = " << ntaxa << "  NCHAR = " << fixed << setprecision(0) << Mlength << ";" << "   FORMAT DATATYPE = ";
	if (model_type == model::Type::aminoacid) {
	    (*results) << "PROTEIN";
	} else {
	    (*results) << "DNA";
	}
	(*results) << "  MISSING = ?  GAP = - ;" << endl << "   MATRIX" << endl;
    }

    // header for phylip format
    if (outputtype == Output::phylip) {
	(*results) << ntaxa << "  " << fixed << setprecision(0) << Mlength << "" << endl;
    }
    int printcount = 0;

    //	int seqINTsize=sequencesINT.size();
    //	(*results)<<"     "<<endl; (*results2)<<"     "<<endl;
    for (int i = 1; i < ntaxa + 1; i++) {
	//Leaf sequences (true sequences in (*results2), and true alignment in results)
	(*results2) << ">";
	if (outputtype == Output::fasta) {
	    (*results) << ">";
	} else if (outputtype == Output::nexus) {
	    (*results) << "   ";
	}
	(*results) << taxaspacenames.at(i);
	(*results2) << taxaspacenames.at(i);
	if (outputtype == Output::fasta) {
	    (*results) << "" << endl;
	}
	(*results2) << "" << endl;

	printcount++;


	for (int j = 0; j < TsequencesINT.size(); j++) {
	    makeprintseqLEAF(j, lastlengthcount, (TinsPOS.at(j)).at(i), (TsequencesINT.at(j)).at(i), (TinsINT.at(j)).at(i), i, (*results), (*results2), 0);
	}
	(*results) << "     " << endl;
	(*results2) << "     " << endl;
    }

//cout<<endl;


    if (ancestralprint) {
	if (ancestralfile) {
	    //ofstream results3(filename3.c_str());
	    if (fileperrep) {
		//	(*results3)=new ofstream;

		(*results3).close();
		(*results3).clear();
		(results3)->open(filename3.c_str());
	    }

	    if (!printonlyroot) {//cout<<endl<<" Printing ancestral sequences:"<<endl;
		// This if loop will print the ancestral sequences into results3
		// but only if there is one partition.
		// This is because a multi-partition simulation could have different underlying guide tree topologies
		for (int i = ntaxa + 1; i < (TsequencesINT.at(0)).size() - 1; i++) {
		    //Ancestral Sequences into results3
		    //if(hf2!=seqINTsize-1)
		    //{
		    if (outputtype == Output::fasta) {
			(*results3) << ">";
		    }
		    (*results3) << "N" << i << "\t";
		    if (outputtype == Output::fasta) {
			(*results3) << "" << endl;
		    }

		    printcount++;

		    for (int j = 0; j < TsequencesINT.size(); j++) {
			makeprintseqINT(j, lastlengthcount, (TinsPOS.at(j)).at(i), (TsequencesINT.at(j)).at(i), (TinsINT.at(j)).at(i), i, (*results3), 1);
		    }
		    (*results3) << "     " << endl;

		    //}
		}
	    }

	    //Root Sequnces

	    if (outputtype == Output::fasta) {
		(*results3) << ">";
	    }
	    (*results3) << "ROOT\t";
	    if (outputtype == Output::fasta) {
		(*results3) << "" << endl;
	    }

	    printcount++;

	    for (int j = 0; j < TsequencesINT.size(); j++) {
		makeprintseqINT(j, lastlengthcount, (TinsPOS.at(j)).at(0), (TsequencesINT.at(j)).at(0), (TinsINT.at(j)).at(0), 0, *results3, 1);
	    }
	    (*results3) << "     " << endl;
	} else {
	    // no seperate ancestral file, ancestral sequences go in main alignment file

	    int dsize = (taxaspacenames.at(1)).size();
	    if (!printonlyroot) {
		//	cout<<endl<<" Printing ancestral sequences:"<<endl;
		// This if loop will print the ancestral sequences into (*results3)
		// but only if there is one partition.
		// This is because a multi-partition simulation could have different underlying guide tree topologies
		for (int i = ntaxa + 1; i < (TsequencesINT.at(0)).size() - 1; i++) {
		    //Ancestral Sequences into (*results3)
		    if (outputtype == Output::fasta) {
			(*results) << ">";
		    }
		    stringstream hg;
		    hg << i;
		    string rf = "N" + hg.str();
		    //results<<"N"<<i<<"\t";
		    (*results) << rf;
		    for (int gv = 0; gv < dsize - rf.size(); gv++) {
			(*results) << " ";
		    }
		    if (outputtype == Output::fasta) {
			(*results) << "" << endl;
		    }

		    printcount++;

		    for (int j = 0; j < TsequencesINT.size(); j++) {
			makeprintseqINT(j, lastlengthcount, (TinsPOS.at(j)).at(i), (TsequencesINT.at(j)).at(i), (TinsINT.at(j)).at(i), i, *results, 1);
		    }
		    (*results) << "     " << endl;
		}
	    }
	    //Root Sequnces

	    if (outputtype == Output::fasta) {
		(*results) << ">";
	    }
	    (*results) << "ROOT";
	    for (int gv = 0; gv < dsize - 4; gv++) {
		(*results) << " ";
	    }

	    if (outputtype == Output::fasta) {
		(*results) << "" << endl;
	    }

	    printcount++;

	    for (int j = 0; j < TsequencesINT.size(); j++) {
		makeprintseqINT(j, lastlengthcount, (TinsPOS.at(j)).at(0), (TsequencesINT.at(j)).at(0), (TinsINT.at(j)).at(0), 0, *results, 1);
	    }
	    (*results) << "     " << endl;
	}
    }


    if (outputtype == Output::nexus) {
	(*results) << "   ;" << endl << "END;" << endl;
    }

//		(*results)<<"     "<<endl;(*results2)<<"     "<<endl;

    if (fileperrep) {
	(*results) << paupend;
    } else {
	(*results) << paupmiddle;
    }

    (*results) << "     " << endl;
    (*results2) << "     " << endl;
}


int printfreqs(vector<int>& seq)
{
    // calculates f3x4 frequencies from a sequence
//	cout<<"WHERE2a"<<endl;
    double mysize = double(seq.size());
//	cout<<"WHERE2b"<<endl;
    double counts[3][4] = { { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 } };
//	cout<<"WHERE2c"<<endl;
    int i, y, stopcount = 0;

//	cout<<"WHERE2d"<<endl;
    vector<int> stops = getstops((*m).geneticcode);

    for (y = 1; y < mysize; y++) {
	vector<string> CDU = myCDUletters;

	int seqpos = seq.at(y);

	string cod = CDU.at(seqpos);

	for (int k = 0; k < stops.size(); k++) {
	    if (seqpos == stops.at(k)) {
		if (stopcount < 10) {
		    cout << "STOP CODON " << seqpos << "  " << cod << endl;
		}
		stopcount++;
	    }
	}

	for (i = 0; i < 3; i++) {
	    if (cod[i] == 'T') {
		(counts[i][0])++;
	    }
	    if (cod[i] == 'C') {
		(counts[i][1])++;
	    }
	    if (cod[i] == 'A') {
		(counts[i][2])++;
	    }
	    if (cod[i] == 'G') {
		(counts[i][3])++;
	    }
	}
    }

    cout << "There were " << stopcount << " stop codons in total" << endl;
    for (i = 0; i < 3; i++) {
	for (y = 0; y < 4; y++) {
	    cout << (counts[i][y]) / mysize << "\t";
	}
	cout << endl;
    }
    return 0;
}


void setuprootseq(RATES& rates, vector<int>& rootseqINT, int partitionnumber, vector<vector<int> >& sequencesINT, vector<vector<int> >& insPOS,
		  vector<vector<insert> >& insINT, SUMS& sums)
{
    int mysize = taxanames.size();


    vector<int> blank1;
    blank1.reserve(rates.rootlength);
    vector<insert> blank2;      //blank2.assign(rates.rootlength, insert());

    sequencesINT.clear();
    sequencesINT.assign(mysize, blank1);
    insPOS.clear();
    insPOS.assign(mysize, blank1);
    insINT.clear();
    insINT.assign(mysize, blank2);


    if (partitionnumber == 0) {
	dashlength = 0;                         // reset the count of how many inserted characters have occurred - treewide
    }
    if (rootseqINT.empty()) {   //cout<<"rates.rootlength "<<rates.rootlength<<endl;
	makeseq(rates.rootlength, sequencesINT.at(0));
    } else {
	(sequencesINT.at(0)).push_back(-1);
	(sequencesINT.at(0)).insert((sequencesINT.at(0)).end(), rootseqINT.begin(), rootseqINT.end());
	rootseqINT.clear();
    }

    vector<site> blank;
    (insINT.at(0)).push_back(insert(blank, 0, 0));

    (insPOS.at(0)).assign(rates.rootlength, -1);

/////////////////////
    (insPOS.at(0)).at(0) = 0;
/////////////////////

    if (oldmethod) {
	buildsumsold(rates, sums, sequencesINT.at(0), insINT.at(0), insPOS.at(0));
    } else {
	buildsumsnew(rates, sums, sequencesINT.at(0), insINT.at(0), insPOS.at(0));
    }
}


////////////////////////////////////////////////////////////////////////////////////

void printinsertinfo()
{
    // if indels have occurred print the true information about them for a particular block

    insertiontotlength /= insertioncount;
    deletiontotlength  /= deletioncount;
    double T1 = insertioncount + deletioncount;
    double T2 = T1 + substitutioncount;
    (*LOG) << endl;

    string tss, css, ass, gss;
    format(tss, css, ass, gss, insertioncount / T1, deletioncount / T1, T1 / T2, substitutioncount / T2);

    (*LOG) << "  Actual Insertion : Deletion Ratio  " << tss << " : " << css << "\t(Total indel rate = 1)" << endl;
    if (oldmethod) {
	(*LOG) << "  Actual Indel : Substitution Ratio  " << ass << " : " << gss << "\t(Total event rate = 1)" << endl;
    }
    (*LOG) << "  Actual average insertion length    " << insertiontotlength << endl;
    (*LOG) << "  Actual average deletion length     " << deletiontotlength << endl;
    (*LOG) << "  Number of insertion events         " << insertioncount << endl;
    (*LOG) << "  Number of deletion events          " << deletioncount << endl;
    if (oldmethod) {
	(*LOG) << "  Number of substitution events      " << substitutioncount << endl;
    }
}


////////////////////////////////////////////////////////////////////////////////////



void printoutrates(int& sitecount, int partition, ofstream& ratesout, vector<int>& inspos, vector<insert>& insertstuff)
{
//	ratesout<<"Site\tRate\tInserted?\tPartition"<<endl;

    if (model_type != model::Type::codon) {
	int i;
	//eternal link
	if (inspos.at(0) != -1) {
	    insert *ii = &(insertstuff.at(inspos.at(0)));

	    vector<site> *s = &((*ii).insertvec);

	    for (int j = 0; j < (*s).size(); j++) {
		sitecount++;
		ratesout << sitecount << "\t";
		if ((*m).continuousgamma) {
		} else {
		    ratesout << ((*s).at(j)).siteclass + 1 << "\t";
		}
		ratesout << ((*s).at(j)).rate << "\t" << partition << "\tINSERTION" << endl;
	    }
	}

	vector<string> ratevec2, sitevec2, insYN;
	int maxsize = 0, maxsize2 = 0, diff2, diff;

	for (i = 1; i < ratevec.size(); i++) {
	    if (i != 0) {
		stringstream ds;
		ds << ratevec.at(i);
		string sd = ds.str();
		diff = sd.size();
		if (diff > maxsize) {
		    maxsize = diff;
		}
		ratevec2.push_back(sd);
		insYN.push_back(" ");
		stringstream dq;
		dq << siteclassvec.at(i) + 1;
		string sq = dq.str();
		diff2 = sq.size();
		if (diff2 > maxsize2) {
		    maxsize2 = diff;
		}
		sitevec2.push_back(sq);
	    }

	    if (inspos.at(i) != -1) {
		insert *ii = &(insertstuff.at(inspos.at(i)));

		//	if((*ii).length==0) continue;

		vector<site> *s = &((*ii).insertvec);

		for (int j = 0; j < (*s).size(); j++) {
		    stringstream ds;
		    ds << ((*s).at(j)).rate;
		    string sd = ds.str();
		    diff = sd.size();
		    if (diff > maxsize) {
			maxsize = diff;
		    }
		    ratevec2.push_back(sd);
		    insYN.push_back("INSERTION");
		    stringstream dq;
		    dq << ((*s).at(j)).siteclass + 1;
		    string sq = dq.str();
		    diff = sq.size();
		    if (diff2 > maxsize2) {
			maxsize2 = diff;
		    }
		    sitevec2.push_back(sq);
		}
	    }
	}


	if ((*m).continuousgamma) {
	    for (i = 0; i < ratevec2.size(); i++) {
		sitecount++;
		string s1 = ratevec2.at(i);
		diff = maxsize - s1.size();
		for (int l = 0; l < diff; l++) {
		    s1 += " ";
		}
		string s2 = sitevec2.at(i);
		diff = maxsize2 - s2.size();
		for (int l2 = 0; l2 < diff; l2++) {
		    s2 += " ";
		}

		ratesout << sitecount << "\t" << s1 << "   \t" << partition << "\t" << insYN.at(i) << endl;
	    }
	} else {
	    for (i = 0; i < ratevec2.size(); i++) {
		sitecount++;
		string s1 = ratevec2.at(i);
		diff = maxsize - s1.size();
		for (int l = 0; l < diff; l++) {
		    s1 += " ";
		}
		string s2 = sitevec2.at(i);
		diff = maxsize2 - s2.size();
		for (int l2 = 0; l2 < diff; l2++) {
		    s2 += " ";
		}

		ratesout << sitecount << "\t" << s2 << "\t" << s1 << "   \t" << partition << "\t" << insYN.at(i) << endl;
	    }
	}
    } else {
	vector<double> myrates = (*m).myomegas;
	//eternal link
	if (inspos.at(0) != -1) {
	    insert *ii = &(insertstuff.at(inspos.at(0)));

	    vector<site> *s = &((*ii).insertvec);

	    //for(int j=0; j<(*s).size(); j++) {sitecount++; ratesout<<sitecount<<"\t"<<myrates.at(((*s).at(j)).siteclass)<<"\tY\t"<<partition<<endl;}
	    for (int j = 0; j < (*s).size(); j++) {
		sitecount++;
		ratesout << sitecount << "\t" << ((*s).at(j)).siteclass << "\t" << partition << "\tINSERTION" << endl;
	    }
	}


	for (int i = 1; i < ratevec.size(); i++) {
	    sitecount++;

	    ratesout << sitecount << "\t" << siteclassvec.at(i) << "\t" << partition << endl;

	    if (inspos.at(i) != -1) {
		insert *ii = &(insertstuff.at(inspos.at(i)));

		//	if((*ii).length==0) continue;

		vector<site> *s = &((*ii).insertvec);

		for (int j = 0; j < (*s).size(); j++) {
		    sitecount++;
		    ratesout << sitecount << "\t" << ((*s).at(j)).siteclass << "\t" << partition << "\tINSERTION" << endl;
		}
	    }
	}
    }
}


void usage(const char *program_name)
{
    cerr << "Usage:  " << program_name << " options [ inputfile ... ]\n";
    cerr << "  -h  --help             Display this usage information.\n";
    cerr << "  -s  --seed <value>     integer that seed should start from.\n";
    cerr << "  -v  --verbose          Print verbose messages.\n";
    cerr << "  -V  --version          Print version number and exit.\n";
}


/////////////////////////////////////   MAIN PROGRAM
int verbose = 0;

int main(int argc, char *argv[])
{
    // parse options
    int c;
    int this_option_optind = optind ? optind : 1;
    int option_index       = 0;
    static struct option long_options[] = {
	{ "help",    no_argument,       0, 'h' },
	{ "seed",    required_argument, 0, 's' },
	{ "verbose", no_argument,       0, 'v' },
	{ "version", no_argument,       0, 'V' },
	{         0,                 0, 0,   0 }
    };

    int seed = 0;
    string masterfilename;

    while ((c = getopt_long(argc, argv, "s:vh",
			    long_options, &option_index)) != -1) {
	switch (c) {
	case 'h': // help message
	    usage(argv[0]);
	    exit(-1);

	case 's':
	    mtrand1.seed(atoi(optarg));
	    break;

	case 'v':
	    verbose++;
	    break;

	case 'V':
	    cout << VersionNumber << endl;
	    exit(0);
	    break;

	case '?':
	    if (optopt == 'c') {
		fprintf(stderr, "Option -%c requires an argument.\n", optopt);
	    } else if (isprint(optopt)) {
		fprintf(stderr, "Unknown option `-%c'.\n", optopt);
	    } else {
		fprintf(stderr,
			"Unknown option character `\\x%x'.\n",
			optopt);
	    }
	    return 1;

	default:
	    abort();
	}
    }

    // optional first and second positional parameters.
    // First argument is the control file name (relative to output directory if you give relative path).
    if (argc > optind) {
	masterfilename = argv[optind];
    }

    // Second argument is the output directory.
    if (argc > optind + 1) {
	assert(chdir(argv[optind + 1]) >= 0);
    }

    ofstream treesout(("trees" + simtime + ".txt").c_str());

    delete LOG;
    LOG = new ofstream;

    (*LOG).close();
    (*LOG).clear();

    (LOG)->open(("LOG" + simtime + ".txt").c_str());


    // print some initial output to screen and (*LOG) files and set up some timing factors

    start = clock();
    time_t starttime, endtime, startofblock, endofblock;
    struct tm *timeinfo;
    time(&starttime);
    timeinfo = localtime(&starttime);
    (*LOG) << "********************************************************************************" << endl;
    (*LOG) << "INDELible V" << VersionNumber << " by Will Fletcher : Simulation began at: " << asctime(timeinfo);
    (*LOG) << "********************************************************************************" << endl;
    if (verbose) cout << "INDELible V" << VersionNumber << " by Will Fletcher: Simulation began at " << asctime(timeinfo) << endl ;
    
    int isthereanerror = parse_control_file(masterfilename); // parse and process control file

    if (isthereanerror == -1) {                              // if parsing and processing of control file returns an error then quits.
	delete results;
	delete results2;
	delete results3;
	delete LOG;
	return -1;
    }

    treesout << "Tree file : INDELible V" << VersionNumber << " (best viewed cut and paste in e.g. excel)" << endl << endl;
    treesout << "N.B. Simulation blocks with no random trees will have their trees printed for the first replicate only." << endl;

    if (ancestralprint) {
	treesout << endl << "     Interior node labels on the trees match ancestral sequence labels in sequence files." << endl;
    }

    treesout << endl << "FILE\tTREE\tNTAXA\tREP\tPART\tLENGTH\tDEPTH\tMAX PAIRWISE DISTANCE\tTREE STRING" << endl;

    int numberofevolveblocks = totalevolve.size();


    for (int blocknumber = 1; blocknumber < numberofevolveblocks + 1; blocknumber++) {

	if (verbose) cout << "Processing partition " <<blocknumber << endl;

	// cycles through every block of simulations defined in the control file
	// get settings for this block
	e            = &(totalevolve.at(blocknumber - 1));
	reps         = (*e).reps;
	filenamestub = (*e).filenamestub;
	p            = &(totalpartitions.at((*e).partitionpos));


	ntaxa = (*p).ntaxa;
	bool therearerandomtrees = (*p).randomtrees;
	vector<int> treeposvec   = (*p).treeposvec;


//		vector<bool> isitsitesv=(*p).isitsites;
//		vector<bool> isitrandomv=(*p).isitrandomvec;
	vector<bool> isitbranchesv = (*p).isitbranches;
	int numberofpartitions     = isitbranchesv.size();


	string filename1 = filenamestub + "_TRUE";              // true alignment
	string filename2 = filenamestub + "";                   // true unaligned sequences
	string filename3 = filenamestub + "_ANCESTRAL";         // Ancestral sequences

	// this sets up output files for the case where all replicates are printed to the same file.
	string endbit = ".";

	// choose file extension depending on desired output type
	if (outputtype == Output::fasta) {
	    endbit += fastaextension;
	} else if (outputtype == Output::phylip) {
	    endbit += phylipextension;
	} else if (outputtype == Output::nexus) {
	    endbit += nexusextension;
	}

	filename1 += endbit;
	filename2 += ".";
	filename2 += fastaextension;                  // unaligned sequences always output in FASTA format
	filename3 += endbit;


	if (!fileperrep) {
	    if (ancestralprint && ancestralfile) {
		//	(*results3)=new ofstream;
		(*results3).close();
		(*results3).clear();
		(results3)->open(filename3.c_str());
	    }


	    //	(*results)=new ofstream;
	    (*results).close();
	    (*results2).close();
	    (*results).clear();
	    (*results2).clear();

	    (results)->open(filename1.c_str());

	    (*results) << paupstart;

	    //(*results2)=new ofstream;
	    (results2)->open(filename2.c_str());
	}



	startofblock = clock();

	// Set everything in place
	TsequencesINT.clear();
	sequencesINT.clear();
	TsequencesINT.assign(numberofpartitions, sequencesINT);
	TinsINT.clear();
	insINT.clear();
	TinsINT.assign(numberofpartitions, insINT);
	TinsPOS.clear();
	insPOS.clear();
	TinsPOS.assign(numberofpartitions, insPOS);



	ofstream ratesout;


	if (printrates) {
	    string ratefilename = filenamestub + "_RATES.txt";
	    ratesout.open(ratefilename.c_str());
	    ratesout << "Rates file : INDELible V" << VersionNumber << " " << endl << endl; // (best viewed cut and paste in e.g. excel)"<<endl<<endl;
	}



	for (int thisrep = 1; thisrep < reps + 1; thisrep++) {
	    if (printrates) {
		ratesout << endl << "Replicate: " << thisrep << endl;
	    }

	    globalthisrep = thisrep;

	    // this for loop performs the required number of replicates.

	    //(*LOG)<<endl<<"  * Block "<<blocknumber<<"\tReplicate "<<thisrep<<"\tFilename: "<<filenamestub<<"_TRUE"; if(fileperrep) (*LOG)<<"_"<<thisrep; (*LOG)<<endbit<<endl<<endl;

	    inspinvcount = corepinvcount = 0;

	    if (thisrep == 1) {
		substitutioncount  = 0;         // counts the number of substitution events in a block
		insertioncount     = 0;         // counts the number of insertion events in a block
		insertiontotlength = 0;         // tracks cumulative insertion length in a block
		deletioncount      = 0;         // counts the number of deletion events in a block
		deletiontotlength  = 0;         // tracks cumulative deletion length in a block
	    }



	    int sitecount = 0;

	    ///////////////////////////
	    // OLD POSITION OF PRINT RATES STUFF

	    // totallength = original root sequence length + dashlength ends up being length of sequences in true alignment
	    dashlength = 0;

	    totallength = 0;


	    if (numberofpartitions > 1) {
		(*LOG) << "  This dataset contains " << numberofpartitions << " partitions as follows:" << endl << endl << "\tLength\tFrom\tTo\tModel\tPartition" << endl;
	    }

	    int lastlength = 0;

	    vector<int> partitionlengths;
	    partitionlengths.push_back(0);

	    mycodes.clear();
	    vector<int> bl;
	    bl.assign(ntaxa * 200, 0);
	    mycodes.assign(numberofpartitions, bl);


	    for (int partitionnumber = 0; partitionnumber < numberofpartitions; partitionnumber++) {
		// this loop cycles through the different partitions

		//	if(partitionnumber>0) partitionboundariesL.push_back(partitionboundariesR.at(partitionnumber-1)+1);

		weareon = false;
		//			isitrandom   = isitrandomv.at(partitionnumber);
		//			isitsites    = isitsitesv.at(partitionnumber);

		// whether this partition contains a branch class or a model
		isitbranches = isitbranchesv.at(partitionnumber);

		// guide tree set to become tree from treeposvec that matches the partition
		treeC = &(totaltrees.at(treeposvec.at(partitionnumber)));

		if (therearerandomtrees) {
		    // if it is random tree then the tree in the class (*treeC) is relaced with new tree

		    if ((*treeC).treerandom) {
			(*treeC).newrandom(thisrep);
		    }
		}

		originaltree = (*treeC).tree;

		// set tree length for percentage done calculations ------> not used any more.
		currenttreelength = treelength = (*treeC).treelength;

		// sorts out guide tree, adding node labels for navigation and storage
		startnodelabel = 0;
		int error = settreeup(partitionnumber, originaltree, nonamestree, origtreewithnodes, nonamestreewithnodes, startnodelabel, taxanames, ntaxa);


		if ((thisrep == 1) && (partitionnumber == 0)) {
		    treesout << endl;
		}

		if (therearerandomtrees) {
		    treesout << filenamestub << "\t" << (*treeC).name << "\t" << ntaxa << "\t" << thisrep << "\t" << partitionnumber + 1 << "\t" << (*treeC).treelength << "\t" << (*treeC).treedepth << "\t" << (*treeC).max_distance << "\t";

		    if (ancestralprint) {
			treesout << origtreewithnodes << endl;
		    } else {
			treesout << (*treeC).tree << endl;
		    }
		}

		if (!therearerandomtrees && (thisrep == 1)) {
		    treesout << filenamestub << "\t" << (*treeC).name << "\t" << ntaxa << "\t" << thisrep << "\t" << partitionnumber + 1 << "\t" << (*treeC).treelength << "\t" << (*treeC).treedepth << "\t" << (*treeC).max_distance << "\t";
		    //					treesout<<filenamestub<<"\t"<<(*treeC).name<<"\t"<<ntaxa<<"\t"<<thisrep<<"\t"<<partitionnumber+1<<"\t"<<(*treeC).treelength<<"\t"<<(*treeC).treedepth<<"\t";
		    //				treesout<<filenamestub<<"\t"<<(*treeC).name<<"\t"<<ntaxa<<"\t"<<thisrep<<"\t"<<partitionnumber+1<<"\t"<<(*treeC).treelength<<"\t-\t";

		    if (ancestralprint) {
			treesout << origtreewithnodes << endl;
		    } else {
			treesout << (*treeC).tree << endl;
		    }
		}



		if (error == -1) {
		    return 0;
		}

		// create instances of SUMS and RATES to pass down and inherit through the recursive algorithms
		SUMS sums;
		RATES rates = RATES();

		// rates.partitionlength is set to and remains as user-specified root length for duration of simulation
		rates.partitionlength = ((*p).rootlengthvec).at(partitionnumber) + 1;

		// these others begin as rates.partitionlength but change during course of simulation
		rates.corelength = rates.rootlength = rates.partitionlength;


// these are used when checking that substitution works in inserted sites

		// this is position of model in totalmodels when not a branch class
		// otherwise it is position of branchclass in total branches
		mypos = ((*p).mbsposvec).at(partitionnumber);


		// need to point to root model so that correct stationary frequencies are used for root sequence creation
		if (isitbranches) {
		    b = &(totalbranches.at(mypos));

		    m = &(totalmodels.at(((*b).rootmodelpos)));


		    (mycodes.at(partitionnumber)).at(0) = (*m).geneticcode;
		} else {
		    m = &(totalmodels.at(mypos));
		}

		mymods.clear();
		mymods.push_back(mypos);
		// this is done once here to set up the definitions of constants used in Zipfian random number generation
		changezipfrandoms();


		// set up site classes for codon site-models and relative rates for gamma/pinv distributions

		// weareon=(*m).codonratestrue;						// one used for branch classes and
		// if(!weareon) weareon=(*m).codonratestrue;		// one used for models - not sure which
		SetSiteRates2(rates, thisrep, rootlength);


		// points to root sequence for this partition if defined by user rather than created in setuprootseq
		vector<int> *myrootseq = &(((*p).rootseqints).at(partitionnumber));


		// create root sequence if necessary and setup
		setuprootseq(rates, *myrootseq, partitionnumber, TsequencesINT.at(partitionnumber), TinsPOS.at(partitionnumber), TinsINT.at(partitionnumber), sums);


#ifdef timewatch
		startofevolution = clock();
		(*LOG) << endl << "  * Set-up completed.    Time taken: " << (double)(startofevolution - startofblock) / CLOCKS_PER_SEC << " seconds." << endl;
#endif

		// not currently used?
		noinsertionsyet = true;
		// do the evolving
		sorttreeout(TinsPOS.at(partitionnumber), rates, TinsINT.at(partitionnumber), startnodelabel, nonamestreewithnodes, thisrep, ntaxa,
			    taxanames, TsequencesINT.at(partitionnumber), blocknumber, thisrep, numberofevolveblocks, reps, partitionnumber,
			    numberofpartitions, sums);


#ifdef timewatch
		startofprinting = clock();
		(*LOG) << endl << "  * Evolution completed. Time taken: " << (double)(startofprinting - startofevolution) / CLOCKS_PER_SEC << " seconds." << endl;
#endif

		// updates total length for each partition
		totallength += rates.rootlength;


		//random_shuffle(partitionpositions.begin(), partitionpositions.end()); /* shuffle elements --------> takes too long, so another method used */

		//Tpartitionpositions.push_back(partitionpositions);

		if (numberofpartitions > 1) {
		    int R = totallength + dashlength - partitionnumber - 1, bit = 1;
		    partitionlengths.push_back(R);

		    //if(printcodonsasDNA && type==3) {R*=3;bit=3;}

		    int L = lastlength + bit;
		    lastlength = R;

		    (*LOG) << "\t" << R - L + 1 << "\t" << L << "\t" << R << "\t";
		    if (isitbranches) {
			(*LOG) << (*b).name;
		    } else {
			(*LOG) << (*m).name;
		    }
		    (*LOG) << "\t" << partitionnumber + 1 << endl;
		}


		if (printrates) {
		    if (model_type != model::Type::codon) {
			if ((*m).continuousgamma && isitbranches) {
			    ratesout << endl << "This partition uses continuous gamma rate variation and a non-homogenous model." << endl;
			    ratesout << "The rates are not currently logged if alpha has changed on the tree. " << endl;
			    ratesout << "The rates given below correspond to alpha=" << (*m).alpha << " for model " << (*m).name << " only" << endl << endl;
			    ratesout << "Site\tRate\tPartition\tInserted?" << endl;

			    //ratesout<<endl<<"Partition "<<partitionnumber+1<<" uses  Rates are not logged for continuous gamma that changes over the tree. "<<endl<<"Please use discrete gamma if you need to know the specific rates for each site in a changing gamma model. "<<endl<<"Continuous gamma rates *are* logged for models where alpha does not change. Other parameters of the model may still change. "<<endl<<"If this is feature is a necessity for your work then use the [printallrates] command to log rates\nfor non-homogenous continuous gamma models, albeit in a rather crude way."<<endl;
			} else {
			    if (isitbranches) {
				ratesout << endl << "  N.B. For discrete gamma models that change over the tree, individual rates on" << endl;
				ratesout << "       different branches can be calculated from the rates table below. " << endl;
				ratesout << "       Just remember that the site *class* does not change over the tree." << endl << endl;

				int iii, jjj;
				m = &(totalmodels.at(mymods.at(0)));
				int gg = ((*m).Rrates).size();

				ratesout << "The columns below are the rates for the different classes under the different models." << endl;
				ratesout << "The columns correspond to the following values of alpha used in the simulation:" << endl << endl;

				for (iii = 0; iii < mymods.size(); iii++) {
				    m = &(totalmodels.at(mymods.at(iii)));
				    ratesout << (*m).alpha << "    ";
				}


				ratesout << endl << endl << "Class" << endl;

				vector<string> vs;
				vector<vector<string> > vvs;
				int maxsize = 0;
				for (jjj = 0; jjj < gg; jjj++) {
				    vs.clear();
				    for (iii = 0; iii < mymods.size(); iii++) {
					m = &(totalmodels.at(mymods.at(iii)));

					vector<double> rrates = (*m).Rrates;
					stringstream jh;
					jh << rrates.at(jjj);
					string fv = jh.str();
					vs.push_back(fv);
					if (maxsize < fv.size()) {
					    maxsize = fv.size();
					}
				    }
				    vvs.push_back(vs);
				}

				for (jjj = 0; jjj < gg; jjj++) {
				    ratesout << jjj + 1 << "\t";
				    vs = vvs.at(jjj);

				    for (iii = 0; iii < mymods.size(); iii++) {
					string za = vs.at(iii);
					ratesout << za;
					for (int dc = 0; dc < maxsize - za.size(); dc++) {
					    ratesout << " ";
					}
					ratesout << "\t";
				    }
				    ratesout << endl;
				}
				m = &(totalmodels.at(((*b).rootmodelpos)));
			    } else {
				ratesout << endl << "alpha was " << (*m).alpha << " in this partition.";
			    }

			    if ((*m).continuousgamma) {
				ratesout << endl << endl << "Site\tRate\tPartition\tInserted?" << endl;
			    } else {
				ratesout << endl << endl << "Site\tClass\tRate\tPartition\tInserted?" << endl;
			    }
			}
		    } else if (partitionnumber == 0) {
			ratesout << "  N.B. Site classes are numbered from lowest to highest values of omega." << endl;
			ratesout << "       The omegas are not given explicitly as they are permitted to change on different" << endl;
			ratesout << "       branches whereas the site classes are not allowed to change." << endl;
			ratesout << "       Please consult your control file to find the corresponding omega values" << endl << endl;
			ratesout << "Site\tClass\tPartition\tInserted?" << endl;
		    }

		    printoutrates(sitecount, partitionnumber + 1, ratesout, (TinsPOS.at(partitionnumber)).at(0), (TinsINT.at(partitionnumber)).at(0));
		}
	    } // end of for loop for partitions in a replicate


	    // print results to output files
	    printresultsout(thisrep, filenamestub, ntaxa, totallength, taxanames, blocknumber, numberofevolveblocks, reps, numberofpartitions);

#ifdef timewatch
	    endofprinting = clock();
	    (*LOG) << "  * Printing completed.  Time taken: " << (double)(endofprinting - startofprinting) / CLOCKS_PER_SEC << " seconds." << endl;
#endif



#ifdef splittree
	    // this compiler option is designed to split a 24 taxa tree into three 8 taxa trees at output to check that different parts of tree can evolve with different models.
	    // e.g. a 24 taxa tree where the three 8 taxa subtrees could be
	    // 1) GTR + continuous gamma with alpha = 0.5 + base frequencies 0.1 0.2 0.3 0.4
	    // 2) HKY + discrete gamma with alpha 2 + base frequencies 0.4 0.3 0.2 0.1
	    // 3) UNREST

	    if (ntaxa != 24) {
		cout << " WARNING: splittree compiler option only works on 24 taxa trees. No split done." << endl << endl;
		continue;
	    }

	    ofstream waste;

	    if (!fileperrep && (reps != 1)) {
		cout << "ERROR: splittree debug compile code is not meant to be used with fileperrep" << endl << "Please change the setting and re-try analysis" << endl << endl;
		return -1;
	    }

	    if (outputtype != 2) {
		cout << "ERROR: splittree debug compile code is only meant to be used with phylip format" << endl << endl;
		return -1;
	    }


	    stringstream kk;
	    kk << thisrep;
	    string gg = kk.str();

	    string filename = filenamestub + "_TRUE_" + gg;

	    string endbit = ".";
	    endbit += phylipextension;

	    vector<string> myletters;

	    if (model_type == model::Type::nucleotide) {
		myletters = myDUletters;
	    } else if (model_type == model::Type::aminoacid) {
		myletters = myAUletters;
	    } else if (model_type == model::Type::codon) {
		myletters = myCDUletters;
	    }


	    for (int a = 1; a < 4; a++) {
		stringstream dw;
		dw << a;
		string dd = dw.str();
		int myextra;

		ofstream iff;
		iff.open((filename + "_subtree" + dd + endbit).c_str());

		iff << 8 << "   " << totallength + dashlength - numberofpartitions << endl;

		if (a == 1) {
		    myextra = 0;
		} else if (a == 2) {
		    myextra = 8;
		} else {
		    myextra = 16;
		}

		for (int nx = 1 + myextra; nx < 9 + myextra; nx++) {
		    iff << taxaspacenames.at(nx);

		    for (int px = 0; px < numberofpartitions; px++) {
			makeprintseqLEAF(0, totallength + dashlength - numberofpartitions, (TinsPOS.at(px)).at(nx), (TsequencesINT.at(px)).at(nx), (TinsINT.at(px)).at(nx), nx, iff, waste, 0);
		    }

		    iff << endl;
		}
	    }
#endif      // end-ifdef splittree


////////////////////////////////////////////////
#ifdef rippartitions
	    // this is designed to rip the partitions from a partitioned simulated dataset and spit them in to separate files
	    // that can be subjected to separate analyses to check that the partition system works ok.  number of partitions unlimited for test.

	    if (numberofpartitions == 1) {
		continue;
	    }

	    ofstream waste;

	    cout << " WARNING: deletion rate > 0 in any model in any partition will probably" << endl << "          cause the rippartitions debug code to crash." << endl << endl;

	    if (!fileperrep) {
		cout << "ERROR: rippartitions debug compile code is not meant to be used with fileperrep" << endl << "Please delete the file paupmiddle.txt and re-try analysis" << endl << endl;
		return -1;
	    }

	    if (outputtype != 2) {
		cout << "ERROR: rippartitions debug compile code is only meant to be used with phylip format" << endl << endl;
		return -1;
	    }


	    stringstream kk;
	    kk << thisrep;
	    string gg = kk.str();

	    string filename = filenamestub + "_TRUE_" + gg;

	    string endbit = ".";
	    endbit += phylipextension;

	    vector<string> myletters;

	    if (model_type == model::Type::nucleotide) {
		myletters = myDUletters;
	    } else if (model_type == model::Type::aminoacid) {
		myletters = myAUletters;
	    } else if (model_type == model::Type::codon) {
		myletters = myCDUletters;
	    }



	    for (int px = 0; px < partitionlengths.size() - 1; px++) {
		stringstream b;
		b << px + 1;
		string bb = b.str();

		ofstream iff;
		iff.open((filename + "_part" + bb + endbit).c_str());

		int L = partitionlengths.at(px), R = partitionlengths.at(px + 1), partlength = R - L;
		if (model_type == model::Type::codon) {
		    partlength *= 3;
		}

		iff << ntaxa << "   " << partlength << endl;

		for (int nx = 1; nx < ntaxa + 1; nx++) {
		    iff << taxaspacenames.at(nx);

		    makeprintseqLEAF(0, partlength, (TinsPOS.at(px)).at(nx), (TsequencesINT.at(px)).at(nx), (TinsINT.at(px)).at(nx), nx, iff, waste, 0);

		    iff << endl;
		}
	    }
#endif      // end-ifdef rippartitions
//////////////////////////////////////////////////////////////////
	} // end of for loop for replicates in a block

	if (printrates) {
	    ratesout.close();
	}


	//prints out paupend if there is one file for all replicates.
	if (!fileperrep) {
	    (*results) << paupend << endl;
	}

	// if indels have occurred print the true information about them for this block
	if (((*m).insertrate != 0) || ((*m).deleterate != 0)) {
	    printinsertinfo();
	}
	//	cout<<clock()<<endl;

	endofblock = clock();
	(*LOG) << endl;
	cout << "  * Block " << blocknumber << " completed.   Time taken: " << (double)(endofblock - startofblock) / CLOCKS_PER_SEC << " seconds." << endl;
	(*LOG) << "  * Block " << blocknumber << " was completed in " << (double)(endofblock - startofblock) / CLOCKS_PER_SEC << " seconds." << endl;


	(*LOG) << endl << "********************************************************************************" << endl << endl;
    } // end of commandblocks for loop, i.e. blocks in a control file

    (*results).close();
    (*results2).close();
    (*results3).close();

    // print out final screen and log output and calculate simulation finish time etc

    finish = clock();
    time(&endtime);
    timeinfo = localtime(&endtime);

    if (numberofevolveblocks > 1) {
	cout << "  * All blocks complete. Time taken: " << (double)(finish - start) / CLOCKS_PER_SEC << " seconds." << endl << endl;
    }
    cout << "\n\n *** SIMULATION COMPLETED - PLEASE CONSULT OUTPUT FILES ***                                                                     " << endl;


    (*LOG) << "  * Simulation completed. Whole batch took: " << (double)(finish - start) / CLOCKS_PER_SEC << " seconds." << endl << endl;

    (*LOG) << "INDELible V" << VersionNumber << " Simulations completed at: " << asctime(timeinfo);

    (*LOG) << endl << "********************************************************************************" << endl << endl;

    cout << endl << endl << " INDELible V" << VersionNumber << " Simulations completed at: " << asctime(timeinfo) << endl << endl;


    (*LOG) << endl << " Original Control File " << endl;
    (*LOG) << endl << "-----------------------" << endl << endl << endl;


    for (int bv = 0; bv < originalcontrol.size() - 1; bv++) {
	(*LOG) << originalcontrol.at(bv);
    }

    delete results;
    delete results2;
    delete results3;
    delete LOG;


    return 0;
}
