#include <vector>
#include <string>
#include <fstream>      // std::ofstream
#include "MersenneTwister.h"
#include "models.h"

class siteclass
{
   // site class is now disabled.  It allows the use of different models in different proportions along a sequence.

public:

   std::vector<int>                  allmodelsused;
   std::vector<bool>                 trueformodel; //true for model, false for branch class, in mbnames, mbpositions
   std::vector<int>                  mbpositions;
   std::vector<std::string>          mbnames;
   std::vector<double>               props; // cumulative props remember!
   std::vector<std::vector<double> > trootbasefreqs;

   int         geneticcode;
   int         error;
   std::string mbtree;
   std::string name;

   bool therearemodels;
   bool therearebranches;

//	siteclasscopy(site copysite)
//	{
//		//??????
//	}

//	createsiteprops()
//	{
//		//????
//	}

   siteclass(std::string sitename, std::vector<std::string>& mynames, std::vector<double>& myprops, std::vector<std::string>& totalmodelnames, std::vector<std::string>& totalbranchnames)
   {
      error            = 1;
      therearemodels   = false;
      therearebranches = false;
      std::vector<std::string> mbtrees;

      name    = sitename;
      props   = myprops;
      mbnames = mynames;

      if (error != -1) { error = sortoutnames(mbtrees, totalmodelnames, totalbranchnames); }
      if (error != -1) { error = sortouttrees(mbtrees); }
      if (error != -1) { error = testgeneticcode(); }
   }

///////////////////////////////////
private:

   int testgeneticcode();

   int sortoutnames(std::vector<std::string>& mbtrees, std::vector<std::string>& totalmodelnames, std::vector<std::string>& totalbranchnames);

   int sortouttrees(std::vector<std::string>& mbtrees);
};  // end class siteclass



class branchclass
{
   // this class provides the framework for simulations that use different models on different branches in the guide tree.

public:
   std::string              treewithrootmodel; // original tree given in control file, has a model at the root.
   std::string              tree;              // tree minus model at root.
   std::string              baretree;          // tree expressed as only a pattern of parentheses and commas
   std::string              name;              // name of branch class
   std::vector<std::string> modelnames;        // names of models used, given in same order as modelpositions
   //std::vector<model> branchmodels;
   std::vector<int>    modelpositions;         // position of models used in totalmodels.
   std::vector<int>    allmodelsused;          // list of all models used in branch class (listed by model number)
   std::vector<double> rootbasefreqs;          // base frequencies of model defined at root of branch class tree.
   double              insertrate;             // insertion rate at root?
   double              deleterate;             // deletion rate at root?
   int                 rootmodelpos;           // position of root model in totalmodels.

   bool geneticcodefixed;                      // whether to allow geneticcode to change on trees.
   int  geneticcode;                           // genetic code used (for codon simulations) must be the same across entire tree
   int  error;                                 // whether branch class successfully set itself up without errors.

   int numcats;                                // number of categories in discrete gamma model or in site/branch-site models
   //    must be the same across the whole tree.

   bool changecats;                        // says whether to rescale root cats and therefore check all for change.

   std::vector<double> catprops;           // relative proportion of categories (must be same in different models).

   branchclass(std::string& ttree, std::string& myname, std::vector<std::string>& totalmodelnames, bool iscodefixed)
   {
      // constructor
      geneticcodefixed = iscodefixed;

      error             = 1;
      treewithrootmodel = ttree;
      for (int fd = 0; fd < treewithrootmodel.size(); fd++) { char c = treewithrootmodel[fd]; if ((c == '(') || (c == ')') || (c == ',') || (c == ';')) { baretree += c; }
      }
      name = myname;

      if (error != -1) { error = getrootfreqs(treewithrootmodel); }
      if (error != -1) { error = buildbranches(tree, myname, modelnames); }

      if (error != -1) {
         error = getbranchmodels(modelnames, totalmodelnames, modelpositions);          //branchmodels);
      }
      if (error != -1 /*&& geneticcodefixed*/) { error = testgeneticcode2(modelpositions); }
   }

private:


   int getrootfreqs(std::string& testtree2);

///////////////////////////

   int testgeneticcode(std::vector<int> modelpositions);

///////////////////////////

   int testgeneticcode2(std::vector<int> modelpositions);

///////////////////////////
   int getbranchmodels(std::vector<std::string>& modelnames, std::vector<std::string>& totalmodelnames, std::vector<int>& modelpositions); //std::vector<model> &branchmodels);


   ///////////////

   int buildbranches(std::string& testtree, std::string& branchesname, std::vector<std::string>& modelnames);
};  // end class branchclass



class partitionclass
{
   // this class stores the information for simulations that use different partitions
   // A partition may be a branch class, site class or model class.


public:
   std::string name;
   int         ntaxa;
   bool        randomtrees;

   std::vector<int>               geneticcodes;
   std::vector<std::vector<int> > rootseqints;
   std::vector<int>               blank;
   std::vector<int>               mbsposvec;
   std::vector<int>               rootmodelpos;
   std::vector<int>               rootlengthvec;
//	std::vector<int> ntaxavec;
   std::vector<int> treeposvec;

   std::vector<bool> isitsites;
   std::vector<bool> isitbranches;
//	std::vector<bool> isitrandomvec;



   partitionclass(std::string myname, std::vector<int> rootlengths, std::vector<std::string>& rootseqtxtvec, std::vector<std::string> rootfilenames,
                  std::vector<int> mmbsposvec, std::vector<std::string> mbstypes, std::vector<std::string> mbnames, std::vector<int> mrootmodelpos,
                  std::vector<int> geneticcodevec, std::vector<int> mytreeposvec, std::vector<std::string> mytreenamevec, int myntaxa, bool isitrandom);


private:

   int checktaxaintrees(std::vector<int> posvec);

   int makerootseqints(std::vector<int>& rootseqint,
		       int rootlength,
		       const std::string& rootseqtxt,
		       int mbspos,
		       const std::string& rootfilename,
		       const std::string& mbstype,
		       const std::string& mbname,
		       int geneticcode);
};  // end-class partitionclass


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class evolve
{
   // real simple class!
   // this is the one that links it all together.
   // Each instance of the class evolve is a standalone simulation of a partition that contains a branch class, or a model or so on....

public:
   int         reps;
   int         partitionpos;
   std::string filenamestub;

   evolve(int myreps, int mypartitionpos, std::string myfilenamestub)
   {
      reps         = myreps;         // number of replicates
      partitionpos = mypartitionpos; // position in totalpartitions of the partition in question
      filenamestub = myfilenamestub; // filenamestub which is used to create output filenames
   }
};                                   // end-class evolve



class Tree
{
   // this class stores a tree as defined by the user or creates random trees based on user defined settings...

public:
   bool treerandom;
//	bool randomeveryrep;
   bool        seedgiven;
   std::string tree;
   std::string name;
   int         ntaxa, seed, guidetreeroot, last;
   double      scaler, scalerd, scalerm, birth, death, sample, mut, max_distance;
   double      treelength, treedepth;
   MTRand      mymtrand;
   std::string doctorlengths;

   Tree(bool myseedgiven, std::string mtree, double mscaler, double mscalerd, int mntaxa, double mbirth, double mdeath, double msample, double mmut, int mseed, int mguidetreeroot, std::string mname, std::string mydoctor, double mscalerm);

   ////////////
   void newrandom(int rep); //, unsigned long &last);

private:
   ////////////


   double correctlastdepth(double maxdepth, std::string originaltree, double depthsofar, std::vector<std::string>& originals, std::vector<std::string>& changes);

   std::string randomise_ultrametric(std::string originaltree);

   double getmaxtreedepth(std::string originaltree);

   double getnextdepth(std::string originaltree, double depthsofar);

   ////////////
   double gettreelength(std::string originaltree);

   std::string rescaledtree(std::string originaltree, double treelength, double newtreelength);
};  // end-class Tree


extern std::ofstream *LOG;               // used to print out information from simulations
extern bool          fixtrueproportions; //  if true then the proportions of codon sites models will be true rather than random.
extern bool          oldmethod;          // false for method 1, true for waiting time method
extern bool          phylipnametruncate; // truncates taxa names to 10 characters to comply with PHYLIP format

extern bool globalseed;                  // whether one seed sets both random trees, and sequence generation.

extern bool insertaslowercase;           // prints inserted sites as lower case and core sites as uppercase

extern bool fixtrueproportions;          //  if true then the proportions of codon sites models will be true rather than random.

extern bool markdeletedinsertions;       //
extern bool fileperrep;                  // prints all replicates in same file if false, otherwise in their own files.

extern bool printrates;                  // whether to print out relative rate information or not.
extern bool printallrates;               // whether to print out relative rate information or not for every branch

extern bool printcodonsasDNA;            // whether to print codons as DNA bases or not (true for DNA/false for amino acids)


extern bool printinsert;            // whether to print indel statistics
extern bool ancestralprint;         // whether to print ancestral sequences to another file
extern bool ancestralfile;          // whether ancestral sequences go in same or different file.  true for different file

extern bool phylipnametruncate;     // truncates taxa names to 10 characters to comply with PHYLIP format
extern bool guidetreebinary;        // whether guide tree contains polytomies
extern bool guidetreerooted;        // whether guide tree is unrooted or rooted

extern std::string phylipextension; // used for output file filename extension
extern std::string fastaextension;  // used for output file filename extension
extern std::string nexusextension;  // used for output file filename extension

enum OutputType {
   OUTPUT_FASTA=1, OUTPUT_PHYLIP=2, OUTPUT_NEXUS=3
};

extern int outputtype;          // 1 for FASTA, 2 for PHYLIP, 3 for NEXUS
extern int guidetreetype;       // used for logfile

extern std::vector<model> totalmodels;

extern std::string paupstart, paupmiddle, paupend;

extern std::string                 simtime;         //disabled time stamp so that output is over-written.
extern std::vector<evolve>         totalevolve;     // storage for the evolve class instances.
extern std::vector<partitionclass> totalpartitions; // storage for partitions
extern std::vector<Tree>           totaltrees;      // storage for trees!
extern std::vector<char>           originalcontrol; // used to store control file used.
extern std::vector<branchclass>    totalbranches;   // storage for branch classes



extern void controlerrorprint2(std::string blocktype, std::string blockname, std::string commandname, std::string instring, std::string myline);
extern int parse_control_file(const std::string& filename);
