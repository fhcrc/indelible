#ifndef _models_h_
#define _models_h_ 1

#include <vector>
#include <string>

class indelmodel
{
	// all indel models conform to these parameters.
	// not all are needed for each model
	// but they must all be present for the pointer to the functions below to work.
public:
	double meansize;
	int type;
	int r;
	double q;
	double a;
	double b;
	int M;

	std::vector<double> usermodel;

	indelmodel()
	{
		//these settings chosen to cause a crash when indel model parameters not specified.
		type=r=M=-1;
		meansize=q=a=b=-1;
	}
};

enum ModelType { nucleotide=1, aminoacid=2, codon=3 };


class model
{
	// this BIG class contains all details of substitution/indel model and all construtors etc
public:

	// parameters for indel model
	double q1,  q2,  Hx0,  s,  Himax;

	double insD;
	int insI;
	std::vector<double> insV;

	double delD;
	int delI;
	std::vector<double> delV;

	double delmeansize;		// used in formula for total deletion rate

	// pointer for indel size generation
	int (*delrandomsize)(int, double, std::vector<double>&) ;
	int (*insrandomsize)(int, double, std::vector<double>&) ;

	// test function pointers
	int (*pf)();
	int (*pf2)(double, int);

	// model name
	std::string name;

	int modelnumber;	//?
	int geneticcode;

	bool copiedmodel;			// whether the model is copied - defunct?
	bool codonratestrue;		// whether model has different codon position specific relative rates

//	std::vector<double> insertrates;
//	std::vector<double> deleterates;

	double insertrate;			// relative instantaneous rate of insertion
	double deleterate;			// relative instantaneous rate of deletion
	double indelrate;			// relative instantaneous rate of insertion and deletion
	double subrate;				// relative instantaneous rate of substitution (=1)

	int modelpos;				// position in totalmodels
	ModelType type;				// nucleotide=1, aminoacid=2, codon=3
	int error;
	int rootlength;				// ?
	double alpha;				// alpha for gamma models
	double pinv;
	int ngamcat;				// number of gamma categories for discrete gamma

	double codonrates[3];		// relative substitution rates for codon positions

	bool continuousgamma;		// whether this model has continuous gamma rate heterogeneity

	int numberofsiteclasses;	// number of classes in codon sites model, or number of gamma categories too..

	int medianORmean;			// 1 = use medians, 0 = use means, to represent categories in discrete gamma rate variation

	std::vector<double> cumfreqs;	// cumulative frequencies for discrete gamma, or codon sites models
	std::vector<double> Rrates;		// relative rates for discrete gamma categories
	std::vector<double> myomegas;	// different omegas for different site classes.

	std::vector<double> rootbasefreqs;	// base frequencies used in root sequence creation if model at root
	std::vector<double> basefreqs;		// stationary frequencies of model
	std::vector<double> insertfreqs;		// frequencies used to generateinserted sequences
	std::vector<double> myrates;			// std::vector of substitution rates "away" from a given state, used in method 2 jump chain

	std::vector<std::vector<double> > myratesvec;		//std::vector of different "myrates" std::vectors for different site classes in codon site models

	std::vector<int> ratevec;			// not used? ---> done in main skeleton file now.

	std::vector<std::vector<double> > Jvec;	// transition matrix of the jump chain
	std::vector<std::vector<double> > Qvec;	// transition probabilities

	std::vector<double> scalefactors;	// used when scaling different Qvec from different site classes

	std::vector<std::vector<std::vector<double> > > Jvecs;	///////   collection of Jvecs from different site classes
	std::vector<std::vector<std::vector<double> > > Qvecs;	///////   collection of Qvecs from different site classes

	model(int mymodelpos, ModelType mytype, const std::string &myname, int mymodelnumber, int mygeneticode,
	      bool mycopiedmodel, double myinsertrate, double mydeleterate, double myalpha, double mypinv,
	      int myngamcat, double mycodonrates[], const std::vector<double> &mybasefreqs, const std::vector<double> &myrootbasefreqs,
	      const std::vector<double> &myinsertfreqs, std::vector<double> &myparams, const std::vector<double> &aamodel, const indelmodel &insertmodel,
	      const indelmodel &deletemodel);
	void changeQandJ(int numcats);

private:

	void testmyfreqs(std::vector<double> &basefreqs, std::string mycommand);
	void d(std::vector<std::vector<double> > &Q, double S);
	std::vector<double> fx(const std::vector<double> &basefreqs, int which);
	void makeequalfreqs(const ModelType type, std::vector<double> &tbasefreqs);
	void getJvec(double S, std::string name, std::vector<double> &myrates, std::vector<std::vector<double> > Qvec, std::vector<std::vector<double> > &Jvec, std::vector<double> &basefreqs);
	std::vector<std::vector<double> > getDNA( std::string name, std::vector<double> nstnums, std::vector<double> &basefreqs, int mymodel);
	std::vector<std::vector<double> > getAA( std::string name, std::vector<double> params, std::vector<double> &basefreqs, int modelnumber, std::vector<double> aamodel);
	std::vector<std::vector<double> > getECMr();
	std::vector<std::vector<double> > getECMu();
	std::vector<std::vector<double> > getCOD(std::string &name, std::vector<double> &basefreqs, int mymodel, double kappa, double omega);

}; // end class model


extern std::vector<int> allowedcodes(int gencode);
extern std::vector<int> getstops(int geneticcode);

extern ModelType model_type;

extern const char GeneticCodeTable[][65];

#endif  // _models__h_ 
