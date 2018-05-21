//-----------------------------------------------------------------------------
// prnnet.c
//
// Network of Pinsky-Rinzel Two-Chamber Neurons
//   Current Options:
//     Pinsky-Rinzel Single Neuron ( -pr )
//     Nearest-neighbor synaptic network
//     Extracellular resistive lattice
//     Applied electric field - more than one field per program execution  (-ef)
//     Applied current pulses - multiple pulses, always 0.5 sec in duration(-ip)
//     Randomization of parameters and initial values ( -rand )
//     Spike propagation analysis ( -p )
//
// Origninal comments (1) moved to end of this file.
//
// Revisions:
// Ver  Date          Description
// ---- ---------- ----------------------------------------------------- ------
// 1.00  2005/02/11 Cloned and modified from gc11.cpp.                   C.Fan
// 2.00  2006/03/11 Implemented -nr option as a chain of neurons with    C.Fan
//                  vertical extracellular resistance.
// 2.01  2006/04/11 Fixed equation for V2.0 and added the reading of     C.Fan
//                  initial conditions from input file.
// 2.02  2006/06/01 Added -oeq option to output last variable values.    C.Fan
//                  Used of logging equilibrium values.
// 2.03  2006/06/07 Fixed applied field onset.                           C.Fan
// 2.04  2006/06/14 Log terminal neuron as a function of applied field.  C.Fan
// 2.05  2006/08/29 Implement capability of 2 sets of neurons with       C.Fan
//                  with distinct applied potentials.
// 2.06 2006/10/15  Allow initializing chain for resistive grid.         C.Fan
//-----------------------------------------------------------------------------

#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

//
// Offset to an element (e) in a structure (s)
//
#define OFFSET(s,e)   ( (int )( &(((s *)0)->e) ) )

//
#define  PI           3.1415926535897932385
#define  TWOPI        6.283185307179586477
#define  SMALL        0.00000001
#define  TINY         (1.0*pow(10.0,-100.0));   // had been ,-100.0 until 11/15/2

//
// Index into matrix dy[num_neurons][] - integration variables, for each neuron
//
#define  GATE_H   0             // h gate
#define  GATE_N   1             // n gate
#define  GATE_S   2             // s gate
#define  GATE_C   3             // c gate
#define  GATE_Q   4             // q gate
#define  CONC_CA  5             // calcium concentration
#define  I_SNMDA  6             // NMDA weighting 
#define  I_WAMPA  7             // AMPA weighting
#define  V_SOMA   8             // somatic potential
#define  V_DEND   9             // dendritic potential
#define  N_VARS  10             // Number of variables

#define  SP_THRESH     40.0     // Spike detection threshold in mV
#define  SP_REFRACT     0.02    // Spike detection "refractory" period in sec.

#define  MAX_IFIELDS    64      // maximum number of incremental fields

#define  MAX_IPULSES    64      // maximum number of current pulses
#define  MAX_IPSPECS    64
#define  INIT_TIME_MAX   2.0    // maximim randomized initiation time
//
// Constants for resistive grid defintion
//
#define  RR_TG     12.0  
#define  RR_NN      0.1  
#define  RR_AA      0.07568260064501945   
#define  RR_BB      2.383944022328329
//
// Type Definitions
//
typedef double  *Vector;
typedef double **Matrix;

typedef enum { DT_NULL = 0,
               DT_FLAG,
               DT_FLOAT,
               DT_INTEGER,
               DT_STRING,
               DT_SPECIAL } DataType;

typedef struct _gate {          // channel kinetics variables
   double m;                    // sodium channel - activation (not integrated)
   double h;                    // sodium channel - deactivation 
   double n;                    // potassium DR   - activation
   double s;                    // calcium - activation
   double c;                    // postassim - calcium activation
   double q;                    // potassium AHP - activation 
   double Ca;                   // calcium concentration
} gate;

typedef struct _conductance {
   double c;                    // difference between dendrite and soma (gc)
   double L;                    // leakage
   double Na;                   // sodium channel
   double KDR;                  // potassium - delayed-rectifier (DR)
   double Ca;                   // calcium
   double KAHP;                 // potassium - slow after-hyperpolarization (AHP)
   double KCa;                  // potassium - calcium activated
   double NMDA;                 // NMDA
   double AMPA;                 // AMPA
} conductance;     

typedef struct _voltage {
   double soma;                 // transmembrane potential at somatic     (vsoma)
   double dend;                 // transmembrane potential at dendritic   (vdena)
   double syn;                  // synaptic potential
   double dso;                  // extracellular dendrite-soma potaential (vdsout)
   double L;                    // leakage  potential
   double Na;                   // reversal potential sodium
   double K;                    // reversal potential - potassium
   double Ca;                   // reversal potential - calcium
   double last;                 // vsoma one timestep earlier - to determine spikes
} voltage;

typedef struct _current {
   double soma;                 // transmembrane current at soma
   double dend;                 // transmembrane current at dendrite
   double nmda;                 // NMDA
   double ampa;                 // AMPA
   double syn;                  // synaptic current = iNMDA + iAMPA 
   double dsi;                  // inreacellular dendrite-to-soma currrent
} current;

typedef struct _currentinj {    // injected current, steady + pulse
   double  dend;                // dendrite, per unit area
   double  soma;                // soma,     per unit area
   double  denda;               // dendrite
   double  somaa;               // soma
} currentinj;

typedef struct _neuron {      // parameters of a neuron 
   double      init;          // initiation time seconds between similation starts  
   double      Area;          // surface area of neuron
   double      P;             // proportion in area for soma
   double      cm;            // membrane capacitance
   gate        k;             // channel gate and other kinetic variables
   conductance g;             // conductances
   voltage     v;             // potentials
   current     i;             // currents
   currentinj  inj;           // injected currents = steady + pulse
   double      SNMDA;         // NMDA synaptic weighing variable
   double      WAMPA;         // AMPA synaptic weighing variable
   char        ready;
   int         tally;         // to suppress activity
   double      Spike;         // spike time
   double      Spike0;        // last spike time
   double      speed;         // propagation speed (elapsed time per neuron)
} neuron;

typedef struct _EField {
  double  applied;              // applied electric field
  double  potential;            // field potential
} EField;

typedef struct _rratio {    // grid resistance ratio
   double tgR;              // top-ti-ground, "vertical"
   double nnR;              // neuron-to-neuron, "horizontal"
   double aaZ;              // terminal - vertcial
   double bbZ;              // terminal - horizonatl
} rratio;

typedef struct _rgrid {         // resistive lattice
   double TD;
   double SG;
   double DS;
   double SS;
   double DD;
   double AA;                  // pre-calc constant
   double BB;                  // pre-calc constant
} rgrid;

typedef struct _kinetics {      // gating kinetic variables
   double alpha;                // 
   double beta;                 // 
   double inf;                  // infinity (steady state)
   double tau;                  // time constant   
} kinetics;

typedef struct _gatevars {      // gating variables
   kinetics m;                  // sodium channel - activation
   kinetics h;                  // sodium channel - deactivation 
   kinetics n;                  // potassium DR   - activation
   kinetics s;                  // calcium - activation
   kinetics c;                  // postassim - calcium activation
   kinetics q;                  // potassium AHP - activation
} gatevars;

typedef struct _currentpulse {  // specification of current pulses
   double magn;                 // magnitude in uA/cm^2 
   double dur;                  // duration in seconds 
   double onset[MAX_IPULSES+1]; // times of pulses
   int    count;                // number of pulses
   int    index;                // for counting the numebr of pulses
   char   site;                 // 's' for soma; 'd' for dendrite
} currentpulse;

typedef struct _appliedfield {  // specification of applied electric field
   double  min;                 // min magnitude in mV/cm 
   double  incr;                // increment     in mV/cm 
   double  max;                 // max magnitude in mV/cm 
   double  onset;               // onset time
} appliedfield;

typedef struct _initfield {  // specification of initiating electric field v2.05
   double  ef;                  // magnitude in mV 
   int     ineurons;            // number of iniiating neurons in initiating set 
} initfield;

typedef struct _ipulsespn {     // temp buffer for parsing -ip option 
   char    neurons[128];
   char    times[128];
   char    site;
   double  magn;
   double  dur;
} ipulsespn;

typedef struct _ipulsespec {    // temp buffer for parsing -ip option 
   int       count;
   ipulsespn spec[MAX_IPSPECS]; // one per -ip option, max 64
} ipulsespec;

typedef struct _ifield {     // temp buffer for parsing -ip option 
   double  magn;
   double  onset;
} ifield;

typedef struct _ifields {    // temp buffer for parsing -if option 
   int       count;
   ifield    ifs[MAX_IFIELDS];  // one per field increment, max 64
} ifields;

typedef struct _propagatespec { // parameters for propagation analysis
   int    spike;  // if set, record spike times and compute instantaneous speed
   int    speed;  // if set, compute speed by fitting spike times to a straight line   
   double tstart; // time in sec to start analysis, should sync with current pulse
   int    nstart; // neuron number to start for speed computation (first neuron=0)
   int    nend;   // neuron number to end for speed computation (default=last neuron)
} propagatespec;

typedef struct _tokenlist {     // for parsing -o option 
  char  *tok;                   // if matches this token. . .
  int   *flag;                  // set flag so that file will be opened later.
} tokenlist;

typedef struct _outputlist {    // for writing to output file; write_output() 
  FILE **fp;                    // ptr to file ptr of output file
  int    data;                  // offset to structure "neuron" of element to write
} outputlist;

typedef struct _filelist {      // Time-series output files - flags and pointers
   int wvs;    FILE  *vs;       // V-soma     (v.soma)
   int wvd;    FILE  *vd;       // V-dendrite (v.dend)
   int wis;    FILE  *is;       // i-soma     (i.soma)
   int wid;    FILE  *id;       // i-dendrite (i.dend)
   int wca;    FILE  *ca;       // [Ca++]
   int wkh;    FILE  *kh;       //
   int wkn;    FILE  *kn;       //
   int wks;    FILE  *ks;       //
   int wkc;    FILE  *kc;       //
   int wkq;    FILE  *kq;       //
   int wsn;    FILE  *sn;       //
   int wwa;    FILE  *wa;       //
// int wef;    FILE  *ef;       // applied electric field
   int wix;    FILE  *ix;       // synaptic current          (i.syn)
   int wii;    FILE  *ii;       // internal d-s current      (i.dsi)
   int wvo;    FILE  *vo;       // external d-s voltage diff (v.dso, vdsout)
   int wsp;    FILE  *sp;       // spike times
} filelist;

typedef struct _outputspec {    // Time-series output file - specifications
   double  start;               // start time, inclusive, in seconds
   double  end;                 // end   time, inclusive, in seconds
   int     interv;              // write one output every "interval" points        
} outputspec;

typedef struct _filetable {
   int   *open;     // 1: file is enabled, 0: file is disabled
   FILE **ptr;      //
   char  *name;     // file name prefix
   char  *descr;
} filetable;

typedef struct _progrparams { // program parameters
   int        argc;
   char     **argv;
   struct tm *datetime;
   filetable *files;
} progrparams;

typedef struct _modelparams { // model parameters
   int            tneurons;   // total number of neurons: init chain + prop chain
   int            nneurons;   // number of neurons in propagation chain
   int            nlattice;
   neuron        *iineuron;   // initial paramaters for initiating neurons   v2.05
   neuron        *ineuron;    // initial neuron parameters
   neuron        *neurons;    // array of neurons
   EField        *ef;         // applied field for each neuron
   gatevars      *gates;
   double         Kout;       // extracellular [K+];
   rratio        *rr;         // resistive grid resistance ratios  
   appliedfield  *af;         // applied field spec
   initfield     *efi;        // initiating field spec   v2.05
   ifields       *ief;        // incremental fields
   currentpulse  *ip;         // applied current pulse info
   propagatespec *pspec;      // propagation analysis spec (flags & values)
   double         field;      // applied electric field : input   to model()
   double         speed;      // propagation speed      : output from model()
   double         speed_sd;   // propagation speed standar deviation
   int            termneuron; // terminal neuron   v2.04
} modelparams;

typedef struct _derivparams { //
   int           nneurons;    //
   neuron       *neurons;
   gatevars     *gatekin;     // gvar - gate kinetics variables
   EField       *ef;          // electric field for each neuron
   Matrix        rg;          // RGMATRIX
   Vector        vds;         // VMatrix
   Vector        vinj;        // CRInj
   Matrix        yrginv;      // YRGinv
   double       *sum;         // TSUM          
   Matrix        syn;         // isynap
   rgrid        *r;           // resistances of grid    v2.00
} derivparams;

typedef struct _initcondval { // an entry of initial values     v2.02
   double   ef;    // key: applied field
   double   vs;    // transmembrane potentiat at the soma
   double   vd;    // transmembrane potentiat at the dendrite
   double   ca;    // [Ca++]
   double   kh;    // h gate
   double   kn;    // n gate
   double   ks;    // s gate
   double   kc;    // c gate
   double   kq;    // q gate
   double   sn;    // synaptic weight - NMDA
   double   wa;    // synaptic weight - AMPA
} initcondval;

typedef struct _initcondlist {     // v2.02
  initcondval    *curr;
  void           *next;
} initcondlist;

//=== GLOBAL VARIABLE DECLARATIONS ============================================

// Command line options:
// 
// Flags:
int   always       = 1;             // always on
int   never        = 0;             // always off
int   trace        = 0;
int   logfile      = 0;             // if on, write message to file instead of stdout
int   autostop     = 0;             // terminate program automatically
int   no_r         = 0;             // disable resistive grid
int   no_gnn       = 0;             // no conductance between neurons; high "horizontal" R
int   no_rnn       = 0;             // no resistance  between neurons; low  "horizontal" R
int   filetimestamp= 0;             // if on, appending timestamp to files
int   nofilesuffix = 0;             // if on, disable appending e-field & timestamp to files
int   rand_vsoma=0, rand_vdend=0,   // randomization
      rand_gc   =0, rand_init =0;
int   pr_model     =  0;
int   opt_ef = 0, opt_if = 0;       // single ef and incremnetal ef options
//
// Model parameter value:
int           tot_neurons = -1;     // total number of neurons
int           num_neurons = -1;     // num_neurons required to be specified v2.05
double        dKOut       = 3.50;   // extracell [K+]. 3.50 spiking, 8.50 bursting
initfield     efi    = {0.0,0};     // initiating E-Field              V2.05
appliedfield  aefld  = {0,0,0,0};   // applied E-Field
currentpulse *ipulse;               // injected current pulse, 0.500
double        threshold = 25.0;     // (mV) threshold to determin spike.
ifields       ifspec = { 0 };       // for storing i-field option
ipulsespec    ipspec = { 0 };       // for storing i-pulse option for 2nd-pass parsing
propagatespec propaspec = {0,0};    // specification for propagation analysis
// Output file parameters
outputspec outspec = { 0.0, -1.0, 0 };      // time-series output range (in time)
char       outdir [128+1]  = "\0";          // output file directory
char       filenote[64+1]  = "\0";          // file note after filename suffix
char       initcond[128+1] = "\0";          // input file for initial conditions v2.01
char       initconi[128+1] = "\0";          // input file for init cond for init chain 2.05
char       outequil[128+1] = "\0";          // output file for equil init values    v2.01
char       outtermn[128+1] = "\0";          // terminal neuron log   v2.04
//
// Time series output files - one flag and one file pointer for each
filelist flist;
//
// Other output files
FILE *f_prop   = (FILE *)0;         // spike propagation
FILE *f_pvel   = (FILE *)0;         // spike propagation velocity, one file per run
FILE *f_term   = (FILE *)0;         // terminal neuron log   v2.04
FILE *f_param  = (FILE *)0;         // simulation parameters
FILE *f_trace  = (FILE *)0;         // show vectors and matrices
FILE *f_log    = (FILE *)0;         // stdout;   // default to standard output
//
char *pname  = (char *)0;           // program name, points to argv[0]
//
neuron     *neur, ineur, iineur;
gatevars    gates;  
// resistance  dR = { 0.0, 0.0, 0.0 };
rratio      drr = { 0.0, 0.0, 0.0, 0.0 };
rgrid      *em_dR;
EField     *efield;
//
Matrix    spikes;                   // in seconds, spikes[num neurons][num of spikes] (future)
Matrix    isynap;
Matrix    dy, dy0, dydt;
Matrix    YRG, YRGinv, yid, dy, dy0;
Matrix    GMatrix, RGMatrix, FMatrix, Finv;
Vector    VMatrix, TSUM;
Vector    CInj, CRInj;
//
double    dTime = 0.0, dTotal = 12.5;

//=== Function Declarations ==================================================
//
// P-R model simulation functions
//
int       model( modelparams *mp, progrparams *pp, derivparams *dp );
double    model_derivatives( double x, Matrix dy, Matrix dydt, void *dp );
double    vr_potassium( double ko );
neuron   *update_neurons  ( neuron *neurons,   Matrix  variables, int nstart, int nend );
Matrix    update_variables( Matrix  variables, neuron *neurons,   int nstart, int nend );
Matrix    connect_synapses( Matrix  isyn,                         int nstart, int nend );
int       compute_speed( modelparams *m );
double    ratio(   rgrid *em_dR, double k_out, rratio *rr, neuron *neur, int nneurons );
double    network( rgrid *em_dR, double **FMatrix,         int nlattice, int nneurons );
kinetics *Gate_m( kinetics *, double volt );
kinetics *Gate_h( kinetics *, double volt );
kinetics *Gate_n( kinetics *, double volt );
kinetics *Gate_s( kinetics *, double volt );
kinetics *Gate_c( kinetics *, double volt );
kinetics *Gate_q( kinetics *, double conc );
double    zeta(double dDel);
//
// Numerical and matrix computation functions (generic)
//
Matrix  rk4m( Matrix y, Matrix dydt, int r, int c, double x, double h,
              double (*derivs)(), void *dp );
double  LUDCMP( double **a, int n, int *indx );
double *LUBKSB( double **a, int n, int *indx, double *b );
Matrix  identityMatrix( double **m,              int n );
Matrix  inverseMatrix( double **inv, double **m, int n );
void fit( Vector x, Vector y, int ndata,
          double *a,    double *b,
          double *siga, double *sigb, double *chi2, double *q);
char *get_time( char *str );
//
// Command line parse functions
//
int  get_command( int argc, char **argv );                 // later, pass progrparams
int  parse_options  ( void *tb, char *args, char *sw );    // later, pass progrparams
int  parse_outrange ( void *tb, char *args, char *sw );
int  parse_e_field  ( void *tb, char *args, char *sw );    // later, pass progrparams
int  parse_iefield  ( void *tb, char *args, char *sw );    // later, pass progrparams  v2.05
int  parse_i_field  ( void *tb, char *args, char *sw );    // later, pass progrparams
int  parse_i_pulse  ( void *tb, char *args, char *sw );    // later, pass progrparams
int  parse_i_pulse2 ( ipulsespec *sp  );                   // later, pass progrparams
int  parse_propagate( void *tb, char *args, char *sw );    // later, pass progrparams
//
// Output processing functions
//
FILE   *open_output_file( FILE **f_ptr,   char *f_prefix, char *dir, double ef,
                          struct tm *bt,  char *f_note,   char *f_descr );
void write_output( outputlist *f, double t, neuron *nu, int nstart, int nend );
void write_speed     ( FILE *fp, modelparams *mp );
void write_parameters( FILE *fp, modelparams *m, progrparams *p );
//
// Data structure management functions
//
void memory_allocate( modelparams * );
void memory_clear   ( modelparams * );
void memory_free    ();
void  *vector_s( int n, size_t size );
Vector vector  ( int n );
Matrix matrix  ( int r, int c );
void   clear_vector( Vector v, int r, int c );
void   clear_matrix( Matrix m, int r, int c );
void   free_vector( void  *v );
void   free_matrix( void **m );
void   show_vector( FILE *fp, Vector v, int n,        char *name );
void   show_matrix( FILE *fp, Matrix m, int r, int c, char *name );

//=============================================================================
//
// Command line options:
//
//
// Options for the "-o" option. Selects output types to be written to file.
// Format: e.g. . . . -o vs,vd,is
// Parsed by parse_option().
//
// When you need a new output file, add an entry to this table.
//
static tokenlist out_options[] = {  // output file for:
  { "sp", &(flist.wsp) }, // spike times
  { "vs", &(flist.wvs) }, // V-soma
  { "vd", &(flist.wvd) }, // V-dendrite
  { "is", &(flist.wis) }, // i-soma
  { "id", &(flist.wid) }, // i-dendrite
  { "ca", &(flist.wca) }, // [Ca++]
  { "kh", &(flist.wkh) }, // h
  { "kn", &(flist.wkn) }, // n
  { "ks", &(flist.wks) }, // s
  { "kc", &(flist.wkc) }, // c
  { "kq", &(flist.wkq) }, // q
  { "sn", &(flist.wsn) }, // S-nmda
  { "wa", &(flist.wwa) }, // W-ampa
  { "ix", &(flist.wix) }, // synaptic current
  { "ii", &(flist.wii) }, // intravellular dendrite->soma current
  { "vo", &(flist.wvo) }, // extracellular dendrite->soma potential
//{ "ef", &(flist.wef) }, // appl. e-field & total inj. current.    
  { (char *)0,(int *)0 }
};
//
// Options for the "-rand" option. Selects parameters and init values to randomize.
// Format: e.g. . . . -rand vs,vd,gc
// Parsed by parse_option().
//
// When you need a new parameter to be randomized, add an entry to this table.
//
static tokenlist rand_options[] = { // randomize parameter:
  { "vs",   &rand_vsoma },          // V-soma
  { "vd",   &rand_vdend },          // V-dendrite
  { "gc",   &rand_gc    },          // g-c
  { "ni",   &rand_init  },          // neuron initiation
  { (char *)0, (int *)0 }
};
//
// Main Command Line Options Table::
//
static struct {
   char      *param;                       // parameter switch ( "-..." )
   void      *value;                       // parameter value  (optional)
   DataType   type;                        // data type
   int        size;                        // max length for type DT_STRING
   int  (*func)( void *, char *, char * ); // conversiuon function
} options[] = {
   //
   // switch       location for         data type  max length  special parsing
   //          parsed/converted value                             function
   //
   { "-nn",    (void *)(&num_neurons),   DT_INTEGER, 0,            0 }, // num neurons
   { "-pr",    (void *)(&pr_model),      DT_FLAG,    0,            0 }, // P-R 1-neuron model
   { "-area",  (void *)(&ineur.Area),    DT_FLOAT,   0,            0 }, // num neurons
   { "-ap",    (void *)(&ineur.P),       DT_FLOAT,   0,            0 }, // proportion
   { "-gc",    (void *)(&ineur.g.c),     DT_FLOAT,   0,            0 }, // g-c
   { "-gnmda", (void *)(&ineur.g.NMDA),  DT_FLOAT,   0,            0 }, // g-NMDA
   { "-gampa", (void *)(&ineur.g.AMPA),  DT_FLOAT,   0,            0 }, // g-AMPA
   { "-vl",    (void *)(&ineur.v.L),     DT_FLOAT,   0,            0 }, // leakage  potential
   { "-ko",    (void *)(&dKOut),         DT_FLOAT,   0,            0 }, // extracellular [K+]
   { "-ca",    (void *)(&ineur.k.Ca),    DT_FLOAT,   0,            0 }, // intracellular [Ca++]
   { "-isoma", (void *)(&ineur.inj.soma),DT_FLOAT,   0,            0 }, // steady i-inj at soma
   { "-idend", (void *)(&ineur.inj.dend),DT_FLOAT,   0,            0 }, // steady i-inj dendrite
   { "-vsoma", (void *)(&ineur.v.soma),  DT_FLOAT,   0,            0 }, // Vm at soma
   { "-vdend", (void *)(&ineur.v.dend),  DT_FLOAT,   0,            0 }, // Vm at dendrite
   { "-snmda", (void *)(&ineur.SNMDA),   DT_FLOAT,   0,            0 }, // synaptic weight NMDA
   { "-wampa", (void *)(&ineur.WAMPA),   DT_FLOAT,   0,            0 }, // synaptic weight AMPA
   { "-ic",    (void *)(initcond),       DT_STRING,128,            0 }, // pathname of init cond
   { "-ici",   (void *)(initconi),       DT_STRING,128,            0 }, // pathname of init cond for init chain
   { "-oe",    (void *)(outequil),       DT_STRING,128,            0 }, // pathname for equil vals
   { "-oterm", (void *)(outtermn),       DT_STRING,128,            0 }, // pathname for term neuron log
   { "-efi",   (void *)(&efi),           DT_SPECIAL, 0,parse_iefield }, // initiating e-field
   { "-ef",    (void *)(&aefld),         DT_SPECIAL, 0,parse_e_field }, // applied e-field
   { "-if",    (void *)(&ifspec),        DT_SPECIAL, 0,parse_i_field }, // incremental field
   { "-ip",    (void *)(&ipspec),        DT_SPECIAL, 0,parse_i_pulse }, // applied current pulse
   { "-rand",  (void *)(&rand_options),  DT_SPECIAL, 0,parse_options }, // randomize specified vars
   { "-ngnn",  (void *)(&no_gnn),        DT_FLAG,    0,            0 }, // no extracell. n-n g
   { "-nrnn",  (void *)(&no_rnn),        DT_FLAG,    0,            0 }, // no extracell. n-n R
   { "-nr",    (void *)(&no_r),          DT_FLAG,    0,            0 }, // no resistive net
   { "-astop", (void *)(&autostop),      DT_FLAG,    0,            0 }, // stop program auto
   { "-o",     (void *)(&out_options),   DT_SPECIAL, 0,parse_options }, // write specified outputs
   { "-or",    (void *)(&outspec),       DT_SPECIAL, 0,parse_outrange}, // output time range
   { "-oi",    (void *)(&outspec.interv),DT_INTEGER, 0,            0 }, // output interval
   { "-p",     (void *)(&propaspec),     DT_SPECIAL, 0,parse_propagate},// analyze propagation
   { "-dur",   (void *)(&dTotal),        DT_FLOAT,   0,            0 }, // duration of run
   { "-trace", (void *)(&trace),         DT_FLAG,    0,            0 }, // writes matrices
   { "-dir",   (void *)(outdir),         DT_STRING,128,            0 }, // output file dir
   { "-log",   (void *)(&logfile),       DT_FLAG,    0,            0 }, // messages to file
   { "-fts",   (void *)(&filetimestamp), DT_FLAG,    0,            0 }, // no file timestamp
   { "-nfs",   (void *)(&nofilesuffix),  DT_FLAG,    0,            0 }, // no file suffix
   { "-fn",    (void *)(filenote),       DT_STRING, 64,            0 }, // file note after suffux
   { (char *)0,(void *)0,                DT_NULL,    0,            0 }
};
//
// List of output files
//
static filetable filetb[] = {
 { &never,    &f_log,   "pr_mslog","Message log file"          },    // always first entry
// { &always,   &f_param, "pr_param","Parameters"                         },  
 { &flist.wvs,&flist.vs,"pr_vsoma","Vm at Soma: t (s) | mV, 1..n neurons"    },  
 { &flist.wvd,&flist.vd,"pr_vdend","Vm at Dendrite: t (s) | mV, 1..n neurons" },
 { &flist.wis,&flist.is,"pr_isoma","Im at Soma: t (s) | uA/cm, 1..n neurons"    },
 { &flist.wid,&flist.id,"pr_idend","Im at Dendrite: t (s) | uA/cm, 1..n neurons"},
 { &flist.wca,&flist.ca,"pr_conca","[Ca++]: t (s) | mM, 1..n neurons" }, 
 { &flist.wkh,&flist.kh,"pr_gateh","h: t (s) | *, 1..n neurons" }, 
 { &flist.wkn,&flist.kn,"pr_gaten","n: t (s) | *, 1..n neurons" }, 
 { &flist.wks,&flist.ks,"pr_gates","s: t (s) | *, 1..n neurons" }, 
 { &flist.wkc,&flist.kc,"pr_gatec","c: t (s) | *, 1..n neurons" }, 
 { &flist.wkq,&flist.kq,"pr_gateq","q: t (s) | *, 1..n neurons" }, 
 { &flist.wsn,&flist.sn,"pr_snmda","S-NMDA: t (s) | *, 1..n neurons" }, 
 { &flist.wwa,&flist.wa,"pr_wampa","W-AMPA: t (s) | *, 1..n neurons" }, 
//{&flist.wef,&flist.ef,"pr_field","Applied E-Field: t (s) | mV/cm^2 1..n neurons" },
 { &flist.wix,&flist.ix,"pr_synap","Synaptic current: t (s) | uA/cm, 1..n neurons" },
 { &flist.wii,&flist.ii,"pr_id2si","Internal Dendrite-to-Soma Current: t (s) | uA/cm, 1..n neurons"},
 { &flist.wvo,&flist.vo,"pr_vd2so","External Dendrite-Soma Potential Diff: t (s) | mV, 1..n neurons"},
 { &flist.wsp,&flist.sp,"pr_spike","Spikes: list of occurence in time (s), 1..n neurons" },
 { &(propaspec.spike),
              &f_prop,  "pr_props","Propagation of a spike for neurons in order: time: instant speed"},
 { &never,    &f_pvel,  "pr_propc","Speed of spike propagation as a function of applied field:"},
 { &trace,    &f_trace, "pr_trace","Intermediate matrices"     },
 { (int *)0,   (FILE **)0, (char *)0, (char *)0 }
};
//
// Table for write_output()
//
static outputlist outfiles[] = {
  { &(flist.vs), OFFSET(neuron,v.soma) }, // V-soma
  { &(flist.vd), OFFSET(neuron,v.dend) }, // V-dendrite
  { &(flist.is), OFFSET(neuron,i.soma) }, // i-soma
  { &(flist.id), OFFSET(neuron,i.dend) }, // i-dendrite
  { &(flist.ca), OFFSET(neuron,k.Ca  ) }, // [Ca++]
  { &(flist.kh), OFFSET(neuron,k.h   ) }, // h
  { &(flist.kn), OFFSET(neuron,k.n   ) }, // n
  { &(flist.ks), OFFSET(neuron,k.s   ) }, // s
  { &(flist.kc), OFFSET(neuron,k.c   ) }, // c
  { &(flist.kq), OFFSET(neuron,k.q   ) }, // q
  { &(flist.sn), OFFSET(neuron,SNMDA ) }, // S-NMDA
  { &(flist.wa), OFFSET(neuron,WAMPA ) }, // W_AMPA
  { &(flist.ix), OFFSET(neuron,i.syn ) }, // synaptic current
  { &(flist.ii), OFFSET(neuron,i.dsi ) }, // intravellular dendrite->soma current
  { &(flist.vo), OFFSET(neuron,v.dso ) }, // extracellular dendrite->soma potential
  { &(flist.sp),  -1                   }, // spike times
//{ &(flist.ef),  -1                   }, // appl. e-field & total inj. current.
  { (FILE **)0,   0                    }
};

//=========================== Program Entry ===================================
//
// function: int main();
// 
int main( int argc, char **argv )
{
   int       i  = 0;
   int       rc = 1;
   char      str[160] = "";
   char      timestr[32] = "";
   double    rr_factor_h = 1.0;
   double    rr_factor_v = 1.0;
   time_t    systime;       // move to main
   struct tm brktime;
   FILE     *fic = (FILE *)0; // initial condition input file
   FILE     *fii = (FILE *)0; // initial condition input file for init chain
   FILE     *foe = (FILE *)0; // variable values output file    v2.02
   double    tmp_ef = 0.0;
   int       ii_found = 0;
   char      *cstr = (char *)0;
   size_t     clen;
   ssize_t    llen;
   progrparams pparams;
   modelparams mparams;
   derivparams dparams;
   initcondlist *iclist = (initcondlist *)0;     // list of init conds as a function of E v2.02
   initcondlist *iccurr = (initcondlist *)0;
   initcondlist *icprev = (initcondlist *)0;
   
   f_log = stdout;    // default is to write messages to standard output
   (void )time( &systime );                        // get time
   (void )localtime_r( &systime, &brktime );
   //
   (void )memset( &ineur, 0, sizeof( neuron ) );
   (void )memset( &flist, 0, sizeof(filelist) );
   //
   // Default Values - Parameters
   //
   ineur.P        =   0.50;
   ineur.Area     =   6.0*pow(10.0,-6.0);  // 0.0000060;
   ineur.cm       =   3.0;
   ineur.v.Na     = 120.0;
   ineur.v.Ca     = 140.0;
   ineur.v.K      = vr_potassium( dKOut ); // -15.0 bursting, -38.56 spiking
   ineur.v.L      =   0.0;    // 0.0 in PR
   ineur.v.syn    =  60.0;    // synaptic potential
   ineur.g.c      =   2.1;    // 3.1
   ineur.g.L      =   0.1;    // conductance - leakage
   ineur.g.Na     =  30.0;    // conductance - Na+ channel
   ineur.g.KDR    =  15.0;    // conductance - K+
   ineur.g.KAHP   =   0.80;   // conductance - 
   ineur.g.KCa    =  15.00;   // conductance - 
   ineur.g.Ca     =  10.00;   // conductance - 
   ineur.g.NMDA   =   0.0300; // 0.005; //0.0300 in E-H's case,
   ineur.g.AMPA   =   0.0045; // 0.025; //0.0045 in E-H's case,
   ineur.i.soma   =   0.000;  // init always 0; to be computed  
   ineur.i.dend   =   0.000;  // init always 0; to be computed
   ineur.inj.soma =   0.000;  // steady inj current. P-R -0.5; 0.00 in E-H's case
   ineur.inj.dend =   0.000;  // steady inj current. P-R  0.0; 0.70 in E-H's case (ODE), 0 in PR_Rev2
   //
   // Default Values - Initial Conditions - Integration Variables
   //
   ineur.k.h      =  0.9990;  // h gate
   ineur.k.n      =  0.001;   // n gate
   ineur.k.s      =  0.009;   // s gate
   ineur.k.c      =  0.007;   // c gate
   ineur.k.q      =  0.010;   // q gate
   ineur.k.Ca     =  0.200;   // [Ca++]
   ineur.SNMDA    =  0.6;     // synaptic weight - NMDA; 0.3 E-H ODE
   ineur.WAMPA    =  0.5;     // synaptic weight - AMPA; 0.2 E-H ODE
   ineur.v.soma   = -4.6;     // transmembrane potentiat at the soma
   ineur.v.dend   = -4.5;     // transmembrane potentiat at the dendrite
   //
   // Other default values
   //
   // ipulse.neuron = 1;         // steady current at neuron[0] (neuron one)
   //
   // Parse command line to get options,
   // some of the parameters may be replaced. 
   //
   if( !get_command( argc, argv ) ) return( 0 );
   //
   // Determine initial conditions
   // Use default, read from file or use preset values for single neuron
   // Single neuron case overrides -ic option.
   //
   ineur.v.K = vr_potassium( dKOut );
   if( pr_model ) {      // P-R single neuron model parameters:
      tot_neurons    =  1;
      num_neurons    =  1;
      no_r           =  1;    // No resistive grid     
      dKOut          =  8.50; // to get complex spikes (bursts)    
      ineur.v.K      = vr_potassium( dKOut ); // -15.0 bursting
      ineur.v.L      =  0.0;  // 0.0 in PR
      ineur.g.NMDA   =  0.0;  // no synapses
      ineur.g.AMPA   =  0.0;  // no synapses
      ineur.inj.soma = -0.05; // steady inj current. P-R -0.5; 0.00 in E-H's case
      ineur.inj.dend =  0.00; // steady inj current. P-R  0.0; 0.70 in E-H's case (ODE), 0 in PR_Rev2
   }
   else if( *initcond ) {     // else if non-null file name read init values v2.1
      if( !(fic = fopen( initcond, "r" ) ) ) {
         printf( "%s: Failed to open initial condition input file.\n", argv[0] );
         exit(0);
      }
      (void )getline( &cstr, &clen, fic );                          // discard 1st line
      i=0;
      while( (llen = getline( &cstr, &clen, fic ) ) != -1 ) {
         iccurr       = malloc( sizeof(initcondlist) );
         if( !i ) {
            iclist       = iccurr;
            i++;
         }
         else {
            icprev->next = iccurr;                    
         }
         iccurr->curr = malloc( sizeof(initcondval) );
         iccurr->next = (initcondlist *)0;
         icprev       = iccurr;         
         sscanf( cstr, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &(iccurr->curr->ef),  // applied field
            &(iccurr->curr->vs),  // transmembrane potential at the soma
            &(iccurr->curr->vd),  // transmembrane potential at the dendrite
            &(iccurr->curr->ca),  // [Ca++]
            &(iccurr->curr->kh),  // h gate
            &(iccurr->curr->kn),  // n gate
            &(iccurr->curr->ks),  // s gate
            &(iccurr->curr->kc),  // c gate
            &(iccurr->curr->kq),  // q gate
            &(iccurr->curr->sn),  // synaptic weight - NMDA
            &(iccurr->curr->wa)   // synaptic weight - AMPA
         );
      }
      if( cstr ) free( cstr );
      fclose( fic );
   }
   if( *outequil ) {
      if( !(foe = fopen( outequil, "w" ) ) ) {
         printf( "%s: Failed to open variable values output file.\n", argv[0] );
         exit(0);
      }
      fprintf( foe, "# EF: Vs Vd ca h n s c q sNMDA wAMPA\n" );
   }
   if( *outtermn ) {
      if( !(f_term = fopen( outtermn, "w" ) ) ) {
         printf( "%s: Failed to open terminal neuron log file.\n", argv[0] );
         exit(0);
      }
      fprintf( f_term, "# EF: terminal neuron\n" );
   }
   //
   // Define parameters based on command options
   //
   // If no_gnn, set high extracellular neuron-to-neuron resistance;
   // otherwise use nominal values.
   //
   rr_factor_v = 1.0;
   rr_factor_h = 1.0;
   if(      no_gnn && !no_rnn ) {
      rr_factor_h = 999999999999.0;
   }
   else if( no_rnn && !no_gnn ) {
      rr_factor_h = 0.000000000001;
   }
   drr.tgR = RR_TG * rr_factor_v;   // Top-to-ground    ("vertical")
   drr.nnR = RR_NN * rr_factor_h;   // neuron-to-neuron ("horizontal" )
   drr.aaZ = RR_AA * rr_factor_v;   // terminal - vertical   ends
   drr.bbZ = RR_BB * rr_factor_h;   // terminal - horizontal ends
   //
   // Memories for data structures are allocated in main();
   // they are initialized in model() for each applied field.
   mparams.tneurons  =  tot_neurons;
   mparams.nneurons  =  num_neurons;
   mparams.nlattice  =  5*(num_neurons+1)+3;   // resistive grid size (2 neurons=>18)
   memory_allocate( &mparams );
   // Program Parameters   
   //
   pparams.argc      =  argc;
   pparams.argv      =  argv;
   pparams.datetime  = &brktime;
   pparams.files     =  filetb;
   //
   // Model Parameters - all 
   //
   mparams.ineuron   = &ineur;
   mparams.iineuron  = &iineur;        // initial values for initiaiting neurons v2.05 
   mparams.neurons   =  neur;
   mparams.gates     = &gates;
   mparams.Kout      =  dKOut;
   mparams.rr        = &drr; 
   mparams.field     =  0.0;            // applied field to be assigned for model() loop.
   mparams.ef        =  efield;         // applied field for each neuron ([])
   mparams.af        = &aefld;          // applied field spec
   mparams.efi       = &efi;            // initiating field spec    v2.05
   mparams.ief       = &ifspec;         // incremental field spec
   mparams.ip        =  ipulse;         // applied current pulse spec - array for neurons
   mparams.pspec     = &propaspec;      // propagation analysis spec 
   mparams.speed     =  0.0;            // propogation velocity: output from model().
   mparams.speed_sd  =  0.0;            // propogation velocity standard deviation.
   //
   // Derivatives Parameters passed to model_derivatives()
   //
   dparams.tneurons  =  tot_neurons;    //
   dparams.nneurons  =  num_neurons;    //
   dparams.neurons   =  neur;
   dparams.gatekin   = &gates;          // gvar - gate kinetics variables
   dparams.ef        =  efield;         // applied field for each neuron ([])
   dparams.rg        =  RGMatrix;
   dparams.vds       =  VMatrix;
   dparams.vinj      =  CRInj;          // potential due to total injected current
   dparams.yrginv    =  YRGinv;
   dparams.sum       =  TSUM;          
   dparams.syn       =  isynap;
   dparams.r         =  em_dR;          // resistances of grid   v2.0 
   // Create directory
   // There is no error checking. If it fails it could mean that the dir exists,
   // otherewise, fopen() would fail and program would exit. 
   //
   if( *outdir ) {
      sprintf( str, "mkdir %s", outdir );     // for Unix only.
      system(  str );
   }
   //
   // Open message log file if specified in command line
   //
   if( logfile ) { 
      for( i=0;  filetb[i].name && strcmp( filetb[i].name, "pr_mslog" ); i++ );
      if( filetb[i].name ) {
         if( !(f_log = open_output_file( &f_log,   filetb[i].name, outdir,
                                  99999, &brktime, filenote, filetb[i].descr ) ) ) {
            return( 0 );
         }
      }
   }
   //
   // Get parameters intial values for initiating neurons   v2.05
   //
   (void )memcpy( &iineur, &ineur, sizeof(neuron) );
   if( efi.ineurons > 0 && *initconi ) {
   // Open file and get one record
      if( !(fii = fopen( initconi, "r" ) ) ) {
         printf( "%s: Failed to open initial chain initial condition input file.\n", argv[0] );
         exit(0);
      }
      (void )getline( &cstr, &clen, fii );                          // discard 1st line
      i=0;
      while( (llen = getline( &cstr, &clen, fic ) ) != -1 ) {
         sscanf( cstr, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &tmp_ef,              // applied field
            &iineur.v.soma,       // transmembrane potentiat at the soma
            &iineur.v.dend,       // transmembrane potentiat at the dendrite
            &iineur.k.Ca,         // [Ca++]
            &iineur.k.h,          // h gate
            &iineur.k.n,          // n gate
            &iineur.k.s,          // s gate
            &iineur.k.c,          // c gate
            &iineur.k.q,          // q gate
            &iineur.SNMDA,        // synaptic weight - NMDA
            &iineur.WAMPA         // synaptic weight - AMPA
         );
         if( temp_ef == efi.ef ) {    // found V_top_init
            fii_found = 1;
            fprintf( f_log, "%f: %f %f .. %f %f %f %f %f %f %f %f\n",
               efi.ef,
               iineur.v.soma,       // transmembrane potential at the soma
               iineur.v.dend,       // transmembrane potential at the dendrite
               iineur.k.Ca,         // [Ca++]
               iineur.k.h,          // h gate
               iineur.k.n,          // n gate
               iineur.k.s,          // s gate
               iineur.k.c,          // c gate
               iineur.k.q,          // q gate
               iineur.SNMDA,        // synaptic weight - NMDA
               iineur.WAMPA         // synaptic weight - AMPA
            );        
            break;
         }                     
      }
      if( !fii_found ) {
         fprintf( f_log, "Error exit: Initial values for initial chain %f mV/cm not found.\n", efi.ef );
      }
   }          
   //
   // Open file that records propagation velocity as a function of applied field
   //
   if( propaspec.speed ) { 
      for( i=0;  filetb[i].name && strcmp( filetb[i].name, "pr_propc" ); i++ );
      if( filetb[i].name ) {
         if( !(f_pvel = open_output_file( &f_pvel, filetb[i].name, outdir,
                                  99999, &brktime, filenote, filetb[i].descr ) ) ) {
            return( 0 );
         }
      }
   }
   //
   // Model simulation as a fuction of applied fields 
   //
   fprintf( f_log,
      "Applied field (mV/cm) incrementing by %f from %f to %f\n",
                                aefld.incr, aefld.min, aefld.max );
   for( mparams.field  = aefld.min;
        mparams.field  < aefld.max + 0.0005;
        mparams.field += aefld.incr ) {
      fprintf( f_log, "Start: %s\n", get_time( timestr ) );
      //
      // If instructed, read init values.    v2.02
      // 
      if( *initcond ) {
         for( iccurr=iclist; iccurr->next &&
                           !(iccurr->curr->ef > mparams.field - 0.00001 &&
                             iccurr->curr->ef < mparams.field + 0.00001);
                             iccurr=iccurr->next );
         if( !(iccurr->curr->ef > mparams.field - 0.00001 &&
               iccurr->curr->ef < mparams.field + 0.00001) ) {  // don't check "iccurr->next"
            fprintf( f_log, "Error exit: Initial values for %f mV/cm not found.\n", mparams.field );
            continue;                // skip field         
         }
         else {
            ineur.v.soma =  iccurr->curr->vs;  // transmembrane potentiat at the soma
            ineur.v.dend =  iccurr->curr->vd;  // transmembrane potentiat at the dendrite
            ineur.k.Ca   =  iccurr->curr->ca;  // [Ca++]
            ineur.k.h    =  iccurr->curr->kh;  // h gate
            ineur.k.n    =  iccurr->curr->kn;  // n gate
            ineur.k.s    =  iccurr->curr->ks;  // s gate
            ineur.k.c    =  iccurr->curr->kc;  // c gate
            ineur.k.q    =  iccurr->curr->kq;  // q gate
            ineur.SNMDA  =  iccurr->curr->sn;  // synaptic weight - NMDA
            ineur.WAMPA  =  iccurr->curr->wa;  // synaptic weight - AMPA
            fprintf( f_log, "%f: %f %f .. %f %f %f %f %f %f %f %f\n",
               mparams.field,
               ineur.v.soma,       // transmembrane potentiat at the soma
               ineur.v.dend,       // transmembrane potentiat at the dendrite
               ineur.k.Ca,         // [Ca++]
               ineur.k.h,          // h gate
               ineur.k.n,          // n gate
               ineur.k.s,          // s gate
               ineur.k.c,          // c gate
               ineur.k.q,          // q gate
               ineur.SNMDA,        // synaptic weight - NMDA
               ineur.WAMPA         // synaptic weight - AMPA
            );
         }          
      }
      mparams.termneuron = 0;
      if ( !(rc = model( &mparams, &pparams, &dparams ) ) ) break;
      //
      //  Write final variable values -- normally for logging equilibrium values.
      //
      if( *outequil ) {      
         fprintf( foe, "%f %24.17f %24.17f %24.17f %24.17f %24.17f %24.17f %24.17f %24.17f %24.17f %24.17f\n",
            mparams.field,
            neur->v.soma, // transmembrane potentiat at the soma
            neur->v.dend, // transmembrane potentiat at the dendrite
            neur->k.Ca,   // [Ca++]
            neur->k.h,    // h gate
            neur->k.n,    // n gate
            neur->k.s,    // s gate
            neur->k.c,    // c gate
            neur->k.q,    // q gate
            neur->SNMDA,  // synaptic weight - NMDA
            neur->WAMPA   // synaptic weight - AMPA                                     
         );
         fflush( foe ); 
      }
      fprintf( f_log, "End  : %s\n", get_time( timestr ) );
      if( f_term ) {                                             // v2.04
         fprintf( f_term, "%f %i\n", mparams.field, mparams.termneuron );
         fflush(  f_term );
      }
      if( !aefld.incr ) break; 
      if( mparams.pspec->speed ) write_speed( f_pvel, &mparams );
   }
   memory_free();
   if( *outequil ) fclose( foe );
   if( *outtermn ) fclose( f_term );
   if( f_pvel ) fclose( f_pvel ); 
   return( rc );
}

//============================================================================
//
// function - int model();
//
int model( modelparams *mp, progrparams *pp, derivparams *dp )
{
   int     i=0, j=0;
   double  dtemp = 0.0;
   int     ifcount  = 0;
// int     ipn      = mp->ip->neuron-1; // index to neuron that receives current pulse
   double  em_dTime = 0.0;
   double  DelTime  = 0.0005;   // 0.005;    // 0.05;// Timestep in msec. changed Nov. 11 from 0.005;
   double  timestep = DelTime/1000.0;   // Timestep in seconds
   double  Fract    = 0.05;             // should not exceed .10 (10%)
   double  interval = outspec.interv * DelTime * 0.001;
   int     iiN  = mp->efi->ineurons;    // number of initiaiting neurons  v2.05
   int     iN   = mp->nneurons;
   int     iLAT = mp->nlattice;
   neuron *neur = mp->neurons;

   memory_clear( mp );
   //
   // Open files that are specified in file list for writing
   //
   for( i=0; filetb[i].name; i++ ) {
      if( *filetb[i].open ) {        // if flag is on, open file
         if( !( *(filetb[i].ptr) = open_output_file( filetb[i].ptr,
                                                     filetb[i].name, outdir,
                                mp->field, pp->datetime, filenote, filetb[i].descr ) ) )
            return( 0 );
      }
   }
   //
   // Make a copy of parameters and initial values for each neuron.
   // Copy for initiating neurons, if any (iiN=0 if none).            v2.05
   //
   for( i=0;   i<iiN; i++ ) (void )memcpy( &((mp->neurons)[i]), mp->iineuron, sizeof(neuron) );   
   for( i=iiN; i<iN;  i++ ) (void )memcpy( &((mp->neurons)[i]), mp->ineuron,  sizeof(neuron) );

   //--------------------------------------------------------------------------
   //            INITIAL CONDITIONS FOR 10 VARIABLES
   //--------------------------------------------------------------------------
   (void )update_variables( dy,  neur, 0, iN-1 );   // copy neurons to variables
   (void )update_variables( dy0, neur, 0, iN-1 );
   //
   // Varies initial value one or more of: v-soma, v-dendrite & dendrite-soma conductance
   //
   srand(7);
   if( rand_vsoma ) {
      for( i=0; i<iN; i++ ) dy[i][V_SOMA] = dy0[0][V_SOMA]*
                               ( 1.0 + Fract*(-0.50 + (fmod(rand()/10000.0,1.0)) ) );
   }
   if( rand_vdend ) {
      for( i=0; i<iN; i++ ) dy[i][+V_DEND] = dy0[0][V_DEND]*
                               ( 1.0 + Fract*(-0.50 + (fmod(rand()/10000.0,1.0)) ) );
   }
   if( rand_gc ) {
      for( i=0; i<iN; i++ ) neur[i].g.c = neur[i].g.c*
                               ( 1.0 + Fract*(-0.50 + (fmod(rand()/10000.0,1.0)) ) );
   }
   if( rand_init ) {                      // 0->1.0 seconds
      for( i=0; i<iN; i++ ) neur[i].init =  INIT_TIME_MAX * fmod( rand()/10000.0, 1.0 );
   }
   (void )update_neurons( neur, dy, 0, iN-1 ); // write initial variable values back to neurons

   //--------------------------------------------------------
   //
   // Compute matrix if resistive lattice enabled
   //
   (void )ratio( em_dR, mp->Kout, mp->rr, neur, iN );   // v2.0
   if( no_r ) {                                         // v2.0
      //
      // No resistive grid
      // Precalculate constants to be used in model_derivatives()  
      //
      for( i=0; i<iN; i++ ) {
         em_dR[i].AA = ( em_dR[i].TD + em_dR[i].DS + em_dR[i].SG ) / em_dR[i].DS;
         em_dR[i].BB = ( em_dR[i].TD + em_dR[i].SG ) * neur[i].g.c * neur[i].Area;
      }
   }
   else {
      //
      // Setup resistive network
      //
      (void )network( em_dR, FMatrix,          iLAT, iN );
      show_matrix( f_trace, FMatrix, iLAT, iLAT, "FMatrix" );
      (void )inverseMatrix( Finv, FMatrix, iLAT );  // take inverse of FMatrix
      show_matrix( f_trace, Finv,    iLAT, iLAT, "Finv" );
      //
      // Original comments (2) moved to end of this file.
      //
      for( i=0; i<iN; i++ ) {
         neur[i].inj.denda = ( neur[i].Area * neur[i].inj.dend );
         neur[i].inj.somaa = ( neur[i].Area * neur[i].inj.soma );
      }
      for( i=0; i<iN; i++ ) {
         for( j=0; j<iN; j++ ) {
            CInj[i] += Finv[(i*5)+ 1][(j*5)+3] * neur[j].inj.denda
                     + Finv[(i*5)+ 1][(j*5)+4] * neur[j].inj.somaa;
         }
      }
      for( i=0; i<iN; i++ ) CRInj[i]=em_dR[0].DS*CInj[i];
      for( i=0; i<iN; i++ ) {
         for( j=0; j<iN; j++ ) {
            GMatrix[i][j] = -Finv[(i*5)+1][(j*5)+3] * neur[j].g.c * neur[j].Area
                           + Finv[(i*5)+1][(j*5)+4] * neur[j].g.c * neur[j].Area;
                     // original comment: unsure about indexing: neur[i] or neur[j] for g.c and Area.
            RGMatrix[i][j]= em_dR[0].DS*GMatrix[i][j];
         }
      }
      show_matrix( f_trace,  GMatrix, iN, iN,  "GMatrix" );
      show_matrix( f_trace, RGMatrix, iN, iN, "RGMatrix" );
      (void )identityMatrix( yid, iN );
      for( i=0; i<iN; i++ ) {
         for( j=0; j<iN; j++ ) {
            YRG[i][j] = yid[i][j] - RGMatrix[i][j];
         }
      }
      show_matrix( f_trace, YRG,    iN,  iN, "YRG" );
      (void )inverseMatrix( YRGinv, YRG, iN );    // take inverse of YRG
      show_matrix( f_trace, YRGinv, iN,  iN, "YRGinv" );
      //
      //--------------------------------------------------------------------------
      //         WE HAVE INVERSE OF FMatrix (WITH REIMPOSED SYMMETRY) !!!
      //--------------------------------------------------------------------------
      fprintf( f_log, "Resistive grid computed.\n" );
   }
   if( f_trace ) fflush( f_trace ); // write to file when done with matrices
   //
   // Randomize neurons' initial state without synapses  (later: before resistive net)
   //
   if( rand_init ) {
      if( f_log ) for( i=0; i<iN; i++ ) fprintf( f_log, "init %i : %f\n", i, neur[i].init );
      for( em_dTime = 0.0, dTime = 0.0; dTime < INIT_TIME_MAX; )  {
         model_derivatives( em_dTime, dy, dydt, (void *)dp );
         rk4m( dy, dydt, iN, N_VARS, em_dTime, DelTime, model_derivatives, (void *)dp );
         em_dTime += DelTime;          // in msec
         dTime = em_dTime/1000.0;      // in seconds
         for( i=0; i<iN; i++ ) {
            if( dTime < neur[i].init ) {                // if not end of init run
               (void )update_neurons( neur, dy, i, i ); // update paramaters of one neuron
            }
         }
      }
      (void )update_variables( dy,  neur, 0, iN-1 );   // copy neurons to variables
      (void )update_variables( dy0, neur, 0, iN-1 );
   }
   (void )connect_synapses( isynap, 0, iN-1 );
   show_matrix( f_trace, isynap, iN, iN, "Synaptic connections" );
   if( f_trace ) fflush( f_trace ); // write to file before simulation loop
//   write_parameters( f_param, mp, pp );
   if( f_param ) fflush( f_param ); // write to file when done with output
   fflush( f_log );
// if( outspec.start == 0.0 ) write_output( outfiles, 0.0, neur, 0.0, iN-1 );

   // debug
   for( i=0; i<iN; i++ ) mp->ip[i].index = 0;

   //
   // Simulation
   //
   for( em_dTime = 0.0, dTime = 0.0; dTime < dTotal; em_dTime += DelTime )  {
      dTime = em_dTime/1000.0;                    // in seconds
      //======================= Impose Electric Field =========================
      //
      if( opt_ef ) {
         if( dTime > (mp->af->onset - 0.10) && dTime < (mp->af->onset + 0.10) ) {  // 5.0 => 4.9, 5.1
            for( i=0; i<iN; i++ ) {
               if( dTime >= mp->af->onset && dTime < (mp->af->onset + timestep) ) {
                  mp->ef[i].applied = (i < iiN) ? mp->efi->ef : mp->field;  // apply initiating field, if any  v2.05
                  if( i==0 )
                     fprintf( f_log, "\n%f sec   applying E field of %f mV/cm",
                                dTime, dp->ef[i].applied );  // make sure fld val pt'r to by dp
               } 
               //=============== TALLYING RESULTS OF SOMATIC SPIKING ============= ?????
               if( (neur[i].v.soma > 25.0) && neur[i].ready ) {
                  neur[i].tally++;
                  if( (fmod(dTime, 0.5000000) == 0.000) && (mp->neurons[i].tally > 0) )  // every 0.5 sec
                     mp->neurons[i].tally = 0;
               }
            }
/***************
            for( i=0; i<iN; i++ ) {
               if( mp->neurons[i].ready ) {
                  if( mp->neurons[i].v.soma > 25.0 ) {   // threshold = 25.0
                     if( flist.ef ) fprintf( flist.ef, "%f %f ", dTime, mp->ef[0].applied );
                     for( i=0; i<iN; i++ ) {
                        if( flist.ef ) fprintf( flist.ef, "%f ", mp->neurons[i].v.soma );
                     }
                     if( flist.ef ) fprintf( flist.ef, "\n" );
                     mp->neurons[i].ready = 0;
                  }
               }
               else if( !mp->neurons[i].ready ) {
                  if( mp->neurons[i].v.soma < 0.0 ) mp->neurons[i].ready = 1;
               }
            }
********************/
         }
      }
      model_derivatives( em_dTime, dy, dydt, (void *)dp );
      rk4m( dy, dydt, iN, N_VARS, em_dTime, DelTime, model_derivatives, (void *)dp ); //???
//    em_dTime += DelTime;                        // in msec
//    dTime = em_dTime/1000.0;                    // in seconds
      (void )update_neurons( neur, dy, 0, iN-1 ); // update paramaters with variables
//      if( fmod( dTime, 0.50000 ) < timestep ) {   // was 0.00005 9/17/2009
//         fprintf( f_log,"\n%f sec", dTime );
//         fflush ( f_log );
//      }
/*
       // temporary: raise membrane potential at dendrite v2.1
       if( mp->pspec->spike && dTime >= mp->pspec->tstart   // at the start of propagation analysis
                            && dTime <  mp->pspec->tstart + timestep ) {
//        dy[0][V_DEND ] = 60;              // raise above threshold
          dy[0][V_SOMA ] = 60;              // raise above threshold
          printf( "raised to 60mV\n" );
       }
*/
      //===================== Process each neuron ==============================
      //
      for( i=0; i<iN; i++ ) {
         neur[i].v.last      = neur[i].v.soma;
         mp->ef[i].potential = neur[i].v.dso - mp->ef[i].applied; // field potential. applied assigned?

         //
         // (temporary patch) - If analyze propagation, turn pulse back on
         //
//       if( mp->pspec->spike && dTime > mp->pspec->tstart) { // if analyze propagation and have started
//          mp->ip[0].index = 0



         //===================== Apply current pulse ===============================
         //
         if( mp->ip[i].index < mp->ip[i].count        &&  // if more pulses to apply
             dTime > mp->ip[i].onset[mp->ip[i].index] &&  // and time for the current pulse
             dTime < mp->ip[i].onset[mp->ip[i].index]+timestep ) {
            if( mp->ip[i].site == 's' ) { 
               neur[i].inj.soma  += mp->ip[i].magn;
               neur[i].inj.somaa  = neur[i].Area * neur[i].inj.soma;
            }
            else {
               neur[i].inj.dend  += mp->ip[i].magn;
               neur[i].inj.denda  = neur[i].Area * neur[i].inj.dend;
            }
            fprintf( f_log, "\n%f sec   applying current pulse of %f uA/cm^2 for neuron %i at %s",
                     dTime, mp->ip[i].magn, i+1, (mp->ip[i].site == 's' ? "soma" : "dendrite" ) );
            for( j=0; j<iN; j++ ) {
               CInj[i]+= Finv[0*5+1][(j*5)+3] * neur[j].inj.denda
                       + Finv[0*5+1][(j*5)+4] * neur[j].inj.somaa;
            }
            CRInj[i] = em_dR[0].DS * CInj[i];
         } 
         //===================== Remove current pulse =============================
         //
         if( mp->ip[i].index < mp->ip[i].count                      &&
             dTime > mp->ip[i].onset[mp->ip[i].index]+mp->ip[i].dur &&
             dTime < mp->ip[i].onset[mp->ip[i].index]+mp->ip[i].dur+timestep ) {
            if( mp->ip[i].site == 's' ) { 
               neur[i].inj.soma  -= mp->ip[i].magn;
               neur[i].inj.somaa  = neur[i].Area * neur[i].inj.soma;
            }
            else {
               neur[i].inj.dend  -= mp->ip[i].magn;
               neur[i].inj.denda  = neur[i].Area * neur[i].inj.dend;
            }
            fprintf( f_log, "\n%f sec   removing current pulse of %f uA/cm^2 for neuron %i at %s",
                     dTime, mp->ip[i].magn, i+1, (mp->ip[i].site == 's' ? "soma" : "dendrite" ) );
            for( j=0; j<iN; j++ ) {
               CInj[i]+= Finv[0*5+1][(j*5)+3] * neur[j].inj.denda
                       + Finv[0*5+1][(j*5)+4] * neur[j].inj.somaa;
            }
            CRInj[i] = em_dR[0].DS * CInj[i];
            mp->ip[i].index++;                             // point to next current pulse
         }
      }
      //======================= Increment Electric Field =========================
      //
      if( opt_if && ifcount < mp->ief->count ) {
         if( dTime >=  mp->ief->ifs[ifcount].onset &&
             dTime <= (mp->ief->ifs[ifcount].onset + timestep) ) {
            fprintf( f_log, "\n%f sec   applying E field of %f mV/cm (increment)",
                            dTime, mp->ief->ifs[ifcount].magn );
            for( i=0; i<iN; i++ ) {
               mp->ef[i].applied = mp->ief->ifs[ifcount].magn;
            } 
            ifcount++;
         }
      }
      //===================== Record spike propagation ========================
      //
      if( flist.sp || mp->pspec->spike ) {         // if record spikes or analyze propagate
         for( i = 0; i<iN; i++ ) {
            if( mp->neurons[i].v.soma > SP_THRESH &&            // if crosses threshold and
                dTime > mp->neurons[i].Spike + SP_REFRACT ) {   // passed refractory period (...if not first spike)

               if( flist.sp ) {                       // write spike time
                  fprintf( flist.sp, "%i %f\n", (i+1)-iiN, dTime );  // v2.05
      //          fprintf( flist.sp, "%i %f\n", (i+1), dTime );
               }
               //
               // If analyze propagation, after start-analysis time and spike not yet recorded,
               // record one spike per neuron and compute instantaneous speed
               //
               if( mp->pspec->spike && dTime > mp->pspec->tstart &&
                  !mp->neurons[i].Spike0 ) {

                  if( i ) {
                     dtemp = dTime - mp->neurons[i-1].Spike;
                     mp->neurons[i].speed  = (dtemp < SMALL)? 99999999.0 : (1 / dtemp);
                  }
                  mp->neurons[i].Spike0 = 1;
                  mp->termneuron = (i+1)-iiN;                // v2.05                        
                  if( f_prop ) {
                     fprintf( f_prop, "%i %f %f\n", mp->termneuron, dTime, mp->neurons[i].speed );
                     fflush ( f_prop );
                     if( autostop && i+1 == mp->tneurons ) { // move code to function later
                        fprintf( f_prop, "\n");
                        fprintf( f_log,  "\nTerminates at last spike\n");
                        // Close opened files 
                        for( i=0; filetb[i].name; i++ ) {   // do not close standard output
                           if( *filetb[i].open && *filetb[i].ptr ) fclose( *(filetb[i].ptr) );
                        }
                        return( 1 );
                     }                  
                  }
               }
               mp->neurons[i].Spike  = dTime;                    // current spike
            }
         }
      }
      if( dTime >= outspec.start && dTime <= outspec.end ) {
         if( !outspec.interv || fmod( dTime, interval ) < timestep ) {
            write_output( outfiles, dTime, neur, 0, iN-1 );
         }
      }
   }                                // end of while loop -- total time is done
   if( mp->pspec->speed ) {         // fit spikes to compute speed
      (void )compute_speed( mp );
      fprintf( f_log, "Propagation speed (neuron %i to %i) = %f +/- %f\n",
                       mp->pspec->nstart+1, mp->pspec->nend+1, mp->speed, mp->speed_sd );
   }
   //
   if( f_prop ) fprintf( f_prop, "\n");
   fprintf( f_log,  "\n");

   // Close opened files 
   for( i=0; filetb[i].name; i++ ) {   // do not close standard output
      if( *filetb[i].open && *filetb[i].ptr ) fclose( *(filetb[i].ptr) );
   }
   return( 1 );
}                                               // end of model() 

//=============================================================================
//
// function - double *vr_postassium();
//
// Compute reversal potential of potassium given extracellular potassisum concentration. 
//
double vr_potassium( double ko )
{
   return( -15.0 + 61.15*log10(ko/8.50) );
}                                               // end of vr_potassium()

//=============================================================================
//
// function - neuron *update_neurons();
//
// nu     - pointer to array of neurons
// var    - pointer to matrix of number of neurons x number of variables
// nstart - start with this neuron, index to array nu to start
// nend   - end   with this neuron, index to array nu to end
//
// For all neurons, nstart = 0, nend = tot_neurons-1
//
neuron *update_neurons( neuron *nu, double **var, int nstart, int nend )
{
   int i = 0;
   for( i = nstart; i <= nend; i++ ) {
      nu[i].k.h      = var[i][GATE_H ];
      nu[i].k.n      = var[i][GATE_N ];
      nu[i].k.s      = var[i][GATE_S ];
      nu[i].k.c      = var[i][GATE_C ];
      nu[i].k.q      = var[i][GATE_Q ];
      nu[i].k.Ca     = var[i][CONC_CA];
      nu[i].SNMDA    = var[i][I_SNMDA];
      nu[i].WAMPA    = var[i][I_WAMPA];
      nu[i].v.soma   = var[i][V_SOMA ];
      nu[i].v.dend   = var[i][V_DEND ];
   }
   return( nu );
}                                               // end of update_neurons()

//=============================================================================
//
// function - Matrix update_variables();
//
// var    - pointer to matrix of number of neurons x number of variables
// nu     - pointer to array of neurons
// nstart - start with this neuron, index to array nu to start
// nend   - end   with this neuron, index to array nu to end
//
// For all neurons, nstart = 0, nend = tot_neurons-1
//
Matrix update_variables( Matrix var, neuron *nu, int nstart, int nend )
{
   int i = 0;
   for( i = nstart; i <= nend; i++ ) {
      //
      // original comment:
      // The mid point of "hgate" is needed to be shifted,
      // because the max of hgate should not be greater than 1;
      // max=1, min=0., since hgate is "inactivation" probability:
      // non-permissive state of gating
      //
      var[i][GATE_H]  = nu[i].k.h;
      var[i][GATE_N ] = nu[i].k.n;   
      var[i][GATE_S ] = nu[i].k.s;
      var[i][GATE_C ] = nu[i].k.c;
      var[i][GATE_Q ] = nu[i].k.q;
      var[i][CONC_CA] = nu[i].k.Ca;
      var[i][I_SNMDA] = nu[i].SNMDA;
      var[i][I_WAMPA] = nu[i].WAMPA;
      var[i][V_SOMA ] = nu[i].v.soma;
      var[i][V_DEND ] = nu[i].v.dend;
   }
   return( var );
}                                               // end of update_variables()

//=============================================================================
//
// function - Matrix *connect_synapses();
//
//
// isyn   - pointer to matrix of number of neurons x number of neurons
// nstart - start with this neuron, index to array nu to start
// nend   - end   with this neuron, index to array nu to end
//
// For all neurons, nstart = 0, nend = tot_neurons-1
//
Matrix connect_synapses( double **isyn, int nstart, int nend )
{
   int i=0, j=0;
   for( i = nstart; i <= nend; i++ ) {
      for( j = nstart; j <= nend; j++ ) {
         if( i == j )
            isyn[i][j] = 0;
         else if( abs(i-j) < 2 )     // was fabs
            isyn[i][j]=1;
         else
            isyn[i][j]=0;
      }
   }
   return( isyn );
}                                               // end of update_synapses()

//=============================================================================
//
// function - int *compute_speed();
//
int compute_speed( modelparams *m )
{
   int nn = m->pspec->nend - m->pspec->nstart + 1;
   double sptime[nn], spnum[nn];
   int i=0, j=0;
   double chi, q, b, sb;     // dummies  

   fprintf( f_log, "\n" );
   for( j=0, i = m->pspec->nstart; i <= m->pspec->nend; j++, i++ ) {
      spnum[j]  = (double )i;
      sptime[j] = m->neurons[i].Spike;         // spike time
      fprintf( f_log, "%i:%i::%f:%f ", j, i, sptime[j], spnum[j] );
   }
   fprintf( f_log, "\n" );
   fit( sptime, spnum, nn, &(m->speed), &b, &(m->speed_sd), &sb, &chi, &q );
   return( nn );
}                                               // end of compute_speed()

//============================================================================
//
// function - double ratio();
//
double ratio( rgrid *em_dR, double k_out, rratio *rr, neuron *neur, int nn )
{
   int i = 0;
   double DS = 0.0, DD = 0.0, SS = 0.0;
   double dRNominal = 0.0;
   //
   // Nominal Values:
   // dR->rDS = 0.1; k_out = 8.5 => DS = 0.15, DD = SS = 0.015
   //           0.1;         3.5 =>      0.10            0.010
   //              
   DS = 0.1 + (k_out-3.5)*0.05/5.0;
   DD = rr->nnR * DS;
   SS = DD;

   if( f_trace ) fprintf( f_trace, "ratio() - entry. DS = %f  DD = %f  SS = %f\n",
                                           DS,  DD,  SS );

   dRNominal = DS/( neur[0].g.c*neur[0].Area );

   for( i=0; i<nn; i++ ) {
      em_dR[i].TD = rr->tgR * dRNominal;
      em_dR[i].DD = (DD/DS) * dRNominal;
      em_dR[i].DS =           dRNominal;
      em_dR[i].SS = (SS/DS) * dRNominal;
      em_dR[i].SG = rr->tgR * dRNominal;
   }
   //--------------------------------------------------------------------------
   //	(2)	Give terminal resistance associated with boundary
   //		TERMINAL: LEFT BOUNDARY AND RIGHT BOUNDARY
   em_dR[nn+0].TD = rr->aaZ * em_dR[ 0].TD;
   em_dR[nn+1].TD =           em_dR[nn].TD;
   em_dR[nn+0].DS = rr->aaZ * em_dR[ 0].DS;
   em_dR[nn+1].DS =           em_dR[nn].DS;
   em_dR[nn+0].SG = rr->aaZ * em_dR[ 0].SG;
   em_dR[nn+1].SG =           em_dR[nn].SG;
   em_dR[nn+0].DD = rr->bbZ * em_dR[ 0].DD;
   em_dR[nn+1].DD =           em_dR[nn].DD;
   em_dR[nn+0].SS = rr->bbZ * em_dR[ 0].SS;
   em_dR[nn+1].SS =           em_dR[nn].SS;
   return( 1 );
}                                               // end of ratio()

//============================================================================
//
// function - network()
//
double network( rgrid *em_dR, double **FMatrix, int nl, int nn )
{
   int i=0, j=0, iOff=0;
  
   if( f_trace ) {
      fprintf( f_trace, "network() - entry. TD, DD, DS, SS, SG\n" );
      for( i=0; i<nn+2; i++ ) fprintf( f_trace, "     %-3i:  %g %g %g %g %g\n", i,
               em_dR[i].TD, em_dR[i].DD, em_dR[i].DS, em_dR[i].SS, em_dR[i].SG );
   }
   //		RESISTIVE NETWORK BODY

   for(i=0; i<nn; i++) {
      iOff=5*i;
      //-- vertical current ----
      j = iOff;
      FMatrix[j][iOff]   =  em_dR[i  ].TD;
      FMatrix[j][iOff+1] =  em_dR[i  ].DS;
      FMatrix[j][iOff+2] =  em_dR[i  ].SG;
      j = iOff + 1;
      FMatrix[j][iOff]   =  em_dR[i  ].TD;
      FMatrix[j][iOff+5] = -em_dR[i+1].TD;
      //---- horizontal current ----
      FMatrix[j][iOff+8] =  em_dR[i+1].DD;
      //----------------------------
      j=iOff + 2;
      FMatrix[j][iOff+2] =  em_dR[i  ].SG;
      FMatrix[j][iOff+7] = -em_dR[i+1].SG;
      //--horizontal current--------
      FMatrix[j][iOff+9] = -em_dR[i+1].SS;
      //----------------------------
      j=iOff + 3;
      FMatrix[j][iOff]   =  1.0;
      FMatrix[j][iOff+1] = -1.0;
      FMatrix[j][iOff+3] =  1.0;
      FMatrix[j][iOff+8] = -1.0;
      //----------------------------
      j=iOff+4;
      FMatrix[j][iOff+1] =  1.0;
      FMatrix[j][iOff+2] = -1.0;
      FMatrix[j][iOff+4] =  1.0;
      FMatrix[j][iOff+9] = -1.0;
   }
   //-------------- RIGHT BOUNDARY at i=nn=N_NEURONS
   i=nn;
   iOff = 5*i;
   j=iOff;
   FMatrix[j][iOff]   = em_dR[i].TD;
   FMatrix[j][iOff+1] = em_dR[i].DS;
   FMatrix[j][iOff+2] = em_dR[i].SG;

   j=iOff+1;
   FMatrix[j][iOff]   =  1.0;
   FMatrix[j][iOff+1] = -1.0;
   FMatrix[j][iOff+3] =  1.0;

   j=iOff+2;
   FMatrix[j][iOff+1] =  1.0;
   FMatrix[j][iOff+2] = -1.0;
   FMatrix[j][iOff+4] =  1.0;
   //
   //----------- LEFT BOUNDARY at i=N_NEURONS effectively at i=-1
   i=nn+1;
   iOff = 5*i;

   j = 5*nn+3;
   FMatrix[j][iOff]   = em_dR[i].TD;
   FMatrix[j][iOff+1] = em_dR[i].DS;
   FMatrix[j][iOff+2] = em_dR[i].SG;

   j = 5*nn+4;
   FMatrix[j][iOff]   = em_dR[i].TD;
   FMatrix[j][0]      =-em_dR[0].TD;
   FMatrix[j][3]      = em_dR[i].DD;

   j = 5*nn+5;
   FMatrix[j][iOff+2] = em_dR[i].SG;
   FMatrix[j][2]      =-em_dR[0].SG;
   FMatrix[j][4]      =-em_dR[i].SS;

   j = 5*nn+6;
   FMatrix[j][iOff]   = 1.0;
   FMatrix[j][iOff+1] =-1.0;
   FMatrix[j][3]      =-1.0;

   j = 5*nn+7;
   FMatrix[j][iOff+1] = 1.0;
   FMatrix[j][iOff+2] =-1.0;
   FMatrix[j][4]      =-1.0;
  return (1);                 //???
}                                           // end of network()

//=============================================================================
//
// function - double model_derivatives();
//
double model_derivatives( double dte, Matrix dy, Matrix dydt, void *vp )
{
   int i=0, j=0;
   neuron      *nu  = (neuron *)0;
   derivparams *p   = (derivparams *)vp;
   int          nn  = p->tneurons;
   gatevars    *gk  = p->gatekin;
   double       H_S_Sum=0.0, H_W_Sum=0.0;
   double       rinf = 0.0, dChiCa = 0.0;

// for( i=0; i<nn; i++ ) p->vds[i] = 0.0;
   for( i=0; i<nn; i++ ) p->neurons[i].v.dso  = 0.0;

   for( i=0; i<nn; i++ ) {
      nu = &(p->neurons[i]);
      // M gate                                 // comupute all gate variables here
      Gate_m( &(gk->m), dy[i][V_SOMA] );
      //
      // Compute total somatic transmembrabe current
      //
      nu->i.soma = nu->P*(
           (nu->g.L                                  *(dy[i][V_SOMA]-nu->v.L ) )
         + (nu->g.Na*pow(gk->m.inf,2.0)*dy[i][GATE_H]*(dy[i][V_SOMA]-nu->v.Na) )
         + (nu->g.KDR                  *dy[i][GATE_N]*(dy[i][V_SOMA]-nu->v.K ) ) );
      //
      // Compute total dendritic transmembrane current
      //
      rinf = 1.0/(1.0 + 0.28*exp(-0.062*( dy[i][V_DEND] - nu->v.syn )));

      if( dy[i][CONC_CA] < 250.0 ) dChiCa = dy[i][CONC_CA]/250.0;
      else                         dChiCa = 1.0;

      nu->i.nmda  = (1.0 - nu->P)*(                             // NMDA current
           (nu->g.NMDA     * dy[i][I_SNMDA]*rinf * (dy[i][V_DEND]-nu->v.syn) ) );

      nu->i.ampa  = (1.0 - nu->P)*(                             // AMPA current
           (nu->g.AMPA     * dy[i][I_WAMPA]      * (dy[i][V_DEND]-nu->v.syn) ) );

      nu->i.syn  = nu->i.nmda + nu->i.ampa;                     // total synatpic current

      nu->i.dend = nu->i.syn + (1.0 - nu->P)*(
         + (nu->g.L                              * (dy[i][V_DEND]-nu->v.L  ) )
         + (nu->g.Ca   * pow(dy[i][GATE_S],2.0)  * (dy[i][V_DEND]-nu->v.Ca ) )
         + (nu->g.KAHP      *dy[i][GATE_Q]       * (dy[i][V_DEND]-nu->v.K  ) )
         + (nu->g.KCa*dChiCa*dy[i][GATE_C]       * (dy[i][V_DEND]-nu->v.K  ) ) );
      //
      // If no resistive grid, compute V-dsout here.   v2.0 v2.1
      // AA and BB are constants pre-calculated when resistances were computed.
      //
      if( no_r ) {
         nu->v.dso = ( ( p->r[i].TD * nu->inj.denda ) - ( p->r[i].SG * nu->inj.somaa )
                     - ( p->r[i].BB * ( dy[i][V_DEND] - dy[i][V_SOMA] ) )
                     +   p->ef[i].applied )
                     / ( p->r[i].AA + p->r[i].BB );
      }
   }
   //------------------------------------------------------------------------
   //		               Compute V-dsout
   // v2.0
   // If resistive grid is included computed V-dsout with matrix
   // otherwise, V-dsout is computed in the above loop. 
   //------------------------------------------------------------------------
   if( !no_r ) {
      for( i=0; i<nn; i++ ) {
         p->vds[i] = 0.0;
         for( j=0; j<nn; j++ ) {
            p->vds[i] += p->rg[i][j]*( dy[j][V_DEND] - dy[j][V_SOMA] );
         }
      }
      for( i=0; i<nn; i++ ) p->sum[i] = p->vinj[i] + p->vds[i] + p->ef[i].applied/25.0;
      for( i=0; i<nn; i++ ) {
         p->neurons[i].v.dso = 0.0;
         for( j=0; j<nn; j++ ) {
            p->neurons[i].v.dso += p->yrginv[i][j]*p->sum[j];
         }
      }
   }
   //--------------------------------------------------------------------------
   //		     Compute the derivatives for each neuron
   //--------------------------------------------------------------------------
   for(i=0; i<nn; i++) {
      nu = &(p->neurons[i]);
      // H gate
      Gate_h( &(gk->h), dy[i][V_SOMA] );
      dydt[i][GATE_H] = (gk->h.inf-dy[i][GATE_H])/gk->h.tau;
      // N gate
      Gate_n( &(gk->n), dy[i][V_SOMA] );
      dydt[i][GATE_N] = (gk->n.inf-dy[i][GATE_N])/gk->n.tau;
      // S gate
      Gate_s( &(gk->s), dy[i][V_DEND] );
      dydt[i][GATE_S] = (gk->s.inf-dy[i][GATE_S])/gk->s.tau;
      // C gate
      Gate_c( &(gk->c), dy[i][V_DEND] );
      dydt[i][GATE_C] = (gk->c.inf-dy[i][GATE_C])/gk->c.tau;
      // Q gate
      Gate_q( &(gk->q), dy[i][CONC_CA] );
      dydt[i][GATE_Q] = (gk->q.inf-dy[i][GATE_Q])/gk->q.tau;
      // Ca concentration
      // Traub (pow(s,5))---------------------------------------------------------
      //    dydt[i][CONC_CA] = -0.075*dy[i][CONC_CA]
      //           -0.13*nu->g.Ca*pow(dy[i][GATE_S],5)*(dy[i][V_DEND]-nu->v.Ca);
      // Pinsky-Rinzel (pow(s,2))-------------------------------------------------
      dydt[i][CONC_CA]= -0.075*dy[i][CONC_CA]
            -0.13*nu->g.Ca*pow(dy[i][GATE_S],2.0)*(dy[i][V_DEND]-nu->v.Ca);
      //
      // S and W synaptic derivs -------------------------------------------------
      //
      H_S_Sum = 0.0;
      H_W_Sum = 0.0;
      for( j=0; j<nn; j++ ) {
         if( p->syn[i][j] != 0 ) {
            if( dy[j][V_SOMA] >= 10.0 ) H_S_Sum = H_S_Sum + 1.0;
            if( dy[j][V_SOMA] >= 20.0 ) H_W_Sum = H_W_Sum + 1.0;
         }
      }
      dydt[i][I_SNMDA] = H_S_Sum - (dy[i][I_SNMDA]/150.0);

      // SATURATION OF AMPA

      if( (dy[i][I_SNMDA] > 125.0) && (dydt[i][I_SNMDA] > 0) )
         dydt[i][I_SNMDA] = 0.0;
      dydt[i][I_WAMPA] = H_W_Sum - (dy[i][I_WAMPA]/2.0);
      //
      // Derivative for somatic and dendritic transmembrane potentials 
      //
      nu->i.dsi = nu->g.c*
         (dy[i][V_DEND] + nu->v.dso - dy[i][V_SOMA]);   // internal dendrite-to-soma current
      dydt[i][V_SOMA] = -(nu->i.soma - nu->i.dsi - nu->inj.soma)/(nu->cm*nu->P );
      dydt[i][V_DEND] = -(nu->i.dend + nu->i.dsi - nu->inj.dend)/(nu->cm*(1.0-nu->P));
   }
   return(1);
}                                               // endof model_derivatives()

//============================================================================
//
// function - kinetics *Gate_m(), *Gate_h(), *Gate_n(), *Gate_s(), *Gate_c(), *Gate_q(), 
//            double   *zeta();
//
// Kinetic functions for computing gate kinetic variable.
//
kinetics *Gate_m( kinetics *var, double volt )
{
   double adDel, bdDel;

   adDel      = (13.10-volt)/4.0;
   bdDel      = (volt-40.10)/5.0;
   var->alpha = 4.0*0.320*zeta(adDel);
   var->beta  = 5.0*0.280*zeta(bdDel);
   var->tau   = 1.0 / (var->alpha + var->beta);
   var->inf   = var->alpha * var->tau;
   return( var );
}                                               // end of Gate_m()
//
kinetics *Gate_h( kinetics *var, double volt )
{
   var->alpha = 0.1280*exp((17.0-volt)/18.0);
   var->beta  = 4.00/(1.0+exp((40.0-volt)/5.0) );
   var->tau   = 1.0 / (var->alpha + var->beta);
   var->inf   = var->alpha * var->tau;
   return( var );
}                                               // end of Gate_h()
//
kinetics *Gate_n( kinetics *var, double volt )
{
   var->alpha = 0.0160 * 5.0*zeta((35.1-volt)/5.0);
   var->beta  = 0.250 * exp(0.50-0.0250*volt);
   var->tau   = 1.0 / (var->alpha + var->beta);
   var->inf   = var->alpha * var->tau;
   return( var );
}                                               // end if Gate_n()
//
kinetics *Gate_s( kinetics *var, double volt )
{
   var->alpha = 1.6000/(1.00+exp(-0.072*(volt-65.0)));
   var->beta  = 0.020*5.00*zeta( (volt-51.1)/5.0 );
   var->tau   = 1.0 / (var->alpha + var->beta);
   var->inf   = var->alpha * var->tau;
   return( var );
}                                               // end of Gate_s()
//
kinetics *Gate_c( kinetics *var, double volt )
{
   double c_inf_const;

   c_inf_const = 0.50*exp(40.0/11.0);
   var->tau    = 0.50*exp((volt-6.50)/27.0);
   if( volt <= 50.0 ) var->inf = 0.50*exp((volt-10.0)/11.0)/c_inf_const;
   else               var->inf = 1.00;
   return( var );
}                                               // end of Gate_c()
//
kinetics *Gate_q( kinetics *var, double conc )
{
   if( conc < 500.0 ) var->alpha = conc/50000.0;
   else               var->alpha = 0.010;
   var->beta = 0.0010;
   var->tau  = 1.00/(var->alpha + var->beta);
   var->inf  = var->alpha * var->tau;
   return( var );
}                                               // end of Gate_q()
//
double zeta( double dDel )
{
   double zeta_var = 0.0;
   double a =  0.010/(exp( 0.010) - 1.0);       // make this a constant
   double b = -0.010/(exp(-0.010) - 1.0);       // make this a constant
   double em_zeta_slope = (a-b)/(2.00*0.010);
   //--new--- adding .d0: 0.0001d0 and 1.0d0-----------------------------
   if((dDel*dDel) > 0.000100) zeta_var = dDel/(exp(dDel) - 1.00); 
   else                       zeta_var = 1.0+em_zeta_slope*dDel;
   return( zeta_var );
}                                               // end of zeta()
//
//======================= end of kinetic functions ===========================

//============================================================================
//
// function - Matrix rk4m();
//
// Fourth Order Runge-Kutta Methos for a Matrix
//
Matrix rk4m( Matrix y, Matrix dydx, int r, int c, double x, double h,
             double (*derivs)(), void *dp )
{
   int i=0, j=0;
   double xh, hh, h6;
   Matrix dym, dyt, yt;
   
   dym = matrix( r, c );
   dyt = matrix( r, c );
   yt  = matrix( r, c );
   hh  = h*0.5;
   h6  = h/6.0;
   xh  = x + hh;      // dTE
   for( i=0; i<r; i++ ) for( j=0; j<c; j++ ) yt[i][j]=y[i][j]+hh*dydx[i][j];
   (*derivs)( xh, yt, dyt, dp );
   for( i=0; i<r; i++ ) for( j=0; j<c; j++ ) yt[i][j]=y[i][j]+hh*dyt[i][j];
   (*derivs)( xh, yt, dym, dp );
   for( i=0; i<r; i++ ) for( j=0; j<c; j++ ) {
      yt [i][j]  =   y[i][j] + h*dym[i][j];
      dym[i][j] += dyt[i][j];
   }
   (*derivs)( x+h, yt, dyt, dp );
   for( i=0; i<r; i++ ) for( j=0; j<c; j++ )
      y[i][j] = y[i][j] + h6*( dydx[i][j] + dyt[i][j] + 2.0*dym[i][j] );  // y=yout	
   free_matrix( (void *) yt );
   free_matrix( (void *)dyt );
   free_matrix( (void *)dym );
   return( y );
}                                               // end of rk4m()

//============================================================================
//
// function - double LUDCMP();
//
double LUDCMP( double **a, int n, int *indx ) {

   int i, imax = 0, j, k;
   double big, dum, sum = 0.0, temp, *vv, d;

   vv = (double *)malloc( sizeof(double)*n );
   d=1.0;
   for (i=0; i<n; i++) {
      big=0.0;
      for (j=0; j<n; j++) {
         temp = fabs( a[i][j] );
         if(temp > big)	big=temp;
      }
      if(big == 0.0)
          fprintf( f_log, "singular matrix in LUDCMP.\n" );
      vv[i]=1.0/big;
   }
   for( j=0; j<n; j++ ) {
      if( j>0 ) {                    // ?
         for( i=0; i<=j-1; i++ ) {
            sum=a[i][j];
            if( i>0 ) {              // ?
               for( k=0; k<i; k++ ) sum -= a[i][k] * a[k][j];
               a[i][j] = sum;
            }
         }
      }
      big = 0.0;
      for( i=j; i<n; i++ ) {
         sum=a[i][j];
         for( k=0; k<(j); k++ ) sum -= a[i][k]*a[k][j];
         a[i][j] = sum;
         dum = vv[i]*fabs(sum);
         if( dum >= big ) {
            big = dum;
            imax = i;
         }
      }
      if( j != imax ) {
         for( k=0; k<n; k++ ) {
            dum=a[imax][k];
            a[imax][k]=a[j][k];
            a[j][k]=dum;
         }
         d = -d;
         vv[imax] = vv[j];
      }
      indx[j] = imax;
      if( j != n-1 ) {                           // ? 
         if( a[j][j] == 0.0 ) a[j][j] = TINY;
         dum = 1.0/a[j][j];
         for( i=j+1; i<n; i++ ) a[i][j] *= dum;
      }
   }
   if( a[n-1][n-1] == 0.0 ) a[n-1][n-1] = TINY;  // ?
   free( vv );
   return d;
// return a;
}                                               // end of LUDCMP();

//============================================================================
//
// function - double LUBKSB();
//
double *LUBKSB( double **a, int n, int *indx, double *b )
{
   int i, ii, j, ip;
   double sum = 0.0;

   ii=-1;
   for( i=0; i<n; i++ ) {
      ip    = indx[i];
      sum   = b[ip];
      b[ip] = b[i];
      if( ii != -1 )
         for( j=ii; j<=(i-1); j++ ) sum -= a[i][j]*b[j];
      else if(sum != 0.0) ii=i;
      b[i]=sum;
   }
   for( i=n-1; i>=0; i-- ) {
      sum=b[i];
      if( i<(n-1) )             // ?
         for( j=i+1; j<n; j++ ) sum -= a[i][j] * b[j];
      b[i] = sum/a[i][i];
   }
   return( b );
}                                               // end of LUBSKB()

//============================================================================
//
// function - double **identityMatrix();
//
Matrix identityMatrix( Matrix m, int n )
{
   int i = 0, j = 0;
   for( i=0; i<n; i++ )
      for( j=0; j < n; j++ ) {
         if( i==j ) m[i][j]=1.0;
         else       m[i][j]=0.0;
   }
   return( m );
}                                               // end of identityMatrix()

//============================================================================
//
// function - double **inverseMatrix();
//
Matrix inverseMatrix( Matrix x, Matrix m, int n )
{                           
   int i = 0, j = 0;
   int    *index  = (int    *)malloc( sizeof(int   )*n );
   double *invcol = (double *)malloc( sizeof(double)*n ); // inverse column 

   (void )memset( index,  0, sizeof(int   )*n );
   (void )memset( invcol, 0, sizeof(double)*n );

   (void )LUDCMP( m, n, index );            // decompose matrix just once
   for( j=0; j<n; j++ ) {                        // find inverse by columns
      for ( i=0; i<n; i++ ) invcol[i] = 0.0;
      invcol[j] = 1.0;
      LUBKSB( m, n, index, invcol );
      for( i=0; i<n; i++ ) x[i][j] = invcol[i];
                                                // from Numerical Recipes in C
   }
   free( invcol );
   free( index  );
   return( x );
}                                               // end of inverseMatrix()

//============================================================================
//
// function - void fit()
//
// Fit to a straight without weight, y = ax + b
//
void fit( double *x, double *y, int ndata,
          double *a,    double *b,
          double *siga, double *sigb, double *chi2, double *q)
{
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
   float sqrarg;

   int i = 0;
   double t, sxoss, sx=0.0, sy=0.0, st2=0.0, sigdat;
   double ss = (double )ndata;

   *b=0.0;
   for (i=0;i<ndata;i++) {
      sx += x[i];
      sy += y[i];
   }
   sxoss=sx/ss;
   for (i=0;i<ndata;i++) {
      t=x[i]-sxoss;
      st2 += t*t;
      *b += t*y[i];
   }
   *b /= st2;
   *a=(sy-sx*(*b))/ss;
   *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
   *sigb=sqrt(1.0/st2);
   *chi2=0.0;
   *q=1.0;
   for (i=0;i<ndata;i++)
      *chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
   sigdat=sqrt((*chi2)/(ndata-2));   // -1, -2?
   *siga *= sigdat;
   *sigb *= sigdat;
}

//============================================================================
//
// function - open_out_files();
//
// File format: <dir>[/]<f_name>[<e-field><timestamp].txt
//
//              e-field:   in 10-microvolts, 6 digits; if <0 neg sign + 5 digits
//              timestamp: YYMMDDmmhhss
//              If <dir> does not end with '/', add '/'.
//              If filetimestamp is on, do not include timstamp
//              If nofilesuffix  is on, do not include e-field and timstamp
//
// future - crate directory if it does not exist
//
FILE *open_output_file( FILE **f_ptr, char *f_name, char *dir, double ef,
                      struct tm *bt,  char *f_note, char *f_descr )
{
   char      f_prefix[128] = "";
   char      f_suffix[ 64] = "";
   char      f_path  [128] = "";
   char      f_time  [ 14] = "";
   int       ef_sign       = 1;

   *f_ptr = (FILE *)0;
   // prepend dir, if it does not terminate with '/', add it
   //
   sprintf( f_prefix, "%s%s%s", dir,
                   ( (!dir[0] || dir[strlen(dir)-1] == '/') ? "" : "/" ), f_name );

   // append e-field value and timestamp
   //
   if( !nofilesuffix ) {
      if( filetimestamp ) {
         sprintf( f_time, "%02i%02i%02i%02i%02i%02i",
                          bt->tm_year-100, bt->tm_mon+1, bt->tm_mday,
                          bt->tm_hour,     bt->tm_min,   bt->tm_sec  );
      }
      if( ef < 0 ) ef_sign = -1;
      sprintf( f_suffix, "%07i%s", (int )((ef*1000.0)+(ef_sign*0.1)), f_time );
   }
   // build output file pathname
   //
   sprintf( f_path, "%s%s%s.txt", f_prefix, f_suffix, f_note );
   //
   // If file opens, write header, otherwise write error message
   //
   if( (*f_ptr = fopen( f_path, "w" ) ) ) {
      fprintf( *f_ptr, "#  nn=%d in=%d ef=%f mV/cm. %s\n",
                       num_neurons, tot_neurons-num_neurons, ef, f_descr );
      fprintf(  f_log, "%s: File %s opened.\n", pname, f_path ); 
   }
   else {
      fprintf( f_log, "%s: Failed to open file %s.\n", pname, f_path );
   }
   return( *f_ptr );
}
//============================================================================
//
// function - get_time();
//
// Parse command line to get command options.
//
char *get_time( char *str )
{
   time_t    systime;
   (void )time( &systime );                        // get time
   return( ctime_r( &systime, str ) );
}

//============================================================================
//
// function - get_command();
//
// Parse command line to get command options.
//
int get_command( int argc, char **argv )
{
   int   i = 0, j = 0;
   int   opt_index = 0;
   int   expect_value = 0;    // init to expect switch
   char *pswitch = (char *)0; // parameter switch 
   
   pname = argv[0];
   for( i=1; i<argc; i++ ) {
      if( !expect_value ) {
	 //
	 // Expect parameter switch, error if it does not start with '-' or
         // if a match not found in options[] table
         //
         if( *(argv[i]) != '-' ) {
            fprintf( f_log, "%s: Command line error at argument number %i <%s>. %s\n",
	            argv[0], i, argv[i], "Switch expected." );
            return( 0 );
         }
         for( j=0;         options[j].param &&
                   strcmp( options[j].param, argv[i] ); j++ );
         if( !options[j].param ) {
            fprintf( f_log, "%s: Command line error at argument number %i <%s>. %s\n",
	            argv[0], i, argv[i], "Invalid option." );
            return( 0 );
         }
         if( options[j].type == DT_FLAG ) {    // It is a flag,
	    *(int *)(options[j].value) = 1;    // so set flag on and
            expect_value = 0;                  // expect switch next
            opt_index    = -1;
         }
         else {                               // Not a flag
	    expect_value = 1;                 // so expect value next and
            opt_index = j;                    // save index to options[]
            pswitch = argv[i];
         }
      }
      else {
	 //
	 // Expect parameter value. If it starts with a '-',
         // make sure that it is not a switch.
         //
         if( *(argv[i]) == '-' ) {
	    for( j=0; options[j].param; j++ ) {
	       if( !strcmp( options[j].param, argv[i] ) ) { 
                  fprintf( f_log, "%s: Command line error at argument number %i <%s>. %s\n",
	                argv[0], i, argv[i], "Value expected." );
                  return( 0 );
               }
            }
         }
	 switch( options[opt_index].type ) {
            case DT_FLOAT:
               if( !sscanf( argv[i], "%lf", (double *)(options[opt_index].value) ) ) {
                  fprintf( f_log, "%s: Command line error at argument number %i <%s>. %s.\n",
                          argv[0], i, argv[i], "Enter floating point value" );
                  return( 0 );
               }
               break;
            case DT_INTEGER:                            // use scanf integer
               if( !sscanf( argv[i], "%i", (int *)(options[opt_index].value) ) ) {
                  fprintf( f_log, "%s: Command line error at argument number %i <%s>. %s.\n",
                           argv[0], i, argv[i], "Enter integer value" );

                  return( 0 );
               }
               break;
	    case DT_STRING:
	       if( strlen(argv[i]) > (int )options[opt_index].size ) {
                  fprintf( f_log, "%s: Command line error at argument number %i <%s>. %s%d.\n",
                          argv[0], i, argv[i], "Value exceeds mximum length of ",
	                  (int )options[opt_index].size-1 );
                  return( 0 );
               }
	       else {
		 strcpy( (char *)options[opt_index].value, argv[i] );
	       }
               break;
	    case DT_SPECIAL:
               if( options[opt_index].func &&
                  !options[opt_index].func( options[opt_index].value, argv[i], pswitch ) )
	          return( 0 );                  // error exit
               break;
	    default:
	       break;
         }
         expect_value = 0;                     // expect switch next;
      }
   }
   if( expect_value ) {
      fprintf( f_log, "%s: Command line error at argument number %i <%s>. %s\n",
	                argv[0], i-1, argv[i-1], "Command ends when value is expected." );
      return( 0 );
   }
   if( num_neurons < 0 ) {
      fprintf( f_log, "%s: -nn <number of neurons> must be specified.\n", argv[0] );
      return( 0 );
   }
   //
   // total number of neurons = initiating neurons + num neurons in expeimental chain
   //
   tot_neurons = num_neurons + efi.ineurons; // calc total # neurons before used in parse_i_pulse2() v2.05
   //
   // Parse current pulses (second pass) for each neuron
   //
   if( !(parse_i_pulse2( &ipspec ) ) ) return( 0 );
   //
   // Sematic error check
   //
   if( outspec.end == -1 ) outspec.end = dTotal;
   if( outspec.start < 0.0 || outspec.end < 0.0 || outspec.end < outspec.start ) {
      fprintf( f_log, "%s: Output range <-or> consistency error. %s%s.\n", argv[0],
               "Both start and end time must be positive and ",
               "end time must be greaer then start time" );
      return( 0 );
   }
   if( opt_ef && opt_if ) {
      fprintf( f_log, "%s: Both <-ef> and <-if> options are specified. %s.\n", argv[0],
               "Only one of the two options is permitted" );
      return( 0 );
   }
   if( no_gnn && no_rnn ) {
      fprintf( f_log, "%s: Both <-ngnn> and <-nrnn> options are specified. %s.\n", argv[0],
               "Only one of the two options is permitted" );
      return( 0 );
   }
   if( *initcond && *outequil ) {           // v2.2
      fprintf( f_log, "%s: Both <-ic> and <-oe> options are specified. %s.\n", argv[0],
               "Only one of the two options is permitted" );
      return( 0 );
   }
   return( 1 );
}
//
// function: parse_options();
//
// Parse argument to get list of options.
// If a keyword (token) in the tokenlist is found, set flag in that entry of tokenlist.
//
// input: tb  -- pointer to tokenlist
//        str -- input string to be parsed
//        sw  -- switch for this set of options
//
int parse_options( void *tb, char *args, char *sw )
{
   int i=0;
   char *token = "";
   char  str[128];
   tokenlist *optb = (tokenlist *)tb;
   
   strcpy( str, args );                    // preserve argument string
   for( token=strtok(str,","); token; token=strtok( (char *)0, "," ) ) {
      for( i=0; optb[i].tok && strcmp( optb[i].tok, token ); i++ );     
      if( !optb[i].tok ) {
         fprintf( f_log, "%s: Invalid option <%s> value <%s>.\n", pname, sw, token );
         return( 0 );
      }
      else {
         *(optb[i].flag) = 1;
      }
   }
   return( 1 );
}
//
// function: parse_outrange();
//
// Parse argument to patameters for applied field and its onset times.
// Format: -or <start>[:<end>]
//
// input: tb  -- table to store parsed values
//        str -- input sring to be parsed
//        sw  -- switch for this set of options
//
int parse_outrange( void *tb, char *args, char *sw )
{
   char *token = "";
   char  str[128];
   outputspec *or = (outputspec *)tb;

   strcpy( str, args );                    // preserve argument string
   if( !(token=strtok( str, ":" ) ) ) {    // start time
      fprintf( f_log, "%s: Output range <%s> parameters missing. %s.\n",
                      pname, sw,
                      "Format is <starttime>[:<end time>]" );
      return( 0 );
   }
   else if( !sscanf( token, "%lf", &(or->start) ) ) {
      fprintf( f_log, "%s: Output range <%s> start time <%s> error. %s.\n",
                      pname, sw, token, "Enter floating point value in seconds" );
      return( 0 );
   }
   //
   if( !(token=strtok( (char *)0, ":" ) ) ) {    // get end time
      or->end = -1;                              // to be set to simulation duration
   }
   else if( !sscanf( token, "%lf", &(or->end) ) ) {
      fprintf( f_log, "%s: Output range <%s> end time <%s> error. %s.\n",
                      pname, sw, token, "Enter floating point value in seconds" );
      return( 0 );
   }
   return( 1 );
}
//
// function: parse_i_pulse();
//
// Parse argument to specifiaction for current pulses - First Pass:
// Format: -ip <neuron spec>,<site>,<magnitude>,<duration>,<onset time spec>
//         where <neuron spec> specifies neurons to which this spec applies
//         where <site> = s or d
//         where <onset time spec> specifies list of onset times for this pulse spec
//         The -ip option may be specified more than once to a max of 63 time.
//         <neuron spec> and <onset time spec> will be parsed by the seconde pass,
//         parse_i_pulse2() after the number of neurons is defined.
//
// input: tb  -- table to store parsed values
//        str -- input sring to be parsed
//        sw  -- switch for this set of options
//
int parse_i_pulse( void *tb, char *args, char *sw )
{
   char  fmtstr[] = "Format is <neurons>,<site>,<magnitide>,<duration>,<onset times>";
   char *token = "";
   char  str[256];
   ipulsespec *ips = (ipulsespec *)tb;

   strcpy( str, args );                    // preserve argument string
   //
   // save neuron list - to be paresed in parse_i_pulse2()
   //
   if ( !( token = strtok( str, "," ) ) ) {
      fprintf( f_log, "%s: Current pulses <%s> neuron or neuron list not specified. %s.\n",
                           pname, sw, fmtstr );
      return( 0 ); 
   }
   else if( strlen( token ) > 64 ) { 
      fprintf( f_log, "%s: Current pulse <%s> neuron list too long. %s.\n",
                           pname, sw, "Neuron list can be at most 63 characters long" );
      return( 0 );      
   }
   else {
      strcpy( ips->spec[ips->count].neurons, token );   // to be parsed in 2nd pass
   }
   //
   // parse injection site
   //
   if ( !( token = strtok( (char *)0, "," ) ) ) {
      fprintf( f_log, "%s: Current pulse <%s> injection site not specified. %s.\n",
                           pname, sw, fmtstr );
      return( 0 ); 
   }
   else if( strcmp( "s", token ) && strcmp( "d", token ) ) {
      fprintf( f_log, "%s: current pulse <%s> injection site <%s> specification error. %s\n",
                           pname, sw, str, "Enter 's' for soma, 'd' for dendrite." );
      return( 0 ); 
   }
   else {
      ips->spec[ips->count].site = token[0];
   }
   //
   // parse injection current pulse magnitude
   //
   if ( !( token = strtok( (char *)0, "," ) ) ) {
      fprintf( f_log, "%s: Current pulse <%s> magnitude not specified. %s.\n",
                           pname, sw, fmtstr );
      return( 0 ); 
   }
   else if( !sscanf( token, "%lf", &(ips->spec[ips->count].magn) ) ) {
      fprintf( f_log, "%s: Current pulse <%s> magnitude <%s> error. %s\n",
                      pname, sw, token, "Enter a floating poiunt value." );
      return( 0 );
   }
   //
   // parse injection current pulse duration
   //
   if ( !( token = strtok( (char *)0, "," ) ) ) {
      fprintf( f_log, "%s: Current pulse <%s> duration not specified. %s.\n",
                           pname, sw, fmtstr );
      return( 0 ); 
   }
   else if( !sscanf( token, "%lf", &(ips->spec[ips->count].dur) ) ) {
      fprintf( f_log, "%s: Current pulse <%s> duration <%s> error. %s\n",
                      pname, sw, token, "Enter a floating poiunt value." );
      return( 0 );
   }
   //
   // save onset time - to be parsed in parse_i_pulse2()
   //
   if ( !( token = strtok( (char *)0, "," ) ) ) {
      fprintf( f_log, "%s: Current pulses <%s> onset times not specified. %s %s.\n",
                      pname, sw, fmtstr, "where <onset times> = <time1>:<time2>. . ." );
      return( 0 ); 
   }
   else if( strlen( token ) > 64 ) { 
      fprintf( f_log, "%s: Current pulses <%s> onset time list too long. %s.\n",
                            pname, sw, "List can be at most 64 characters long" );
      return( 0 );      
   }
   else {
      strcpy( ips->spec[ips->count].times, token );   // to be parsed in 2nd pass
   }
   ips->count++;
   ips->spec[ips->count].site = '\0';    // mark end
   return( 1 );
}

//
// function: parse_i_pulse2();   (second pass)
//
// Parse neuron list and onset times or -ip option.
// Format: <neuron spec> = all | <n1>:<n2> | <n3>. . .
//                         where <n#> specifies 1 neuron or
//                                    specifies a range <nA>-<nB>
// Format: <onset times> = <time1>:<time2>:<time3>. . .
//
// input: tb  -- table to store parsed values
//        str -- input sring to be parsed
//        sw  -- switch for this set of options
//
int parse_i_pulse2( ipulsespec *sp )
{
   int    i = 0, j = 0;
   char   sw[] = "-ip";
   currentpulse  ip;
   int    nstart, nend;         // neuron range
   char  *token = "", *nu = "";
   char   str[128];
   char  *s1, *s2;

   ipulse = (currentpulse *)(vector_s( tot_neurons, sizeof(currentpulse) ) );
   (void )memset( (void *)ipulse, 0, sizeof(currentpulse)*tot_neurons );

   for( i=0; i < sp->count; i++ ) {  // loop through the number of "-ip" specs
      //
      // get a list of onset time
      //
      strcpy( str, sp->spec[i].times );             // preserve argument string
      for( token = strtok(     str,   ":" ), ip.count=0;
           token;
           token = strtok( (char *)0, ":" ), ip.count++ ) {

         if( ip.count >= MAX_IPULSES ) {
	    fprintf( f_log, "%s:  Number of current pulses <%s> execeeds maximum of %d.\n",
                            pname, sw, MAX_IPULSES );
            return( 0 );
         }
         else if( !sscanf( token, "%lf", &(ip.onset[ip.count]) ) ) {
            fprintf( f_log, "%s: Current pulses <%s> onset time <%s> error. %s\n",
                            pname, sw, token, "Enter floating point value in seconds." );
            return( 0 );
         }
      }
      ip.site = sp->spec[i].site;
      ip.magn = sp->spec[i].magn;
      ip.dur  = sp->spec[i].dur;
      //
      // get neuron numbers
      //
//    s1 = malloc( sizeof(char)*128 );  
//    s2 = malloc( sizeof(char)*128 );  
      strcpy( str, sp->spec[i].neurons );           // preserve argument string
      for( token=strtok_r(str,":", &s1); token; token=strtok_r( (char *)0, ":", &s1) ) {
         if( !strcmp( "all", token ) ) { 
            nstart = 1;
            nend   = tot_neurons;
         }
         else if( !( nu = strtok_r( token, "-", &s2 ) ) ) {
            fprintf( f_log, "%s: Current pulses <%s> neuron specification <%s> error. %s\n",
                            pname, sw, token, "Enter 'all', an neuron number or a range." );
            return( 0 );
                    
         }
         else if( !sscanf( nu, "%i", &nstart ) ) {              // range start
            fprintf( f_log, "%s: Current pulses <%s> neuron specification <%s> error. %s\n",
                            pname, sw, token, "Enter 'all', an neuron number or a range." );
            return( 0 );
         }            
         else if( !( nu = strtok_r( (char *)0, "-", &s2 ) ) ) { // one neruon num
            nend = nstart;
         }
         else if( !sscanf( nu, "%i", &nend ) ) {                  // range end
            fprintf( f_log, "%s: Current pulses <%s> neuron range-end <%s> error. %s\n",
                            pname, sw, token, "Enter neuron number as an positive integer." );
            return( 0 );
         }
         if( nstart < 1 || nstart > tot_neurons || nend < 1 || nend > tot_neurons || nstart > nend ) {
            fprintf( f_log, "%s: Current pulses <%s> neuron specifications error <%s> error. %s\n.",
                            pname, sw, token,
                            "Enter neuron number between 1 and number of neurons" );
            return( 0 );           
         }
         nstart--;
         nend--;
         for( j = nstart; j <= nend  ; j++ ) {    // make a copy for each specified neuron
            (void )memcpy( &(ipulse[j]), &ip, sizeof(currentpulse) );
         }
      }
   }
   return( 1 );
}
//
// function: parse_i_field();
//
// Parse argument to parameters for incremntal applied field.
// Format: -if <time>:<magnitude>,<time>:<magnitude>. . . 
//
// input: tb  -- table to store parsed values
//        str -- input sring to be parsed
//        sw  -- switch for this set of options
//
int parse_i_field( void *tb, char *args, char *sw )
{
   int   i = 0;
   char *token = "", *tok = "";
   char  str[128], *s1, *s2;
   ifields *iftb = (ifields *)tb;

   opt_if = 1;
   strcpy( str, args );                    // preserve argument string
   for( i = 0,             token = strtok_r(str,",", &s1 );
        i < MAX_IFIELDS && token;
        i++,               token = strtok_r( (char *)0, ",", &s1 ) ) {
      if( !(tok = strtok_r( token, ":", &s2 ) ) ) {           // onset time
         fprintf( f_log, "%s: Incemental field <%s> onset time missing. %s.\n",
                         pname, sw,
                         "Enter floating point value in seconds" );
         return( 0 );
      }
      else if( !sscanf( tok, "%lf", &(iftb->ifs[i].onset) ) ) {
         fprintf( f_log, "%s: Incremental field <%s> onset time <%s> error. %s.\n",
                      pname, sw, tok, "Enter floating point value in seconds" );
         return( 0 );
      }
      if( !(tok = strtok_r( (char *)0, ":", &s2 ) ) ) {   // field magnitude   (':' not necessary)
         fprintf( f_log, "%s: Incremental field <%s> magnitude missing. %s.\n",
                         pname, sw,
                         "Enter floating point value in mV/cm" );
         return( 0 );
      }
      else if( !sscanf( tok, "%lf", &(iftb->ifs[i].magn) ) ) {
         fprintf( f_log, "%s: Incremental field <%s> magnitude <%s> error. %s.\n",
                      pname, sw, tok, "Enter floating point value in mV/cm" );
         return( 0 );
      }
   }
   if( i > MAX_IFIELDS ) {                // Too many increments
         fprintf( f_log, "%s: Incremental fields <%s> exceeds maximum of %i. %s.\n",
                      pname, sw, MAX_IFIELDS, "Reduce the number of field increments" );
      return( 0 );
   }
   iftb->count = i + 1;
   return( 1 );
}
//
// function: parse_e_field();
//
// Parse argument to parameters for applied field and its onset times.
// Format: -ef <onset time>,<min value>:<increment>:<max value>
//
// input: tb  -- table to store parsed values
//        str -- input sring to be parsed
//        sw  -- switch for this set of options
//
int parse_e_field( void *tb, char *args, char *sw )
{
   char *token = "";
   char  str[128];
   appliedfield *eftb = (appliedfield *)tb;

   opt_ef = 1;
   strcpy( str, args );                    // preserve argument string
   if( !(token=strtok( str, "," ) ) ) {          // onset time
      fprintf( f_log, "%s: Applied field <%s> parameters missing. %s.\n",
                      pname, sw,
                      "Format is <onset time>,<min value>:<increment>:<max value>" );
      return( 0 );
   }
   else if( !sscanf( token, "%lf", &(eftb->onset) ) ) {
      fprintf( f_log, "%s: Applied field <%s> onset time <%s> error. %s\n",
                      pname, sw, token, "Enter floating point value in mV/cm." );
      return( 0 );
   }
   //
   if( !(token=strtok( (char *)0, ":" ) ) ) {    // get min magnitude
      fprintf( f_log, "%s: Applied field <%s> minimum value missing. %s\n",
                      pname, sw, "Enter floating point value in mV/cm." );
      return( 0 );
   }
   else if( !sscanf( token, "%lf", &(eftb->min) ) ) {
      fprintf( f_log, "%s: Applied field <%s> minimum value <%s> error. %s\n",
                      pname, sw, token, "Enter floating point value in mV/cm." );
      return( 0 );
   }
   //
   if( !(token=strtok( (char *)0, ":" ) ) ) {   // get magnitude increment
      eftb->incr = 0;                           // default to one iteration
      eftb->max  = eftb->min;
      return( 1 );
   }
   else if( !sscanf( token, "%lf", &(eftb->incr) ) ) {
      fprintf( f_log, "%s: Applied field <%s> incremental value <%s> error. %s\n",
                      pname, sw, token, "Enter floating point value in mV/cm." );
      return( 0 );
   }
   //
   if( !(token=strtok( (char *)0, ":" ) ) ) {   // get max magnitude
      eftb->max = eftb->min + eftb->incr;       // default to 2 iterations
      return( 1);
   } 
   else if( !sscanf( token, "%lf", &(eftb->max) ) ) {
      fprintf( f_log, "%s: Applied field <%s> maximum value <%s> error. %s\n",
                      pname, sw, token, "Enter floating point value in mV/cm," );
      return( 0 );
   }
   //
   if( eftb->incr < 0 ) {
      fprintf( f_log, "%s: Applied field <%s> error. Increment <%f> must be positive.\n",
                      pname, sw, eftb->incr );
      return( 0 );
   }
   else if( eftb->min > eftb->max ) {
      fprintf( f_log, "%s: Applied field <%s> error. Minimum <%f> must be less than maximum <%f>.\n",
                      pname, sw, (eftb->min), (eftb->max) );
      return( 0 );
   }
   return( 1 );
}
//
// function: parse_iefield();
//
// Parse argument to patameters for applied field and its onset times.
// Format: -efi <magnitude>,<num neurons>
//
// input: tb  -- table to store parsed values
//        str -- input sring to be parsed
//        sw  -- switch for this set of options
//
int parse_iefield( void *tb, char *args, char *sw )
{
   char *token = "";
   char  str[128];
   initfield *efi = (initfield *)tb;

   strcpy( str, args );                    // preserve argument string
   if( !(token=strtok( str, "," ) ) ) {    // magnitude
      fprintf( f_log, "%s: Initiating field magnitude <%s> missing. %s.\n",
                      pname, sw,
                      "Format is <magnitude>,<number of neurons>" );
      return( 0 );
   }
   else if( !sscanf( token, "%lf", &(efi->ef) ) ) {
      fprintf( f_log, "%s: Initiating field <%s> start time <%s> error. %s.\n",
                      pname, sw, token, "Enter floating point value in mV." );
      return( 0 );
   }
   //
   if( !(token=strtok( (char *)0, "," ) ) ) {    // get end time
      fprintf( f_log, "%s: Initiating field number of neurons <%s> missing. %s.\n",
                      pname, sw,
                      "Format is <magnitude>,<number of neurons>" );
      return( 0 );

   }
   else if( !sscanf( token, "%i", &(efi->ineurons) ) ) {
      fprintf( f_log, "%s: Initiating field <%s> number of neurons <%s> error. %s.\n",
                      pname, sw, token, "Enter initeger value" );
      return( 0 );
   }
   return( 1 );
}

//
// function: parse_propagate();
//
// Parse argument to parameters of spike propagation analysis.
// Format: -p <time start>,<neuron start>:<neuron end>    neuron start/end = 1...n
//
// input: tb  -- table to store parsed values
//        str -- input string to be parsed
//        sw  -- switch for this set of options
//
int parse_propagate( void *tb, char *args, char *sw )
{
   char *token = "";
   char  str[128];
   propagatespec *pspec = (propagatespec *)tb;

   strcpy( str, args );                    // preserve argument string
   if( !(token=strtok( str, "," ) ) ) {
      fprintf( f_log, "%s: Propagation specification <%s> values missing. %s.\n",
                      pname, sw,
                      "Format is <time start>,<neuron start>:<neuron end>" );
      return( 0 );
   }
   else if( !sscanf( token, "%lf", &(pspec->tstart) ) ) {
      fprintf( f_log, "%s: Propagation specification <%s> start time <%s> error. %s\n",
                      pname, sw, token, "Enter floating point value in second." );
      return( 0 );
   }
   pspec->spike = 1;
   //
   if( !(token=strtok( (char *)0, ":" ) ) ) {    // get neuron start
      return( 1 );                               // do not compute fit
   }
   else if( !sscanf( token, "%d", &(pspec->nstart) ) ) {
      fprintf( f_log, "%s: Propagation specification <%s> start neuron <%s> error. %s\n",
                      pname, sw, token, "Enter neuron number in integer." );
      return( 0 );
   }
   pspec->speed = 1;                            // compute speed
   //
   if( !(token=strtok( (char *)0, ":" ) ) ) {   // get neuron end
      pspec->nend = num_neurons + efi.ineurons; // default to last neuron  v2.05
   }
   else if( !sscanf( token, "%d", &(pspec->nend) ) ) {
      fprintf( f_log, "%s: Propagation specification <%s> end neuron <%s> error. %s\n",
                      pname, sw, token, "Enter neuron number in integer." );
      return( 0 );
   }
   pspec->nstart--;   // index to neuron array starts with 0
   pspec->nend--;
   return( 1 );
}

//============================================================================
//
// function - write_output();
//
// Write the neurons' time-series output to file.
//
// Go through outputlist, for each output type write result for the range of neurons
// specified by nstart and nend (inclusive). The element "data" is the offset of the
// element to data struct "neurons" containing data for the output type.
//
void write_output( outputlist *f, double t, neuron *nu, int nstart, int nend )
{
   int i = 0, j = 0;

   for( j = 0; f[j].fp; j++ ) {
      if( *(f[j].fp) && f[j].data >= 0 ) {
         fprintf( *(f[j].fp), "%f", t );
         for( i = nstart; i <= nend; i++ ) {
            fprintf( *(f[j].fp), " %f", *((double *)( (int )(&(nu[i])) + f[j].data)) );
         } 
         fprintf( *(f[j].fp), "\n" );
      } 
   }
   return;
}

//============================================================================
//
// function - write_speed();
//
// Write once record of spike propagation velocity as a function of applied fiels.
//
void write_speed( FILE *fp, modelparams *mp )
{
   fprintf( fp, "%f %f %f\n", mp->field, mp->speed, mp->speed_sd );
   return;
}

//============================================================================
//
// function - write_parameters();
//
// Write parameters to file.
//
void write_parameters( FILE *fp, modelparams *m, progrparams *p )
{
   int i=0;
   neuron *nu = m->neurons;

   fprintf( fp, "%s  - %s\n", p->argv[0], asctime( p->datetime ) );
   fprintf( fp, "\nCommand line options:\n" );
   for( i=1; i < p->argc; i++ ) fprintf( fp, "%s ", p->argv[i] );
   fprintf( fp, "\n" );
   
   fprintf( fp, "Pinsky-Rinzel 2-chamber neurons - nearest-neighbor network\n\n" ); 

   fprintf( fp, "applied field         %f %s\n",  m->field,       "mV/cm"   );
   fprintf( fp, "at elapsed time       %f %s\n",  m->af->onset,   "sec"     );
   fprintf( fp, "\n" );
//   fprintf( fp, "applied current pulse %f %s\n",  m->ip->magn,    "uA/cm^2"   );
//   fprintf( fp, "at compartment        %s   \n", (m->ip->site=='s'? "soma":"dendrite") );
//   fprintf( fp, "at elapsed times     " );
//   for( i=0; i<m->ip->count; i++ ) fprintf( fp, " %f", m->ip->onset[i] );
//   fprintf( fp, " sec\n" );
   fprintf( fp, "\n" );
   fprintf( fp, "number of neurons     %i %s\n",  m->nneurons,   ""        );
   fprintf( fp, "resistive lattice     %s   \n",  ( no_r ? "No" : "Yes" )  );
   fprintf( fp, "area                  %f %s\n",  nu[0].Area,    "cm^2"    );
   fprintf( fp, "soma proportion       %f %s\n",  nu[0].P,       ""        );
   fprintf( fp, "c_m                   %f %s\n",  nu[0].cm,      "uF/cm^2" );
   fprintf( fp, "extracell. [K+]       %f %s\n",  m->Kout,       "mM"      );
   fprintf( fp, "intracell. [Ca++]     %f %s\n",  nu[0].v.Ca,    "mM"      );
   fprintf( fp, "gc(0)                 %f %s\n",  nu[0].g.c,     "mS/cm^2" );  
   if( m->nneurons > 1 )
   fprintf( fp, "gc(1)                 %f %s\n",  nu[1].g.c,     "mS/cm^2" );
   fprintf( fp, "gNa                   %f %s\n",  nu[0].g.Na,    "mS/cm^2" );
   fprintf( fp, "gKDR                  %f %s\n",  nu[0].g.KDR,   "mS/cm^2" );
   fprintf( fp, "gCa(0)                %f %s\n",  nu[0].g.Ca,    "mS/cm^2" );
   if( m->nneurons > 1 )
   fprintf( fp, "gCa(1)                %f %s\n",  nu[1].g.Ca,    "mS/cm^2" );
   fprintf( fp, "gKAHP                 %f %s\n",  nu[0].g.KAHP,  "mS/cm^2" );
   fprintf( fp, "gKCa                  %f %s\n",  nu[0].g.KCa,   "mS/cm^2" );
   fprintf( fp, "gNMDA                 %f %s\n",  nu[0].g.NMDA,  "mS/cm^2" );
   fprintf( fp, "gAMPA                 %f %s\n",  nu[0].g.AMPA,  "mS/cm^2" );
   fprintf( fp, "vL                    %f %s\n",  nu[0].v.L,     "mV"      );
   fprintf( fp, "vK (from [K+])        %f %s\n",  nu[0].v.K,     "mV"      );
   fprintf( fp, "vNa                   %f %s\n",  nu[0].v.Na,    "mV"      );
   fprintf( fp, "vsyn                  %f %s\n",  nu[0].v.syn,   "mV"      );
   fprintf( fp, "injected soma cur     %f %s\n",  nu[0].i.soma,  "uA/cm^2" );
   fprintf( fp, "injected dend cur     %f %s\n",  nu[0].i.dend,  "uA/cm^2" );

   // write randomized conditons - later

   //
   // Scan list of files in file table (filetb[]), display those that are enabled.
   //
   for( i=0; p->files[i].name; i++ ) {
      if(  *(p->files[i].open) ) {
         fprintf( fp, "%s  %s\n", p->files[i].name, p->files[i].descr );
      }
   }
   return;
}
//
void memory_allocate( modelparams *mp )
{
   int tn = mp->tneurons;
   int nn = mp->nneurons;
   int nl = mp->nlattice;

   neur     = (neuron *)vector_s( nn,   sizeof(neuron)     );

   em_dR    = (rgrid  *)vector_s( nn+2, sizeof(rgrid )     );
   efield   = (EField *)vector_s( nn,   sizeof(EField)     );

   CInj     = vector( 2*nn );
   CRInj    = vector( 2*nn );    

   isynap   = matrix( tn, nn );

   FMatrix  = matrix( nl, nl );
   Finv     = matrix( nl, nl );
   YRG      = matrix( nn, nn );
   YRGinv   = matrix( nn, nn );
   yid      = matrix( nn, nn );

   dydt     = matrix( tn, N_VARS );
   dy       = matrix( tn, N_VARS );
   dy0      = matrix( tn, N_VARS );  // init cond midpoint

   GMatrix  = matrix( nn, nn );
   RGMatrix = matrix( nn, nn );

   VMatrix  = vector( nn );
   TSUM     = vector( nn );

   return;
}
//
void memory_clear( modelparams *mp )
{
   int nn = mp->nneurons;
   int nl = mp->nlattice;

   (void )memset( (void *)(&gates),   0, sizeof(gatevars) );
   
   (void )memset( (void *)neur,       0, sizeof(neuron)*nn     );
   (void )memset( (void *)em_dR,      0, sizeof(rgrid) *(nn+2) );
   (void )memset( (void *)efield,     0, sizeof(EField)*nn     );

   (void )memset( (void *)CInj,       0, sizeof(double)*(2*nn) );
   (void )memset( (void *)CRInj,      0, sizeof(double)*(2*nn) );

   clear_matrix( isynap,   nn, nn );
   
   clear_matrix( FMatrix,  nl, nl );   
   clear_matrix( Finv,     nl, nl );   
   clear_matrix( YRG,      nn, nn );   
   clear_matrix( YRGinv,   nn, nn );   
   clear_matrix( yid,      nn, nn );

   clear_matrix( dydt,     nn, N_VARS );   
   clear_matrix( dy,       nn, N_VARS );   
   clear_matrix( dy0,      nn, N_VARS );
   
   clear_matrix( GMatrix,  nn, nn);   
   clear_matrix( RGMatrix, nn, nn );   

   return;
}
//
void memory_free()
{
   free( neur     );
   free( em_dR    );
   free_matrix( (void **)isynap   );
   free_matrix( (void **)FMatrix  );
   free_matrix( (void **)Finv     );
   free_matrix( (void **)dydt     );
   free_matrix( (void **)dy       );
   free_matrix( (void **)dy0      );
   free_matrix( (void **)GMatrix  );
   free_matrix( (void **)RGMatrix );
   free_matrix( (void **)YRG      );
   free_matrix( (void **)YRGinv   );
   free_matrix( (void **)yid      );
   free_vector( (void **)TSUM     );
   free( CInj     );
   free( CRInj    );   
   free( VMatrix  );
   return;
}
//
//
//
Vector vector( int n )
{
   Vector v;
   v = (double *)malloc( sizeof(double)*n );
   if( !v ) fprintf( f_log, "vector_d(): allocation failure in double vector.\n");
   return( v );
}

void *vector_s( int n, size_t size )
{
   void *v;
   v = malloc( size*n );
   if( !v ) fprintf( f_log, "vector_d(): allocation failure in structural vector.\n");
   return( v );
}

Matrix matrix( int r, int c )
{
   int    i = 0;
   Matrix m;

   // allocate pointers to rows
   m = (Matrix )malloc( sizeof(double *)*r );
   if( !m )    fprintf( f_log, "matrix_d(): allocation failure in double matrix rows.\n");

   // allocate matrix at first row and set pointers for each row
   m[0]=(double *)malloc( sizeof(double)*r*c );
   if( !m[0] ) fprintf( f_log, "matrix_d(): allocation failure in double matrix columns.\n");

   for( i=1; i< r; i++ ) m[i]=m[i-1]+c;

   return( m );
}

void clear_matrix( Matrix m, int r, int c )
{
   int i = 0;
   for( i=0; i<r; i++ ) memset( m[i], 0, sizeof(double)*c );
}

void free_vector( void *v )
{
   free( v);
}

void free_matrix( void **m )
{
	free( m[0] );
	free( m    );
}
//============================================================================
//
// function - show_vector(), show_matrix();
//
// If trace option is on, write vector or matrix to file
//
void show_vector( FILE *fp, double *v, int n, char *name )
{
   int i=0;
   if (!fp) return;
   fprintf( fp, name );
   fprintf( fp, ":\n" );
   for( i=0; i<n; i++ ) {
      fprintf( fp, "%-4g ", v[i] ); 
   }
   fprintf( fp, "\n" );
   return;
}

void show_matrix( FILE *fp, Matrix m, int r, int c, char *name )
{
   int i=0, j=0;
   if (!fp) return;
   fprintf( fp, name );
   fprintf( fp, ":\n" );
   for( i=0; i<r; i++ ) {
      fprintf( fp, "%-3i: ", i ); 
      for( j=0; j<c; j++ ) {
         fprintf( fp, "%-4g ", m[i][j] ); 
      }
      fprintf( fp, "\n" );
   }
   return;
}

// Origninal comments (1):
//    INCLUDING
//       1) TERMINAL RESISTANCE
//       2) CORRECTION OF KVL
//       3) gCa VALUE IS FIXED (AS 10)
//       4) PHYSICAL DIMENSIONS(G) INCLUDEd(LUdmp and LUBKSB)
//       5) INVERSE MATRIX (TO SOLVE LINEAR EQ: A x = b; CALCULATE CURRENTS)
//       6) RE-IMPOSING SYMMETRY TO SOME ELEMENTS THAT ARE INVOLVED IN TWO CURRENTS
//          AS OUTPUT CURRENTS AMONG TOTAL 18 CURRENTS (COUPLED TWO NEURONS CASE ONLY)
//          INVERSE MATRIX : Finv(i,j)
//             reimposing: Finv(2,1)=Finv(7,6)
//             reimposing: Finv(2,6)=Finv(7,1)
//             reimposing: Finv(2,11)=Finv(7,14)
//             reimposing: Finv(2,14)=Finv(7,11)
//
//             reimposing: Finv(2,4)=Finv(7,9)
//             reimposing: Finv(2,5)=Finv(7,10)
//             reimposing: Finv(2,9)=Finv(7,4)
//             reimposing: Finv(2,10)=Finv(7,5)
//
//             reimposing: Finv(2,4)=-Finv(2,5)
//             reimposing: Finv(2,9)=-Finv(2,10)
//
// DATED: 2002, DEC. 18

// Original comments (2):
// Note that FORTRAN stores two-dimensional matrice by column,
// so Finv[1,j] the address of the jth column of Finv
// Inverse of Matrix is done !!
//----------------------------------
//        Finv[2, 4] =Finv[7, 9]=Finv[12, 4] =a5;
//        Finv[2, 5] =Finv[7,10]=Finv[12, 5] =a6;
//        Finv[2, 9] =Finv[7, 4]=Finv[12, 9] =a7;
//        Finv[2,10] =Finv[7, 5]=Finv[12,10] =a8;
//                reimposing: Finv(2,1)=Finv(7,6)
//                reimposing: Finv(2,6)=Finv(7,1)
//                reimposing: Finv(2,11)=Finv(7,14)
//                reimposing: Finv(2,14)=Finv(7,11)
//----------------------------------
//   double a5= -0.10679859827759630164;
//   double a6=  0.10679859827759630164;
//   double a7= -0.06792345226484586274;
//   double a8=  0.06792345226484586274;
//------Area multiplication---------
