#include <cmath>
#include "funs.h"
#include <map>
#ifdef MPIBART
#include "mpi.h"
#endif

//--------------------------------------------------
// normal density N(x, mean, variance)
double pn(double x, double m, double v)
{
   double dif = x-m;
   return exp(-.5*dif*dif/v)/sqrt(2*PI*v);
}
double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
 
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
 
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    return 0.5*(1.0 + sign*y);
}
//--------------------------------------------------
// draw from discrete distributin given by p, return index
int rdisc(double *p, RNG& gen)
{

   double sum;
   double u = gen.uniform();

    int i=0;
    sum=p[0];
    while(sum<u) {
       i += 1;
       sum += p[i];
    }
    return i;
}
//--------------------------------------------------
//evalute tree tr on grid given by xi and write to os
void grm(tree& tr, xinfo& xi, std::ostream& os) 
{
   size_t p = xi.size();
   if(p!=2) {
      cout << "error in grm, p !=2\n";
      return;
   }
   size_t n1 = xi[0].size();
   size_t n2 = xi[1].size();
   tree::tree_cp bp; //pointer to bottom node
   double *x = new double[2];
   for(size_t i=0;i!=n1;i++) {
      for(size_t j=0;j!=n2;j++) {
         x[0] = xi[0][i]; 
         x[1] = xi[1][j]; 
         bp = tr.bn(x,xi);
         os << x[0] << " " << x[1] << " " << bp->getm() << " " << bp->nid() << endl;
      }
   }
   delete[] x;
}
//--------------------------------------------------
//does this bottom node n have any variables it can split on.
//CHANGE!!!
bool cansplit(tree::tree_p n, xinfo& xi,int* vartype)
{
   int L,U;
   bool v_found = false; //have you found a variable you can split on
   size_t v=0;
   while(!v_found && (v < xi.size())) { //invar: splitvar not found, vars left
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U,vartype);
      if(U>=L) v_found=true;
      v++;
   }
   return v_found;
}
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
//CHANGE!!!
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots,int* vartype)
{
   double pb;  //prob of birth to be returned
   tree::npv bnv; //all the bottom nodes
   t.getbots(bnv);
   for(size_t i=0;i!=bnv.size();i++) 
     if(cansplit(bnv[i],xi,vartype)) goodbots.push_back(bnv[i]);
   if(goodbots.size()==0) { //are there any bottom nodes you can split on?
      pb=0.0;
   } else { 
      if(t.treesize()==1) pb=1.0; //is there just one node?
      else pb=pi.pb;
   }
   return pb;
}
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
//CHANGE!!!
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars,int* vartype)
{
   int L,U;
   for(size_t v=0;v!=xi.size();v++) {//try each variable
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U,vartype);
      if(U>=L) goodvars.push_back(v);
   }
}
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
//CHANGE!!!
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi,int* vartype)
{
  if(cansplit(n,xi,vartype)) {
      return pi.alpha/pow(1.0+n->depth(),pi.beta);
   } else {
      return 0.0;
   }
}
//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv)
{
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x
   double y;          //current y

   bnv.clear();
   x.getbots(bnv);

   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();
   sv.resize(nb);

   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      y=di.y[i];

      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];

      ++(sv[ni].n);
      sv[ni].sy += y;
      sv[ni].sy2 += y*y;
   }
}
#ifdef MPIBART
void MPImasterallsuff(tree& x, tree::npv& bnv, std::vector<sinfo>& sv, size_t numslaves)
{
	bnv.clear();
	x.getbots(bnv);

   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();
   sv.resize(nb);

	int bufsz=((int)nb)*sizeof(int)+2*nb*sizeof(double);
	unsigned int *n = new unsigned int[nb];
	double *sy = new double[nb];
	double *sy2 = new double[nb];
	char *buffer = new char[bufsz];
	int position;
	MPI_Status status;

	for(size_t i=1;i<=numslaves;i++)
	{
		position=0;
		MPI_Recv(buffer,bufsz,MPI_PACKED,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
		MPI_Unpack(buffer,bufsz,&position,n,(int)nb,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,bufsz,&position,sy,(int)nb,MPI_DOUBLE,MPI_COMM_WORLD);
		MPI_Unpack(buffer,bufsz,&position,sy2,(int)nb,MPI_DOUBLE,MPI_COMM_WORLD);
		for(size_t j=0;j<nb;j++)
		{
			sv[j].n += (size_t)n[j];
			sv[j].sy += sy[j];
			sv[j].sy2 += sy2[j];
		}		
	}

	delete[] buffer;
	delete[] n;
	delete[] sy;
	delete[] sy2;
}
void MPIslaveallsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv)
{
	tree::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x
   double y;          //current y

	bnv.clear();
	x.getbots(bnv);

   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();

	int bufsz=((int)nb)*sizeof(int)+2*nb*sizeof(double);
	unsigned int *n = new unsigned int[nb];
	double *sy = new double[nb];
	double *sy2 = new double[nb];
	char *buffer = new char[bufsz];
	int position=0;

	std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++){
		bnmap[bnv[i]]=i;
		n[i]=0;
		sy[i]=0.0;
		sy2[i]=0.0;
	}

	for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      y=di.y[i];

      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];

		++n[ni];
		sy[ni]+=y;
		sy2[ni]+=y*y;
   }
	MPI_Pack(n,(int)nb,MPI_UNSIGNED,buffer,bufsz,&position,MPI_COMM_WORLD);
	MPI_Pack(sy,(int)nb,MPI_DOUBLE,buffer,bufsz,&position,MPI_COMM_WORLD);
	MPI_Pack(sy2,(int)nb,MPI_DOUBLE,buffer,bufsz,&position,MPI_COMM_WORLD);
	MPI_Send(buffer,bufsz,MPI_PACKED,0,0,MPI_COMM_WORLD);

	delete[] buffer;
	delete[] sy2;
	delete[] sy;
	delete[] n;
}
#endif
//--------------------------------------------------
//get sufficient stats for children (v,c) of node nx in tree x
void getsuff(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
   double *xx;//current x
   double y;  //current y
   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
         y = di.y[i];
         if(xx[v] < xi[v][c]) {
               sl.n++;
               sl.sy += y;
               sl.sy2 += y*y;
          } else {
               sr.n++;
               sr.sy += y;
               sr.sy2 += y*y;
          }
      }
   }
}
//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
void getsuff(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
   double *xx;//current x
   double y;  //current y
   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);
      if(bn==nl) {
         y = di.y[i];
         sl.n++;
         sl.sy += y;
         sl.sy2 += y*y;
      }
      if(bn==nr) {
         y = di.y[i];
         sr.n++;
         sr.sy += y;
         sr.sy2 += y*y;
      }
   }
}
#ifdef MPIBART
//MPI version of get sufficient stats - this is the master code
void MPImastergetsuff(tree::tree_cp nl, tree::tree_cp nr, sinfo &sl, sinfo &sr, size_t numslaves)
{
	sl.n=0;sl.sy=0.0;sl.sy2=0.0;
	sr.n=0;sr.sy=0.0;sr.sy2=0.0;
	char buffer[48];
	int position=0;
	sinfo slavel,slaver;
	MPI_Status status;
	MPI_Request *request=new MPI_Request[numslaves];
	const int tag=0; //tag=0 means it's not a v,c type get sufficient stats.
	unsigned int nlid,nrid,ln,rn;

	nlid=(unsigned int)nl->nid();
	nrid=(unsigned int)nr->nid();

	// Pack and send info to the slaves
	MPI_Pack(&nlid,1,MPI_UNSIGNED,buffer,48,&position,MPI_COMM_WORLD);
	MPI_Pack(&nrid,1,MPI_UNSIGNED,buffer,48,&position,MPI_COMM_WORLD);
	for(size_t i=1; i<=numslaves; i++) {   
		MPI_Isend(buffer,48,MPI_PACKED,i,tag,MPI_COMM_WORLD,&request[i-1]);
	}
	MPI_Waitall(numslaves,request,MPI_STATUSES_IGNORE);

	// MPI receive all the answers from the slaves
	for(size_t i=1; i<=numslaves; i++) {
		position=0;
		MPI_Recv(buffer,48,MPI_PACKED,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
		MPI_Unpack(buffer,48,&position,&ln,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&slavel.sy,1,MPI_DOUBLE,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&slavel.sy2,1,MPI_DOUBLE,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&rn,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&slaver.sy,1,MPI_DOUBLE,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&slaver.sy2,1,MPI_DOUBLE,MPI_COMM_WORLD);
		slavel.n=(size_t)ln;
		slaver.n=(size_t)rn;
		sl.n+=slavel.n;
		sl.sy+=slavel.sy;
		sl.sy2+=slavel.sy2;
		sr.n+=slaver.n;
		sr.sy+=slaver.sy;
		sr.sy2+=slaver.sy2;
	}
	
	delete[] request;
}
//--------------------------------------------------
//MPI version of get sufficient stats for children (v,c) of node nx in tree x - this is the master code
void MPImastergetsuffvc(tree::tree_cp nx, size_t v, size_t c, xinfo& xi, sinfo& sl, sinfo& sr, size_t numslaves)
{
   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;
	char buffer[48];
	int position=0;
	sinfo slavel,slaver;
	MPI_Status status;
	MPI_Request *request = new MPI_Request[numslaves];
	const int tag=1; //tag=1 means it is a v,c type get sufficient stats
	unsigned int vv,cc,nxid,ln,rn;

	vv=(unsigned int)v;
	cc=(unsigned int)c;
	nxid=(unsigned int)nx->nid();

	// Pack and send info to the slaves
	MPI_Pack(&nxid,1,MPI_UNSIGNED,buffer,48,&position,MPI_COMM_WORLD);
	MPI_Pack(&vv,1,MPI_UNSIGNED,buffer,48,&position,MPI_COMM_WORLD);
	MPI_Pack(&cc,1,MPI_UNSIGNED,buffer,48,&position,MPI_COMM_WORLD);
	for(size_t i=1; i<=numslaves; i++) {
		MPI_Isend(buffer,48,MPI_PACKED,i,tag,MPI_COMM_WORLD,&request[i-1]);
	}
	MPI_Waitall(numslaves,request,MPI_STATUSES_IGNORE);

	// MPI receive all the answers from the slaves
	for(size_t i=1; i<=numslaves; i++) {
		position=0;
		MPI_Recv(buffer,48,MPI_PACKED,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
		MPI_Unpack(buffer,48,&position,&ln,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&slavel.sy,1,MPI_DOUBLE,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&slavel.sy2,1,MPI_DOUBLE,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&rn,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&slaver.sy,1,MPI_DOUBLE,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&slaver.sy2,1,MPI_DOUBLE,MPI_COMM_WORLD);
		slavel.n=(size_t)ln;
		slaver.n=(size_t)rn;
		sl.n+=slavel.n;
		sl.sy+=slavel.sy;
		sl.sy2+=slavel.sy2;
		sr.n+=slaver.n;
		sr.sy+=slaver.sy;
		sr.sy2+=slaver.sy2;
	}

	delete[] request;
}
//--------------------------------------------------
//MPI version of add birth to a tree on the compute nodes
void MPImastersendbirth(tree::tree_p nx, size_t v, size_t c, double mul, double mur, size_t numslaves)
{
	char buffer[40];
	int position=0;
	const int tag=1; //tag=1 for birth
	MPI_Request *request = new MPI_Request[numslaves];
	unsigned int nxid, vv, cc;

	vv=(unsigned int)v;
	cc=(unsigned int)c;
	nxid=(unsigned int)nx->nid();

	MPI_Pack(&nxid,1,MPI_UNSIGNED,buffer,40,&position,MPI_COMM_WORLD);
	MPI_Pack(&vv,1,MPI_UNSIGNED,buffer,40,&position,MPI_COMM_WORLD);
	MPI_Pack(&cc,1,MPI_UNSIGNED,buffer,40,&position,MPI_COMM_WORLD);
	MPI_Pack(&mul,1,MPI_DOUBLE,buffer,40,&position,MPI_COMM_WORLD);
	MPI_Pack(&mur,1,MPI_DOUBLE,buffer,40,&position,MPI_COMM_WORLD);
	for(size_t i=1; i<=numslaves; i++) {
		MPI_Isend(buffer,40,MPI_PACKED,i,tag,MPI_COMM_WORLD,&request[i-1]);
	}
	MPI_Waitall(numslaves,request,MPI_STATUSES_IGNORE);
	delete[] request;
}
//--------------------------------------------------
//MPI version of add death to a tree on the compute nodes
void MPImastersenddeath(tree::tree_p nx, double mu, size_t numslaves)
{
	char buffer[40];
	int position=0;
	const int tag=0; //tag=0 for death
	MPI_Request *request = new MPI_Request[numslaves];
	unsigned int nxid;

	nxid=(unsigned int)nx->nid();
	MPI_Pack(&nxid,1,MPI_UNSIGNED,buffer,40,&position,MPI_COMM_WORLD);	
	MPI_Pack(&mu,1,MPI_DOUBLE,buffer,40,&position,MPI_COMM_WORLD);
	for(size_t i=1; i<=numslaves; i++) {
		MPI_Isend(buffer,40,MPI_PACKED,i,tag,MPI_COMM_WORLD,&request[i-1]);
	}
	MPI_Waitall(numslaves,request,MPI_STATUSES_IGNORE);
	delete[] request;
}
//--------------------------------------------------
//MPI master code to send no birth/death to slave, ie the mcmc rejected the birth/death step
void MPImastersendnobirthdeath(size_t numslaves)
{
	const int tag=2; //tag=2 when there is no birth and no death at current step of the mcmc
	MPI_Request *request = new MPI_Request[numslaves];

	for(size_t i=1; i<=numslaves; i++) {
		MPI_Isend(0,0,MPI_PACKED,i,tag,MPI_COMM_WORLD,&request[i-1]);
	}
	MPI_Waitall(numslaves,request,MPI_STATUSES_IGNORE);
	delete[] request;
}
//--------------------------------------------------
//MPI version of add birth or death to a tree on the compute nodes - slave code.
void MPIslaveupdatebirthdeath(tree& x)
{
	double mul, mur, mu;
	tree::npv nv;
	MPI_Status status;
	char buffer[40];
	int position=0;
	unsigned int v,c,nxid;

	MPI_Recv(buffer,40,MPI_PACKED,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	if(status.MPI_TAG==0) //death
	{
		MPI_Unpack(buffer,40,&position,&nxid,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,40,&position,&mu,1,MPI_DOUBLE,MPI_COMM_WORLD);
		x.death((size_t)nxid,mu);
	}
	else if(status.MPI_TAG==1) //birth
	{
		MPI_Unpack(buffer,40,&position,&nxid,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,40,&position,&v,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,40,&position,&c,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,40,&position,&mul,1,MPI_DOUBLE,MPI_COMM_WORLD);
		MPI_Unpack(buffer,40,&position,&mur,1,MPI_DOUBLE,MPI_COMM_WORLD);
		x.birth((size_t)nxid,(size_t)v,(size_t)c,mul,mur);
	}
	//else, no birth death so do nothing.
}
//--------------------------------------------------
//MPI version of get sufficient stats - this is the slave code
void MPIslavegetsuff(tree& x, xinfo& xi, dinfo& di)
{
	sinfo sl, sr;  // what we will send back to the master.
	unsigned int nxid,nlid,nrid,v,c,ln,rn;
	tree::tree_cp nl,nr;
	tree::tree_p nx;
	tree::npv bnv, tnv;
	char buffer[48];
	int position=0;
	MPI_Status status;

	//	MPI receive the nlid and nrid.
	MPI_Recv(buffer,48,MPI_PACKED,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	if(status.MPI_TAG==0)
	{
		MPI_Unpack(buffer,48,&position,&nlid,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&nrid,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		position=0;

		nl=x.getptr((size_t)nlid);
		nr=x.getptr((size_t)nrid);
		getsuff(x,nl,nr,xi,di,sl,sr);
	
		// Pack the result and MPI send it.
		ln=(unsigned int)sl.n;
		rn=(unsigned int)sr.n;
		MPI_Pack(&ln,1,MPI_UNSIGNED,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Pack(&sl.sy,1,MPI_DOUBLE,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Pack(&sl.sy2,1,MPI_DOUBLE,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Pack(&rn,1,MPI_UNSIGNED,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Pack(&sr.sy,1,MPI_DOUBLE,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Pack(&sr.sy2,1,MPI_DOUBLE,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Send(buffer,48,MPI_PACKED,0,0,MPI_COMM_WORLD);
	}
	else //tag==1 => get suff stats using children (v,c) of node nx in tree x.
	{
		MPI_Unpack(buffer,48,&position,&nxid,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&v,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		MPI_Unpack(buffer,48,&position,&c,1,MPI_UNSIGNED,MPI_COMM_WORLD);
		position=0;

		nx=x.getptr((size_t)nxid);
		getsuff(x,nx,(size_t)v,(size_t)c,xi,di,sl,sr);

		ln=(unsigned int)sl.n;
		rn=(unsigned int)sr.n;
		MPI_Pack(&ln,1,MPI_UNSIGNED,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Pack(&sl.sy,1,MPI_DOUBLE,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Pack(&sl.sy2,1,MPI_DOUBLE,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Pack(&rn,1,MPI_UNSIGNED,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Pack(&sr.sy,1,MPI_DOUBLE,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Pack(&sr.sy2,1,MPI_DOUBLE,buffer,48,&position,MPI_COMM_WORLD);
		MPI_Send(buffer,48,MPI_PACKED,0,0,MPI_COMM_WORLD);
	}
}
#endif
//--------------------------------------------------
//log of the integrated likelihood
double lil(size_t n, double sy, double sy2, double sigma, double tau)
{
   double yb,yb2,S,sig2,d;
   double sum, rv;

   yb = sy/n;
   yb2 = yb*yb;
   S = sy2 - (n*yb2);
   sig2 = sigma*sigma;
   d = n*tau*tau + sig2;
   sum = S/sig2 + (n*yb2)/d;
   rv = -(n*LTPI/2.0) - (n-1)*log(sigma) -log(d)/2.0;
   rv = rv -sum/2.0;
   return rv;
}
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, std::vector<double>& fv)
{
   double *xx;
   tree::tree_cp bn;
   fv.resize(di.n);
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);
      fv[i] = bn->getm();
   }
}
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, double* fv)
{
   double *xx;
   tree::tree_cp bn;
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);
      fv[i] = bn->getm();
   }
}
//--------------------------------------------------
//partition
void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv)
{
   double *xx;
   tree::tree_cp bn;
   pv.resize(di.n);
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);
      pv[i] = bn->nid();
   }
}
//--------------------------------------------------
// draw all the bottom node mu's

void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen)
{
   tree::npv bnv;
   std::vector<sinfo> sv;
   allsuff(t,xi,di,bnv,sv);

   double a = 1.0/(pi.tau * pi.tau);
   double sig2 = pi.sigma * pi.sigma;
   double b,ybar;

   for(tree::npv::size_type i=0;i!=bnv.size();i++) {
      b = sv[i].n/sig2;
      if(sv[i].n==0){ybar=0.0;}
      else{
	ybar = sv[i].sy/sv[i].n;}
      bnv[i]->setm(b*ybar/(a+b) + gen.normal()/sqrt(a+b));
   }
}

#ifdef MPIBART
//-----------------------------------------------------
// Draw all the bottom node mu's -- slave code for the MPI version
void MPIslavedrmu(tree& t, xinfo& xi, dinfo& di)
{
	tree::npv bnv;
	int position=0;
	char *buffer;
	double temp;
	MPI_Status status;

	//slave code to support the master computing all the sufficient statistics	
	MPIslaveallsuff(t,xi,di,bnv);

	//now sync the new draws made on the master back to all the slaves
	int bufsz=((int)bnv.size())*sizeof(double);
	buffer = new char[bufsz];
	MPI_Recv(buffer,bufsz,MPI_PACKED,0,0,MPI_COMM_WORLD,&status);
	for(tree::npv::size_type i=0;i!=bnv.size();i++) {
		MPI_Unpack(buffer,bufsz,&position,&temp,1,MPI_DOUBLE,MPI_COMM_WORLD);
		bnv[i]->setm(temp);
	}
	delete[] buffer;
}

//---------------------------------------------------
// Draw all the bottom node mu's -- master code for MPI version
void MPImasterdrmu(tree& t, xinfo& xi, pinfo& pi, RNG& gen, size_t numslaves)
{
   tree::npv bnv;
   std::vector<sinfo> sv;
	MPI_Request *request=new MPI_Request[numslaves];
	char *buffer;
	double temp;
	int position=0;

	MPImasterallsuff(t,bnv,sv,numslaves);
	int bufsz=((int)bnv.size())*sizeof(double);
	buffer = new char[bufsz];

   double a = 1.0/(pi.tau * pi.tau);
   double sig2 = pi.sigma * pi.sigma;
   double b,ybar;

   for(tree::npv::size_type i=0;i!=bnv.size();i++) {
      b = sv[i].n/sig2;
      ybar = sv[i].sy/sv[i].n;
      bnv[i]->setm(b*ybar/(a+b) + gen.normal()/sqrt(a+b));
		temp=bnv[i]->getm();
		MPI_Pack(&temp,1,MPI_DOUBLE,buffer,bufsz,&position,MPI_COMM_WORLD);
   }

	for(size_t i=1;i<=numslaves;i++)
		MPI_Isend(buffer,bufsz,MPI_PACKED,i,0,MPI_COMM_WORLD,&request[i-1]);
	MPI_Waitall(numslaves,request,MPI_STATUSES_IGNORE);

	delete[] request;
	delete[] buffer;
}
#endif
//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi)
{
   cout << "xinfo: \n";
   for(size_t v=0;v!=xi.size();v++) {
      cout << "v: " << v << endl;
      for(size_t j=0;j!=xi[v].size();j++) cout << "j,xi[v][j]: " << j << ", " << xi[v][j] << endl;
   }
   cout << "\n\n";
}
//--------------------------------------------------
//make xinfo = cutpoints
// NOTE!!!!! add categorical data
// catnum is the number of categories-1
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc,int *vartype)
{
   double xinc;

   //compute min and max for each x
   std::vector<double> minx(p,INFINITY);
   std::vector<double> maxx(p,-INFINITY);
   double xx;
   for(size_t i=0;i<p;i++) {
     if(vartype[i]==0){
      for(size_t j=0;j<n;j++) {
         xx = *(x+p*j+i);
         if(xx < minx[i]) minx[i]=xx;
         if(xx > maxx[i]) maxx[i]=xx;
      }}
   }
   //make grid of nc cutpoints between min and max for each x.
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
     if(vartype[i]==0){
      xinc = (maxx[i]-minx[i])/(nc+1.0);
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
     } else {xi[i].resize(1.0);
             xi[i][0]=0.5;}
   }
}
// get min/max needed to make cutpoints
void makeminmax(size_t p, size_t n, double *x, std::vector<double> &minx, std::vector<double> &maxx)
{
   double xx;

   for(size_t i=0;i<p;i++) {
      for(size_t j=0;j<n;j++) {
         xx = *(x+p*j+i);
         if(xx < minx[i]) minx[i]=xx;
         if(xx > maxx[i]) maxx[i]=xx;
      }
   }
}
//make xinfo = cutpoints give the minx and maxx vectors
void makexinfominmax(size_t p, xinfo& xi, size_t nc, std::vector<double> &minx, std::vector<double> &maxx,int* vartype)
{
	double xinc;
   //make grid of nc cutpoints between min and max for each x.
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
     if(vartype[i]==0){
      xinc = (maxx[i]-minx[i])/(nc+1.0);
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
     } else {xi[i].resize(1.0);
             xi[i][0]=0.5;}
   }
}

#ifdef MPIBART
// construct prediction of f given a single tree, overwrites existing values in ppredmean if there are any.
void makepred(dinfo dip, xinfo &xi, std::vector<tree> &t, double *ppredmean)
{
	double* fpredtemp=0; //temporary fit vector to compute prediction
	size_t m;

	fpredtemp = new double[dip.n];
	m=t.size();
    for(size_t i=0;i<dip.n;i++) ppredmean[i]=0.0;
	for(size_t j=0;j<m;j++) {
		fit(t[j],xi,dip,fpredtemp);
		for(size_t k=0;k<dip.n;k++) ppredmean[k] += fpredtemp[k];
	}
	delete[] fpredtemp;
}

// construct prediction of y given a single tree, overwrites existing values in ppredmean if there are any.
void makeypred(dinfo dip, xinfo &xi, std::vector<tree> &t, double sigma, double *ppredmean)
{
	uint seed;
	seed=(unsigned int)time(NULL);
	RNG rnrm(seed);
	double* fpredtemp=0; //temporary fit vector to compute prediction
	size_t m;

	fpredtemp = new double[dip.n];
	m=t.size();
    for(size_t i=0;i<dip.n;i++) ppredmean[i]=rnrm.normal(0.0,sigma);
	for(size_t j=0;j<m;j++) {
		fit(t[j],xi,dip,fpredtemp);
		for(size_t k=0;k<dip.n;k++) ppredmean[k] += fpredtemp[k];
	}
	delete[] fpredtemp;

}

// construct E[f] over sample of draws from the posterior, overwrites exisiting values in postp
// if there are any.
void makepostpred(dinfo dip, xinfo &xi, std::vector< std::vector<tree> > &t, double *postp)
{
	double* fpredtemp=0; //temporary fit vector to compute prediction
	double* ppredmean=0; //temporary fit vector for mean from 1 tree
	size_t m,ndraws;

	fpredtemp = new double[dip.n];
	ppredmean = new double[dip.n];
	ndraws=t.size();
	m=t[0].size();
	for(size_t i=0;i<dip.n;i++) { 
		ppredmean[i]=0.0;
		postp[i]=0.0;
	}

	for(size_t i=0;i<ndraws;i++) {
		for(size_t j=0;j<m;j++) {
			fit(t[i][j],xi,dip,fpredtemp);
			for(size_t k=0;k<dip.n;k++) ppredmean[k] += fpredtemp[k];
		}
		if(i>0)
			for(size_t k=0;k<dip.n;k++)
			{
				postp[k] *= (int)i;
				postp[k] += ppredmean[k];
				postp[k] /= (int)(i+1);
				ppredmean[k] = 0.0;
			}
		else
			for(size_t k=0;k<dip.n;k++)
			{
				postp[k] = ppredmean[k];
				ppredmean[k] = 0.0;
			}
	}
}

// construct E[f] and E[f^2] over sample of draws from the posterior, overwrites exisiting values in postp,
// postp2 if there are any.
void makepostpred2(dinfo dip, xinfo &xi, std::vector< std::vector<tree> > &t, double *postp, double *postp2)
{
	double* fpredtemp=0; //temporary fit vector to compute prediction
	double* ppredmean=0; //temporary fit vector for mean from 1 tree
	size_t m,ndraws;

	fpredtemp = new double[dip.n];
	ppredmean = new double[dip.n];
	ndraws=t.size();
	m=t[0].size();
	for(size_t i=0;i<dip.n;i++) { 
		ppredmean[i]=0.0;
		postp[i]=0.0;
		postp2[i]=0.0;
	}

	for(size_t i=0;i<ndraws;i++) {
		for(size_t j=0;j<m;j++) {
			fit(t[i][j],xi,dip,fpredtemp);
			for(size_t k=0;k<dip.n;k++) ppredmean[k] += fpredtemp[k];
		}
		if(i>0)
			for(size_t k=0;k<dip.n;k++)
			{
				postp[k] *= (int)i;
				postp[k] += ppredmean[k];
				postp[k] /= (int)(i+1);
				postp2[k] *= (int)i;
				postp2[k] += ppredmean[k]*ppredmean[k];
				postp2[k] /= (int)(i+1);
				ppredmean[k] = 0.0;
			}
		else
			for(size_t k=0;k<dip.n;k++)
			{
				postp[k] = ppredmean[k];
				postp2[k] = ppredmean[k]*ppredmean[k];
				ppredmean[k] = 0.0;
			}
	}
}

double logcalp(std::vector<double> &theta, dinfo dip, xinfo &xi, std::vector<tree> &t, double sigmae, double sigma, size_t pth, size_t myrank)
{
	double logptemp=0.0, logp=0.0;
	double *ppredmean;
	
	if(myrank==0)
	{
		for(size_t j=0;j<(dip.p-pth);j++)
			if(theta[j]<0.0 || theta[j]>1.0) 
				logptemp+=-INFINITY;
	}
	else //myrank>0
	{
		ppredmean=new double[dip.n];
		makeypred(dip,xi,t,sigma,ppredmean);
		for(size_t i=0;i<dip.n;i++)
			logptemp+=(dip.y[i]-ppredmean[i])*(dip.y[i]-ppredmean[i]);
		logptemp/=sigmae*sigmae;
		logptemp/=-2;
		delete[] ppredmean;
	}
	MPI_Allreduce(&logptemp,&logp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	return logp;
}
#endif
