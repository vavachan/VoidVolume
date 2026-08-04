/*
 * ============================================================================
 *  voidVolume3D -- exact void/free-volume analysis for sphere packings
 * ============================================================================
 *
 *  Computes the exact geometric volume of the void (unoccupied) space in a
 *  3-D packing of possibly-overlapping, possibly-polydisperse spheres, by
 *  constructing the (weighted / radical) Delaunay tessellation of the sphere
 *  centers and integrating the empty space within each Delaunay tetrahedron
 *  that is not covered by any sphere.
 *
 *  Algorithm reference:
 *    Sastry, S., Corti, D. S., Debenedetti, P. G., & Stillinger, F. H. (1997).
 *    "Statistical geometry of particle packings. I. Algorithm for exact
 *    determination of connectivity, volume, and surface areas of void space
 *    in monodisperse and polydisperse sphere packings."
 *    Physical Review E, 56(5), 5524.
 *
 *  Build:
 *    g++ voidVolume3D.cpp -std=c++11 -O3 -o voidVolume3D
 *    (or simply: make)
 *
 *  Usage:
 *    ./voidVolume3D <input_file> <probe_radius> [output_dir] [sample_atom_index]
 *
 *  See README.md for a full description of the input file format, all
 *  command-line arguments, every output file that is produced, and known
 *  limitations (in particular the periodic-boundary-condition caveat about
 *  particles tessellating with their own periodic image).
 * ============================================================================
 */

#include<iostream>
#include<math.h>
#include<fstream>
#include<ostream>
#include<iomanip>
#include<vector>
#include <set>
#include <limits>
#include <functional>
#include <numeric>
#include <algorithm>
#include <deque>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

// ---------------------------------------------------------------------------
// Global simulation state.
//
// The original implementation threads state (box geometry, periodicity
// flags, running totals, and work queues for the breadth-first Delaunay
// construction) through global variables rather than a context object or
// class members. That architecture has been preserved here to avoid
// touching the numerically delicate geometry routines below; the globals
// are documented in one place instead of being scattered as implicit
// parameters through dozens of function signatures.
// ---------------------------------------------------------------------------
long double boxx, boxy, boxz;        // Simulation box lengths (x, y, z).
int PBCx = 1;                        // Periodic boundary conditions are
int PBCy = 1;                        // applied along all three axes; see
int PBCz = 1;                        // README.md for the PBC caveat.
long double r_cut = 0.0;             // Probe radius, added to every particle
				     // radius before the geometric analysis.

long double convex_vol = 0.;         // Running sum of Delaunay tetrahedron
				     // volumes (diagnostic: should equal
				     // boxx*boxy*boxz for a valid PBC
				     // tessellation).
int nAtoms = 0;                      // Number of particles read from input.
std::vector<int> incompleteAtoms;    // Particles whose local Delaunay
				     // neighborhood is not yet fully closed.
std::vector<int> allAtoms;           // Particles not yet visited by the
				     // breadth-first traversal in main().
std::deque<int> atomsToAnalyze;              // BFS work queue.
std::unordered_set<int> processedOrQueuedAtoms;  // BFS "already seen" set.

// Opened in main() once the output directory has been resolved from the
// command line; receives a VMD Tcl script drawing every edge of every
// Delaunay tetrahedron (useful for visually auditing the tessellation).
std::ofstream outputFile;

struct face;
struct vertice;
struct vect
{
	long double x=0;
	long double y=0;
	long double z=0;
};
//#include <unordered_map>

struct pair_hash {
	std::size_t operator()(const std::pair<int, int>& p) const {
		return std::hash<int>()(p.first) ^ std::hash<int>()(p.second);
	}
};
struct faceHash{
	std::size_t operator()(const std::array<int,3> p) const {
		return (std::hash<int>()(p[0]) ^ std::hash<int>()(p[1]))^ std::hash<int>()(p[2]);
	}
};
long double distancesq(vect p,vect q)
{
	long double dist=0.;
	dist=(p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y)+(p.z-q.z)*(p.z-q.z);
	return dist;
}
//struct vect
//{
//	 long double x=0;
//	 long double y=0;
//	 long double z=0;
//};

vect vectDiff(vect *a, vect *b)
{
	vect dr;
	dr.x=a->x-b->x;
	dr.y=a->y-b->y;
	dr.z=a->z-b->z;

	dr.x=(dr.x-(boxx*PBCx*lroundl(dr.x/boxx)));
	dr.y=(dr.y-(boxy*PBCy*lroundl(dr.y/boxy)));
	dr.z=(dr.z-(boxz*PBCz*lroundl(dr.z/boxz)));
	return dr;
}

long double innerProduct(vect* a, vect* b)
{
	return a->x*b->x+a->y*b->y+a->z*b->z;
}
vect cross_product(vect* a, vect* b)
{
	vect result;
	result.x = a->y * b->z - a->z * b->y;
	result.y = a->z * b->x - a->x * b->z;
	result.z = a->x * b->y - a->y * b->x;
	return result;
}

long double magnitudeSq(vect *v)
{
	return v->x*v->x+v->y*v->y+v->z*v->z;
}
struct face
{
	vect A1,A2,A3;
	vect B;
	face *next=nullptr;
};

struct delunay
{
	int AT[4]= {0};
	vect MID[4][4][4];
	vect MIDP[4][4];
	int EDGE[4][4];
	int FACE[4][4][4];
	std::unordered_map<std::array<int,3>, delunay*, faceHash> neighDel;
	std::unordered_map<std::array<int,3>, bool, faceHash> isVoid;
	int numNeighbors=0;
	vertice *v=nullptr;
	vect circumCenter;
	//long double circum_x=0.;
	//long double circum_y=0.;
	//long double circum_z=0.;
	delunay *next=nullptr;
	delunay *prev=nullptr;
	int hull=0;
	int clusterIndex=-1;
	//int solid=0;

	bool checkIfNeighbor(delunay *D)
	{
		for(auto nD:neighDel)
		{
			if(nD.second==D)
				return true;
		}
		return false;
	}
	delunay* findNeighbor(std::array<int,3> face)
	{
		std::sort(face.begin(),face.end());
		if(neighDel.find(face)!=neighDel.end())
		{
			return neighDel[face];
		}
		else
		{
			return nullptr;
			//std::cerr<<"Face ["<<face[0]<<"\t"<<face[1]<<"\t"<<face[2]<<"] does not exist in tetrahedron ["<<AT[0]<<"\t"<<AT[1]<<"\t"<<AT[2]<<"\t"<<AT[3]<<"]\n";
			//exit(EXIT_FAILURE);
		}
	}
	void addNeighbor(delunay* D, std::array<int,3> face, bool bondInVoid)
	{
		std::sort(face.begin(),face.end());
		if(!findNeighbor(face))
		{
			neighDel[face]=D;
			isVoid[face]=bondInVoid;
			numNeighbors++;
		}
		else
		{
			if(findNeighbor(face) != D)
				std::cerr<<"error\n";
		}
	}
	void delNeighbor(delunay* D)
	{
		std::vector<std::pair<std::array<int,3>, delunay*>> entries(neighDel.begin(), neighDel.end());
		for(auto nD:entries)
		{
			if(nD.second==D)
				neighDel.erase(nD.first);
		}
	}
};
struct DelanuayTessellaton
{
	std::vector<delunay*> tetrahedrons;
	void addTetrahedron(delunay* D) 
	{
		tetrahedrons.push_back(D);
	}
} fullSet;
struct vertice
{
	struct vect *p=nullptr;
	struct half_edge *leaving=nullptr;
	struct vertice *next=nullptr;
	struct vertice *prev=nullptr;
	int dangling=0;
	int vivectd=0;
	delunay *D=nullptr;
	int is_void=0;
	//int cluster_index=-1;
	long double r;
	long double ball_r=0.;
	//vert_list *V;
	vertice *neib_vert[10]= {nullptr};
	int neib_ed[10];
	int v_neigh_count=0;
};

class convex_hull
{
	public : 
		face *initial=nullptr;
		face *end=nullptr;
		void insert_face(face *f);
}CH,solid_wall;

void convex_hull :: insert_face(face *f1)
{
	if(!initial)
	{
		initial=f1;
		end=f1;
	}
	else
	{
		end->next=f1;
		end=end->next;
	}
}
struct atom
{
	vect p;
	//std::vector<int> neighlist;
	std::vector<int> contiguous;
	//std::unordered_map<int,int> conti_index;
	std::unordered_map<std::pair<int, int>, int, pair_hash> part_c;
	std::unordered_map<int,int> edge_index;
	std::unordered_map<std::pair<int, int>, vect ,pair_hash> MIDP;
	std::unordered_map<int,vect> RMID;
	std::unordered_map<std::pair<int, int>, std::pair<delunay*,delunay*>,pair_hash> DelaunayTessellations;
	int bor=0;
	std::vector<delunay*> delunayTetrahedrons;
	int numFaces=0;
	struct face *F=nullptr;
	long double radius=1.;
	int ignore=0;
	int index;
	long double vor_vol;

	bool check_if_contiguous(int I)
	{
		bool exists = (std::find(contiguous.begin(), contiguous.end(), I) != contiguous.end());
		return exists;
	}
	int get_part_c(int I, int J)
	{
		if(not(check_if_contiguous(I)))// and check_if_contiguous(J)))
		{
			std::cerr << "Error: "<<I<<" not part\n";
			exit(EXIT_FAILURE);
		}
		if(not(check_if_contiguous(J)))// and check_if_contiguous(J)))
		{
			std::cerr << "Error: "<<J<<" not part\n";
			exit(EXIT_FAILURE);
		}
		if (I < J) 
		{
			std::cerr << "Error: part_c[" << I << "," << J << "] wrong insertion \n";
			exit(EXIT_FAILURE);
		}
		std::pair<int, int> key = std::make_pair(I, J);
		return part_c[key];
	}
	void increment_part_c(int I, int J) 
	{
		if (I < J) 
		{
			std::cerr << "Error: part_c[" << I << "," << J << "] wrong insertion \n";
			exit(EXIT_FAILURE);
		}

		std::pair<int, int> key = std::make_pair(I, J);

		// Check if it exists and increment
		int &counter = part_c[key];
		counter++;

		if (counter > 2) {
			std::cerr << "Error: part_c[" << I << "," << J << "] exceeded 2\nNumber of contiguous\t"<<contiguous.size()<<"\n";;
			exit(EXIT_FAILURE);
		}
	}
	bool checkIfIncompleteFace()
	{
		bool yes=false;
		for (const auto& element : part_c) 
		{
			if(element.second==1)
			{
				yes=true;
				break;
			}
		}
		return yes;
	}
};
// (unused debug helpers 'print_face' and 'print_delunay_solid' removed;
// they printed VMD draw commands to stdout for a single face/tetrahedron
// and were never called anywhere in the pipeline.)
void print_delunay(std::ostream& outFile,delunay *D,atom Atoms[],int TYPE=0)
{
	outFile<<"##\t";
	for(int i=0; i<4; i++)
	{
		outFile<<D->AT[i]<<"\t";
	}
	//for(int i=0;i<4;i++)
	//	for(int j=0;j<4;j++)
	//		for(int k=0;k<4;k++)
	//			for(uu
	outFile<<"\n";
	//return ;
	atom *ATOM;
	ATOM=&(Atoms[D->AT[0]]);
	long double Sx,Sy,Sz;
	long double Px,Py,Pz;
	Sx=ATOM->p.x-Atoms[D->AT[1]].p.x;
	Sy=ATOM->p.y-Atoms[D->AT[1]].p.y;
	Sz=ATOM->p.z-Atoms[D->AT[1]].p.z;
	Sx=(Sx-(boxx*PBCx*lroundl(Sx/boxx)));
	Sy=(Sy-(boxy*PBCy*lroundl(Sy/boxy)));
	Sz=(Sz-(boxz*PBCz*lroundl(Sz/boxz)));
	outFile<<"draw line\t{";
	outFile<<ATOM->p.x<<"\t"<<ATOM->p.y<<"\t"<<ATOM->p.z<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	4\n";
	Sx=ATOM->p.x-Atoms[D->AT[2]].p.x;
	Sy=ATOM->p.y-Atoms[D->AT[2]].p.y;
	Sz=ATOM->p.z-Atoms[D->AT[2]].p.z;
	Sx=Sx-boxx*PBCx*lroundl(Sx/boxx);
	Sy=Sy-boxy*PBCy*lroundl(Sy/boxy);
	Sz=Sz-boxz*PBCz*lroundl(Sz/boxz);

	outFile<<"draw line\t{";
	outFile<<ATOM->p.x<<"\t"<<ATOM->p.y<<"\t"<<ATOM->p.z<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	4\n";
	Sx=ATOM->p.x-Atoms[D->AT[3]].p.x;
	Sy=ATOM->p.y-Atoms[D->AT[3]].p.y;
	Sz=ATOM->p.z-Atoms[D->AT[3]].p.z;
	Sx=Sx-boxx*PBCx*lroundl(Sx/boxx);
	Sy=Sy-boxy*PBCy*lroundl(Sy/boxy);
	Sz=Sz-boxz*PBCz*lroundl(Sz/boxz);

	outFile<<"draw line\t{";
	outFile<<ATOM->p.x<<"\t"<<ATOM->p.y<<"\t"<<ATOM->p.z<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	4\n";
	Sx=ATOM->p.x-Atoms[D->AT[2]].p.x;
	Sy=ATOM->p.y-Atoms[D->AT[2]].p.y;
	Sz=ATOM->p.z-Atoms[D->AT[2]].p.z;
	Sx=Sx-boxx*PBCx*lroundl(Sx/boxx);
	Sy=Sy-boxy*PBCy*lroundl(Sy/boxy);
	Sz=Sz-boxz*PBCz*lroundl(Sz/boxz);

	Px=ATOM->p.x-Atoms[D->AT[1]].p.x;
	Py=ATOM->p.y-Atoms[D->AT[1]].p.y;
	Pz=ATOM->p.z-Atoms[D->AT[1]].p.z;
	Px=Px-boxx*PBCx*lroundl(Px/boxx);
	Py=Py-boxy*PBCy*lroundl(Py/boxy);
	Pz=Pz-boxz*PBCz*lroundl(Pz/boxz);

	outFile<<"draw line\t{";
	outFile<<ATOM->p.x-Px<<"\t"<<ATOM->p.y-Py<<"\t"<<ATOM->p.z-Pz<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	4\n";
	Sx=ATOM->p.x-Atoms[D->AT[3]].p.x;
	Sy=ATOM->p.y-Atoms[D->AT[3]].p.y;
	Sz=ATOM->p.z-Atoms[D->AT[3]].p.z;
	Sx=Sx-boxx*PBCx*lroundl(Sx/boxx);
	Sy=Sy-boxy*PBCy*lroundl(Sy/boxy);
	Sz=Sz-boxz*PBCz*lroundl(Sz/boxz);

	Px=ATOM->p.x-Atoms[D->AT[1]].p.x;
	Py=ATOM->p.y-Atoms[D->AT[1]].p.y;
	Pz=ATOM->p.z-Atoms[D->AT[1]].p.z;
	Px=Px-boxx*PBCx*lroundl(Px/boxx);
	Py=Py-boxy*PBCy*lroundl(Py/boxy);
	Pz=Pz-boxz*PBCz*lroundl(Pz/boxz);

	outFile<<"draw line\t{";
	outFile<<ATOM->p.x-Px<<"\t"<<ATOM->p.y-Py<<"\t"<<ATOM->p.z-Pz<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	4\n";
	Sx=ATOM->p.x-Atoms[D->AT[2]].p.x;
	Sy=ATOM->p.y-Atoms[D->AT[2]].p.y;
	Sz=ATOM->p.z-Atoms[D->AT[2]].p.z;
	Sx=Sx-boxx*PBCx*lroundl(Sx/boxx);
	Sy=Sy-boxy*PBCy*lroundl(Sy/boxy);
	Sz=Sz-boxz*PBCz*lroundl(Sz/boxz);

	Px=ATOM->p.x-Atoms[D->AT[3]].p.x;
	Py=ATOM->p.y-Atoms[D->AT[3]].p.y;
	Pz=ATOM->p.z-Atoms[D->AT[3]].p.z;
	Px=Px-boxx*PBCx*lroundl(Px/boxx);
	Py=Py-boxy*PBCy*lroundl(Py/boxy);
	Pz=Pz-boxz*PBCz*lroundl(Pz/boxz);

	outFile<<"draw line\t{";
	outFile<<ATOM->p.x-Px<<"\t"<<ATOM->p.y-Py<<"\t"<<ATOM->p.z-Pz<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	4\n";
	//outFile<<"draw sphere\t{";
	//outFile<<Atoms[D->AT[0]].p.x<<"\t"<<Atoms[D->AT[0]].p.y<<"\t"<<Atoms[D->AT[0]].p.z<<"}\tradius\t"<<Atoms[D->AT[0]].radius<<"\n";
	//outFile<<"draw sphere\t{";
	//outFile<<Atoms[D->AT[1]].p.x<<"\t"<<Atoms[D->AT[1]].p.y<<"\t"<<Atoms[D->AT[1]].p.z<<"}\tradius\t"<<Atoms[D->AT[1]].radius<<"\n";
	//outFile<<"draw sphere\t{";
	//outFile<<Atoms[D->AT[2]].p.x<<"\t"<<Atoms[D->AT[2]].p.y<<"\t"<<Atoms[D->AT[2]].p.z<<"}\tradius\t"<<Atoms[D->AT[2]].radius<<"\n";
	//outFile<<"draw sphere\t{";
	//outFile<<Atoms[D->AT[3]].p.x<<"\t"<<Atoms[D->AT[3]].p.y<<"\t"<<Atoms[D->AT[3]].p.z<<"}\tradius\t"<<Atoms[D->AT[3]].radius<<"\n";
}
vect cross_product(vect a1,vect a2)
{
	vect cp;
	cp.x=a2.z*a1.y-a1.z*a2.y;
	cp.y=a1.z*a2.x-a1.x*a2.z;
	cp.z=a2.y*a1.x-a1.y*a2.x;
	return cp;
}
/** 3x3 matrix determinant (used to solve the small linear systems that
 *  locate radical-plane intersection points and circumcenters). */
long double determinant( long double a[3][3])
{   
	return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
}
/** Solves the 3x3 linear system a*x = b by Cramer's rule. */
vect cramer( long double a[3][3], long double b[3],int debug=0)
{   
	long double x[3];
	long double det_a=determinant(a);
	long double det_a1=0.;
	long double a_1[3][3];
	if(debug)
	{
		cout<<"[";
		for(int i=0;i<3;i++)
		{
			cout<<"[";
			for(int j=0;j<3;j++)
			{
				cout<<a[i][j]<<",";
				//a_1[i][j]=a[i][j];
			}
			cout<<"]";
			//cout<<"\n";
		}
		cout<<"]\n";
		cout<<"[";
		for(int i=0;i<3;i++)
		{
			cout<<b[i]<<",";
		}
		cout<<"]\n";
		cout<<det_a<<"\n";
	}

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
				a_1[j][k]=a[j][k];
			}
		}       
		for(int j=0;j<3;j++)
		{
			a_1[j][i]=b[j]/det_a;
		}
		//for(int i=0;i<3;i++)
		//{
		//    for(int j=0;j<3;j++)
		//    {
		//        cout<<a_1[i][j]<<"\t";
		//    }
		//    cout<<"\n";
		//}
		//cout<<"\n\n";
		det_a1=determinant(a_1);
		//cout<<det_a1<<"\n";
		x[i]=det_a1;
		////for(int j=0;j<3;j++)
		////{
		////    a_1[j][i]=a[j][i];
		////}       
	}
	vect p;
	p.x=x[0];
	p.y=x[1];
	p.z=x[2];
	return p;
}
/** Radical-plane center of three (possibly differently sized) spheres
 *  centered at A1, A2, A3 with radii rS, rA, rB: the point equidistant, in
 *  the power-distance sense, from all three -- the polydisperse analogue of
 *  a circumcenter used to build the weighted Delaunay tessellation. */
vect center_of_triangle(vect A1,vect A2,vect A3, long double rS, long double rA, long double rB,int debug=0)
{
	vect v21,v31;

	vect center;
	v21=vectDiff(&A2,&A1);
	v31=vectDiff(&A3,&A1);
	// long double rA,rS,rB;
	long double DISA;
	long double DISB;
	DISA=sqrtl(magnitudeSq(&v21));
	DISB=sqrtl(magnitudeSq(&v31));
	vect v21crossV31=cross_product(&v21,&v31);
	//a=ZB*YA-ZA*YB;
	//b=ZA*XB-XA*ZB;
	//c=YB*XA-YA*XB;
	long double B[3],A[3][3];

	B[0]=(DISA*DISA+rS*rS-rA*rA)/2.;
	B[1]=(DISB*DISB+rS*rS-rB*rB)/2.;
	B[2]=0.;

	A[0][0]=v21.x;
	A[0][1]=v21.y;
	A[0][2]=v21.z;

	A[1][0]=v31.x;
	A[1][1]=v31.y;
	A[1][2]=v31.z;

	A[2][0]=v21crossV31.x;
	A[2][1]=v21crossV31.y;
	A[2][2]=v21crossV31.z;

	center=cramer(A,B,debug);
	center.x=center.x+A1.x;
	center.y=center.y+A1.y;
	center.z=center.z+A1.z;
	return center;
}
/** Signed-then-absolute volume of one Delaunay tetrahedron, applying the
 *  minimum-image convention so tetrahedra spanning a periodic boundary are
 *  measured correctly. */
long double volume_delunay(delunay *D,atom Atoms[])
{
	long double a[3][3];
	a[0][0]=Atoms[D->AT[1]].p.x-Atoms[D->AT[0]].p.x;
	a[0][1]=Atoms[D->AT[1]].p.y-Atoms[D->AT[0]].p.y;
	a[0][2]=Atoms[D->AT[1]].p.z-Atoms[D->AT[0]].p.z;

	a[0][0]=a[0][0]-boxx*PBCx*lroundl(a[0][0]/boxx);
	a[0][1]=a[0][1]-boxy*PBCy*lroundl(a[0][1]/boxy);
	a[0][2]=a[0][2]-boxz*PBCz*lroundl(a[0][2]/boxz);


	a[1][0]=Atoms[D->AT[2]].p.x-Atoms[D->AT[0]].p.x;
	a[1][1]=Atoms[D->AT[2]].p.y-Atoms[D->AT[0]].p.y;
	a[1][2]=Atoms[D->AT[2]].p.z-Atoms[D->AT[0]].p.z;

	a[1][0]=a[1][0]-boxx*PBCx*lroundl(a[1][0]/boxx);
	a[1][1]=a[1][1]-boxy*PBCy*lroundl(a[1][1]/boxy);
	a[1][2]=a[1][2]-boxz*PBCz*lroundl(a[1][2]/boxz);


	a[2][0]=Atoms[D->AT[3]].p.x-Atoms[D->AT[0]].p.x;
	a[2][1]=Atoms[D->AT[3]].p.y-Atoms[D->AT[0]].p.y;
	a[2][2]=Atoms[D->AT[3]].p.z-Atoms[D->AT[0]].p.z;

	a[2][0]=a[2][0]-boxx*PBCx*lroundl(a[2][0]/boxx);
	a[2][1]=a[2][1]-boxy*PBCy*lroundl(a[2][1]/boxy);
	a[2][2]=a[2][2]-boxz*PBCz*lroundl(a[2][2]/boxz);


	long double det=determinant(a);
	return abs(det/6.);
}
/** Records particles i and j as Delaunay-adjacent ("contiguous"), unless
 *  the bond is already recorded from either side. */
void connectAtoms(int i, int j, atom Atoms [])
{
	//std::cout<<i<<"\t"<<j<<"\n";
	//if(Atoms[i].conti_index.find(j)==Atoms[i].conti_index.end() and Atoms[j].conti_index.find(i)==Atoms[j].conti_index.end())
	//bool existsIJ = (std::find(Atoms[i].contiguous.begin(), Atoms[i].contiguous.end(), j) != Atoms[i].contiguous.end());
	//bool existsJI = (std::find(Atoms[j].contiguous.begin(), Atoms[j].contiguous.end(), i) != Atoms[j].contiguous.end());
	if(!Atoms[i].check_if_contiguous(j) and !Atoms[j].check_if_contiguous(i))
	{
		atom *Ai = &Atoms[i];
		atom *Aj = &Atoms[j];

		//contiguous[i][j]=1;
		//contiguous[j][i]=1;

		//Ai->conti_index[j]=Ai->conti;
		//Aj->conti_index[i]=Aj->conti;

		//Ai->conti_index[j]=Ai->contiguous.size();
		//Aj->conti_index[i]=Aj->contiguous.size();

		Ai->contiguous.push_back(j);
		Aj->contiguous.push_back(i);

		//Ai->conti++;
		//Aj->conti++;
	}
}
/** Links two Delaunay tetrahedra (D_ONE, D_TWO) that share the triangular
 *  face {ATOM, I, J} as Voronoi-dual neighbors, and classifies the shared
 *  Voronoi facet/edge as lying inside a sphere or in the void: it is a void
 *  edge when the segment between the two circumcenters exits the covering
 *  sphere of every incident particle. */
void connectDelunayTessellations(delunay* D_ONE, delunay* D_TWO, atom* ATOM, atom Atoms[], int I, int J)
{
	if(not(ATOM->check_if_contiguous(I) and ATOM->check_if_contiguous(J)))
	{
		std::cerr<<"wrong delunay in connect atoms\n";
		exit(EXIT_FAILURE);
	}

	vertice *temp_vert_o=nullptr;
	vertice *temp_vert_d=nullptr;

	temp_vert_o=D_ONE->v;
	temp_vert_d=D_TWO->v;

	vect vI = vectDiff(&Atoms[I].p,&ATOM->p);
	vect vJ = vectDiff(&Atoms[J].p,&ATOM->p);

	vect vert1 = vectDiff(&D_ONE->circumCenter,&ATOM->p);
	vect vert2 = vectDiff(&D_TWO->circumCenter,&ATOM->p);
	long double rS=ATOM->radius;

	vect vJcrossI = cross_product(&vJ,&vI);//a=(Y1*Z2-Z1*Y2);

	long double overlap1,overlap2;
	int sign1,sign2;
	long double disV1,disV2,disA;

	disV1=sqrtl(magnitudeSq(&vert1));
	disV2=sqrtl(magnitudeSq(&vert2));
	if((disV1*disV1-rS*rS)>0.)
	{
		temp_vert_o->is_void=1;
		temp_vert_o->ball_r=disV1-rS;
	}
	//disV2=sqrtl((V2x*V2x+V2y*V2y+V2z*V2z));
	if((disV2*disV2-rS*rS)>0.)
	{
		temp_vert_d->is_void=1;
		temp_vert_d->ball_r=disV2-rS;//sqrtl(disV2*disV2-rS*rS);
	}

	disA=sqrtl(magnitudeSq(&vJcrossI));
	overlap1=innerProduct(&vJcrossI,&vert1);
	overlap2=innerProduct(&vJcrossI,&vert2);
	if(overlap1<0.)
		sign1=1;
	else
		sign1=-1;
	//overlap2=a*V2x+b*V2y+c*V2z;
	if(overlap2<0.)
		sign2=1;
	else
		sign2=-1;
	bool isBondInVoid=false;
	if(sign1==sign2)
	{

		if(overlap1/(disV1*disA)<overlap2/(disV2*disA))
		{
			if(temp_vert_o->is_void==1)
			{
				//ATOM->D3bondinvoid[{I,J}]=1;
				//ATOM->D3bondinvoid[{J,I}]=1;
				isBondInVoid = true;
			}
		}
		else
		{
			if(temp_vert_d->is_void==1)
			{
				//ATOM->D3bondinvoid[{I,J}]=1;
				//ATOM->D3bondinvoid[{J,I}]=1;
				isBondInVoid = true;
			}
		}
	}
	else
	{
		//long double midx,midy,midz;
		//int p,q,r;
		//cout<<ATOM->MIDP[i][j].x<<"\t"<<ATOM->MIDP[i][j].y<<"\t"<<ATOM->MIDP[i][j].z<<"\n";
		vect mid = vectDiff(&ATOM->MIDP[{I,J}],&ATOM->p);
		//midx=ATOM->MIDP[{i,j}].x-ATOM->p.x;
		//midy=ATOM->MIDP[{i,j}].y-ATOM->p.y;
		//midz=ATOM->MIDP[{i,j}].z-ATOM->p.z;
		//midx=midx-boxx*PBCx*lroundl(midx/boxx);
		//midy=midy-boxy*PBCy*lroundl(midy/boxy);
		//midz=midz-boxz*PBCz*lroundl(midz/boxz);


		long double dismsq=magnitudeSq(&mid);
		if((dismsq-rS*rS)>0.)
		{
			//ATOM->D3bondinvoid[{I,J}]=1;
			//ATOM->D3bondinvoid[{J,I}]=1;
			isBondInVoid = true;
		}
	}
	D_ONE->addNeighbor(D_TWO,std::array<int,3>{ATOM->index,I,J},isBondInVoid);
	D_TWO->addNeighbor(D_ONE,std::array<int,3>{ATOM->index,I,J},isBondInVoid);
}
/** Finalizes a newly constructed Delaunay tetrahedron D: computes its
 *  circumcenter/Voronoi vertex, accumulates its volume into convex_vol, and
 *  registers it against the four particles it is built from. */
void create_delunay(atom Atoms[],delunay *D)
{
	//cout<<"ghere\n";
	convex_vol=convex_vol+volume_delunay(D,Atoms);
	vertice *temp_v=nullptr;
	temp_v=new vertice;
	temp_v->p=new vect;
	temp_v->p = &D->circumCenter;
	//temp_v->p->x=D->circum_x;
	//temp_v->p->y=D->circum_y;
	//temp_v->p->z=D->circum_z;
	temp_v->D=D;
	D->v=temp_v;
	for(int a=0; a<3; a++)
	{
		for(int b=a+1; b<4; b++)
		{
			vect a1,a2,a3,center;
			long double r1,r2,r3;
			a1=Atoms[D->AT[a]].p;
			r1=Atoms[D->AT[a]].radius;
			a2=Atoms[D->AT[b]].p;
			r2=Atoms[D->AT[b]].radius;
			//D->MIDP[a][b]=
			long double DISA,l;
			vect a21 = vectDiff(&a2,&a1);

			DISA=sqrtl(magnitudeSq(&a21));

			l=0.5*(DISA+(r1*r1-r2*r2)/DISA);

			D->MIDP[a][b].x=l/DISA*a21.x+a1.x;
			D->MIDP[a][b].y=l/DISA*a21.y+a1.y;
			D->MIDP[a][b].z=l/DISA*a21.z+a1.z;

			D->MIDP[b][a].x=l/DISA*a21.x+a1.x;
			D->MIDP[b][a].y=l/DISA*a21.y+a1.y;
			D->MIDP[b][a].z=l/DISA*a21.z+a1.z;
			////{
			////	cout<<"B\t"<<a<<"\t"<<b<<"\n";;
			////	display_SITE(&D->MIDP[a][b]);
			////}
			//cout<<a<<"\t"<<b<<"\n";
			for(int c=b+1; c<4; c++)
			{
				a3=Atoms[D->AT[c]].p;

				r3=Atoms[D->AT[c]].radius;

				//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";
				////display_SITE(&a1);
				////display_SITE(&a2);
				////display_SITE(&a3);

				center=center_of_triangle(a1,a2,a3,r1,r2,r3);

				D->MID[a][b][c]=center;
				D->MID[a][c][b]=center;

				D->MID[b][a][c]=center;
				D->MID[b][c][a]=center;

				D->MID[c][a][b]=center;
				D->MID[c][b][a]=center;

			}	
		}
	}
	for(int i=0; i<3; i++)
	{
		for(int j=i+1; j<4;j++)
		{
			connectAtoms(D->AT[i],D->AT[j],Atoms);
		}
	}
	for(int k=0; k<4; k++)
	{
		atom *ATOM;
		atom *Ai;
		atom *Aj;
		atom *Ak;
		ATOM=&(Atoms[D->AT[k]]);
		Ai=&(Atoms[D->AT[(k+1)%4]]);
		Aj=&(Atoms[D->AT[(k+2)%4]]);
		Ak=&(Atoms[D->AT[(k+3)%4]]);

		int a1=D->AT[(k+1)%4];
		int a2=D->AT[(k+2)%4];
		int a3=D->AT[(k+3)%4];

		//if(ATOM->part_c[a1][a2]==0)
		//if(ATOM->part_c[a1][a3]==0)
		//if(ATOM->part_c[a2][a3]==0)

		for (const std::pair<int, int>& pair : { std::make_pair(a1, a2), std::make_pair(a1, a3), std::make_pair(a2, a3) }) 
		{
			int I,J;
			if(pair.first > pair.second)
			{
				I=pair.first;
				J=pair.second;
			}
			else
			{
				I=pair.second;
				J=pair.first;
			}
			//if(ATOM->part_c.find({I,J}) == ATOM->part_c.end())
			//{
			//	ATOM->edge_index[I]++;
			//	ATOM->edge_index[J]++;
			//}
			//ATOM->part_c[{I,J}]++;
			if(ATOM->part_c.find({I,J})!=ATOM->part_c.end())
			{
				if(ATOM->get_part_c(I,J)==2)
				{
					std::cerr << "\nError: part_c[" << I << "," << J << "] exceeded 2\n";
					std::cerr << ATOM->index<<"\t"<<Ai->index<<"\t"<<Aj->index<<"\t"<<Ak->index<<"\n";
					std::cerr << a1 <<"\t"<<a2<<"\t"<<a3<<"\n";
					std::cerr << ATOM->contiguous.size()<<"\n";
					exit(EXIT_FAILURE);
				}
			}
			ATOM->increment_part_c(I,J);
			ATOM->numFaces++;
			ATOM->MIDP[{I,J}]=D->MID[k][(k+1)%4][(k+2)%4];
			if(!ATOM->checkIfIncompleteFace())
			{
				/* This atom is completed */
				incompleteAtoms.erase(std::remove(incompleteAtoms.begin(), incompleteAtoms.end(), ATOM->index), incompleteAtoms.end());
			}
			if(ATOM->DelaunayTessellations.find({I,J}) == ATOM->DelaunayTessellations.end())
				ATOM->DelaunayTessellations[{I,J}] = std::make_pair(D,nullptr);
			else
			{
				ATOM->DelaunayTessellations[{I,J}].second=D;
				connectDelunayTessellations(ATOM->DelaunayTessellations[{I,J}].first,D,ATOM,Atoms,I,J);
			}
		}
		// Use pair.first and pair.second

		//ATOM->RMID[a1]=D->MIDP[k][(k+1)%4];
		//ATOM->RMID[a2]=D->MIDP[k][(k+2)%4];
		//ATOM->RMID[a3]=D->MIDP[k][(k+3)%4];
		ATOM->delunayTetrahedrons.push_back(D);
	}
}
/** Builds the very first Delaunay tetrahedron incident to ATOM (its three
 *  nearest, non-collinear/non-coplanar neighbors plus ATOM itself), seeding
 *  the local tessellation that completeDelunayTessellation() then expands
 *  outward from by rotating around each edge. */
void first_delunay(atom *ATOM,atom Atoms[])
{
	long double min=1e10;
	int A1,A2,A3,A4;
	A1=ATOM->index;
	atom a1=Atoms[A1];
	//cout<<ATOM->index<<"\t"<<ATOM->neighbours<<"\n";
	//cout<<ATOM->radius<<"\n";
	for(int i=0; i<nAtoms; i++)
	{
		long double rA,rS,DIS,l;
		vect neighImAtom = vectDiff(&Atoms[i].p,&ATOM->p);

		DIS=sqrtl(magnitudeSq(&neighImAtom));
		rA=Atoms[i].radius;
		rS=ATOM->radius;
		l=0.5*(DIS+(rS*rS-rA*rA)/DIS); // This is the distance between the atom center and radical plane (a the line connecting the centers of the atoms)

		long double squaredTangent =l*l-rS*rS; // This is the squared tangented distance 
		if(squaredTangent<min)
		{
			//std::cout<<squaredTangent<<"\t"<<i<<"\n";
			min=squaredTangent;
			A2=i;
		}
	}
	atom a2=Atoms[A2];
	min = 1e10;
	for(int i=0; i<nAtoms; i++)
	{
		if(i!=A2)
		{
			atom a3=Atoms[i];
			vect center;
			center=center_of_triangle(a1.p,a2.p,a3.p,a1.radius,a2.radius,a3.radius);
			long double squaredTangent = distancesq(center,a1.p)-a1.radius*a1.radius;	
			if(squaredTangent<min)
			{
				min=squaredTangent;
				A3=i;
			}
		}
	}
	atom a3=Atoms[A3];
	min=1e10;
	long double circx_s;//=circx;
	long double circy_s;//=circy;
	long double circz_s;//=circz;
			    //	cout<<DIS_atom<<"\n";
			    //atom a1=Atoms[A1];
			    //atom a2=Atoms[A2];
			    //atom a3=Atoms[A3];
	vect center;
	for(int i=0; i<nAtoms; i++)
	{
		if(i!=A2 && i!=A3)
		{
			atom a4=Atoms[i];

			vect a21 = vectDiff(&a2.p,&a1.p);
			vect a31 = vectDiff(&a3.p,&a1.p);
			vect a41 = vectDiff(&a4.p,&a1.p);

			long double B[3],A[3][3];

			long double r1=a1.radius;
			long double r2=a2.radius;
			long double r3=a3.radius;
			long double r4=a4.radius;

			long double DIS2=sqrtl(magnitudeSq(&a21));
			long double DIS3=sqrtl(magnitudeSq(&a31));
			long double DIS4=sqrtl(magnitudeSq(&a41));

			B[0]=(DIS2*DIS2+r1*r1-r2*r2)/2.;
			B[1]=(DIS3*DIS3+r1*r1-r3*r3)/2.;
			B[2]=(DIS4*DIS4+r1*r1-r4*r4)/2.;

			A[0][0]=a21.x;
			A[0][1]=a21.y;
			A[0][2]=a21.z;
			A[1][0]=a31.x;
			A[1][1]=a31.y;
			A[1][2]=a31.z;
			A[2][0]=a41.x;
			A[2][1]=a41.y;
			A[2][2]=a41.z;

			center=cramer(A,B);

			long double squaredTangent =center.x*center.x+center.y*center.y+center.z*center.z-r1*r1;

			if(squaredTangent < min)
			{
				min = squaredTangent;
				A4  = i;
				circx_s=center.x+a1.p.x;
				circy_s=center.y+a1.p.y;
				circz_s=center.z+a1.p.z;

			}
		}
	}

	std::vector<int> EV1 {A1,A2,A3,A4};
	//std::sort(EV1.begin(),EV1.end());
	delunay *D;

	D=new delunay;
	fullSet.addTetrahedron(D);


	D->AT[0]=EV1[0];
	D->AT[1]=EV1[1];
	D->AT[2]=EV1[2];
	D->AT[3]=EV1[3];

	if(D->AT[0]!=ATOM->index)
	{
		//THis shouldn't happen because if this happens it means in first delunay you made a delunay which involves atoms which were analyzed p
		//previously. If this delunay was acceptable then this atoms delunay should have been found before and we wouldn't be here
		//I have to do this to get this thing working. A good progrmmer would find some other way
		// All that is left in front of me this ad hoc trick.
		//cout<<"error\n";
		return;
	}

	D->circumCenter.x=circx_s;
	D->circumCenter.y=circy_s;
	D->circumCenter.z=circz_s;
	create_delunay(Atoms,D);

}
/** Given an existing Delaunay edge (I, J) around ATOM and one known face
 *  neighbor K, walks the incremental-flip / "gift wrapping" construction to
 *  find the fourth vertex that completes the next Delaunay tetrahedron on
 *  the far side of that edge, and creates it via create_delunay(). */
delunay* constr_del(atom *ATOM,atom Atoms[],vect vI,vect vJ,int sign,int I,int J,int K, long double rI, long double rJ,delunay *D,int debug=0)
{
	long double rS=ATOM->radius;
	long double Y_MIN=1e10;
	long double circx,circy,circz;
	int A1,A2,A3,A4;
	A1=ATOM->index;
	A2=I;
	A3=J;
	A4=-1;
	//if(not (I==11027 and J==5430 and K==3940))
	//{
	//	debug=0;
	//}
	if(not(ATOM->check_if_contiguous(I) and ATOM->check_if_contiguous(J) and ATOM->check_if_contiguous(K)))
	{
		std::cerr<<"wrong delunay in constr del\n";
		std::cerr<<ATOM->index<<"\t"<<I<<"\t"<<J<<"\t"<<K<<"\n";
		exit(EXIT_FAILURE);
	}
	if(debug)
	{
		std::cout<<std::setprecision(16);
	}
	for(int atom=0; atom<nAtoms; atom++)
		//for(auto atom : incompleteAtoms)
	{
		int L =atom;
		if(L==K || L ==I || L==J)
			continue;

		/* if L is to be successful candidate, then 
		   Then we will be adding L,I,J L,I,ATOM, L,J,ATOM, I,J,ATOM 
		   faces to the system. 

		   We need to make sure that adding so does not create 
		 */

		int sign_N;
		vect vL = vectDiff(&Atoms[L].p,&ATOM->p);
		long double rL=Atoms[L].radius;

		vect vJIcross = cross_product(&vJ,&vI);
		long double overlap=innerProduct(&vL,&vJIcross);
		if(overlap<0.)
			sign_N=1;
		else
			sign_N=-1;
		if(sign!=sign_N) // The atom J lies on the +ve z axis 
		{
			long double DISIsq=magnitudeSq(&vI);
			long double DISJsq=magnitudeSq(&vJ);
			long double DISLsq=magnitudeSq(&vL);
			long double B[3],A[3][3];
			B[0]=(DISIsq+rS*rS-rI*rI)/2.;
			B[1]=(DISJsq+rS*rS-rJ*rJ)/2.;
			B[2]=(DISLsq+rS*rS-rL*rL)/2.;
			A[0][0]=vI.x;
			A[0][1]=vI.y;
			A[0][2]=vI.z;

			A[1][0]=vJ.x;
			A[1][1]=vJ.y;
			A[1][2]=vJ.z;

			A[2][0]=vL.x;
			A[2][1]=vL.y;
			A[2][2]=vL.z;
			vect center;
			center=cramer(A,B);
			long double norm=sqrtl(innerProduct(&vJIcross,&vJIcross));
			vect X = center;
			X.x=center.x;
			X.y=center.y;
			X.z=center.z;
			long double Y_AXIS;//=sqrtl(powl(X-X2,2)+powl(Y-Y2,2)+powl(Z-Z2,2));
					   ////cout<<Y_AXIS<<" 1\n";
			Y_AXIS=1./norm*(innerProduct(&X,&vJIcross));
			Y_AXIS=abs(Y_AXIS);
			int sign_C;
			long double overlap=(innerProduct(&X,&vJIcross)); 
			if(overlap<0.) // The CMSTD is negative 
				sign_C=1;
			else
				sign_C=-1;
			if(sign_C!=sign_N)
			{
				Y_AXIS=-1.*Y_AXIS;
			}
			//cout<<Y_AXIS<<"\n";
			X.x=X.x+ATOM->p.x;
			X.y=X.y+ATOM->p.y;
			X.z=X.z+ATOM->p.z;

			if(Y_AXIS<Y_MIN)
			{
				Y_MIN=Y_AXIS;
				A4=L;
				circx=X.x;
				circy=X.y;
				circz=X.z;
			}
		}
		//}//j loop
}
if(A4!= -1)
{
	std::vector<int> EV {A1,A2,A3,A4};
	//std::sort(EV.begin(),EV.end());
	D=new delunay;
	fullSet.addTetrahedron(D);
	D->AT[0]=EV[0];
	D->AT[1]=EV[1];
	D->AT[2]=EV[2];
	D->AT[3]=EV[3];
	D->circumCenter.x=circx;
	D->circumCenter.y=circy;
	D->circumCenter.z=circz;
	create_delunay(Atoms,D);
	return D;
}
else
{
	return D;
}
}
/** Marks Delaunay tetrahedron D (built from ATOM, I, J, K) as lying on the
 *  outer convex hull of ATOM's local neighbor set when no fourth-vertex
 *  candidate was found on one side of a face, i.e. that face bounds empty
 *  space rather than another tetrahedron. */
void convexHull(delunay* D, atom* ATOM, atom Atoms[], int I, int J, int K)
{
	atom* Ai;
	atom* Aj;
	atom* Ak;

	Ai=&(Atoms[I]);
	Aj=&(Atoms[J]);
	Ak=&(Atoms[K]);

	if(I>J)
	{
		ATOM->increment_part_c(I,J);
	}
	else
	{
		ATOM->increment_part_c(J,I);
	}

	if(ATOM->index > Aj->index)
		Ai->increment_part_c(ATOM->index,Aj->index);
	else
		Ai->increment_part_c(Aj->index,ATOM->index);

	if(ATOM->index>Ai->index)
		Aj->increment_part_c(ATOM->index,Ai->index);
	else
		Aj->increment_part_c(Ai->index,ATOM->index);
	if(!ATOM->checkIfIncompleteFace())
	{
		/* This atom is completed */
		incompleteAtoms.erase(std::remove(incompleteAtoms.begin(), incompleteAtoms.end(), ATOM->index), incompleteAtoms.end());
	}
	if(!Ai->checkIfIncompleteFace())
	{
		/* This atom is completed */
		incompleteAtoms.erase(std::remove(incompleteAtoms.begin(), incompleteAtoms.end(), Ai->index), incompleteAtoms.end());
	}
	if(!Aj->checkIfIncompleteFace())
	{
		/* This atom is completed */
		incompleteAtoms.erase(std::remove(incompleteAtoms.begin(), incompleteAtoms.end(), Aj->index), incompleteAtoms.end());
	}
	D->hull=1;
}
/** Expands the local Delaunay tessellation around ATOM until every face
 *  incident to it is shared by exactly two tetrahedra (or flagged as a
 *  convex-hull face), i.e. until ATOM's full local neighborhood is closed.
 *  This is the main per-particle driver called from the BFS loop in
 *  main(). */
void completeDelunayTessellation(atom *ATOM,atom Atoms[],int nAtoms)
{
	//long double MIN = 1e10;
	bool newElementAdded=1;
	/* first we have to go through existing delaney tetrahedrons this 
	   atom has and find an incomplete side, i.e a face that has the 
	   current atom but is not part of TWO tetrahedrons. */
	//if(!ATOM->checkIfIncompleteFace())
	//{
	//	return;
	//}	
	while(newElementAdded)
	{
		std::cout<<std::flush;
		newElementAdded=false;
		std::vector<std::pair<std::pair<int,int>, int>> entries(ATOM->part_c.begin(), ATOM->part_c.end());
		//std::cout<<ATOM->contiguous.size()<<"\n";
		for(const std::pair<std::pair<int,int>,int>& pairIJ: entries)
		{
			if(ATOM->get_part_c(pairIJ.first.first,pairIJ.first.second)==1) // incomplete face
			{
				int I = pairIJ.first.first;
				int J = pairIJ.first.second;
				int K;
				for(int k=0; k<4; k++)
				{
					if(ATOM->index != ATOM->DelaunayTessellations[pairIJ.first].first->AT[k]) // So that K is not this atom 
					{
						if(ATOM->DelaunayTessellations[pairIJ.first].first->AT[k]!=I and ATOM->DelaunayTessellations[pairIJ.first].first->AT[k]!=J)
						{
							K = ATOM->DelaunayTessellations[pairIJ.first].first->AT[k];
							break;
						}
					}
				}
				vect vI = vectDiff(&Atoms[I].p,&ATOM->p);
				vect vJ = vectDiff(&Atoms[J].p,&ATOM->p);
				vect vK = vectDiff(&Atoms[K].p,&ATOM->p);

				long double rI=Atoms[I].radius;
				long double rJ=Atoms[J].radius;
				int sign;
				vect vJIcross = cross_product(&vJ,&vI);
				long double overlap=innerProduct(&vK,&vJIcross);
				if(overlap<0.)
					sign=1;
				else
					sign=-1;
				delunay *D_TWO = nullptr;
				D_TWO=constr_del(ATOM,Atoms,vI,vJ,sign,I,J,K,rI,rJ,D_TWO);
				//if(D_TWO)
				//{
				//	print_delunay(outputFile,D_TWO,Atoms);
				//	outputFile<<std::flush;
				//}
				if(!D_TWO)
				{
					/* part of convex hull */
					convexHull(ATOM->DelaunayTessellations[pairIJ.first].first,ATOM,Atoms,I,J,K);
				}
				newElementAdded=true;
			}
		}
	}
}

long double volume_tetrahedron( long double Ax, long double Ay, long double Az, long double Ex, long double Ey, long double Ez, long double Bx, long double By, long double Bz, long double vx, long double vy, long double vz, long double r,int compliment=0)
{
	long double ABx,ABy,ABz,EBx,EBy,EBz,VEx,VEy,VEz,DISB,DISBE,DISVE;
	//cout<<"draw sphere\t{";
	//cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\tradius 0.3\n";
	//cout<<"draw sphere\t{";
	//cout<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\tradius 0.3\n";
	//cout<<"draw sphere\t{";
	//cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\tradius 0.3\n";
	//cout<<"draw sphere\t{";
	//cout<<vx<<"\t"<<vy<<"\t"<<vz<<"}\tradius 0.3\n";
	ABx=Ax-Bx;
	ABy=Ay-By;
	ABz=Az-Bz;
	EBx=Ex-Bx;
	EBy=Ey-By;
	EBz=Ez-Bz;
	VEx=Ex-vx;
	VEy=Ey-vy;
	VEz=Ez-vz;
	ABx=ABx-boxx*PBCx*lroundl(ABx/boxx);
	ABy=ABy-boxy*PBCy*lroundl(ABy/boxy);
	ABz=ABz-boxz*PBCz*lroundl(ABz/boxz);

	EBx=EBx-boxx*PBCx*lroundl(EBx/boxx);
	EBy=EBy-boxy*PBCy*lroundl(EBy/boxy);
	EBz=EBz-boxz*PBCz*lroundl(EBz/boxz);

	VEx=VEx-boxx*PBCx*lroundl(VEx/boxx);
	VEy=VEy-boxy*PBCy*lroundl(VEy/boxy);
	VEz=VEz-boxz*PBCz*lroundl(VEz/boxz);


	DISB=ABx*ABx+ABy*ABy+ABz*ABz;
	DISBE=EBx*EBx+EBy*EBy+EBz*EBz;
	DISVE=VEx*VEx+VEy*VEy+VEz*VEz;
	Ax=0.;
	Ay=0.;
	Az=0.;
	Bx=sqrtl(DISB);
	By=0.;
	Bz=0.;
	Ex=sqrtl(DISB);
	Ey=sqrtl(DISBE);
	Ez=0.;
	vx=sqrtl(DISB);
	vy=sqrtl(DISBE);
	vz=sqrtl(DISVE);
	long double x0,y0,z0;
	long double x0sq,y0sq,z0sq;
	x0=Bx;
	y0=Ey;
	z0=vz;
	x0sq=DISB;
	y0sq=DISBE;
	z0sq=DISVE;
	long double rB,rV,rE;
	rV=sqrtl(vx*vx+vy*vy+vz*vz);
	rE=sqrtl(Ex*Ex+Ey*Ey+Ez*Ez);
	rB=sqrtl(DISB);
	long double theta,x2,y2;
	x2=r*x0/rE;
	y2=r*y0/rE;
	theta=atanl(z0/y0);
	long double Vc,Vt;
	Vt=(x0*y0*z0)/6.;
	Vc=0.;
	if(r<rB)
	{
		if(abs((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq)))>1)
			Vc=(r*r*r/6.)*(2*theta);
		else
			Vc=(r*r*r/6.)*(2*theta-M_PI/2.-asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq))));
		//cout<<Vc<<"\n";
	}
	if(rB<r && r<rE)
	{
		if(abs((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq)))>1)
			Vc=theta/2.*(r*r*x0-x0*x0*x0/3.);
		else
			Vc=theta/2.*(r*r*x0-x0*x0*x0/3.)-(r*r*r/6.)*(M_PI/2.+asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq))));
		//cout<<Vc<<"\n";
	}
	if(rE<r && r<rV)
	{
		Vc=0.5*(theta-M_PI/2.+asin(y0/sqrtl(r*r-x0sq)))*(r*r*x0-x0*x0*x0/3.)+x0*y0/6.*sqrtl(r*r-rE*rE)+r*r*r/6.*asinl((x2*x2-y2*y2-x0sq)/(r*r-x0sq))-r*r*r/6.*asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq)));
		//cout<<Vc<<"\n";
	}
	if(Vt<Vc)
	{
		//cout<<"help \t"<<Vt<<"\t"<<Vc<<"\t"<<Vt-Vc<<"\n";
	}
	if(!compliment)
		return Vt-Vc;
	else
	{
		//cout<<"###\t"<<Vt<<"\t"<<Vc<<"\t"<<r<<"\t"<<rB<<"\t"<<rE<<"\t"<<rV<<"\n";
		if(Vc==0.)
			return Vt;
		else
			return Vc;
	}
}
/** Orientation predicate: sign of the volume of the tetrahedron
 *  (A1, A2, A3, V), used to test which side of the plane through A1,A2,A3
 *  the point V lies on. */
int sign_aaav( long double A1x, long double A1y, long double A1z, long double A2x, long double A2y, long double A2z, long double A3x, long double A3y, long double A3z, long double A4x, long double A4y, long double A4z, long double Vx, long double Vy, long double Vz)
{
	long double a1x,a1y,a1z;
	long double a2x,a2y,a2z;
	long double a3x,a3y,a3z;
	long double vx,vy,vz;
	a1x=A2x-A1x;
	a1y=A2y-A1y;
	a1z=A2z-A1z;
	a2x=A3x-A1x;
	a2y=A3y-A1y;
	a2z=A3z-A1z;
	a3x=A4x-A1x;
	a3y=A4y-A1y;
	a3z=A4z-A1z;
	vx=Vx-A1x;
	vy=Vy-A1y;
	vz=Vz-A1z;
	a3x=a3x-boxx*PBCx*lroundl(a3x/boxx);
	a3y=a3y-boxy*PBCy*lroundl(a3y/boxy);
	a3z=a3z-boxz*PBCz*lroundl(a3z/boxz);

	vx=vx-boxx*PBCx*lroundl(vx/boxx);
	vy=vy-boxy*PBCy*lroundl(vy/boxy);
	vz=vz-boxz*PBCz*lroundl(vz/boxz);

	a1x=a1x-boxx*PBCx*lroundl(a1x/boxx);
	a1y=a1y-boxy*PBCy*lroundl(a1y/boxy);
	a1z=a1z-boxz*PBCz*lroundl(a1z/boxz);

	a2x=a2x-boxx*PBCx*lroundl(a2x/boxx);
	a2y=a2y-boxy*PBCy*lroundl(a2y/boxy);
	a2z=a2z-boxz*PBCz*lroundl(a2z/boxz);

	long double a,b,c;
	long double overlap1,overlap2;
	int sign1,sign2;
	a=a2z*a1y-a1z*a2y;
	b=a1z*a2x-a1x*a2z;
	c=a2y*a1x-a1y*a2x;
	overlap1=a*a3x+b*a3y+c*a3z;
	if(overlap1<0.)
		sign1=1;
	else
		sign1=-1;
	overlap2=a*vx+b*vy+c*vz;
	if(overlap2<0.)
		sign2=1;
	else
		sign2=-1;
	if(sign1==sign2)
	{
		return 1; 
	}
	else
	{
		return -1;
	}
}
/** Orientation predicate: sign of the projection of edge (A1..E) relative
 *  to the plane through A1, A2, A3, used to test which side of a Delaunay
 *  face an edge direction points toward. */
int sign_aaae( long double A1x, long double A1y, long double A1z, long double A2x, long double A2y, long double A2z, long double A3x, long double A3y, long double A3z, long double Ex, long double Ey, long double Ez)
{
	long double a1x,a1y,a1z;
	long double a2x,a2y,a2z;
	long double ex,ey,ez;
	a1x=A1x-A2x;
	a1y=A1y-A2y;
	a1z=A1z-A2z;
	a2x=A3x-A2x;
	a2y=A3y-A2y;
	a2z=A3z-A2z;
	ex=Ex-A2x;
	ey=Ey-A2y;
	ez=Ez-A2z;
	a1x=a1x-boxx*PBCx*lroundl(a1x/boxx);
	a1y=a1y-boxy*PBCy*lroundl(a1y/boxy);
	a1z=a1z-boxz*PBCz*lroundl(a1z/boxz);

	a2x=a2x-boxx*PBCx*lroundl(a2x/boxx);
	a2y=a2y-boxy*PBCy*lroundl(a2y/boxy);
	a2z=a2z-boxz*PBCz*lroundl(a2z/boxz);

	ex=ex-boxx*PBCx*lroundl(ex/boxx);
	ey=ey-boxy*PBCy*lroundl(ey/boxy);
	ez=ez-boxz*PBCz*lroundl(ez/boxz);

	long double a,b,c;
	long double a1,b1,c1;
	a=a2z*a1y-a1z*a2y;
	b=a1z*a2x-a1x*a2z;
	c=a2y*a1x-a1y*a2x;
	a1=c*a1y-a1z*b;
	b1=a1z*a-a1x*c;
	c1=b*a1x-a1y*a;
	long double overlap1,overlap2;
	int sign1,sign2;
	overlap1=a1*ex+b1*ey+c1*ez;
	if(overlap1<0.)
	{
		sign1=1;
	}
	else
	{
		sign1=-1;
	}
	overlap2=a1*a2x+b1*a2y+c1*a2z;
	if(overlap2<0.)
	{
		sign2=1;
	}
	else
	{
		sign2=-1;
	}
	if(sign1==sign2)
	{
		return 1;
	}
	else
		return -1;
}
/** Tests whether Voronoi vertex v (the circumcenter of tetrahedron D) truly
 *  lies inside D's own tetrahedron. Used to peel away convex-hull
 *  tetrahedra whose circumcenter falls outside their own volume (which
 *  would otherwise be wrongly counted as void space at the packing's
 *  boundary). */
bool inside_delunay(vertice *v,delunay *D,atom Atoms[],int nAtoms)
{
	int sign1,sign2;
	long double overlap1,overlap2;	
	atom *a1=&(Atoms[D->AT[0]]);
	atom *a2=&(Atoms[D->AT[1]]);;
	atom *a3=&(Atoms[D->AT[2]]);;
	atom *a4=&(Atoms[D->AT[3]]);;
	vect v2 = vectDiff(&a2->p,&a1->p);
	vect v3 = vectDiff(&a3->p,&a1->p);
	vect v4 = vectDiff(&a4->p,&a1->p);
	vect V = vectDiff(v->p,&a1->p);
	//V1.x=a2->p.x-a1->p.x;
	//V1.y=a2->p.y-a1->p.y;
	//V1.z=a2->p.z-a1->p.z;

	//V2.x=a3->p.x-a1->p.x;
	//V2.y=a3->p.y-a1->p.y;
	//V2.z=a3->p.z-a1->p.z;

	//V.x=v->p->x-a1->p.x;
	//V.y=v->p->y-a1->p.y;
	//V.z=v->p->z-a1->p.z;

	//A.x=a4->p.x-a1->p.x;
	//A.y=a4->p.y-a1->p.y;
	//A.z=a4->p.z-a1->p.z;

	vect cross;
	cross=cross_product(&v2,&v3);
	overlap1=innerProduct(&cross,&v4);//cross_12.x*A.x+cross_12.y*A.y+cross_12.z*A.z;
	overlap2=innerProduct(&cross,&V);//cross_12.x*V.x+cross_12.y*V.y+cross_12.z*V.z;
	if(overlap1<0.)
		sign1=1;
	else
		sign1=-1;

	if(overlap2<0.)
		sign2=1;
	else
		sign2=-1;

	if(sign1==sign2)
	{
		cross=cross_product(&v2,&v4);
		overlap1=innerProduct(&cross,&v3);//cross_12.x*A.x+cross_12.y*A.y+cross_12.z*A.z;
		overlap2=innerProduct(&cross,&V);//cross_12.x*V.x+cross_12.y*V.y+cross_12.z*V.z;
		if(overlap1<0.)
			sign1=1;
		else
			sign1=-1;

		if(overlap2<0.)
			sign2=1;
		else
			sign2=-1;
		if(sign1==sign2)
		{
			cross=cross_product(&v3,&v4);
			overlap1=innerProduct(&cross,&v2);//cross_12.x*A.x+cross_12.y*A.y+cross_12.z*A.z;
			overlap2=innerProduct(&cross,&V);//cross_12.x*V.x+cross_12.y*V.y+cross_12.z*V.z;
			if(overlap1<0.)
				sign1=1;
			else
				sign1=-1;

			if(overlap2<0.)
				sign2=1;
			else
				sign2=-1;
			if(sign1==sign2)
			{
				vect v2 = vectDiff(&a1->p,&a2->p);
				vect v3 = vectDiff(&a3->p,&a2->p);
				vect v4 = vectDiff(&a4->p,&a2->p);
				vect V = vectDiff(v->p,&a2->p);
				cross=cross_product(&v3,&v4);
				overlap1=innerProduct(&cross,&v2);//cross_12.x*A.x+cross_12.y*A.y+cross_12.z*A.z;
				overlap2=innerProduct(&cross,&V);//cross_12.x*V.x+cross_12.y*V.y+cross_12.z*V.z;
				if(overlap1<0.)
					sign1=1;
				else
					sign1=-1;

				if(overlap2<0.)
					sign2=1;
				else
					sign2=-1;
				if(sign1==sign2)
					return true;
				else 
					return false;
			}
			else
				return false;
		}
		else
			return false;
	}	
	else
		return false;
}
/** Computes the exact void-space volume attributable to Voronoi vertex v,
 *  i.e. the piece of the empty-space region around this vertex that is not
 *  covered by any of the (dilated) particle spheres. Summed over every
 *  void vertex in a cluster, this gives that cavity's total void volume. */
long double void_vol(vertice *v,atom Atoms[])
{
	long double vol=0.;
	vect A1,A2,A3,A4;
	long double r1,r2,r3,r4;
	int S123,S124,S234,S134;
	int A1A2E123,A1A3E123,A2A3E123;
	int A1A2E124,A1A4E124,A2A4E124;
	int A1A3E134,A1A4E134,A3A4E134;
	int A2A3E234,A2A4E234,A3A4E234;
	// long double Vx,Vy,Vz;
	long double Vx=v->p->x;
	long double Vy=v->p->y;
	long double Vz=v->p->z;


	A1=Atoms[v->D->AT[0]].p;
	A2=Atoms[v->D->AT[1]].p;
	A3=Atoms[v->D->AT[2]].p;
	A4=Atoms[v->D->AT[3]].p;

	r1=Atoms[v->D->AT[0]].radius;
	r2=Atoms[v->D->AT[1]].radius;
	r3=Atoms[v->D->AT[2]].radius;
	r4=Atoms[v->D->AT[3]].radius;

	vect E123,E124,E134,E234;
	vect B12,B13,B14,B23,B34,B24;

	S123=sign_aaav(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,Vx,Vy,Vz);
	S124=sign_aaav(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,Vx,Vy,Vz);
	S234=sign_aaav(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,Vx,Vy,Vz);
	S134=sign_aaav(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,Vx,Vy,Vz);

	E123=v->D->MID[0][1][2];
	E124=v->D->MID[0][1][3];
	E134=v->D->MID[0][2][3];
	E234=v->D->MID[1][2][3];

	B12=v->D->MIDP[0][1];
	B13=v->D->MIDP[0][2];
	B14=v->D->MIDP[0][3];
	B23=v->D->MIDP[1][2];
	B24=v->D->MIDP[1][3];
	B34=v->D->MIDP[2][3];

	A1A2E123=sign_aaae(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,E123.x,E123.y,E123.z);
	A1A3E123=sign_aaae(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A2.x,A2.y,A2.z,E123.x,E123.y,E123.z);
	A2A3E123=sign_aaae(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A1.x,A1.y,A1.z,E123.x,E123.y,E123.z);

	A1A2E124=sign_aaae(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,E124.x,E124.y,E124.z);
	A1A4E124=sign_aaae(A1.x,A1.y,A1.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,E124.x,E124.y,E124.z);
	A2A4E124=sign_aaae(A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,E124.x,E124.y,E124.z);

	A1A3E134=sign_aaae(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,E134.x,E134.y,E134.z);
	A1A4E134=sign_aaae(A1.x,A1.y,A1.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,E134.x,E134.y,E134.z);
	A3A4E134=sign_aaae(A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,E134.x,E134.y,E134.z);

	A2A3E234=sign_aaae(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,E234.x,E234.y,E234.z);
	A2A4E234=sign_aaae(A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,E234.x,E234.y,E234.z);
	A3A4E234=sign_aaae(A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,E234.x,E234.y,E234.z);

	vol=vol+S123*A1A2E123*volume_tetrahedron(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1);
	vol=vol+S123*A1A3E123*volume_tetrahedron(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1);
	vol=vol+S123*A1A2E123*volume_tetrahedron(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2);
	vol=vol+S123*A2A3E123*volume_tetrahedron(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2);
	vol=vol+S123*A1A3E123*volume_tetrahedron(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3);
	vol=vol+S123*A2A3E123*volume_tetrahedron(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3);

	vol=vol+S124*A1A2E124*volume_tetrahedron(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1);
	vol=vol+S124*A1A4E124*volume_tetrahedron(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1);
	vol=vol+S124*A1A2E124*volume_tetrahedron(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2);
	vol=vol+S124*A2A4E124*volume_tetrahedron(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2);
	vol=vol+S124*A1A4E124*volume_tetrahedron(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4);
	vol=vol+S124*A2A4E124*volume_tetrahedron(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4);

	vol=vol+S134*A1A3E134*volume_tetrahedron(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1);
	vol=vol+S134*A1A4E134*volume_tetrahedron(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1);
	vol=vol+S134*A1A3E134*volume_tetrahedron(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3);
	vol=vol+S134*A3A4E134*volume_tetrahedron(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3);
	vol=vol+S134*A1A4E134*volume_tetrahedron(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4);
	vol=vol+S134*A3A4E134*volume_tetrahedron(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4);

	vol=vol+S234*A2A3E234*volume_tetrahedron(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2);
	vol=vol+S234*A2A4E234*volume_tetrahedron(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2);
	vol=vol+S234*A2A3E234*volume_tetrahedron(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3);
	vol=vol+S234*A3A4E234*volume_tetrahedron(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3);
	vol=vol+S234*A2A4E234*volume_tetrahedron(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4);
	vol=vol+S234*A3A4E234*volume_tetrahedron(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4);

	return vol;
}
/*
 * Prints command-line usage information to stderr.
 */
static void printUsage(const char *progName)
{
	std::cerr
		<< "Usage: " << progName << " <input_file> <probe_radius> [output_dir] [sample_atom_index]\n\n"
		<< "  input_file          Path to the packing file (see README.md for the exact format).\n"
		<< "  probe_radius        Non-negative probe radius (r_cut), added to every particle's\n"
		<< "                      radius before the void-volume analysis. Use 0 for the exact\n"
		<< "                      void volume of the packing as given.\n"
		<< "  output_dir          Optional. Directory where all output files are written; created\n"
		<< "                      if it does not already exist. Defaults to the current directory.\n"
		<< "  sample_atom_index   Optional. 0-based index of one particle whose local Voronoi cell\n"
		<< "                      is additionally written to 'demo.txt' (and appended to\n"
		<< "                      'voronoiEdges.txt') for a quick, small-scale visual check.\n"
		<< "                      Defaults to 0. Must satisfy 0 <= index < number of particles.\n";
}

int main( int argc, char * argv[] )
{
	if(argc < 3 || argc > 5)
	{
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	// --- probe radius -------------------------------------------------
	long double requestedProbeRadius;
	try
	{
		requestedProbeRadius = std::stold(argv[2]);
	}
	catch(const std::exception &e)
	{
		std::cerr<<"Error: could not parse probe_radius '"<<argv[2]<<"': "<<e.what()<<"\n";
		return EXIT_FAILURE;
	}
	if(requestedProbeRadius < 0.0L)
	{
		std::cerr<<"Error: probe_radius must be >= 0 (got "<<argv[2]<<")\n";
		return EXIT_FAILURE;
	}
	r_cut = requestedProbeRadius;

	// --- output directory ----------------------------------------------
	std::string outDir = (argc >= 4) ? argv[3] : ".";
	while(outDir.size() > 1 && outDir.back() == '/')
		outDir.pop_back();
	if(mkdir(outDir.c_str(), 0755) != 0 && errno != EEXIST)
	{
		std::cerr<<"Error: could not create output directory '"<<outDir<<"': "<<std::strerror(errno)<<"\n";
		return EXIT_FAILURE;
	}
	const std::string outPrefix = outDir + "/";

	// --- input file ------------------------------------------------------
	std::ifstream infile(argv[1]);
	if(!infile)
	{
		std::cerr<<"Error: could not open input file '"<<argv[1]<<"'\n";
		return EXIT_FAILURE;
	}
	infile>>nAtoms;
	if(!infile || nAtoms <= 0)
	{
		std::cerr<<"Error: input file '"<<argv[1]<<"' must begin with a positive particle count.\n";
		return EXIT_FAILURE;
	}

	// --- optional sample-atom index for the small demo/QC export -------
	int sampleAtom = 0;
	if(argc >= 5)
	{
		try
		{
			sampleAtom = std::stoi(argv[4]);
		}
		catch(const std::exception &e)
		{
			std::cerr<<"Error: could not parse sample_atom_index '"<<argv[4]<<"': "<<e.what()<<"\n";
			return EXIT_FAILURE;
		}
	}
	if(sampleAtom < 0 || sampleAtom >= nAtoms)
	{
		std::cerr<<"Error: sample_atom_index ("<<sampleAtom<<") must be in [0, "<<nAtoms<<").\n";
		return EXIT_FAILURE;
	}

	atom *Atoms = new (nothrow) atom[nAtoms];
	if(!Atoms)
	{
		std::cerr<<"Error: failed to allocate memory for "<<nAtoms<<" particles.\n";
		return EXIT_FAILURE;
	}

	outputFile.open(outPrefix + "delaunay_edges.txt");

	long double b,c,d,e;
	// SIGMA is a leftover unit-rescaling hook (always 1 here); kept so the
	// input length scale can be trivially reduced/expanded in the future
	// without touching the geometry routines below.
	long double SIGMA = 1;
	long double volumeBox=0.;
	long double volumeAtoms=0.;
	int counter=0;
	infile>>boxx>>boxy>>boxz;
	if(!infile || boxx <= 0.0L || boxy <= 0.0L || boxz <= 0.0L)
	{
		std::cerr<<"Error: input file '"<<argv[1]<<"' must specify a positive box size (boxx boxy boxz) on the second line.\n";
		delete [] Atoms;
		return EXIT_FAILURE;
	}
	boxx = boxx/SIGMA;
	boxy = boxy/SIGMA;
	boxz = boxz/SIGMA;
	volumeBox=boxx*boxy*boxz;
	while(infile>>b>>c>>d>>e)
	{
		Atoms[counter].p.x=b/SIGMA;
		Atoms[counter].p.y=c/SIGMA;
		Atoms[counter].p.z=d/SIGMA;
		Atoms[counter].radius=(e/SIGMA)+r_cut;
		Atoms[counter].index=counter;
		volumeAtoms=volumeAtoms+Atoms[counter].radius*Atoms[counter].radius*Atoms[counter].radius*4./3.*M_PI;
		incompleteAtoms.push_back(counter);
		allAtoms.push_back(counter);
		counter++;
		if(counter==nAtoms)
		{
			break;
		}
	}
	if(counter != nAtoms)
	{
		std::cerr<<"Error: input file '"<<argv[1]<<"' declares "<<nAtoms
			<<" particles but only "<<counter<<" coordinate lines were found.\n";
		delete [] Atoms;
		return EXIT_FAILURE;
	}

	// Breadth-first traversal of the Delaunay neighbor graph, starting from
	// particle 0. This assumes the packing forms a single connected Delaunay
	// complex reachable from particle 0, which holds for any ordinary sphere
	// packing; see README.md for the periodic-boundary-condition caveat.
	atomsToAnalyze.push_back(0);
	processedOrQueuedAtoms.insert(0);
	//for(int i=0; i< nAtoms; i++)
	int index=0;
	while (!atomsToAnalyze.empty()) 
	{
		// Get the next atom from the front of the queue
		int iAtom = atomsToAnalyze.front();

		atomsToAnalyze.pop_front();
		std::cout << "Processing atom: " << std::left << std::setw(8) << iAtom << "\t"<<index<<"\t" << std::fixed << std::setprecision(2) << std::setw(6) << (index* 100.0 / nAtoms) << "%  \r" << std::flush;
		index++;
		if(Atoms[iAtom].delunayTetrahedrons.size()==0)
		{
			first_delunay(&(Atoms[iAtom]),Atoms);
		}
		if(Atoms[iAtom].delunayTetrahedrons.size())
			completeDelunayTessellation(&(Atoms[iAtom]),Atoms,nAtoms);
		//std::cout<<iAtom<<"\t"<<Atoms[iAtom].contiguous.size()<<"\n";
		allAtoms.erase(std::remove(allAtoms.begin(), allAtoms.end(), iAtom), allAtoms.end());
		//for (const auto& element : allAtoms)
		//{
		//	std::cout << element << " ";
		//}
		//cout<<"\n";
		for(int neighborAtom:Atoms[iAtom].contiguous)
		{
			if (processedOrQueuedAtoms.count(neighborAtom) == 0) {
				// If it's a new atom, add it to both the "to-do" list
				// and the "has-been-seen" list.
				processedOrQueuedAtoms.insert(neighborAtom);
				atomsToAnalyze.push_back(neighborAtom);
			}
		}
	}
	std::cout<<"processing done\n";
	std::cout<<std::setprecision(16);
	// Diagnostic: the sum of all Delaunay tetrahedron volumes should equal
	// the box volume for a geometrically consistent periodic tessellation.
	std::cout<<"Sum of Delaunay tetrahedron volumes vs. box volume:\n";
	std::cout<<convex_vol<<"\t"<<boxx*boxy*boxz<<"\n";
	std::ofstream atomInfoVMDOutputFile(outPrefix + "atoms.txt");
	std::ofstream atomInfoOutputFile(outPrefix + "atomsPOVRAY.txt");
	atomInfoOutputFile<<std::setprecision(16);
	for(int iAtom=0; iAtom< nAtoms; iAtom++)
	{
		atomInfoVMDOutputFile<<"draw sphere\t{";
		atomInfoVMDOutputFile<<Atoms[iAtom].p.x<<"\t"<<Atoms[iAtom].p.y<<"\t"<<Atoms[iAtom].p.z<<"}\tradius\t"<<Atoms[iAtom].radius<<"\t"<<"resolution\t10\n";
		atomInfoOutputFile<<"1,\t";
		atomInfoOutputFile<<Atoms[iAtom].p.x<<",\t"<<Atoms[iAtom].p.z<<",\t"<<Atoms[iAtom].p.y<<","<<Atoms[iAtom].radius;
		atomInfoOutputFile<<",\n";
	}
	outputFile<<"draw color 11\n";
	for(auto D : fullSet.tetrahedrons)
	{
		print_delunay(outputFile,D,Atoms);
	}
	std::ofstream voroInfoOutputFile(outPrefix + "voronoiEdges.txt");
	voroInfoOutputFile<<std::setprecision(16);
	voroInfoOutputFile<<"draw color 11\n";
	for(int atom =0; atom < nAtoms; atom++)
	{
		// Voronoi (radical-plane) edges for every atom, in the same
		// POV-Ray-style CSV format as the atom spheres above: shape code 2
		// followed by the two circumcenter endpoints and a display radius.
		for(auto delFace: Atoms[atom].DelaunayTessellations)
		{
			auto v1 = delFace.second.first->circumCenter;
			auto v2 = delFace.second.second->circumCenter;
			vect V1 = vectDiff(&v1, &Atoms[atom].p);
			vect V2 = vectDiff(&v2, &Atoms[atom].p);
			atomInfoOutputFile<<"2,\t";
			atomInfoOutputFile<<V1.x+Atoms[atom].p.x<<",\t"<<V1.z+Atoms[atom].p.z<<",\t"<<V1.y+Atoms[atom].p.y<<",\t";
			atomInfoOutputFile<<V2.x+Atoms[atom].p.x<<",\t"<<V2.z+Atoms[atom].p.z<<",\t"<<V2.y+Atoms[atom].p.y<<", 5e-4,\n";
		}
	}

	// Small single-particle export ('demo.txt' + a few extra VMD draw
	// commands appended to voronoiEdges.txt) for quickly sanity-checking
	// one particle's local geometry without rendering the whole system.
	// Rendering every particle's Voronoi cell as individual VMD 'draw'
	// commands would be extremely slow, so this is intentionally scoped to
	// a single, user-selectable particle (--sample_atom_index, default 0).
	{
		std::ofstream demoFile(outPrefix + "demo.txt");
		demoFile<<std::setprecision(16);
		const int atom = sampleAtom;

		voroInfoOutputFile<<"# \t"<<atom<<"\n";
		voroInfoOutputFile<<"draw sphere\t{";
		voroInfoOutputFile<<Atoms[atom].p.x<<"\t"<<Atoms[atom].p.z<<"\t"<<Atoms[atom].p.y<<"}\tradius\t"<<Atoms[atom].radius<<"\t"<<"resolution\t10\n";
		demoFile<<"3,\t";
		demoFile<<Atoms[atom].p.x<<",\t"<<Atoms[atom].p.z<<",\t"<<Atoms[atom].p.y<<",\t"<<Atoms[atom].radius<<",\n";
		for(auto neiA : Atoms[atom].contiguous)
		{
			vect V = vectDiff(&Atoms[atom].p,&Atoms[neiA].p);
			demoFile<<"1,\t";
			demoFile<<Atoms[atom].p.x-V.x<<",\t"<<Atoms[atom].p.z-V.z<<",\t"<<Atoms[atom].p.y-V.y<<",\t"<<Atoms[neiA].radius<<",\n";
			voroInfoOutputFile<<"draw sphere\t{";
			voroInfoOutputFile<<Atoms[atom].p.x-V.x<<"\t"<<Atoms[atom].p.z-V.z<<"\t"<<Atoms[atom].p.y-V.y<<"}\tradius\t"<<Atoms[neiA].radius<<"\t"<<"resolution\t10\n";
		}
		for(auto element:Atoms[atom].part_c)
		{
			int I=element.first.first;
			int J=element.first.second;
			vect V1 = vectDiff(&Atoms[atom].p,&Atoms[I].p);
			vect V2 = vectDiff(&Atoms[atom].p,&Atoms[J].p);
			demoFile<<"4,\t";
			demoFile<<Atoms[atom].p.x-V1.x<<",\t"<<Atoms[atom].p.z-V1.z<<",\t"<<Atoms[atom].p.y-V1.y<<",\t";
			demoFile<<Atoms[atom].p.x-V2.x<<",\t"<<Atoms[atom].p.z-V2.z<<",\t"<<Atoms[atom].p.y-V2.y<<", 1e-6,\n";
			demoFile<<"4,\t";
			demoFile<<Atoms[atom].p.x<<",\t"<<Atoms[atom].p.z<<",\t"<<Atoms[atom].p.y<<",\t";
			demoFile<<Atoms[atom].p.x-V2.x<<",\t"<<Atoms[atom].p.z-V2.z<<",\t"<<Atoms[atom].p.y-V2.y<<", 1e-6,\n";
			demoFile<<"4,\t";
			demoFile<<Atoms[atom].p.x<<",\t"<<Atoms[atom].p.z<<",\t"<<Atoms[atom].p.y<<",\t";
			demoFile<<Atoms[atom].p.x-V1.x<<",\t"<<Atoms[atom].p.z-V1.z<<",\t"<<Atoms[atom].p.y-V1.y<<", 1e-6,\n";
		}
		for(auto delFace: Atoms[atom].DelaunayTessellations)
		{
			auto v1 = delFace.second.first->circumCenter;
			auto v2 = delFace.second.second->circumCenter;
			vect V1 = vectDiff(&Atoms[atom].p,&v1);
			vect V2 = vectDiff(&Atoms[atom].p,&v2);
			demoFile<<"2,\t";
			demoFile<<Atoms[atom].p.x-V1.x<<",\t"<<Atoms[atom].p.z-V1.z<<",\t"<<Atoms[atom].p.y-V1.y<<",\t";
			demoFile<<Atoms[atom].p.x-V2.x<<",\t"<<Atoms[atom].p.z-V2.z<<",\t"<<Atoms[atom].p.y-V2.y<<", 1e-6,\n";

			voroInfoOutputFile<<"draw line \t{";
			voroInfoOutputFile<<Atoms[atom].p.x-V1.x<<"\t"<<Atoms[atom].p.z-V1.z<<"\t"<<Atoms[atom].p.y-V1.y<<"}\t{";
			voroInfoOutputFile<<Atoms[atom].p.x-V2.x<<"\t"<<Atoms[atom].p.z-V2.z<<"\t"<<Atoms[atom].p.y-V2.y<<"} width 4\n";
		}
	}
	// Trim away Delaunay tetrahedra on the convex hull whose Voronoi vertex
	// (circumcenter) is not actually contained within the tetrahedron; this
	// iteratively peels the hull until only genuine void-space vertices
	// remain.
	std::vector<delunay *> voidDelToAnalyze =  fullSet.tetrahedrons;
	bool vertexDeleted=true;
	while(vertexDeleted)
	{
		vertexDeleted=false;
		for(auto D : voidDelToAnalyze)
		{
			if(D->hull and D->v->is_void)
			{
				if(!inside_delunay(D->v,D,Atoms,nAtoms))
				{
					for(auto neigh:D->neighDel)
					{
						neigh.second->hull=1;
						neigh.second->delNeighbor(D);
					}
					voidDelToAnalyze.erase(std::remove(voidDelToAnalyze.begin(), voidDelToAnalyze.end(), D), voidDelToAnalyze.end());
					vertexDeleted=true;
				}
			}
		}
	}
	int voidVertCount=0;
	std::ofstream voroVertexOutputFile(outPrefix + "voidVoronoiVertices.txt");
	for(auto D : voidDelToAnalyze)
	{
		if(D->neighDel.size() != 4 and !D->hull)
			std::cout<<"error\n";
		if(D->v->is_void)
		{
			voroVertexOutputFile<<"draw sphere\t{"<<D->circumCenter.x<<"\t"<<D->circumCenter.y<<"\t"<<D->circumCenter.z<<"}\tradius 1e-6\tresolution 10\n";
			voidVertCount++;
		}
	}
	delunay **cavity_list=nullptr;
	cavity_list = new (nothrow) delunay*[voidVertCount];
	int i=0;
	for(auto D : voidDelToAnalyze)
	{
		if(D->v->is_void)
		{
			cavity_list[i]=D;
			D->clusterIndex=i;
			i++;
		}
	}
	int change=1;
	int *old_label;
	old_label=new (nothrow) int[voidVertCount];
	int min;
	while(change)
	{
		for(int i=0; i<voidVertCount; i++)
		{
			old_label[i]=cavity_list[i]->clusterIndex;
		}
		for(int i=0; i<voidVertCount; i++)
		{
			min=cavity_list[i]->clusterIndex;
			for(auto neighB : cavity_list[i]->neighDel)
			{
				if(neighB.second->clusterIndex!=-1 && cavity_list[i]->isVoid[neighB.first])
				{
					if(min > neighB.second->clusterIndex)
					{
						min = neighB.second->clusterIndex;
					}
				}
			}
			for(auto neighB : cavity_list[i]->neighDel)
			{
				if(neighB.second->clusterIndex!=-1 && cavity_list[i]->isVoid[neighB.first])
				{
					neighB.second->clusterIndex=min;
				}
			}
		}
		change=0;
		int flag=0;
		for(int i=0; i<voidVertCount; i++)
		{
			flag=(old_label[i]!=cavity_list[i]->clusterIndex);
			if(flag)
			{
				change=1;
				break;
			}
		}
	}
	delete[] old_label;
	ofstream cav;
	cav.open(outPrefix + "cav"+std::to_string(r_cut)+".dat");
	ofstream cavitiesPOVRAY;
	cavitiesPOVRAY.open(outPrefix + "cavities"+std::to_string(r_cut)+".dat");
	cavitiesPOVRAY<<std::setprecision(16);
	for(int i=0; i<voidVertCount; i++)
	{
		int color=1;
		cav<<"#"<<i<<"\n\n";
		{
			for(int j=0; j<voidVertCount; j++)
			{
				if(cavity_list[j]->clusterIndex==i)
				{
					if(color)
					{
						cav<<"draw color "<<(i+1)%16<<"\n";;
						color=0;
					}
					cav<<"draw sphere\t{";
					cav<<cavity_list[j]->circumCenter.x<<"\t"<<cavity_list[j]->circumCenter.y<<"\t"<<cavity_list[j]->circumCenter.z<<"}\tradius\t"<<1e-3<<"\t"<<"resolution\t10\n";
					cavitiesPOVRAY<<"1,\t"<<cavity_list[j]->circumCenter.x<<",\t"<<cavity_list[j]->circumCenter.y<<",\t"<<cavity_list[j]->circumCenter.z<<",\n";
					for(auto neighB:cavity_list[j]->neighDel)
					{
						if(cavity_list[j]->isVoid[neighB.first] and neighB.second->v->is_void)
						{
							vect v12 = vectDiff(&cavity_list[j]->circumCenter,&neighB.second->circumCenter);
							cav<<"draw line { ";
							cav<<cavity_list[j]->circumCenter.x<<"\t"<<cavity_list[j]->circumCenter.y<<"\t"<<cavity_list[j]->circumCenter.z<<"}\t{";
							cav<<cavity_list[j]->circumCenter.x-v12.x<<"\t"<<cavity_list[j]->circumCenter.y-v12.y<<"\t"<<cavity_list[j]->circumCenter.z-v12.z<<"}\twidth 4\n";
							cavitiesPOVRAY<<"2,\t";
							cavitiesPOVRAY<<cavity_list[j]->circumCenter.x<<",\t"<<cavity_list[j]->circumCenter.y<<",\t"<<cavity_list[j]->circumCenter.z<<",\t";
							cavitiesPOVRAY<<cavity_list[j]->circumCenter.x-v12.x<<",\t"<<cavity_list[j]->circumCenter.y-v12.y<<",\t"<<cavity_list[j]->circumCenter.z-v12.z<<",\n";
						}
					}
				}
			}
		}
	}
	std::cout<<voidVertCount<<" number of vertices in void\n";
	long double *cav_vol = new (nothrow) long double[voidVertCount];
	long double *cav_area = new (nothrow) long double[voidVertCount];
	if(!cav_vol || !cav_area)
	{
		std::cerr<<"Error: failed to allocate memory for "<<voidVertCount<<" void vertices.\n";
		delete [] Atoms;
		return EXIT_FAILURE;
	}
	for(int i=0; i<voidVertCount; i++)
	{
		cav_vol[i]=0;
		cav_area[i]=0;
		for(int j=0; j<voidVertCount; j++)
		{
			if(cavity_list[j]->clusterIndex==i)
			{
				cav_vol[i]=cav_vol[i]+void_vol(cavity_list[j]->v,Atoms);
			}
		}
	}
	long double cav_tot=0.;
	long double ca_per_tot=0.;
	// cluster_index <TAB> void_volume for every non-empty cavity/pocket.
	std::ofstream cavityVolumeOutput(outPrefix + "cavityVolume.txt");
	for(int i=0; i<voidVertCount; i++)
	{
		if(cav_vol[i])
			cavityVolumeOutput<<i<<"\t"<<cav_vol[i]<<"\n";
		cav_tot=cav_tot+cav_vol[i];
		ca_per_tot=ca_per_tot+cav_area[i];
	}
	std::cout<<cav_tot<<"\t"<<cav_tot/(boxx*boxy*boxz)<<"\n";

	// Appended (not overwritten) so this file can accumulate a
	// probe-radius-vs-void-volume curve across multiple runs with
	// different probe_radius values on the same packing.
	std::ofstream outfile(outPrefix + "r_cut_vs_pore_vol.txt", std::ios::app);
	outfile<<std::setprecision(16)<<r_cut<<"\t"<<cav_tot<<"\n";

	std::cout<<std::setprecision(16);
	std::cout<<"probe_radius\tbox_volume\tatoms_volume+void_volume\tvoid_volume\n";
	std::cout<<r_cut<<"\t"<<volumeBox<<"\t"<<volumeAtoms+cav_tot<<"\t"<<cav_tot<<"\n";

	delete [] cav_vol;
	delete [] cav_area;
	delete [] cavity_list;
	delete [] Atoms;
	return EXIT_SUCCESS;
}
