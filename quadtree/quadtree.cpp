#ifdef __APPLE__
#include <GLUT/glut.h>  // GLUT, include glu.h and gl.h
#include <stdlib.h>
#endif
#ifdef __linux__
#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#endif
#include <bitmap.h>

#include <assert.h>
#include <algorithm>
#include <ctime>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

// convergence stuff

bool doneVCycle = false;

int MAX_RELAX_COUNT = 20;
int OPTIMIZED_RELAX_COUNT = 4;

double eps = .0001;

// options
int windowid;
int windowWidth = 640;
int windowHeight = 640;

void statScreenshotDir() {
	 // meh
}


int frameNumber = 0;

void saveScreen() {
    unsigned char * imageBuffer = new unsigned char[3 * windowWidth * windowHeight];
    glReadPixels( 0, 0, windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_BYTE, imageBuffer );

    char * fname = new char[80];
    snprintf(fname, 80, "screenshots/image%d.bmp", frameNumber);

    writeBMP(fname, windowWidth, windowHeight, imageBuffer);

    delete imageBuffer;
    delete fname;

    printf("saved screenshot of frame %d\n", frameNumber);
}

bool headless; // run without ui for profiling
int levels;
int levelToDisplay; // level to draw
bool water;
bool drawVelocity;
bool drawCells;
bool screenshot;
int numToRun = 0;
int warmupRuns = 0;
bool pressureInterp = false;

// profiling
int numLeaves, numNodes;
int* leafLevels;
timeval start, end;


void resetProfiling() {
	if (leafLevels) delete leafLevels;
	leafLevels = new int[levels];
	memset(leafLevels, 0, levels * sizeof(int));
	numNodes = 0;
	numLeaves = 0;
}

void printProfiling() {
	printf("number of leaves at each level: \n");
	for (int i = 0; i < levels; i++) {
		printf("%d: %d\n", i, leafLevels[i]);
	}
	printf("total number of nodes: %d\n", numNodes);
	printf("total number of leaves: %d\n", numLeaves);
	printf("total number of non-max leaves: %d\n", numLeaves - leafLevels[levels - 1]);
}

void startTime() {
	gettimeofday(&start, NULL);
}
double endTime(char* task) {
	gettimeofday(&end, NULL);
	double msdif = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
	printf("%s took time(ms): %.3f\n", task, msdif);
	return msdif;
}



// options for particle rendering
enum ParticleAlgorithm { EULER, PARTICLENONE };
ParticleAlgorithm particleAlgorithm = PARTICLENONE;

struct Particle {
	double x, y;
};

int numParticles;
Particle* particles;

// start state stuff
enum StartState { LEFT, SINK, SRCSINK, POISSONTEST, ADAPTTEST , ERRORTEST, PROJECTTEST };
StartState startState;
enum AdaptScheme { ADAPTNONE, CURL, PGRAD, VNORM };
AdaptScheme adaptScheme = ADAPTNONE;
enum AdaptTestFunc { ADAPTFUNCNONE, ADAPTCIRCLE, ADAPTWAVE };
AdaptTestFunc adaptTestFunc = ADAPTFUNCNONE;
enum PoissonTestFunc { POISSONFUNCNONE, POISSONXY, POISSONCOS, POISSONCOSKL };
PoissonTestFunc poissonTestFunc = POISSONFUNCNONE;
enum ErrorTestFunc { ERRORFUNCNONE, ERRORXY, ERRORCOS, ERRORCOSKL };
ErrorTestFunc errorTestFunc = ERRORFUNCNONE;
enum ProjectTestFunc { PROJECTFUNCNONE, PROJECTXY, PROJECTSIN/*, TODO real h-h function*/ };
ProjectTestFunc projectTestFunc = PROJECTFUNCNONE;

double dt = .03;

// drawing stuff
double minP;
double maxP;
double maxMag;


// navigating 
double deltas[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

// returns velocity with each term bilinearly interpolated from nodal values
/*
 a------b
 |  |   |
 |  |di |
 |--    |
 | dj   |
 |      |
 c------d
 */
double bilinearInterpolation(double a, double b, double c, double d, double di, double dj) {
	return a * (1-di) * (1-dj) + b * (1-di) * dj + c * di * (1-dj) + d * di * dj;
}


int nodeId = 0;
enum NodeValue { P, VX, VY, PHI, DP, DIVV };
class qNode {
	public:
		// tree stuffs
		int id;
		bool leaf, shouldExpand, shouldContract;
		qNode* parent;
		qNode* children[4]; // 2D array in 1D array
		qNode* neighbors[4];
		int level, i, j; // level coordinates

		// math stuffs
		double p, dp, R, phi, divV, temp; // pressure ,x velocity, y velocity, pressure correction, residual, level set value
		double vx, vy, vx2, vy2;
		bool cvxValid[4];
		double cvx[4];
		bool cvyValid[4];
		double cvy[4];

		qNode(qNode *p, int i, int j): parent(p), i(i), j(j) {
			if (parent == NULL) {
				level = 0;
			} else {
				level = 1 + parent->level;
			}
			id = nodeId++;
			leaf = true;
			for (int i = 0; i < 4; i++) {
				children[i] = NULL;
				neighbors[i] = NULL;
			}
		}
		~qNode() {
			// set your neighbors' neighbor pointer at you to be null
			
			for (int i = 0; i < 2; i++) {
				// delete x neighbor
				if (neighbors[3-i] != NULL) {
					neighbors[3-i]->neighbors[2+i] = NULL;
				}
				// delete y neighbor
				if (neighbors[1-i] != NULL) {
					neighbors[1-i]->neighbors[i] = NULL;
				}
			}
			if (!leaf) {
				for (int k = 0; k < 4; k++) {
					delete children[k];
				}
			}
		}

		// moving around tree functions
		void setChildrenNeighbors() {
			// set both ways because the first time a node is expanded the other one may not have been
			// children(i, j)
			// (0, 0) (0, 1)
			// (1, 0) (1, 1)
			// deltas[k]= = {{1, 0}^, {-1, 0}v, {0, 1}>, {0, -1}<}
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					qNode* c = children[i * 2 + j];
					// outside neighbor in x
					c->neighbors[3-j] = NULL;
					if (this->neighbors[3-j] != NULL
							&& !this->neighbors[3-j]->leaf) {
						c->neighbors[3-j] = this->neighbors[3-j]->children[2*i+(1-j)]; 
						this->neighbors[3-j]->children[2*i+(1-j)]->neighbors[2+j] = c; 
					}
					// inside neighbor in x
					c->neighbors[2+j] = this->children[2*i+(1-j)];
					this->children[2*i+(1-j)]->neighbors[3-j] = c;

					// outside neighbor in y
					c->neighbors[1-i] = NULL;
					if (this->neighbors[1-i] != NULL && !this->neighbors[1-i]->leaf) {
						c->neighbors[1-i] = this->neighbors[1-i]->children[(1-i)*2+j]; 
						this->neighbors[1-i]->children[(1-i)*2+j]->neighbors[i] = c; 
					}
					// inside neighbor in y
					c->neighbors[i] = this->children[(1-i)*2+j];
					this->children[(1-i)*2+j]->neighbors[1-i] = c;
				}
			}
		}

		double getVal(NodeValue v) {
			if (v == DP) return dp;
			else if (v == P) return p;
			else if (v == VX) return vx;
			else if (v == VY) return vy;
			else if (v == DIVV) return divV;
			assert(false);
		}

		double getValueInterpCorrection(qNode* original, int ml, int k, NodeValue v) {
			int newk;
			double dpos;
			if (k < 2) {
				// previously computing gradient in y direction, need x direction for interpolation
				double origX = (original->j + 0.5)/(1<<original->level);
				double x = (j + 0.5)/(1<<level);
				dpos = origX - x;
				newk = dpos > 0 ? 2 : 3;
			} else {
				double origY = (original->i + 0.5)/(1<<original->level); 
				double y = (i + 0.5)/(1<<level);
				dpos = origY - y;
				newk = dpos > 0 ? 0 : 1;
			}
			double faceGrad = getFaceGradient(ml, newk, v);
			return fabs(dpos) * faceGrad;
		}

		double getFaceGradient(int ml, int k, NodeValue v) {
			qNode* n = this;
			int oppositeK = (k < 2) ? 1-k : 3-(k%2);//1 - (k%2);
			while (n != NULL && n->neighbors[k] == NULL) n = n->parent;
			if (n != NULL) {
				return n->neighbors[k]->addFaceToGradient(this, ml, oppositeK, v);
			} else {
				return 0;
			}
		}

		double addFaceToGradient(qNode* original, int ml, int k, NodeValue v) {
			if (leaf || level == ml) {
				int origSize = 1<<original->level;
				int size = 1<<level;
				double dside = fmin(1.0/origSize, 1.0/size);
				double h = 0.5/origSize + 0.5/size;
				double d = (dside * origSize)/h;
				double val = getVal(v);
				if (pressureInterp && original->level > level) {
					val += getValueInterpCorrection(original, ml, k, v);
				}
				
				return (val - original->getVal(v)) * d;
			} else {
				double total = 0.0;
				if (k < 2) { // y dir
					int targetI = 1 - k;
					for (int j = 0; j < 2; j++) {
						total += children[targetI*2+j]->addFaceToGradient(original, ml, k, v);
					}
				} else {
					int targetJ = 3 - k;
					for (int i = 0; i < 2; i++) {
						total += children[2*i+targetJ]->addFaceToGradient(original, ml, k, v);
					}
				}
				return total;
			}

		}

		void getLaplacian(int ml, double *aSum, double* bSum, NodeValue v) {
			for (int k = 0; k < 4; k++) {
				int oppositeK = (k < 2) ? 1-k : 3-(k%2); //1 - (k%2);
				qNode* n = this;
				while (n != NULL && n->neighbors[k] == NULL) n = n->parent;
				if (n != NULL) {
					n->neighbors[k]->addFaceToLaplacian(this, ml, oppositeK, aSum, bSum, v);
				} else {
					*aSum -= 1;
					*bSum += getVal(v);
				}
			}
		}

		// note: k is opposite of way used to get there, for example if k in original was (1, 0), this k is (-1, 0)
		void addFaceToLaplacian(qNode* original, int ml, int k, double* aSum, double* bSum, NodeValue v) {
			if (leaf || level == ml) {
				int origSize = 1<<original->level;
				int size = 1<<level;
				double dside = fmin(1.0/origSize, 1.0/size);
				double h = 0.5/origSize + 0.5/size;
				double d = dside/h;
				double val = getVal(v);
				if (pressureInterp && original->level > level) {
					val += getValueInterpCorrection(original, ml, k, v);
				}


				*aSum -= d;
				*bSum += d*val;
			} else {
				if (k < 2) { // y dir
					int targetI = 1 - k;
					for (int j = 0; j < 2; j++) {
						children[targetI*2+j]->addFaceToLaplacian(original, ml, k, aSum, bSum, v);
					}
				} else {
					int targetJ = 3 - k;
					for (int i = 0; i < 2; i++) {
						children[2*i+targetJ]->addFaceToLaplacian(original, ml, k, aSum, bSum, v);
					}
				}
			}
		}

		// utility functions
		std::pair<double, double> getValueGradient(NodeValue v) {
			
			int ml = leaf ? levels - 1 : level;
			std::pair<double, double> newGrad =  std::make_pair((getFaceGradient(ml, 2, v) - getFaceGradient(ml, 3, v))/2.0, (getFaceGradient(ml, 0, v) - getFaceGradient(ml, 1, v))/2.0);

			return newGrad;
		}
		
		// curl(F) = d(Fy)/dx - d(Fx)/dy
		double getCurl() {
			assert (false);
			 //TODO fix
			std::pair<double, double> xgrad = getValueGradient(VX);
			std::pair<double, double> ygrad = getValueGradient(VY);
			return ygrad.first - xgrad.second;
		}

		void computeVelocityDivergence() {
			//std::pair<double, double> grad = getVelocityGradient();
			//divV = grad.first + grad.second;
			if (leaf) {
				int size = 1<<level;
				double a, b;

				// x
				if (j == 0)
					a = 0.0;
				else
					a = vx;
				if (j == size - 1)
					b = 0.0;
				else {
					b = vx2;
				}
				divV = b-a;

				// y
				if (i == 0)
					a = 0.0;
				else
					a = vy;
				if (i == size-1)
					b = 0.0;
				else {
					b = vy2;
				}
				divV += b-a;

				divV *= size;
			} else {
				divV = 0.0;
				for (int k = 0; k < 4; k++) {
					children[k]->computeVelocityDivergence();
				}
			}
		}

		std::pair<double, double> getVelocityAt(double x, double y) {
			int size = 1<<level;
			double minX = ((float)j)/size;
			double minY = ((float)i)/size;
			assert (!(x < minX || y < minY || x > minX + 1.0/size || y > minY + 1.0/size));
			double dj = (x*size)-j;
			double di = (y*size)-i;

			double newvx = bilinearInterpolation(cvx[0], cvx[1], cvx[2], cvx[3], di, dj);
			double newvy = bilinearInterpolation(cvy[0], cvy[1], cvy[2], cvy[3], di, dj);	
			return std::make_pair(newvx, newvy);

		}
		
		double getValueAt(double x, double y, NodeValue v) {
			double val;
			if (v == P)
				val = p;
			else if (v == PHI)
				val = phi;
			else if (v == DP)
				val = dp;
    		std::pair<double, double> grad = getValueGradient(v);
			int size = 1<<level;
			double dx = x - (j + 0.5)/size;
			double dy = y - (i + 0.5)/size;
			if (fabs(dx) > 0.5/size || fabs(dy) > 0.5/size) {
				printf("dx: %f, dy: %f\n", dx, dy);
			}
			return val + dx * grad.first + dy * grad.second;
		}

		// adaptivity
		void expand(bool calculateNewVals) {
			assert(leaf);
    		int size = 1<<level;
    		for (int k = 0; k < 4; k++) {
				this->children[k] = new qNode(this, 2*i+(k/2), 2*j+(k%2));
			}
			if (calculateNewVals) {

				// p will be changed later anyway


				for (int c = 0; c < 4; c++) {
					children[c]->cvx[c] = cvx[c];
					children[c]->cvy[c] = cvy[c];
					children[c]->p = p;
				}
				
				double newv;
				// neg x node, compute y
				newv = (cvy[0] + cvy[2])/2.0;
				children[0]->cvx[2] = vx;
				children[0]->cvy[2] = newv;
				children[2]->cvx[0] = vx;
				children[2]->cvy[0] = newv;

				// pos x node, compute y
				newv = (cvy[1] + cvy[3])/2.0;
				children[1]->cvx[3] = vx2;
				children[1]->cvy[3] = newv;
				children[3]->cvx[1] = vx2;
				children[3]->cvy[1] = newv;

				// neg y node, compute x
				newv = (cvx[0] + cvx[1])/2.0;
				children[0]->cvx[1] = newv;
				children[0]->cvy[1] = vy;
				children[1]->cvx[0] = newv;
				children[1]->cvy[0] = vy;

				// neg y node, compute x
				newv = (cvx[2] + cvx[3])/2.0;
				children[2]->cvx[3] = newv;
				children[2]->cvy[3] = vy2;
				children[3]->cvx[2] = newv;
				children[3]->cvy[2] = vy2;

				// center
				newv = (vx + vx2)/2.0;
				children[0]->cvx[3] = newv;
				children[1]->cvx[2] = newv;
				children[2]->cvx[1] = newv;
				children[3]->cvx[0] = newv;

				newv = (vy + vy2)/2.0;
				children[0]->cvy[3] = newv;
				children[1]->cvy[2] = newv;
				children[2]->cvy[1] = newv;
				children[3]->cvy[0] = newv;

				// compute new face velocities from nodes
				for (int k = 0; k < 4; k++) {
					children[k]->vx = (cvx[0] + cvx[2])/2.0;
					children[k]->vx2 = (cvx[1] + cvx[3])/2.0;
					children[k]->vy = (cvy[0] + cvy[1])/2.0;
					children[k]->vy2 = (cvy[2] + cvy[3])/2.0;
				}

			}
    		
			setChildrenNeighbors();
			leaf = false;
		}

		// x is which axis edge is on (so opposite velocity is needed)
		void getNodalSide(int ni, int nj, int* retL, double* retV, bool x, bool left, int targetLevel) {
			int leveldif = targetLevel - level;
			int middist = 1<<(leveldif-1);
			int midi = i * (1<<leveldif) + middist;
			int midj = j * (1<<leveldif) + middist;
			
			if (leaf) {
				if (*retL > level) return;
				if (!((ni == midi - middist || ni == midi + middist) && (nj == midj - middist || nj == midj + middist))) return;
				if (x) {
					if (nj < midj && left) return;
					if (nj > midj && !left) return;
					*retL = level;
					if (ni < midi)
						*retV = vy;
					else
						*retV = vy2;
				} else {
					if (ni < midi && left) return;
					if (ni > midi && !left) return;
					*retL = level;
					if (nj < midj)
						*retV = vx;
					else
						*retV = vx2;
				}
				return;
			}
			if (ni != midi && nj != midj) {
				int newi = (ni < midi) ? 0 : 1;
				int newj = (nj < midj) ? 0 : 1;
				return children[2*newi+newj]->getNodalSide(ni, nj, retL, retV, x, left, targetLevel);
			} else if (ni != midi) {
				// j is on mid
				assert (!x);
				if (left) {
					for (int k = 2; k < 4; k++)
						children[k]->getNodalSide(ni, nj, retL, retV, x, left, targetLevel);	
				} else {
					for (int k = 0; k < 2; k++)
						children[k]->getNodalSide(ni, nj, retL, retV, x, left, targetLevel);		
				}
			} else if (nj != midj) {
				assert (x);
				if (left) {
					for (int k = 1; k < 4; k+=2)
						children[k]->getNodalSide(ni, nj, retL, retV, x, left, targetLevel);	
				} else {
					for (int k = 0; k < 4; k+=2)
						children[k]->getNodalSide(ni, nj, retL, retV, x, left, targetLevel);		
				}

			} else {
				for (int k = 0; k < 4; k++)
					children[k]->getNodalSide(ni, nj, retL, retV, x, left, targetLevel);
			}
		}

		bool getArtificialValue(int ni, int nj, int nl, bool x, bool left, double* retV) {
			if (leaf) {
				// compute value here

				int leveldif = nl - level;
				double delta = .5/(1<<leveldif);
				double va, vb, vc, vd;
				double di, dj;
				di = (((double)i)/level - ((double)ni)/nl) * (1<<level);
				dj = (((double)j)/level - ((double)nj)/nl) * (1<<level);
				if (x) {
					if (left) {
						assert (dj > .5);
						// right side of cell
						if (!cvyValid[1] || !cvyValid[3]) return false;
						va = vy;
						vb = cvy[1];
						vc = vy2;
						vd = cvy[3];
						dj -= delta;
						dj -= .5;
					} else {
						assert (dj < .5);
						// left side of cell
						if (!cvyValid[0] || !cvyValid[2]) return false;
						va = cvy[0];
						vb = vy;
						vc = cvy[2];
						vd = vy2;
						dj += delta;
					}
					dj *= 2.0;
				} else {
					if (left) {
						assert (di > .5);
						// bottom side of cell
						if (!cvxValid[2] || !cvxValid[3]) return false;
						va = vx;
						vb = vx2;
						vc = cvx[2];
						vd = cvx[3];
						di -= delta;
						di -= .5;
					} else {
						// top side of cell
						assert (di < .5);
						if (!cvxValid[0] || !cvxValid[1]) return false;
						va = cvx[0];
						vb = cvx[1];
						vc = vx;
						vd = vx2;
						di += delta;
					}
					di *= 2.0;
				}
				*retV = bilinearInterpolation(va, vb, vc, vd, di, dj);
			} else {
				int leveldif = nl - level;
				int midi = i * (1<<leveldif) + (1<<(leveldif-1));
				int midj = j * (1<<leveldif) + (1<<(leveldif-1));
				assert (ni != midi && nj != midj); // otherwise we wouldn't be artificially discretizing
				int nexti = (i < midi) ? 0 : 1;
				int nextj = (j < midj) ? 0 : 1;
				return children[nexti*2+nextj]->getArtificialValue(ni, nj, nl, x, left, retV);
			}
		}

		// attempts to compute the given nodal value. If successful, it returns true with the value in ret, otherwise it returns false.
		bool attemptComputeNodalValue(qNode* r, int c, bool x, double* ret) {
			int ci = (c/2);
			int cj = (c%2);
			int ni = i + ci;
			int nj = j + cj;
			int size = 1<<level;
			// handle edge cases differently
			if (x) {
				if (nj == 0 || nj == size) {
					*ret = 0;
					return true;
				}
				if (ni == 0 || ni == size) {
					*ret = (cj == 0) ? vx : vx2;
					return true;	
				}
			} else {
				if (ni == 0 || ni == size) {
					*ret == 0;
					return true;	
				}
				if (nj == 0 || nj == size) {
					*ret = (ci == 0) ? vy : vy2;
					return true;
				}

			}
			int levelL, levelR;
			levelL = -1;
			levelR = -1;
			double valL, valR;
			bool haveL, haveR;
			int leveldif = levels - 1 - level;
			r->getNodalSide(ni * (1<<leveldif), nj * (1<<leveldif), &levelL, &valL, x, true, levels - 1);
			haveL = levelL != -1;
			r->getNodalSide(ni * (1<<leveldif), nj * (1<<leveldif), &levelR, &valR, x, false, levels - 1);
			haveR = levelR != -1;
			bool left = false;
			if (haveR && haveL) {
				double lenL = .5/(1<<levelL);
				double lenR = .5/(1<<levelR);
				double totalLen = lenL + lenR;
				*ret = lenR*valL/totalLen + lenL*valR/totalLen;
				return true;
			} else if (haveR) {
				valL = valR;
				levelL = valR;
				left = true;
				haveL = true;
			}
			assert (haveL);
			assert (levelL >= level);
			if (r->getArtificialValue(ni, nj, level, x, left, &valR)) {
				*ret = (valR + valL)/2.0;
				return true;
			}
			return false;
		}

		// returns the number of nodal values left to compute
		int attemptComputeNodalVals(qNode* r) {
			int total = 0;
			if (leaf) {
				// do stuff

				double v;
				for (int c = 0; c < 4; c++) {
					if (!cvxValid[c] && attemptComputeNodalValue(r, c, true, &v)) {
						cvxValid[c] = true;
						cvx[c] = v;
					}
					if (!cvyValid[c] && attemptComputeNodalValue(r, c, false, &v)) {
						cvyValid[c] = true;
						cvy[c] = v;
					}
				}

				// return stuff				
				for (int c = 0; c < 4; c++) {
					if (!cvxValid[c]) total++;
					if (!cvyValid[c]) total++;
				}
				return total;
			} else {
				for (int k = 0; k < 4; k++) {
					total += children[k]->attemptComputeNodalVals(r);
					if (!cvxValid[k] && children[k]->cvxValid[k]) {
						cvx[k] = children[k]->cvx[k];
						cvxValid[k] = true;
					}
					if (!cvyValid[k] && children[k]->cvyValid[k]) {
						cvy[k] = children[k]->cvy[k];
						cvyValid[k] = true;
					}
				}
				return total;
			}
		}

		void computeVelFromChildren() {
			assert (!leaf);
			vx = (children[0]->vx + children[2]->vx)/2.0;
			vx2 = (children[1]->vx2 + children[3]->vx2)/2.0;
			vy = (children[0]->vy + children[1]->vy)/2.0;
			vy2 = (children[2]->vy2 + children[3]->vy2)/2.0;
		}

		// don't need to average values since adapt function takes care of that
		void contract(bool calculateNewVals) {
			assert(!leaf);
			for (int k = 0; k < 4; k++) {
				assert(children[k]->leaf);
			}
			if (calculateNewVals) {
				computeVelFromChildren();
				p = 0;
				for (int k = 0; k < 4; k++) {
					p += children[k]->p;
				}
				p /= 4.0;
			}
		    //printf("contracting cell %d (%d, %d)\n", level, i, j);
			for (int k = 0; k < 4; k++) {
				delete children[k];
			}
			leaf = true;
		}
		void profile() {
			numNodes++;
			if (leaf) {
				numLeaves++;
				leafLevels[level]++;
			} else {
				for (int k = 0; k < 4; k++) {
					children[k]->profile();
				}
			}
		}

		// current node, target d/i/j
		
		qNode* get(int nd, int ni, int nj) {
			if (leaf) return this;
			int leveldif = nd - level;
			int midsize = 1<<(leveldif-1);
			int midi = i*(1<<leveldif) + midsize;
			int midj = j*(1<<leveldif) + midsize;
			int newi = (ni < midi) ? 0 : 1;
			int newj = (nj < midj) ? 0 : 1;
			return children[newi*2+newj]->get(nd, ni, nj);
		}	
};


qNode* root;
qNode* oldRoot;


// debugging

qNode* getLeaf(qNode* node, double x, double y, int ml) {
	int size = 1<<node->level;
	assert (x >= (node->j + 0.0)/size && x <= (node->j + 1.0)/size);
	assert (y >= (node->i + 0.0)/size && y <= (node->i + 1.0)/size);
	/*if (x*size < node->j || x*size > node->j + 1 || y*size < node->i || y*size > node->i + 1) {
		printf("out of range!\n");
	}*/
	if (node->leaf || node->level == ml) {
		return node;
	}
	double midx = (node->j + 0.5)/size;
	double midy = (node->i + 0.5)/size;
	int newj = x < midx ? 0 : 1;
	int newi = y < midy ? 0 : 1;
	return getLeaf(node->children[2*newi + newj], x, y, ml);
}

qNode*  getSemiLagrangianLookback(qNode* r, double* x, double* y, int steps, int ml, double vx, double vy) {
	double newdt = dt / steps;
	qNode* cur;
	while (steps--) {	
		*x -= vx * newdt;
		*y -= vy * newdt;
		*x = fmin(1.0, fmax(0, *x));
		*y = fmin(1.0, fmax(0, *y));
		cur = getLeaf(r, *x, *y, ml);
		std::pair<double, double> vel = cur->getVelocityAt(*x, *y);
		vx = vel.first;
		vy = vel.second;
	}
	return cur;
}

// not fast but who cares
void printValue(NodeValue v) {
	int size = 1<<(levels - 1);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			qNode* n = root->get(levels - 1, i, j);
			if (n != NULL) {
				printf(" %.4f", n->getVal(v));
			} else {
				printf("      ");
			}
		}
		printf("\n");
	}
}

/* Initialize OpenGL Graphics */
void initGL() {
   // Set "clearing" or background color
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
}

void drawMultilevel(qNode* node, int ml) {
	assert (node->level <= ml);
	int size = 1 << node->level;
	if (node->level == ml || node->leaf) {
		// draw it
		double x = -1.0 + (2.0 * node->j)/size;
		double y = -1.0 + (2.0 * node->i)/size;

		if (water) {
			if (node->phi >= 0) {
				double val = fmin(1.0, node->phi);
				glColor3f(0.0, 0.0, val);
			}
		} else {
			double redPercent = 0.0;
			double bluePercent = 0.0;
			if (node->p > 0 && maxP > 0) {
				redPercent = node->p/maxP;
			}
			if (node->p < 0 && minP < 0) {
				bluePercent = node->p/minP;
			}
			
			glColor3f(redPercent, 0, bluePercent);
		}

		glRectf(x, y, x + 2.0/size, y + 2.0/size);

		if (drawVelocity && node->level <= 8) {
			if (maxMag >= eps) {
				glColor3f(0.0, 1.0, 0.0);
				double vx = node->vx;
				double vy = node->vy;
				double mag = sqrt(vx*vx + vy*vy);
				if (mag >= eps) {
					glBegin(GL_LINES);

					double scale  = maxMag  * size;
					// vx
					glVertex2f(x, y + 1.0/size); 
					glVertex2f(x + vx / scale, y + 1.0/size);

					//vy
					glVertex2f(x + 1.0/size, y);
					glVertex2f(x + 1.0/size, y + vy / scale);

					glEnd();

				}
			}
		}
		// draw cell
		if (drawCells && node->level <= 8) {
			glColor3f(1.0, 1.0, 1.0);
			glBegin(GL_LINE_LOOP);
			glVertex2f(x, y);
			glVertex2f(x, y + 2.0/size);
			glVertex2f(x + 2.0/size, y + 2.0/size);
			glVertex2f(x + 2.0/size, y);
			glEnd();
		}
	} else {
		// draw children
		for (int k = 0; k < 4; k++) {
			drawMultilevel(node->children[k], ml);
		}
	}
}

void computeMax(qNode* node, int ml) {
	if (node->leaf || node->level == ml) {
		if (!water) {
			minP = std::min(minP, node->p);
			maxP = std::max(maxP, node->p);
		}

		if (drawVelocity && ml <= 5) {
			maxMag = std::max(maxMag, sqrt(node->vx*node->vx + node->vy*node->vy));
		}
	} else {
		for (int k = 0; k < 4; k++) {
			computeMax(node->children[k], ml);
		}
	}
}

void drawParticles() {
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4f(0.0, 1.0, 0.0, 0.5);
	glPointSize(10);
	glBegin(GL_POINTS);
	for (int i = 0; i < numParticles; i++) {
		glVertex2f(-1.0 + 2*particles[i].x, -1.0 + 2*particles[i].y);
	}
	glEnd();
}

double adaptScale = 0.6;
double wavePeriod = 2;
double offsetY = 0.0;

bool pointInCircle(double circleCenterX, double circleCenterY, double circleRadius, double x, double y) {
	return (circleCenterX-x)*(circleCenterX-x)+(circleCenterY-y)*(circleCenterY-y) <= (circleRadius * circleRadius);
}

bool lineIntersectCircle(double circleCenterX, double circleCenterY, double circleRadius, double x1, double y1, double x2, double y2) {
	if (pointInCircle(circleCenterX, circleCenterY, circleRadius, x1, y1) && pointInCircle(circleCenterX, circleCenterY, circleRadius, x2, y2)) return false;
	if (pointInCircle(circleCenterX, circleCenterY, circleRadius, x1, y1) ^ pointInCircle(circleCenterX, circleCenterY, circleRadius, x2, y2)) return true;
	
	// build vector a from point to radius, project it onto current vector b
	double bx = x2-x1;
	double by = y2-y1;
	double ax = circleCenterX - x1;
	double ay = circleCenterY - y1;
	
	double adotb = ax*bx + ay*by;
	double magA = sqrt(ax*ax + ay*ay);
	double magB = sqrt(bx*bx + by*by);
	double a1 = adotb/magB;
	
	if (a1 >= 0 && a1 <= magB) {
		double ax1 = x1 + a1*bx/magB;
		double ay1 = y1 + a1*by/magB;
	
		double cx = ax1 - circleCenterX;
		double cy = ay1 - circleCenterY;
		
		return (cx*cx + cy*cy) <= (circleRadius * circleRadius);
	} else {
		bx = circleCenterX - x2;
		by = circleCenterY - y2;
		return (ax * ax + ay*ay) <= (circleRadius*circleRadius) || (bx*bx+by*by) <= (circleRadius * circleRadius);
	}
}

bool rectIntersectFuncCircle(double x1, double y1, double x2, double y2) {
	double circleCenterX = 0.5;
	double circleCenterY = 0.5 + offsetY;
	double circleRadius = .5*adaptScale;

	double xmin = fmin(x1, x2);
	double xmax = fmax(x1, x2);
	double ymin = fmin(y1, y2);
	double ymax = fmax(y1, y2);

	
	if(lineIntersectCircle(circleCenterX, circleCenterY, circleRadius, xmin, ymin, xmax, ymin)) return true;
	if(lineIntersectCircle(circleCenterX, circleCenterY, circleRadius, xmin, ymin, xmin, ymax)) return true;
	if(lineIntersectCircle(circleCenterX, circleCenterY, circleRadius, xmax, ymax, xmax, ymin)) return true;
	if(lineIntersectCircle(circleCenterX, circleCenterY, circleRadius, xmax, ymax, xmin, ymax)) return true;
	return false;
}

double getWaveFunc(double x) {
	return .5 + adaptScale*.5*cos(wavePeriod*2*M_PI*x) + offsetY;
}

bool waveIntersectHoriz(double x1, double x2, double y) {
	// break into parts by period
	int numParts = wavePeriod * 2;
	double periodWidth = 1.0/numParts;
	for (int i = 0; i < numParts; i++) {
		double periodStartX = i*periodWidth;
		double periodEndX = periodStartX + periodWidth;
		if (periodEndX < x1 || periodStartX > x2) continue;
		periodStartX = fmax(periodStartX, x1);
		periodEndX = fmin(periodEndX, x2);
	
		double periodStartY = getWaveFunc(periodStartX);
		double periodEndY = getWaveFunc(periodEndX);

		double ymin = fmin(periodStartY, periodEndY);
		double ymax = fmax(periodStartY, periodEndY);
		if (ymin <= y && y <= ymax) return true;
	}
	return false;
}

bool rectIntersectFuncWave(double x1, double y1, double x2, double y2) {
	double ymin = fmin(y1, y2);
	double ymax = fmax(y1, y2);

	double fy = getWaveFunc(x1);
	if (ymin <= fy && fy <= ymax) return true;
	fy = getWaveFunc(x2);
	if (ymin <= fy && fy <= ymax) return true;

	if (waveIntersectHoriz(x1, x2, ymin)) return true;
	if (waveIntersectHoriz(x1, x2, ymax)) return true;
}

bool rectIntersectFunc(double x1, double y1, double x2, double y2) {
	if (adaptTestFunc == ADAPTCIRCLE) {
		return rectIntersectFuncCircle(x1, y1, x2, y2);
	} else if (adaptTestFunc == ADAPTWAVE) {
		return rectIntersectFuncWave(x1, y1, x2, y2);
	}
	assert (false);
}

void testIntersect() {
	assert (rectIntersectFuncCircle(0.0, 0.0, .25, .25) == false);
	assert (rectIntersectFuncWave(0, .5, .25, .75));
	//rectIntersectFuncWave(.5, 0, 1, .5);
}

void drawWaveSubdivide(double x1, double x2) {
	double y1 = getWaveFunc(x1);
	double y2 = getWaveFunc(x2);
	double slope = (y2-y1)/(x2-x1);
	if (x2-x1 < .01) {
		glBegin(GL_LINES);
		glVertex2f(-1+2*x1, -1+2*y1);
		glVertex2f(-1+2*x2, -1+2*y2);
		glEnd();
	} else {
		double median = (x1+x2)/2.0;
		drawWaveSubdivide(x1, median);
		drawWaveSubdivide(median, x2);
	}
}

void drawAdaptFunction() {
	glColor3f(0.0, 1.0, 0.0);
	if (adaptTestFunc == ADAPTCIRCLE) {
		// draw circle
		double radius = adaptScale;
		glBegin(GL_LINE_LOOP);
		for (int i = 0; i < 360; i++) {
			glVertex2f(radius * cos(i*M_PI/180), radius * sin(i*M_PI/180) +2*offsetY);
		}
		glEnd();
	} else {
		drawWaveSubdivide(0.0, 1.0);
	}
}

/* Handler for window-repaint event. Call back when the window first appears and
   whenever the window needs to be re-painted. */
void display() {
	glClear(GL_COLOR_BUFFER_BIT);   // Clear the color buffer with current clearing color

	printf("display\n");

    if (particleAlgorithm == PARTICLENONE) {
		minP = 100000000.0;
		maxP = -100000000.0;
		maxMag = 0.0;

		computeMax(root, levelToDisplay);

		drawMultilevel(root, levelToDisplay);
    } else {
		drawParticles();
	}

	if (adaptTestFunc != ADAPTFUNCNONE) {
		drawAdaptFunction();
	}
	
	glFlush();  // Render now

	// capture screen if necessary
	if (screenshot) {
		saveScreen();
	}
}

//tests
bool assertNeighbors(qNode* n1, qNode* n2, qNode* realN1, qNode* realN2) {
	printf("asserting neighbors: \n");
	printf("n1: [%d][%d][%d], ", n1->level, n1->i, n1->j);
	if (n2 == NULL) printf("n2: NULL\n");
	else printf("n2: [%d][%d][%d]\n", n2->level, n2->i, n2->j);
	printf("realN1: [%d][%d][%d], ", realN1->level, realN1->i, realN1->j);
	if (realN2 == NULL) printf("realN2: NULL\n");
	else printf("realN2: [%d][%d][%d]\n", realN2->level, realN2->i, realN2->j);
	if (n2 == NULL) {
		return realN2 == NULL && n1 == realN1;
	} else {
		return (n1 == realN1 && n2 == realN2) || (n1 == realN2 && n2 == realN1);
	}
}

// tests neighbor pointers
void testNeighbors() {
	qNode* testRoot = new qNode(NULL, 0, 0);
	testRoot->expand(false);

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			qNode* c = testRoot->children[i * 2 + j];
			c->expand(false);
			if (c->neighbors[2+j] != testRoot->children[i*2+(1-j)]) {
				printf("WRONG INNER X NEIGHBOR (%d, %d)\n", i, j);
				assert(false);
			}
			if (c->neighbors[i] != testRoot->children[(1-i)*2+j]) {
				printf("WRONG INNER Y NEIGHBOR (%d,%d)\n", i, j);
				assert(false);
			}
		}
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			qNode* c = testRoot->children[i * 2 + j];
			for (int ii = 0; ii < 2; ii++) {
				for (int jj = 0; jj < 2; jj++) {
					qNode* cc = c->children[ii* 2 + jj];
					if (cc->neighbors[ii] != c->children[(1-ii) * 2 + jj]) {
						printf("WRONG INNER Y NEIGHBOR of %d    - (%d, %d): (%d, %d)\n", cc->id, i, j,  ii, jj);
						printf("expected: %d, actual: %d\n", cc->neighbors[ii]->id, c->children[(1-ii) * 2 + jj]->id);
						assert(false);
					}			
					if (cc->neighbors[2+jj] != c->children[ii * 2 + (1-jj)]) {
						printf("WRONG INNER X NEIGHBOR of %d    - (%d, %d): (%d, %d)\n", cc->id, i, j,  ii, jj);
						printf("expected: %d, actual: %d\n", cc->neighbors[2+jj]->id, c->children[ii * 2 + (1-jj)]->id);
						assert(false);
					}

					if (i == ii) {
						if (cc->neighbors[1-i] != NULL) {
							printf("null y neighbor incorrect: (%d,%d)(%d,%d)\n", i, j, ii, jj);
							assert(false);
						}
					} else {
						qNode* shouldNeighbor = testRoot->children[(1-i) * 2 + j]->children[(1-ii) * 2 + jj];
						if (cc->neighbors[1-ii] != shouldNeighbor) {
							printf("cross-node y neighbor of %d    incorrect: (%d,%d)(%d,%d)\n", cc->id, i, j, ii, jj);
							printf("expected: %d, actual: %d\n", cc->neighbors[1-ii]->id, shouldNeighbor->id);
							assert(false);
						}
					}
					if (j == jj) {
						if (cc->neighbors[3-j] != NULL) {
							printf("null x neighbor incorrect: (%d,%d)(%d,%d)", i, j, ii, jj);
							assert(false);
						}
					} else {
						qNode* shouldNeighbor = testRoot->children[i * 2 + (1-j)]->children[ii * 2 + (1-jj)];
						if (cc->neighbors[3-jj] != shouldNeighbor) {
							printf("cross-node x neighbor of %d    incorrect: (%d,%d)(%d,%d)\n", cc->id, i, j, ii, jj);
							printf("expected: %d, actual: %d\n", cc->neighbors[3-jj]->id, shouldNeighbor->id);
							assert(false);
						}
					}
				}
			}
		}
	}
}

void assertf (int id, double a, double b) {
	//printf("test for id: %d\n", id);
	if (fabs(a-b) > .00001) {
		printf("id %d incorrect!, a: %f, b: %f\n", id, a, b);
		assert (false);
	}
}

void testPressureInterp() {
	int oldLevels = levels;	
	bool oldPressureInterp = pressureInterp;
	levels = 3;
	pressureInterp = true;

	qNode* testRoot = new qNode(NULL, 0, 0);
	testRoot->expand(false);
	testRoot->children[1]->expand(false);

	testRoot->children[0]->p = 1;
	testRoot->children[2]->p = 6;
	testRoot->children[3]->p = 7;
	testRoot->children[1]->children[0]->p = 2;
	testRoot->children[1]->children[1]->p = 3;
	testRoot->children[1]->children[2]->p = 4;
	testRoot->children[1]->children[3]->p = 5;

	assertf(0, testRoot->children[0]->getFaceGradient(2, 0, P), 10.0);
	assertf(1, testRoot->children[0]->getFaceGradient(2, 1, P), 0);
	assertf(2, testRoot->children[0]->getFaceGradient(2, 2, P), 16.0/3.0);
	assertf(3, testRoot->children[0]->getFaceGradient(2, 3, P), 0);

	assertf(4, testRoot->children[2]->getFaceGradient(2, 0, P), 0);
	assertf(5, testRoot->children[2]->getFaceGradient(2, 1, P), -10.0);
	assertf(6, testRoot->children[2]->getFaceGradient(2, 2, P), 2.0);
	assertf(7, testRoot->children[2]->getFaceGradient(2, 3, P), 0);

	assertf(8, testRoot->children[3]->getFaceGradient(2, 0, P), 0);
	assertf(9, testRoot->children[3]->getFaceGradient(2, 1, P), -20.0/3.0);
	assertf(10, testRoot->children[3]->getFaceGradient(2, 2, P), 0.0);
	assertf(11, testRoot->children[3]->getFaceGradient(2, 3, P), -2.0);

	assertf(12, testRoot->children[1]->children[0]->getFaceGradient(2, 0, P), 8.0);
	assertf(13, testRoot->children[1]->children[0]->getFaceGradient(2, 1, P), 0.0);
	assertf(14, testRoot->children[1]->children[0]->getFaceGradient(2, 2, P), 4.0);
	assertf(15, testRoot->children[1]->children[0]->getFaceGradient(2, 3, P), -8.0/3.0);

	assertf(16, testRoot->children[1]->children[1]->getFaceGradient(2, 0, P), 8.0);
	assertf(17, testRoot->children[1]->children[1]->getFaceGradient(2, 1, P), 0.0);
	assertf(18, testRoot->children[1]->children[1]->getFaceGradient(2, 2, P), 0.0);
	assertf(19, testRoot->children[1]->children[1]->getFaceGradient(2, 3, P), -4.0);

	assertf(20, testRoot->children[1]->children[2]->getFaceGradient(2, 0, P), 22.0/3.0);
	assertf(21, testRoot->children[1]->children[2]->getFaceGradient(2, 1, P), -8.0);
	assertf(22, testRoot->children[1]->children[2]->getFaceGradient(2, 2, P), 4.0);
	assertf(23, testRoot->children[1]->children[2]->getFaceGradient(2, 3, P), -14.0/3.0);

	assertf(24, testRoot->children[1]->children[3]->getFaceGradient(2, 0, P), 16.0/3);
	assertf(25, testRoot->children[1]->children[3]->getFaceGradient(2, 1, P), -8.0);
	assertf(26, testRoot->children[1]->children[3]->getFaceGradient(2, 2, P), 0.0);
	assertf(27, testRoot->children[1]->children[3]->getFaceGradient(2, 3, P), -4.0);

	levels = oldLevels;
	pressureInterp = oldPressureInterp;
	delete testRoot;
}

void computeNodalVelocity(qNode* root) {
	// loop to try to figure out if done
	int oldTotalLeft;
	int totalLeft = 10000000;
	while (totalLeft > 0) {
		oldTotalLeft = totalLeft;
		totalLeft = root->attemptComputeNodalVals(root);
		assert (totalLeft < oldTotalLeft);
	}
}

void testNodalVelocity() {
	int oldLevels = levels;

	levels = 4;

	root = new qNode(NULL, 0, 0);
	root->expand(false);
	root->children[0]->expand(false);
	root->children[0]->children[3]->expand(false);
	root->children[3]->expand(false);

	root->children[0]->vx = 0.0;
	root->children[0]->vx2 = 6;
	root->children[0]->vy = 0.0;
	root->children[0]->vy2 = 13;

	root->children[1]->vx = 6;
	root->children[1]->vx2 = 0.0;
	root->children[1]->vy = 0.0;
	root->children[1]->vy2 = 16.5;

	root->children[2]->vx = 0.0;
	root->children[2]->vx2 = 20;
	root->children[2]->vy = 13;
	root->children[2]->vy2 = 0.0;

	root->children[3]->vx = 20;
	root->children[3]->vx2 = 0.0;
	root->children[3]->vy = 16.5;
	root->children[3]->vy2 = 0.0;

	root->children[0]->children[0]->vx = 0.0;
	root->children[0]->children[0]->vx2 = 1.0;
	root->children[0]->children[0]->vy = 0.0;
	root->children[0]->children[0]->vy2 = 3.0;

	root->children[0]->children[1]->vx = 1.0;
	root->children[0]->children[1]->vx2 = 2.0;
	root->children[0]->children[1]->vy = 0.0;
	root->children[0]->children[1]->vy2 = 4.5;

	root->children[0]->children[2]->vx = 0.0;
	root->children[0]->children[2]->vx2 = 8.5;
	root->children[0]->children[2]->vy = 3;
	root->children[0]->children[2]->vy2 = 11;

	root->children[0]->children[3]->vx = 8.5;
	root->children[0]->children[3]->vx2 = 10.5;
	root->children[0]->children[3]->vy = 4.5;
	root->children[0]->children[3]->vy2 = 14.5;

	root->children[3]->children[0]->vx = 18;
	root->children[3]->children[0]->vx2 = 19;
	root->children[3]->children[0]->vy = 16;
	root->children[3]->children[0]->vy2 = 20;

	root->children[3]->children[1]->vx = 19;
	root->children[3]->children[1]->vx2 = 0.0;
	root->children[3]->children[1]->vy = 17;
	root->children[3]->children[1]->vy2 = 21;

	root->children[3]->children[2]->vx = 22;
	root->children[3]->children[2]->vx2 = 23;
	root->children[3]->children[2]->vy = 20;
	root->children[3]->children[2]->vy2 = 0;

	root->children[3]->children[3]->vx = 23;
	root->children[3]->children[3]->vx2 = 0;
	root->children[3]->children[3]->vy = 21;
	root->children[3]->children[3]->vy2 = 0;

	root->children[0]->children[3]->children[0]->vx = 6;
	root->children[0]->children[3]->children[0]->vx2 = 7;
	root->children[0]->children[3]->children[0]->vy = 4;
	root->children[0]->children[3]->children[0]->vy2 = 9;

	root->children[0]->children[3]->children[1]->vx = 7;
	root->children[0]->children[3]->children[1]->vx2 = 8;
	root->children[0]->children[3]->children[1]->vy = 5;
	root->children[0]->children[3]->children[1]->vy2 = 10;

	root->children[0]->children[3]->children[2]->vx = 11;
	root->children[0]->children[3]->children[2]->vx2 = 12;
	root->children[0]->children[3]->children[2]->vy = 9;
	root->children[0]->children[3]->children[2]->vy2 = 14;

	root->children[0]->children[3]->children[3]->vx = 12;
	root->children[0]->children[3]->children[3]->vx2 = 13;
	root->children[0]->children[3]->children[3]->vy = 10;
	root->children[0]->children[3]->children[3]->vy2 = 15;

	computeNodalVelocity(root);
	
	assertf(0, root->children[0]->cvx[1], 2);
	assertf(1, root->children[0]->children[1]->cvx[1], 2);
	assertf(2, root->children[1]->cvx[0], 2);
	assertf(3, root->children[0]->cvy[1], 0.0);
	assertf(4, root->children[0]->children[1]->cvy[1], 0.0);
	assertf(5, root->children[1]->cvy[0], 0.0);

	assertf(6, root->children[0]->cvy[2], 11);
	assertf(7, root->children[0]->children[2]->cvy[2], 11);
	assertf(8, root->children[2]->cvy[0], 11);
	assertf(9, root->children[0]->cvx[2], 0.0);
	assertf(10, root->children[0]->children[2]->cvx[2], 0.0);
	assertf(11, root->children[2]->cvx[0], 0.0);

	assertf(12, root->children[1]->cvy[3], 17);
	assertf(13, root->children[3]->cvy[1], 17);
	assertf(14, root->children[3]->children[1]->cvy[1], 17);
	assertf(15, root->children[1]->cvx[3], 0.0);
	assertf(16, root->children[3]->cvx[1], 0.0);
	assertf(17, root->children[3]->children[1]->cvx[1], 0.0);

	assertf(18, root->children[3]->cvx[2], 22);
	assertf(19, root->children[3]->children[2]->cvx[2], 22);
	assertf(20, root->children[2]->cvx[3], 22);
	assertf(21, root->children[3]->cvy[2], 0.0);
	assertf(22, root->children[3]->children[2]->cvy[2], 0.0);
	assertf(23, root->children[2]->cvy[3], 0.0);

	// middle corner
	assertf(24, root->children[0]->cvx[3], 15.5);
	assertf(25, root->children[1]->cvx[2], 15.5);
	assertf(26, root->children[2]->cvx[1], 15.5);
	assertf(27, root->children[3]->cvx[0], 15.5);
	assertf(28, root->children[0]->children[3]->cvx[3], 15.5);
	assertf(29, root->children[0]->children[3]->children[3]->cvx[3], 15.5);
	assertf(30, root->children[3]->children[0]->cvx[0], 15.5);
	assertf(31, root->children[0]->cvy[3], 15.5);
	assertf(32, root->children[1]->cvy[2], 15.5);
	assertf(33, root->children[2]->cvy[1], 15.5);
	assertf(34, root->children[3]->cvy[0], 15.5);
	assertf(35, root->children[0]->children[3]->cvy[3], 15.5);
	assertf(36, root->children[0]->children[3]->children[3]->cvy[3], 15.5);
	assertf(37, root->children[3]->children[0]->cvy[0], 15.5);

	//
	assertf(38, root->children[0]->children[0]->cvx[1], 1);
	assertf(39, root->children[0]->children[1]->cvx[0], 1);
	assertf(40, root->children[0]->children[0]->cvy[1], 0.0);
	assertf(41, root->children[0]->children[1]->cvy[0], 0.0);

	assertf(42, root->children[0]->children[0]->cvy[2], 3);
	assertf(43, root->children[0]->children[2]->cvy[0], 3);
	assertf(44, root->children[0]->children[0]->cvx[2], 0.0);
	assertf(45, root->children[0]->children[2]->cvx[0], 0.0);

	assertf(46, root->children[0]->children[1]->cvx[3], 6);
	assertf(47, root->children[0]->children[3]->cvx[1], 6);
	assertf(48, root->children[0]->children[3]->children[1]->cvx[1], 6);
	assertf(49, root->children[0]->children[1]->cvy[3], 103.0/16.0);
	assertf(50, root->children[0]->children[3]->cvy[1], 103.0/16.0);
	assertf(51, root->children[0]->children[3]->children[1]->cvy[1], 103.0/16.0);

	assertf(52, root->children[0]->children[3]->cvx[2], 309.0/32.0);
	assertf(53, root->children[0]->children[3]->children[2]->cvx[2], 309.0/32.0);
	assertf(54, root->children[0]->children[2]->cvx[3], 309.0/32.0);
	assertf(55, root->children[0]->children[3]->cvy[2], 13);
	assertf(56, root->children[0]->children[3]->children[2]->cvy[2], 13);
	assertf(57, root->children[0]->children[2]->cvy[3], 13);

	
	assertf(58, root->children[0]->children[0]->cvx[3], 13.0/3.0);
	assertf(59, root->children[0]->children[0]->cvy[3], 11.0/3.0);
	assertf(60, root->children[0]->children[1]->cvx[2], 13.0/3.0);
	assertf(61, root->children[0]->children[1]->cvy[2], 11.0/3.0);
	assertf(62, root->children[0]->children[2]->cvx[1], 13.0/3.0);
	assertf(63, root->children[0]->children[2]->cvy[1], 11.0/3.0);
	assertf(64, root->children[0]->children[3]->cvx[0], 13.0/3.0);
	assertf(65, root->children[0]->children[3]->cvy[0], 11.0/3.0);
	assertf(66, root->children[0]->children[3]->children[0]->cvx[0], 13.0/3.0);	
	assertf(67, root->children[0]->children[3]->children[0]->cvy[0], 11.0/3.0);

	///
	assertf(68, root->children[3]->children[0]->cvx[1], 275.0/16.0);
	assertf(69, root->children[3]->children[1]->cvx[0], 275.0/16.0);
	assertf(70, root->children[3]->children[0]->cvy[1], 16.5);
	assertf(71, root->children[3]->children[1]->cvy[0], 16.5);

	assertf(72, root->children[3]->children[0]->cvy[2], 57.0/8.0);
	assertf(73, root->children[3]->children[2]->cvy[0], 57.0/8.0);
	assertf(74, root->children[3]->children[0]->cvx[2], 20);
	assertf(75, root->children[3]->children[2]->cvx[0], 20);

	assertf(76, root->children[3]->children[1]->cvx[3], 0.0);
	assertf(77, root->children[3]->children[3]->cvx[1], 0.0);
	assertf(78, root->children[3]->children[1]->cvy[3], 21);
	assertf(79, root->children[3]->children[3]->cvy[1], 21);

	assertf(80, root->children[3]->children[3]->cvx[2], 0.0);
	assertf(81, root->children[3]->children[2]->cvx[3], 0.0);
	assertf(82, root->children[3]->children[3]->cvy[2], 23);
	assertf(83, root->children[3]->children[2]->cvy[3], 23);

	assertf(84, root->children[3]->children[0]->cvx[3], 21);
	assertf(85, root->children[3]->children[0]->cvy[3], 20.5);	
	assertf(86, root->children[3]->children[1]->cvx[2], 21);	
	assertf(87, root->children[3]->children[1]->cvy[2], 20.5);	
	assertf(88, root->children[3]->children[2]->cvx[1], 21);
	assertf(89, root->children[3]->children[2]->cvy[1], 20.5);	
	assertf(90, root->children[3]->children[3]->cvx[0], 21);	
	assertf(91, root->children[3]->children[3]->cvy[0], 20.5);	

	///
	assertf(92, root->children[0]->children[3]->children[0]->cvx[1], 10.0/3.0);
	assertf(93, root->children[0]->children[3]->children[1]->cvx[0], 10.0/3.0);
	assertf(94, root->children[0]->children[3]->children[0]->cvy[1], 4.5);
	assertf(95, root->children[0]->children[3]->children[1]->cvy[0], 4.5);

	assertf(96, root->children[0]->children[3]->children[0]->cvy[2], 23.0/3.0);
	assertf(97, root->children[0]->children[3]->children[2]->cvy[0], 23.0/3.0);
	assertf(98, root->children[0]->children[3]->children[0]->cvx[2], 8.5);
	assertf(99, root->children[0]->children[3]->children[2]->cvx[0], 8.5);

	assertf(100, root->children[0]->children[3]->children[1]->cvx[3], 10.5);
	assertf(101, root->children[0]->children[3]->children[3]->cvx[1], 10.5);
	assertf(102, root->children[0]->children[3]->children[1]->cvy[3], 349.0/32.0);
	assertf(103, root->children[0]->children[3]->children[3]->cvy[1], 349.0/32.0);

	assertf(104, root->children[0]->children[3]->children[3]->cvx[2], 535.0/64.0);
	assertf(105, root->children[0]->children[3]->children[2]->cvx[3], 535.0/64.0);
	assertf(106, root->children[0]->children[3]->children[3]->cvy[2], 14.5);
	assertf(107, root->children[0]->children[3]->children[2]->cvy[3], 14.5);
	
	assertf(108, root->children[0]->children[3]->children[0]->cvx[3], 9.5);
	assertf(109, root->children[0]->children[3]->children[0]->cvy[3], 9.5);	
	assertf(110, root->children[0]->children[3]->children[1]->cvx[2], 9.5);	
	assertf(111, root->children[0]->children[3]->children[1]->cvy[2], 9.5);	
	assertf(112, root->children[0]->children[3]->children[2]->cvx[1], 9.5);
	assertf(113, root->children[0]->children[3]->children[2]->cvy[1], 9.5);	
	assertf(114, root->children[0]->children[3]->children[3]->cvx[0], 9.5);	
	assertf(115, root->children[0]->children[3]->children[3]->cvy[0], 9.5);	

	delete root;
	levels = oldLevels;
}

// poisson solver functions

// returns max(abs(R))
double computeResidual(qNode* node) {
	int size = 1<<node->level;
	if (node->leaf) { // if leaf cell, compute residual
		// compute it: R = div(vx, vy) - 1/(ha)*sum of (s * grad) for each face
		double aSum = 0.0;
		double bSum = 0.0;
		node->getLaplacian(levels - 1, &aSum, &bSum, P);
		double laplacian = (aSum * node->p + bSum) * size * size;
		
		double R = node->divV - laplacian;
		
		node->R = R;
		//printf("at [%d][%d], divV: %f, flux: %f, R: %f\n", node->i, node->j, node->divV, flux, R);
		// if a*R > e, not done
		if (fabs(R) > eps) {
			doneVCycle = false;
			//printf("more work to do at [%d][%d], %f is bigger than epsilon of %f\n", i, j, grid[d][i*size+j].R, eps);
		} else {
			//printf("done with this cell already, %f is smaller than epsilon of %f\n", fabs(grid[d][i*size+j].R), eps);
		}
		return fabs(R);
	} else {
		double maxR = 0.0;
        node->R = 0.0;
		for (int k = 0; k < 4; k++) {
			maxR = std::max(maxR, computeResidual(node->children[k]));
			node->R += node->children[k]->R;
		}
		node->R /= 4.0;
		return maxR;
	}
}


// relax the node at the given multilevel
// puts the new value of dp in temp
bool relaxRecursive(qNode* node, int ml) {
	assert (node->level <= ml);
	int size = 1<<node->level;
	if (node->level == ml || node->leaf) {
		double aSum = 0.0;
        double bSum = 0.0;
        double dp;
        double h;
		double oldDif = node->dp;

		node->getLaplacian(ml, &aSum, &bSum, DP);
		
        // A*R = bSum - aSum*dp, dp = (bSum - A*R)/asum, A = h*h = 1/size^2
		node->temp = -(bSum - node->R/size/size)/aSum;

		double diff = oldDif - node->temp;

        //printf("relaxing[%d][%d][%d]: aSum: %f, bSumL %f, R: %f, result: %f, diff from old: %f\n", d, i, j, aSum, bSum, grid[d][i*size+j].R, grid[d][i*size+j].temp, diff);
		if (fabs(diff) > eps) {
			return false;
		}
		return true;
	} else {
		// relax children
		bool done = true;
		for (int k = 0; k < 4; k++) {
			done &= relaxRecursive(node->children[k], ml);
		}
		return done;
	}
}

void recursiveInject(qNode* node, int ml) {
	if (node->level == ml) {
		node->dp = node->parent->dp;
	}
	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			recursiveInject(node->children[k], ml);
		}
	}
}

void recursiveUpdate(qNode* node, int ml) {	
	node->dp = node->temp;
	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			recursiveUpdate(node->children[k], ml);
		}
	}
}

void relax(int d, int r) {
	assert (d < levels);

	//printf("relaxing level %d\n", d);
	// get initial gress from previous level, if possible
	recursiveInject(root, d);

	bool done = false;
    int totalCycles = 0;
	while (r-- > 0/* && !done*/) {
        totalCycles++;
		done = relaxRecursive(root, d);
		recursiveUpdate(root, d);
	}
	//printf("dp matrix with %d cycles left: \n", r);
    /*printf("relaxing took %d cycles\n", totalCycles);
    for (d = 0; d <= maxlevel; d++) {
        int size = 1<<d;
        //printf("level %d\n", d);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (grid[d][i*size+j].used) {
                    //printf(" %.4f", grid[d][i*size+j].dp);
                }
            }
            //printf("\n");
        }
    }*/
}

// adaptivity

// returns true if leaf should be expanded, false if it should not

double thresh = 0.1;

double curlAdaptFunction(qNode* node) {
	double curl = node->getCurl();
	return fabs(curl) > thresh;
}

double vnormAdaptFunction(qNode* node) {
	int size = 1<<node->level;
	double vnorm = fabs(node->vx2-node->vx)/size + fabs(node->vy2-node->vy)/size;
	return vnorm > thresh;
}

bool pGradAdaptFunction(qNode* node) {
    //std::pair<double, double> pgrad = node->getPressureGradient();
	std::pair<double, double> pgrad = node->getValueGradient(P);
    return fabs(pgrad.first + pgrad.second) > thresh/(1<<node->level);
}

bool adaptFunction(qNode* node) {
	assert(adaptScheme != ADAPTNONE);
	if (adaptScheme == PGRAD) {
    	return pGradAdaptFunction(node);
	} else if (adaptScheme == CURL) {
		return curlAdaptFunction(node);
	} else if (adaptScheme == VNORM) {
		return vnormAdaptFunction(node);
	}
	assert (false);
}

// trees should have the same structure
// use old tree for adapt calculations
// copies new pressure over
// returns true if any nodes were changed, false otherwise
bool recursiveAdaptAndCopy(qNode* node, qNode* oldNode) {
	// TODO implement properly with nodal/face velocities
	node->p = oldNode->p;
	if (node->level == levels - 1) return false;
	assert (node->leaf == oldNode->leaf);
	if (node->leaf) {
		if (adaptFunction(oldNode)) {

			oldNode->expand(true);
			node->expand(false);
			for (int k = 0; k < 4; k++) {
				node->children[k]->p = oldNode->children[k]->p;
				node->children[k]->vx = oldNode->children[k]->vx;
				node->children[k]->vy =  oldNode->children[k]->vy;
				node->children[k]->vx2 = oldNode->children[k]->vx2;
				node->children[k]->vy2 = oldNode->children[k]->vy2;
				for (int c = 0; c < 4; c++) {
					node->cvx[c] = oldNode->cvx[c];
					node->cvy[c] = oldNode->cvy[c];
				}
			}
			// reset old node to old state so it is the same as before, so it can be used in other calculations
			oldNode->contract(false);

			// now node is adapted, with values from oldnode's calculations
			return true;
		}
		return false;
	}
	bool allChildrenLeaves = true;
    for (int k = 0; k < 4; k++) {
        if (!node->children[k]->leaf) {
            allChildrenLeaves = false;
            break;
        }
    }
	if (allChildrenLeaves && !adaptFunction(oldNode)) {
        // this wouldn't be expanded if it was a leaf, so it shouldn't have leaf children
		node->contract(true);
		//node->p = oldNode->p;
		return true;
    } else {
		bool anyChanged = false;
        for (int k = 0; k < 4; k++) {
            anyChanged |= recursiveAdaptAndCopy(node->children[k], oldNode->children[k]);
        }
		return anyChanged;
    }
	
}

void copy(qNode* node, qNode* oldNode) {
	// need to create/delete new nodes if it adapted
	if (node->leaf && !oldNode->leaf) {
		// old node expanded, expand new node
		node->expand(false);
	} else if (!node->leaf && oldNode->leaf) {
		// old node contracted
		node->contract(false); // this doesn't handle multiple contractions per step, but rn only do single contract/expand per step
	}
	
	node->p = oldNode->p;

	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			copy(node->children[k], oldNode->children[k]);
		}
	}
}

// advect face velocities, reading from old tree and writing to new tree
void setNewAdvect(qNode* node, qNode* oldNode) {
	if (node->leaf) {
		int size = 1<<node->level;
		double x, y;
		qNode* last;
		std::pair<double, double> newvel;
		// vx
		if (node->j == 0) {
			node->vx = 0.0;
		} else {
			x = (double(node->j))/size;
			y = (node->i + 0.5)/size;
			last = getSemiLagrangianLookback(oldRoot, &x, &y, 1, levels - 1, oldNode->vx, (oldNode->cvy[0] + oldNode->cvy[2])/2.0);
			// TODO implement
			newvel = last->getVelocityAt(x, y);
			node->vx = newvel.first;
		}
		// vx2
		if (node->j == size - 1) {
			node->vx2 = 0.0;
		} else {
			x = (node->j + 1.0)/size;
			y = (node->i + 0.5)/size;
			last = getSemiLagrangianLookback(oldRoot, &x, &y, 1, levels - 1, oldNode->vx2, (oldNode->cvy[1] + oldNode->cvy[3])/2.0);
			// TODO implement
			newvel = last->getVelocityAt(x, y);
			node->vx2 = newvel.first;		
		}
		// vy
		if (node->i == 0) {
			node->vy = 0.0;
		} else {
			x = (node->j + 0.5)/size;
			y = (double(node->i))/size;
			last = getSemiLagrangianLookback(oldRoot, &x, &y, 1, levels - 1, (oldNode->cvx[0] + oldNode->cvx[1])/2.0, oldNode->vy);
			// TODO implement
			newvel = last->getVelocityAt(x, y);
			node->vy = newvel.second;
		}
		// vx2
		if (node->i == size - 1) {
			node->vy2 = 0.0;
		} else {
			x = (node->j + 0.5)/size;
			y = (node->i + 1.0)/size;
			last = getSemiLagrangianLookback(oldRoot, &x, &y, 1, levels - 1, (oldNode->cvx[2] + oldNode->cvx[3])/2.0, oldNode->vy2);
			// TODO implement
			newvel = last->getVelocityAt(x, y);
			node->vy2 = newvel.second;
		}
	} else {

		// average child velocities to this node
		for (int k = 0; k < 4; k++) {
			setNewAdvect(node->children[k], oldNode->children[k]);
		}
		node->computeVelFromChildren();
	}
}

// assumes nodal velocities are already set on old tree
void advectAndCopy() {

	copy(root, oldRoot);

	// assume nodal values have been computed on old tree
	// set nodal values on new tree to advected nodal values on old tree
	setNewAdvect(root, oldRoot);

	// computeNodalVelocities(root);
}

void correctPressure(qNode* node) {
	//node->p += node->dp;
	if (node->leaf) {
		//printf("at node %d (%d, %d) correcting pressure %f by %f\n", node->level, node->i, node->j, node->p, node->dp);
		node->p += node->dp;
	} else {
		for (int k = 0; k < 4; k++) {
			correctPressure(node->children[k]);
		}
	}
}

void project(qNode* node) {
	// correct velocity with updated pressure field to make non-divergent
	if (node->leaf) {
		//if (node->i == 0 || node->j == 0) return;
		double pgradX, pgradY;
		int size = 1<<node->level;
		if (node->j > 0) {
			pgradX = -node->getFaceGradient(levels - 1, 3, P);
			node->vx -= pgradX;
		}
		if (node->i > 0) {
			pgradY = -node->getFaceGradient(levels - 1, 1, P);
			node->vy -= pgradY;
		}
		if (node->j < size - 1) {
			pgradX = node->getFaceGradient(levels - 1, 2, P);
			node->vx2 -= pgradX;
		}
		if (node->i < size - 1) {
			pgradY = node->getFaceGradient(levels - 1, 0, P);
			node->vy2 -= pgradY;
		}
		
	} else {
		for (int k = 0; k < 4; k++) {
			project(node->children[k]);
		}

		node->vx = (node->children[0]->vx + node->children[2]->vx)/2.0;
		node->vy = (node->children[0]->vy + node->children[1]->vy)/2.0;
		node->vx2 = (node->children[2]->vx2 + node->children[3]->vx2)/2.0;
		node->vy2 = (node->children[1]->vy2 + node->children[3]->vy2)/2.0;
	}
}

void clampVelocity(qNode* node) {
	// clamp velocity
	int size = 1<<node->level;
	if (node->i == 0 || node->i == size-1) {
		node->vy = 0.0;
	}
	if (node->j == 0 || node->j == size-1) {
		node->vx = 0.0;
	}
	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			clampVelocity(node->children[k]);
		}	
	}
}

void eulerAdvectParticle(Particle& p, std::pair<double, double> v) {
	p.x += v.first * dt;
	p.y += v.second * dt;
}

void advectParticles() {
	for (int i = 0; i < numParticles; i++) {
		std::pair<double, double> vGrad = getLeaf(root, particles[i].x, particles[i].y, levels-1)->getVelocityAt(particles[i].x, particles[i].y);
		if (particleAlgorithm == EULER) {
			eulerAdvectParticle(particles[i], vGrad);
		}
	}	
}

double runVCycle() {
	//relax(0, MAX_RELAX_COUNT);
	root->dp = 0;
	// do not relax level 0 since it makes no sense to do so...

	for (int d = 1; d < levels; d++) {
		relax(d, (d==1) ? MAX_RELAX_COUNT : OPTIMIZED_RELAX_COUNT);
	}

	//printf("correcting pressure\n");
	// use corrections to improve pressure
	correctPressure(root);

	// TODO REMOVE
	/*printf("printing pressure after corrections\n");
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			printf("pressure at [%d][%d][%d]: %f\n", k, i, j, grid[k][i*size+j].p);
		}
	}
	*/

	
	// re-compute residual
	doneVCycle = true;
	return computeResidual(root);

	// TODO REMOVE
	// doneVCycle = true;
	//if (numCycles == 1) break;

	//printf("end v cycle\n");

}

void poissonAverageR(qNode* node, double* total) {
	if (node->leaf) {
		int size = 1<<node->level;
		*total += node->R / size / size;
	} else {
		for (int k = 0; k < 4; k++) {
			poissonAverageR(node->children[k], total);
		}
	}
}

void poissonCorrectR(qNode* node, double K) {
	node->R -= K;
	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			poissonCorrectR(node->children[k], K);
		}
	}
}

double getMaxR(qNode* node) {
	if (node->leaf)	return fabs(node->R);
	double max = 0.0;
	for (int k = 0; k < 4; k++)
		max = fmax(max, getMaxR(node->children[k]));
	return max;
}

double fixAndCheckResidual() {
	double avgR = 0.0;
	poissonAverageR(root, &avgR);
	poissonCorrectR(root, avgR);

	double newMaxR = getMaxR(root);
	if (newMaxR < eps) {
		doneVCycle = true;
	}
	return newMaxR;
}

void projectionCheck(qNode* node, double* max, double* avg) {
	if (node->leaf) {
		//printf("divV after projection: %f\n", node->divV);
		*max = fmax(*max, fabs(node->divV));
		int size = 1<<node->level;
		*avg += fabs(node->divV)/size/size;
	} else {
		for (int k = 0; k < 4; k++) {
			projectionCheck(node->children[k], max, avg);
		}
	}
}

void runStep() {
	printf("running step %d\n", frameNumber + 1);
	/*for (int d = 0; d < levels; d++) {
		printf("level %d\n", d);
		int size = 1<<d;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				printf(" %.3f", grid[d][i*size+j].p);
			}
			printf("\n");
		}
	}*/
	std::swap(root, oldRoot);

	double totalTime  = 0.0;

	// compute velocity divergence for advection calculation
	// TODO fix?
	//oldRoot->computeVelocityDivergence();
	
	// advect velocity, copy old stuff
	printf("advecting velocity\n");
	startTime();
	advectAndCopy();
	totalTime += endTime("advecting and copying");
	
	printf("computing divV\n");

	// grid's vx and vy are now provisional velocity
	// compute velocity divergence of new grid to use in residual calcuation
	root->computeVelocityDivergence();

	printf("starting poisson solver\n");

	startTime();

	// compute all residuals, so that the whole bottom multilevel is computed
	doneVCycle = true;
	computeResidual(root);

	fixAndCheckResidual();

	int numCycles = 0;
	
	while (!doneVCycle) {
		numCycles++;
		//printf("start v cycle %d\n", numCycles);

		double residual = runVCycle();	
		double newR = fixAndCheckResidual();
		printf("(corrected) residual after %d cycles: %f\n", numCycles, newR);
	}
	printf("end poisson solver, took %d cycles\n", numCycles);
	totalTime += endTime("poisson solver");
	

	printf ("doing adaptivity\n");
	startTime();
    // given new state, do adaptivity
	if (adaptScheme != ADAPTNONE) {
		std::swap(root, oldRoot);
		recursiveAdaptAndCopy(root, oldRoot);
	}
	totalTime += endTime("adapting");

	startTime();
	double max = 0.0;
	double avg = 0.0;

	root->computeVelocityDivergence();
	projectionCheck(root, &max, &avg);
	printf("max divV before project: %f, avg divV before project: %f\n", max, avg);

	project(root);

	max = 0.0;
	avg = 0.0;

	root->computeVelocityDivergence();
	projectionCheck(root, &max, &avg);
	printf("max divV afer project: %f, avg divV after project: %f\n", max, avg);

	clampVelocity(root);
	totalTime += endTime("projecting and clamping velocity");

	printf("total run step time: %.3f\n", totalTime);
	
	resetProfiling();
	root->profile();
	printProfiling();

	if (particleAlgorithm != PARTICLENONE) {
		advectParticles();
	}
	//printValue(P);

	// increase frame number here instead of in display so that you get 1 of each frame, not dependent on glut's redrawing, like on alt-tabbing or anything
	frameNumber++;
}

void initNodeLeftUniform(qNode* node) {
	int size = 1<<node->level;
	//node->p = 1.0/size/size;
	node->p = 0.0;
	node->vx = 0.0;
	node->vy = 0.0;
	node->phi = 0.0;
	if (node->i == size/2 || node->i == size/2-1) {
		node->vx = -1.0;
	}	
}

void initNodeSink(qNode* node) {
	int size = 1<<node->level;
	node->p = 1.0/size/size;
	node->phi = 0.0;

	double center = 0.5;
	node->vx = center - (node->j + 0.5)/size;
	node->vy = center - (node->i + 0.5)/size;
	double len = sqrt(node->vx * node->vx + node->vy * node->vy);
	node->vx /= len;
	node->vy /= len;
}

void initNodeSrcSink(qNode* node) {
	int size = 1<<node->level;
	node->p = 1.0/size/size;
	//node->p = 1.0;
	node->phi = 0.0;
	
	double src = .125;
	double sink = .875;

	double vxsrc = (node->j + 0.5)/size - src;
	double vysrc = (node->i + 0.5)/size - src;
	double lensrc = sqrt(vxsrc * vxsrc + vysrc * vysrc);

	double vxsink = sink - (node->j + 0.5)/size;
	double vysink = sink - (node->i + 0.5)/size;
	double lensink = sqrt(vxsink * vxsink + vysink * vysink);
	
	node->vx = vxsrc/lensrc + vxsink/lensink;
	node->vy = vysrc/lensrc + vysink/lensink;
}

void initPoissonTest(qNode* node) {
	
	return;
	node->phi = 0;
	node->vx = 0;
	node->vy = 0;
	
	int size = 1<<node->level;
	double x = (node->j + 0.5)/size;
	double y = (node->i + 0.5)/size;

	node->p = 0;
	
}

void initProjectTest(qNode* node) {
	int size = 1<<node->level;
	node->p = 0.0;
	node->dp = 0.0;

	assert (projectTestFunc != PROJECTFUNCNONE);

	double x = ((float)node->j)/size;
	double y = ((float)node->i)/size;

	if (projectTestFunc == PROJECTXY) {
		node->vx = (node->j+0.5)*(node->j+0.5)/size/size;
		node->vy = 1-(node->i+0.5)*(node->i+0.5)/size/size;
	
		node->vx2 = (node->j+1.5)*(node->j+1.5)/size/size;
		node->vy2 = 1-(node->i+1.5)*(node->i+1.5)/size/size;
	} else if (projectTestFunc == PROJECTSIN) {
		// vx
		double x2 = x;
		double y2 = y + 0.5/size;
		node->vx = 2*M_PI*sin(2*M_PI*x2) - 2*M_PI*sin(2*M_PI*x2)*cos(2*M_PI*y2);

		x2 = x+1.0/size;
		node->vx2 = 2*M_PI*sin(2*M_PI*x2) - 2*M_PI*sin(2*M_PI*x2)*cos(2*M_PI*y2);

		//vy
		x2 = x + 0.5/size;
		y2 = y;
		node->vy = 2*M_PI*sin(2*M_PI*y2) - 2*M_PI*cos(2*M_PI*x2)*sin(2*M_PI*y2);

		y2 = y + 1.0/size;
		node->vy2 = 2*M_PI*sin(2*M_PI*y2) - 2*M_PI*cos(2*M_PI*x2)*sin(2*M_PI*y2);	
	}
}

void initNodeFunction(qNode* node) {
	if (startState == LEFT || startState == ERRORTEST)
		initNodeLeftUniform(node);
	else if (startState == SINK)
		initNodeSink(node);
	else if (startState == SRCSINK)
		initNodeSrcSink(node);
	else if (startState == POISSONTEST || startState == ADAPTTEST)
		initPoissonTest(node);
	else if (startState == PROJECTTEST)
		initProjectTest(node);
	else
		printf("invalid start state\n");
}

// inits leaves according to test and averages non-leaves
void initRecursive(qNode* node, int d) {
	if (node->level == d) {
		initNodeFunction(node);
	} else {
		node->expand(false);

		node->p = 0.0;
        node->vx = 0.0;
        node->vy = 0.0;
        //node->phi = 0.0;

        for (int k = 0; k < 4; k++) {
			initRecursive(node->children[k], d);
			node->p += node->children[k]->p;
			node->vx += node->children[k]->vx;
			node->vy += node->children[k]->vy;
			//node->phi += node->children[k]->phi;
        }

        node->p /= 4.0;
        node->vx /= 4.0;
        node->vy /= 4.0;
        //node->phi /= 4.0;
	}
}

int poissonk = 2;
int poissonl = 2;

double getPoissonVX(double x, double y){
	if (poissonTestFunc == POISSONXY) {
		double newy = 2*y-1;
		return 6*x*(x-1) * (newy*newy*newy*newy - 2*newy*newy + 1);
	} else if (poissonTestFunc == POISSONCOS) {
		return M_PI*2*sin(2*M_PI*x)*(1-cos(M_PI*2*y));
	} else if (poissonTestFunc == POISSONCOSKL) {
		return -M_PI*poissonk*sin(M_PI*poissonk*x)*cos(M_PI*poissonl*y);
	}
	assert (false);
}

double getPoissonVY(double x, double y){
	if (poissonTestFunc == POISSONXY) {
		return (2*x*x*x - 3*x*x+ 1) * 32*y*(2*y*y-3*y+1);
	} else if (poissonTestFunc == POISSONCOS) {
		return M_PI*2*sin(2*M_PI*y)*(1-cos(M_PI*2*x));
	} else if (poissonTestFunc == POISSONCOSKL) {
		return -M_PI*poissonl*cos(M_PI*poissonk*x)*sin(M_PI*poissonl*y);
	}
	assert (false);
}



void poissonReset(qNode* node) {
	node->p = 0;
	int size = 1<<node->level;
	double x = (node->j + 0.0)/size;
	double y = (node->i + 0.0)/size;

	// vel
	double x2 = x;
	double y2 = y + 0.5/size;
	node->vx = getPoissonVX(x2, y2);
	x2 = x + 1.0/size;
	node->vx2 = getPoissonVX(x2, y2);

	x2 = x + 0.5/size;
	y2 = y;
	node->vy = getPoissonVY(x2, y2);
	y2 = y + 1.0/size;
	node->vy2 = getPoissonVY(x2, y2);
	
	// set pressure to avg of function so it comes out with right constant
	if (poissonTestFunc == POISSONXY)
		node->p = .2667;
	else if (poissonTestFunc == POISSONCOS)
		node->p = 1.0;
	else if (poissonTestFunc == POISSONCOSKL)
		node->p = 0.0;
	else
		assert(false);
	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			poissonReset(node->children[k]);
		}
	}
}

void computePoissonError(qNode* node, double* total) {
	if (node->leaf) {
		int size = 1<<node->level;
		double x = (node->j + 0.5)/size;
		double y = (node->i + 0.5)/size;
		double correct  = 0;
		if (poissonTestFunc == POISSONXY) {
			double newy = 2*y-1;
			correct = (2*x*x*x - 3*x*x + 1)  *  (newy*newy*newy*newy - 2*newy*newy + 1);
		} else if (poissonTestFunc == POISSONCOS) {
			correct = (1-cos(M_PI * 2 * x)) * (1-cos(M_PI * 2 * y));
		} else if (poissonTestFunc == POISSONCOSKL) {
			correct = cos(M_PI * poissonk * x) * cos(M_PI * poissonl * y);
		} else {
			assert(false);
		}
		*total += fabs(node->p - correct)/size/size;
	} else {
		for (int k = 0; k < 4; k++) {
			computePoissonError(node->children[k], total);
		}
	}
}

void poissonAverage(qNode* node, double* total) {
	if (node->leaf) {
		int size = 1<<node->level;
		*total += node->p / size / size;
	} else {
		for (int k = 0; k < 4; k++) {
			poissonAverage(node->children[k], total);
		}
	}
}

void expandRadius(qNode* node, double radius) {
	if (node->leaf) {
		int size = 1<<node->level;
		double x = (node->j + 0.5)/size;
		double y = (node->i + 0.5)/size;
		double dist = sqrt((0.5-x)*(0.5-x) + (0.5-y)*(0.5-y));
		if (dist < radius) {
			node->expand(false);
		}
	} else {
		for (int k = 0; k < 4; k++) {
			expandRadius(node->children[k], radius);
		}
	}
}

double runPoissonTest(bool print) {
	// set all pressure to 0 again
	assert(poissonTestFunc != POISSONFUNCNONE);

	poissonReset(root);
	
	double initialResidual = computeResidual(root);
	//printf("initial residual: %f\n", initialResidual);

	// try manually fixing residual
	double avgR = 0.0;
	poissonAverageR(root, &avgR);
	//printf("new average residual: %f\n", avgR);
	poissonCorrectR(root, avgR);
	
	// run V cycles, compute stuff
	int i = 0;

	double newR, oldR;
	newR = initialResidual;
	double avgError;

	if (print) {
		startTime();
	}
	while (!doneVCycle) {
		i++;
		oldR = newR;
		newR = runVCycle();
		double total = 0.0;
		//double K = 0.0;
		//poissonAverage(root, &K);
		avgError = 0.0;
		computePoissonError(root,  &avgError);
		//printf("average error after %d vcycles: %f\n", i, avgError);

		// try manually fixing residual
		avgR = 0.0;
		poissonAverageR(root, &avgR);
		poissonCorrectR(root, avgR);
		avgR = 0.0;
		poissonAverageR(root, &avgR);
		//doneVCycle |= fabs(newR-oldR) < eps/100;
		newR = getMaxR(root);
		doneVCycle = newR < eps;

		//printf("residual after %d vcycles: %f\n", i, newR);
	}
	double time = 0.0;
	if (print) {
		time = endTime("poisson test");	
		printf("poisson test took %d vcycles\n", i);
		printf("average error: %f\n", avgError);
	}
	return time;

	//printPressure();
}

// computes error in gradient discretization
void checkPoissonErrorRecursive(qNode* node, double* max, double* avg) {
	if (node->leaf) {
		int size = 1<<node->level;
		for (int k = 0; k < 4; k++) {
			double vel;
			if (k == 0) vel = node->vy2;
			else if (k == 1) vel = -node->vy;
			else if (k == 2) vel = node->vx2;
			else vel = -node->vx;
			double grad = node->getFaceGradient(levels - 1, k, P);
			double error = fabs(grad - vel);
			//printf("d: %d, velocity: %f, pressure gradient :%f, error: %f\n", k, vel, grad, error);
			*max = fmax(*max, error);
			*avg += error/size/size/4;
		}
	} else {
		for (int k = 0; k < 4; k++) {
			checkPoissonErrorRecursive(node->children[k], max, avg);
		}
	}
}

void checkPoissonError() {
	// TODO pass values to sum
	double max = 0.0;
	double avg = 0.0;
	checkPoissonErrorRecursive(root, &max, &avg);	

	printf("gradient error. max: %f, avg: %f\n", max, avg);
}




void setErrorPressure(qNode* node) {
	int size = 1<<node->level;
	double x = (node->j + 0.5)/size;
	double y = (node->i + 0.5)/size;

	if (errorTestFunc == ERRORXY) {
		double newy = 2*y+1;
		node->p = (2*x*x*x - 3*x*x + 1)  *  (newy*newy*newy*newy - 2*newy*newy + 1);
	} else if (errorTestFunc == ERRORCOS) {
		node->p = (1-cos(M_PI * 2 * x)) * (1-cos(M_PI * 2 * y));	
	} else if (errorTestFunc == ERRORCOSKL) {
		int k = 3;
		int l = 3;
		node->p = cos(M_PI * k * x) * cos(M_PI * l * y);
	} else {
		assert(false);
	}
		
	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			setErrorPressure(node->children[k]);
		}
	}
}

double getRealDerivative(double x, double y, int k) {
	if (errorTestFunc == ERRORXY) {
		if (k < 2) // y partial
			return 0.0;
		else // x partial
			return 0.0;
	} else if (errorTestFunc == ERRORCOS) {
		if (k < 2) // y partial
			return 2*M_PI*sin(2*M_PI*y) - 2*M_PI*cos(2*M_PI*x)*sin(2*M_PI*y);
		else // x partial
			return 2*M_PI*sin(2*M_PI*x) - 2*M_PI*sin(2*M_PI*x)*cos(2*M_PI*y);
	} else if (errorTestFunc == ERRORCOSKL) {
		int kk = 3;
		int l = 3;
		if (k < 2) // y partial
			return -M_PI*l*cos(M_PI*kk*x)*sin(M_PI*l*y);
		else //x partial
			return -M_PI*kk*sin(M_PI*kk*x)*cos(M_PI*l*y);
	}
	assert(false);
}

double calculateError(qNode* node, double* avgError) {
	if (node->leaf) {
		int size = 1<<node->level;
		double max = 0.0;
		for (int k = 0; k < 4; k++) {
			// calculate gradient in direction k
			int ml = levels - 1;
			double calc = node->getFaceGradient(ml, k, P);
			
			double x = (node->j + 0.5 + 0.5*deltas[k][1])/size;
			double y = (node->i + 0.5 + 0.5*deltas[k][0])/size;
			double real = getRealDerivative(x, y, k);
			if (k % 2 == 1) {
				// switch directions because other way
				real = -real;
			}
			double error = fabs(real - calc);
			printf("Error for node %d: (%d, %d) in dir %d, ", node->level, node->i, node->j, k);
			if (ml == node->level) {
				printf(" at the same level.");
			} else if (ml < node->level) {
				printf(" up %d levels.", node->level - ml);
			} else {
				printf(" down %d levels.", ml - node->level);
			}
			printf("real: %f, calc: %f, error: %f\n", real, calc, error);
			*avgError += error/size/size;
			max = fmax(max, error);
		}
		return max;
	} else {
		double maxError = 0.0;
		for (int k = 0; k < 4; k++) {
			maxError = fmax(maxError, calculateError(node->children[k], avgError));
		}
		return maxError;
	}
}

void runErrorTest() {
	//adapt
	if (adaptScheme != ADAPTNONE) {
		expandRadius(root, .3);
		expandRadius(root, .2);
	}

	// set initial pressure values
	setErrorPressure(root);
	
	double avgError = 0.0;
	// do calculations
	double maxError = calculateError(root, &avgError);
	printf("maximum error: %f, average error: %f\n", maxError, avgError);
}

void adaptTestRecursive(qNode* node) {
	node->p = 0.0;
	if (node->level == levels - 1) {
		return;
	}
	if (node->leaf) {
		int size = 1<<node->level;
		double x1, y1, x2, y2;
		x1 = ((float)node->j)/size;
		y1 = ((float)node->i)/size;
		x2 = x1 + 1.0/size;
		y2 = y1 + 1.0/size;
		if (rectIntersectFunc(x1, y1, x2, y2)) {
			node->expand(false);
			for (int k = 0; k < 4; k++) {
				adaptTestRecursive(node->children[k]);
			}
		}
	} else {
		for (int k = 0; k < 4; k++) {
			adaptTestRecursive(node->children[k]);
		}
	}
}

void runAdaptTest() {	
	bool done = false;
	int numCycles = -1;

	//root = new qNode(NULL, 0, 0);

	adaptTestRecursive(root);

	resetProfiling();
	root->profile();
	printProfiling();	

}

void runProjectTest() {
	if (adaptScheme != ADAPTNONE) {
		expandRadius(root, .3);
		expandRadius(root, .2);
	}
	
	initProjectTest(root);	

	root->computeVelocityDivergence();

	computeResidual(root);

	// try manually fixing residual
	double avgR = 0.0;
	poissonAverageR(root, &avgR);
	poissonCorrectR(root, avgR);
	
	// run V cycles, compute stuff

	while (!doneVCycle) {
		runVCycle();
		// try manually fixing residual
		avgR = 0.0;
		poissonAverageR(root, &avgR);
		poissonCorrectR(root, avgR);

		//doneVCycle |= fabs(newR-oldR) < eps/100;
		double newR = getMaxR(root);

		doneVCycle = newR < eps;
	}
	

	double avg = 0.0;
	double max = 0.0;

	printf("divV before projection:\n");
	printValue(DIVV);
	
	projectionCheck(root, &max, &avg);
	printf("max: %f, avg: %f\n", max, avg);
	
	project(root);


	root->computeVelocityDivergence();

	printf("divV after projection:\n");
	avg = 0.0;
	max = 0.0;

	printValue(DIVV);
	projectionCheck(root, &max, &avg);

	printf("max: %f, avg: %f\n", max, avg);
	
}

void computeResultStats(int n, double results[]) {
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += results[i];
	}
	double avg = sum / n;
	double var = 0.0;
	for (int i = 0; i < n; i++) {
		var += (avg - results[i]) * (avg - results[i]);
	}
	var /= n-1;
	double stdev = sqrt(var);
	printf("mean: %f, SD: %f\n", avg, stdev);
}

void initSim() {
	printf("init sim\n");

	root = new qNode(NULL, 0, 0);
	oldRoot = new qNode(NULL, 0, 0);

	int startLevel = levels - 2;
	if (startState == POISSONTEST || startState == ERRORTEST || startState == PROJECTTEST) {
		if (adaptScheme == ADAPTNONE)
			startLevel = levels - 1;
		else
			startLevel = 3;
	} else if (startState == ADAPTTEST) {
		//assert(levels > 3);
		//startLevel = std::max(levels - 4, levels/2 + 1);
		//startLevel = levels - 3;
		startLevel = 1;
	}


	initRecursive(root, startLevel);

	if (startState == POISSONTEST) {

		poissonReset(root);
		// adapt
		if (adaptScheme != ADAPTNONE) {
			//
			bool done = false;

			printf("doing poisson adapt first\n");
			
			int totalAdapts = 0;
			while (!done) {
				std::swap(root, oldRoot);
				copy(root, oldRoot);
				poissonReset(root);
				std::swap(root, oldRoot);
				done = !recursiveAdaptAndCopy(root, oldRoot);
				poissonReset(root);
				totalAdapts++;
			}
			printf("took %d adapts\n", totalAdapts);

			resetProfiling();
			root->profile();
			printProfiling();
		}
		// TODO remove

		root->computeVelocityDivergence();

		for (int i = 0; i < warmupRuns; i++) {
			runPoissonTest(false);
			printf("prerun %d\n", i+1);
		}

		if (numToRun == 0) numToRun++;
		double results[numToRun];
		for (int i = 0; i < numToRun; i++) {
		 results[i] = runPoissonTest(true);
		}
		computeResultStats(numToRun, results);
		numToRun = 0;
		checkPoissonError();
		return;
	} else if (startState == ADAPTTEST) {
		runAdaptTest();
		return;
	} else if (startState == ERRORTEST) {
		runErrorTest();
		return;
	} else if (startState == PROJECTTEST) {
		runProjectTest();
		return;
	}

	computeNodalVelocity(root);

	if (particleAlgorithm != PARTICLENONE) {

		particles = new Particle[numParticles];
		std::default_random_engine generator;
		double distmin = 0.0;
		double distmax = 1.0;

		if (startState == LEFT) {
    		int size = 1<<levelToDisplay;
			double h = 1.0 / size;
			distmin = 0.5 - 2 * h;
			distmax = 0.5 + 2 * h;
		}

		std::uniform_real_distribution<double> distribution(distmin, distmax);

		for (int i = 0; i < numParticles; i++) {

				particles[i].x = distribution(generator);
				particles[i].y = distribution(generator);
		}
	}

	// TODO MAKE MULTILEVEL
}


// input
void keyboard(unsigned char key, int x, int y) {
	switch (key) {
		case 32:     // space key
			runStep();
			if (!headless) {
				glutPostRedisplay();
			}
			break;
		case 'w':
			levelToDisplay = std::min(levels - 1, levelToDisplay + 1);
			if (!headless) {
				glutPostRedisplay();
			}
			break;
		case 's':
			levelToDisplay = std::max(0, levelToDisplay - 1);
			if (!headless) {
				glutPostRedisplay();
			}
			break;
		case 27:
			glutDestroyWindow(windowid);
			exit(0);
			break;
	}
}

void parseAdapt(char** argv, int i) {
	
}
 
/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {
	glutInit(&argc, argv);          // Initialize GLUT
	headless = false;
	levels = 6;
	water = false;
	drawCells = false;
	screenshot = false;
	particleAlgorithm = PARTICLENONE;
	adaptScheme = ADAPTNONE;
	adaptTestFunc = ADAPTFUNCNONE;
	poissonTestFunc = POISSONFUNCNONE;
	projectTestFunc = PROJECTFUNCNONE;
	startState = LEFT;
	for (int i = 1; i < argc; i++) {
		char* arg = argv[i];
		if (!strcmp("--headless", arg)) {
			headless = true;
		} else if (!strcmp("-levels", arg)) {
			levels= atoi(argv[++i]);
		} else if (!strcmp("-adapt", arg)) {
			char* scheme = argv[++i];
			if (!strcmp("curl", scheme)) {
				adaptScheme = CURL;
			} else if (!strcmp("pgrad", scheme)) {
				adaptScheme = PGRAD;
			} else if (!strcmp("vnorm", scheme)) {
				adaptScheme = VNORM;
			} else {
				printf("invalid adapt scheme %s\n", scheme);
			}
		} else if (!strcmp("--water", arg)) {
			water = true;
		} else if (!strcmp("--vel", arg)) {
			drawVelocity = true;
		} else if (!strcmp("--screenshot", arg)) {
			screenshot = true;
			statScreenshotDir();
		} else if (!strcmp("--cells", arg)) {
			drawCells = true;
		} else if (!strcmp("-run", arg)) {
            numToRun = atoi(argv[++i]);
        } else if (!strcmp("-particles", arg)) {
			char* alg = argv[++i];
            numParticles = atoi(argv[++i]);
			if (!strcmp("euler", alg)) {
				particleAlgorithm = EULER;
			} else {
				printf("invalid particle algorithm %s\n", alg);
				particleAlgorithm = PARTICLENONE;
				numParticles = 0;
			}
		// start states
		} else if (!strcmp("left", arg)) {
			startState = LEFT;
		} else if (!strcmp("sink", arg)) {
			startState = SINK;
		} else if (!strcmp("srcsink", arg)) {
			startState = SRCSINK;
		} else if (!strcmp("poissontest", arg)) {
			startState = POISSONTEST;
			char* func = argv[++i];
			if (!strcmp("xy", func)) {
				poissonTestFunc = POISSONXY;
			} else if (!strcmp("cos", func)) {
				poissonTestFunc = POISSONCOS;
			} else if (!strcmp("coskl", func)) {
				poissonTestFunc = POISSONCOSKL;
			} else {
				return 1;
			}
		} else if (!strcmp("adapttest", arg)) {
			startState = ADAPTTEST;
			char* func = argv[++i];
			if (!strcmp("circle", func)) {
				adaptTestFunc = ADAPTCIRCLE;
			} else if (!strcmp("wave", func)) {
				adaptTestFunc = ADAPTWAVE;
			} else {
				return 1;
			}
			adaptScheme = CURL;
		} else if (!strcmp("errortest", arg)) {
			startState = ERRORTEST;
			char* func = argv[++i];
			if (!strcmp("xy", func)) {
				errorTestFunc = ERRORXY;
			} else if (!strcmp("cos", func)) {
				errorTestFunc = ERRORCOS;
			} else if (!strcmp("coskl", func)) {
				errorTestFunc = ERRORCOSKL;
			} else {
				return 1;
			}
		} else if (!strcmp("projecttest", arg)) {
			startState = PROJECTTEST;
			char* func = argv[++i];
			if (!strcmp("xy", func)) {
				projectTestFunc = PROJECTXY;
			} else if (!strcmp("sin", func)) {
				projectTestFunc = PROJECTSIN;
			} else {
				return 1;
			}
		} else if (!strcmp("-pre", arg)) {
			warmupRuns = atoi(argv[++i]);
		} else if (!strcmp("-adaptscale", arg)) {
			adaptScale = atof(argv[++i]);
		} else if (!strcmp("-wavePeriod",arg)) {
			wavePeriod = atof(argv[++i]);
		} else if (!strcmp("-offsetY", arg)) {
			offsetY = atof(argv[++i]);
		} else if (!strcmp("-thresh", arg)) {
			thresh = atof(argv[++i]);
		} else if (!strcmp("-k", arg)) {
			poissonk = atoi(argv[++i]);
		} else if (!strcmp("-l", arg)) {	
			poissonl = atoi(argv[++i]);
		} else if (!strcmp("--pressureinterp", arg)) {
			pressureInterp = true;
		}
	}
	//levelToDisplay = levels/2;
    levelToDisplay = levels - 1;
	printf("headless: %d, levels: %d\n", headless, levels);

	// run tests
	//testNeighbors();
	//testMultilevelNeighbors();
	//testIntersect();
	//testPressureInterp();
	testNodalVelocity();

	initSim();
	
    printf("pre-running %d steps.\n", numToRun);
    for (int i = 0; i < numToRun; i++) {
        runStep();
    }

	// TODO don't do if headless
	glutInitWindowSize(windowWidth, windowHeight);   // Set the window's initial width & height
	glutInitWindowPosition(50, 50);
	windowid = glutCreateWindow("Quad Tree");  // Create window with the given title
	glutDisplayFunc(display);       // Register callback handler for window re-paint event
	glutKeyboardFunc(keyboard);
	initGL();                       // Our own OpenGL initialization
	glutMainLoop();                 // Enter the event-processing loop
	return 0;
}
