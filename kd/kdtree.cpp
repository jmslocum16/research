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
int totalLevels;
int levelToDisplay; // level to draw
bool water;
bool drawVelocity;
bool drawCells;
bool screenshot;
int numToRun = 0;
int warmupRuns = 0;
int MAX_LEVEL_DIF = 2;
bool pressureInterp = false;

// profiling
int numLeaves, numNodes;
int* leafLevels;
timeval start, end;


void resetProfiling() {
	if (leafLevels) delete leafLevels;
	leafLevels = new int[totalLevels+1];
	memset(leafLevels, 0, (totalLevels+1) * sizeof(int));
	numNodes = 0;
	numLeaves = 0;
}

void printProfiling() {
	printf("number of leaves at each level: \n");
	for (int i = 0; i <= totalLevels; i++) {
		printf("%d: %d\n", i, leafLevels[i]);
	}
	printf("total number of nodes: %d\n", numNodes);
	printf("total number of leaves: %d\n", numLeaves);
	printf("total number of leaves at non-max level: %d\n", numLeaves - leafLevels[totalLevels]);
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
// the dimension that is cut in half by the split
// vertical is split x, horizontal is split y
enum SplitDir { SPLIT_NONE, SPLIT_X, SPLIT_Y };
class kdNode {
	public:
		// tree stuffs
		int id;
		bool leaf;
		kdNode* parent;
		SplitDir splitDir;
		kdNode* neighbors[4];
		kdNode* children[2];
		int level_i, level_j, i, j; // level coordinates

		// math stuffs
		double p, dp, R, phi, divV, temp; // pressure ,x velocity, y velocity, pressure correction, residual, level set value
		double vx, vy, vx2, vy2;
		bool cvxValid[4];
		double cvx[4];
		bool cvyValid[4];
		double cvy[4];

		kdNode(kdNode *p, int li, int lj, int i, int j): parent(p), level_i(li), level_j(lj), i(i), j(j) {
			id = nodeId++;
			leaf = true;
			for (int i = 0; i < 4; i++) {
				neighbors[i] = NULL;
			}
			children[0] = NULL;
			children[1] = NULL;
			splitDir = SPLIT_NONE;
		}
		~kdNode() {
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
				for (int k = 0; k < 2; k++) {
					delete children[k];
				}
			}
		}

		// moving around tree functions
		void setChildrenNeighbors() {
			// set both ways because the first time a node is expanded the other one may not have been
			// deltas[k]= = {{1, 0}^, {-1, 0}v, {0, 1}>, {0, -1}<}
			assert (splitDir != SPLIT_NONE);

			kdNode* n = NULL;

			int size_i = 1<<level_i;
			int size_j = 1<<level_j;
			if (splitDir == SPLIT_X) {
				// set inside neighbors in x to each other
				children[0]->neighbors[2] = children[1];
				children[1]->neighbors[3] = children[0];

				// child 0 <
				if (j != 0) {
					n = this;
					while (n->neighbors[3] == NULL) n = n->parent;
					assert(n != NULL && n->neighbors[3] != NULL);
					children[0]->neighbors[3] = n->neighbors[3]->get(children[0]->level_i, children[0]->level_j, children[0]->i, children[0]->j-1, false);
					if (children[0]->neighbors[3] != NULL)
						children[0]->neighbors[3]->neighbors[2] = children[0];
				}
				// child 1 >
				if (j != size_j - 1) {
					n = this;
					while (n->neighbors[2] == NULL) n = n->parent;
					assert(n != NULL && n->neighbors[2] != NULL);

					children[1]->neighbors[2] = n->neighbors[2]->get(children[1]->level_i, children[1]->level_j, children[1]->i, children[1]->j+1, false);
					if (children[1]->neighbors[2] != NULL)
						children[1]->neighbors[2]->neighbors[3] = children[1];
				}

				// both v
				if (i != 0) {
					n = this;
					while (n->neighbors[1] == NULL) n = n->parent;
					assert (n != NULL && n->neighbors[1] != NULL);

					children[0]->neighbors[1] = n->neighbors[1]->get(children[0]->level_i, children[0]->level_j, children[0]->i-1, children[0]->j, false);
					if (children[0]->neighbors[1] != NULL)
						children[0]->neighbors[1]->neighbors[0] = children[0];

					children[1]->neighbors[1] = n->neighbors[1]->get(children[1]->level_i, children[1]->level_j, children[1]->i-1, children[1]->j, false);
					if (children[1]->neighbors[1] != NULL)
						children[1]->neighbors[1]->neighbors[0] = children[1];
				}

				// both ^
				if (i != size_i - 1) {
					n = this;
					while (n->neighbors[0] == NULL) n = n->parent;
					assert (n != NULL && n->neighbors[0] != NULL);

					children[0]->neighbors[0] = n->neighbors[0]->get(children[0]->level_i, children[0]->level_j, children[0]->i+1, children[0]->j, false);
					if (children[0]->neighbors[0] != NULL)
						children[0]->neighbors[0]->neighbors[1] = children[0];

					children[1]->neighbors[0] = n->neighbors[0]->get(children[1]->level_i, children[1]->level_j, children[1]->i+1, children[1]->j, false);
					if (children[1]->neighbors[0] != NULL)
						children[1]->neighbors[0]->neighbors[1] = children[1];

				}
			} else {
				// set inside neighbors in y to each other
				children[0]->neighbors[0] = children[1];
				children[1]->neighbors[1] = children[0];

				// child 0 v
				if (i != 0) {
					n = this;
					while (n->neighbors[1] == NULL) n = n->parent;
					assert (n != NULL && n->neighbors[1] != NULL);

					children[0]->neighbors[1] = n->neighbors[1]->get(children[0]->level_i, children[0]->level_j, children[0]->i-1, children[0]->j, false);
					if (children[0]->neighbors[1] != NULL)
						children[0]->neighbors[1]->neighbors[0] = children[0];
				}
				// child 1 ^
				if (i != size_i - 1) {
					n = this;
					while (n->neighbors[0] == NULL) n = n->parent;
					assert (n != NULL && n->neighbors[0] != NULL);

					children[1]->neighbors[0] = n->neighbors[0]->get(children[1]->level_i, children[1]->level_j, children[1]->i+1, children[1]->j, false);
					if (children[1]->neighbors[0] != NULL)
						children[1]->neighbors[0]->neighbors[1] = children[1];
				}

				// both <
				if (j != 0) {
					n = this;
					while (n->neighbors[3] == NULL) n = n->parent;
					assert (n != NULL && n->neighbors[3] != NULL);

					children[0]->neighbors[3] = n->neighbors[3]->get(children[0]->level_i, children[0]->level_j, children[0]->i, children[0]->j-1, false);
					if (children[0]->neighbors[3] != NULL)
						children[0]->neighbors[3]->neighbors[2] = children[0];

					children[1]->neighbors[3] = n->neighbors[3]->get(children[1]->level_i, children[1]->level_j, children[1]->i, children[1]->j-1, false);
					if (children[1]->neighbors[3] != NULL)
						children[1]->neighbors[3]->neighbors[2] = children[1];
				}

				// both >
				if (j != size_j - 1) {
					n = this;
					while (n->neighbors[2] == NULL) n = n->parent;
					assert (n != NULL && n->neighbors[2] != NULL);

					children[0]->neighbors[2] = n->neighbors[2]->get(children[0]->level_i, children[0]->level_j, children[0]->i, children[0]->j-1, false);
					if (children[0]->neighbors[2] != NULL)
						children[0]->neighbors[2]->neighbors[3] = children[0];

					children[1]->neighbors[2] = n->neighbors[2]->get(children[1]->level_i, children[1]->level_j, children[1]->i, children[1]->j-1, false);
					if (children[1]->neighbors[2] != NULL)
						children[1]->neighbors[2]->neighbors[3] = children[1];

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

		/*double getFaceGradient(int ml, int k, NodeValue v) {
			kdNode* n = this;
			int oppositeK = (k < 2) ? 1-k : 3-(k%2);//1 - (k%2);
			while (n != NULL && n->neighbors[k] == NULL) n = n->parent;
			if (n != NULL) {
				return n->neighbors[k]->addFaceToGradient(this, ml, oppositeK, v);
			} else {
				return 0;
			}
		}*/

		double getValueInterpCorrection(kdNode* original, int ml, int k, NodeValue v, int left) {
			int newk;
			double dpos;
			if (k < 2) {
				// previously computing gradient in y direction, need x direction for interpolation
				double origX = (original->j + 0.5)/(1<<original->level_j);
				double x = (j + 0.5)/(1<<level_j);
				dpos = origX - x;
				newk = dpos > 0 ? 2 : 3;
			} else {
				double origY = (original->i + 0.5)/(1<<original->level_i); 
				double y = (i + 0.5)/(1<<level_i);
				dpos = origY - y;
				newk = dpos > 0 ? 0 : 1;
			}
			double faceGrad = getFaceGradient(ml, newk, v, left-1);
			return fabs(dpos) * faceGrad;
		}

		kdNode* getChildInDir(int targetCoord, int targetLevel, int k) {
			if (k < 2) {
				// in y, going off x level + coord
				if (level_j == targetLevel) {
					assert (targetCoord ==  j);
					return this;
				}
				assert (level_j < targetLevel);
				if (leaf) {
					return this;
				} else if (splitDir == SPLIT_X) {
					// choose correct one
					int leveldif = targetLevel - level_j;
					int midcoord = (1<<leveldif)*j + (1<<(leveldif - 1));
					int index = (targetCoord < midcoord) ? 0 : 1;
					return children[index]->getChildInDir(targetCoord, targetLevel, k);
				} else {
					return children[1-(k%2)]->getChildInDir(targetCoord, targetLevel, k);
				}
			} else {
				// in x, going off y level + coord
				if (level_i == targetLevel) {
					assert (targetCoord ==  i);
					return this;
				}
				assert (level_i < targetLevel);
				if (leaf) {
					return this;
				} else if (splitDir == SPLIT_Y) {
					// choose correct one
					int leveldif = targetLevel - level_i;
					int midcoord = (1<<leveldif)*i + (1<<(leveldif - 1));
					int index = (targetCoord < midcoord) ? 0 : 1;
					return children[index]->getChildInDir(targetCoord, targetLevel, k);
				} else {
					return children[1-(k%2)]->getChildInDir(targetCoord, targetLevel, k);
				}

			}
		}

		kdNode* getNeighborInDir(int targetCoord, int targetLevel, int k) {
			kdNode* n = this;
			int oppositeK = (k < 2) ? 1-k : 3-(k%2);//1 - (k%2);
			while (n != NULL && n->neighbors[k] == NULL) n = n->parent;
			if (n != NULL) {
				int newi = i + deltas[k][0];
				int newj = j + deltas[k][1];
				n = n->neighbors[k]->getChildInDir(targetCoord, targetLevel, oppositeK);
			}
			return n;
		}

		double getFaceGradient(int ml, int k, NodeValue v) {
			return getFaceGradient(ml, k, v, levels);
		}

		double getFaceGradient(int ml, int k, NodeValue v, int left) {
			int oppositeK = (k < 2) ? 1-k : 3-(k%2);//1 - (k%2);
			int targetLevel = (k < 2 ? level_j : level_i);
			int targetCoord = (k < 2 ? j : i);
			kdNode* n = getNeighborInDir(targetCoord, targetLevel, k);
			if (n == NULL) {
				return 0;
			} else {
				return n->addFaceToGradient(this, ml, oppositeK, v, left);
			}

		}

		double addFaceToGradient(kdNode* original, int ml, int k, NodeValue v, int left) {
			if (leaf || (level_i+level_j) == ml) {
				int horig = 1<< (k < 2 ? original->level_i : original->level_j);
				int hn = 1<< (k < 2 ? level_i : level_j);
				double h = 0.5/horig + 0.5/hn;

				//printf("Adding node with pressure %f to face gradient of node with pressure %f\n", getVal(v), original->getVal(v));

				int origLevel = (k < 2 ? original->level_j : original->level_i);
				int thisLevel = (k < 2 ? level_j : level_i);
				double val = getVal(v);
				if (pressureInterp && origLevel > thisLevel && left > 0) {
					val += getValueInterpCorrection(original, ml, k, v, left);
				}

				return (val - original->getVal(v))/h;
			} else {
				double total = 0.0;
				if ((k < 2 && splitDir == SPLIT_X) || (k >= 2 && splitDir == SPLIT_Y)) {
					total += children[0]->addFaceToGradient(original, ml, k, v, left);
					total += children[1]->addFaceToGradient(original, ml, k, v, left);
					total /= 2.0;
				} else {
					total += children[1-(k%2)]->addFaceToGradient(original, ml, k, v, left);
				}
				return total;
			}
		}

		void getLaplacian(int ml, double *aSum, double* bSum, NodeValue v) {
			for (int k = 0; k < 4; k++) {
				int oppositeK = (k < 2) ? 1-k : 3-(k%2); //1 - (k%2);
				int targetLevel = (k < 2 ? level_j : level_i);
				int targetCoord = (k < 2 ? j : i);
				kdNode* n = getNeighborInDir(targetCoord, targetLevel, k);
			
				if (n != NULL) {
					
					int newi = i + deltas[k][0];
					int newj = j + deltas[k][1];
					n->addFaceToLaplacian(this, ml, oppositeK, aSum, bSum, v);
				} else {
					*aSum -= 1;
					*bSum += getVal(v);
				}
			}
		}

		// note: k is opposite of way used to get there, for example if k in original was (1, 0), this k is (-1, 0)
		void addFaceToLaplacian(kdNode* original, int ml, int k, double* aSum, double* bSum, NodeValue v) {
			if (leaf || (level_i + level_j) == ml) {
				int horig = 1<< (k < 2 ? original->level_i : original->level_j);
				int hn = 1<< (k < 2 ? level_i : level_j);
				double h = 0.5/horig + 0.5/hn;

				double worig = 1 << (k < 2 ? original->level_j : original->level_i);
				double wn = 1 << (k < 2 ? level_j : level_i);
				double dside = fmin(1.0/worig, 1.0/wn);

				double d = dside/h;

				int origLevel = (k < 2 ? original->level_j : original->level_i);
				int thisLevel = (k < 2 ? level_j : level_i);
				double val = getVal(v);
				if (pressureInterp && origLevel > thisLevel) {
					val += getValueInterpCorrection(original, ml, k, v, levels);
				}


				*aSum -= d;
				*bSum += d*val;
				//printf("adding node with pressure %f to laplacian of node with pressure %f\n", getVal(v), original->getVal(v));
			} else {
				if ((k < 2 && splitDir == SPLIT_X) || (k >= 2 && splitDir == SPLIT_Y)) {
					children[0]->addFaceToLaplacian(original, ml, k, aSum, bSum, v);
					children[1]->addFaceToLaplacian(original, ml, k, aSum, bSum, v);
				} else {
					children[1-(k%2)]->addFaceToLaplacian(original, ml, k, aSum, bSum, v);
				}
			}
		}

		// utility functions
		std::pair<double, double> getValueGradient(NodeValue v) {
			
			int ml = leaf ? totalLevels : level_i + level_j;
			std::pair<double, double> newGrad =  std::make_pair((getFaceGradient(ml, 2, v) - getFaceGradient(ml, 3, v))/2.0, (getFaceGradient(ml, 0, v) - getFaceGradient(ml, 1, v))/2.0);

			return newGrad;
		}
		
		// curl(F) = d(Fy)/dx - d(Fx)/dy
		double getCurl() {
			double y1 = (cvy[0] + cvy[2])/2.0;
			double y2 = (cvy[1] + cvy[3])/2.0;
			double x1 = (cvx[0] + cvx[1])/2.0;
			double x2 = (cvx[2] + cvx[3])/2.0;
			int size_i = 1<<level_i;
			int size_j = 1<<level_j;
			return  -(x2-x1)/size_j + (y2-y1)/size_i;	
		}

		void computeVelocityDivergence() {
			//std::pair<double, double> grad = getVelocityGradient();
			//divV = grad.first + grad.second;
			if (leaf) {
				int size_i = 1<<level_i;
				int size_j = 1<<level_j;
				double a, b;

				// x
				if (j == 0)
					a = 0.0;
				else
					a = vx;
				if (j == size_j - 1)
					b = 0.0;
				else {
					b = vx2;
				}
				divV = (b-a) * size_j;

				// y
				if (i == 0)
					a = 0.0;
				else
					a = vy;
				if (i == size_i-1)
					b = 0.0;
				else {
					b = vy2;
				}
				divV += (b-a)*size_i;
			} else {
				divV = 0.0;
				for (int k = 0; k < 2; k++) {
					children[k]->computeVelocityDivergence();
				}
			}
		}

		std::pair<double, double> getVelocityAt(kdNode* r, double x, double y) {
			int size_i = 1<<level_i;
			int size_j = 1<<level_j;
			double minX = ((float)j)/size_j;
			double minY = ((float)i)/size_i;
			assert (!(x < minX || y < minY || x > minX + 1.0/size_j || y > minY + 1.0/size_i));

			double dj = (x*size_j)-j;
			double di = (y*size_i)-i;
			
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
			int size_i = 1<<level_i;
			int size_j = 1<<level_j;
			double dx = x - (j + 0.5)/size_j;
			double dy = y - (i + 0.5)/size_i;
			if (fabs(dx) > 0.5/size_j || fabs(dy) > 0.5/size_i) {
				printf("dx: %f, dy: %f\n", dx, dy);
			}
			return val + dx * grad.first + dy * grad.second;
		}

		// adaptivity
		void expand(bool calculateNewVals, SplitDir dir) {
			assert (leaf);
			assert (dir != SPLIT_NONE);
	
			leaf = false;
			splitDir = dir;

			// TODO
    		//printf("expanding cell %d, (%d, %d)\n", level, i, j);
    		int size_i = 1<<level_i;
    		int size_j = 1<<level_j;
    		// std::pair<double	, double> levelSetGrad;
			if (dir == SPLIT_X) {
	    		for (int k = 0; k < 2; k++) {
					children[k] = new kdNode(this, level_i, level_j + 1, i, 2*j + k);
    			}
				if (calculateNewVals) {
					children[0]->cvx[0] = cvx[0];
					children[0]->cvy[0] = cvy[0];
					children[0]->cvx[2] = cvx[2];
					children[0]->cvy[2] = cvy[2];
					children[1]->cvx[1] = cvx[1];
					children[1]->cvy[1] = cvy[1];
					children[1]->cvx[3] = cvx[3];
					children[1]->cvy[3] = cvy[3];
					double newv = (cvx[0] + cvx[1])/2.0;
					children[0]->cvx[1] = newv;
					children[1]->cvx[0] = newv;
					newv = (cvy[0] + cvy[1])/2.0;
					children[0]->cvy[1] = newv;
					children[1]->cvy[0] = newv;
					newv = (cvx[2] + cvx[3])/2.0;
					children[0]->cvx[3] = newv;
					children[1]->cvx[2] = newv;
					newv = (cvy[2] + cvy[3])/2.0;
					children[0]->cvy[3] = newv;
					children[1]->cvy[2] = newv;
				}
			} else {
				for (int k = 0; k < 2; k++) {
					children[k] = new kdNode(this, level_i + 1, level_j, 2*i + k, j);	
				}
				if (calculateNewVals) {
					children[0]->cvx[0] = cvx[0];
					children[0]->cvy[0] = cvy[0];
					children[0]->cvx[1] = cvx[1];
					children[0]->cvy[1] = cvy[1];
					children[1]->cvx[2] = cvx[2];
					children[1]->cvy[2] = cvy[2];
					children[1]->cvx[3] = cvx[3];
					children[1]->cvy[3] = cvy[3];
					double newv = (cvx[0] + cvx[2])/2.0;
					children[0]->cvx[2] = newv;
					children[1]->cvx[0] = newv;
					newv = (cvy[0] + cvy[2])/2.0;
					children[0]->cvy[2] = newv;
					children[1]->cvy[0] = newv;
					newv = (cvx[1] + cvx[3])/2.0;
					children[0]->cvx[3] = newv;
					children[1]->cvx[1] = newv;
					newv = (cvy[1] + cvy[3])/2.0;
					children[0]->cvy[3] = newv;
					children[1]->cvy[1] = newv;
				}
			}
			if (calculateNewVals) {
				for (int k = 0; k < 2; k++) {
					children[k]->p = p;
					children[k]->vx = (cvx[0] + cvx[2])/2.0;
					children[k]->vx2 = (cvx[1] + cvx[3])/2.0;
					children[k]->vy = (cvy[0] + cvy[1])/2.0;
					children[k]->vy2 = (cvy[2] + cvy[3])/2.0;
				}
			}
			setChildrenNeighbors();
		}

		// don't need to average values since adapt function takes care of that
		void contract() {
			assert (!leaf);
			assert (splitDir != SPLIT_NONE);
			for (int k = 0; k < 2; k++) {
				assert(children[k]->leaf);
			}
		    //printf("contracting cell %d (%d, %d)\n", level, i, j);
			for (int k = 0; k < 2; k++) {
				delete children[k];
			}
			leaf = true;
			splitDir = SPLIT_NONE;
		}
		void profile() {
			numNodes++;
			if (leaf) {
				numLeaves++;
				leafLevels[level_i+level_j]++;
			} else {
				for (int k = 0; k < 2; k++) {
					children[k]->profile();
				}
			}
		}

		// current node, target li/lj/ni/nj
		
		kdNode* get(int li, int lj, int ni, int nj, bool allowEarly) {
			if (li == level_i && lj == level_j && ni == i && nj == j) {
				return this;
			}
			if (leaf) {
				return allowEarly ? this : NULL;
			}
			assert (splitDir != SPLIT_NONE);
			if (splitDir == SPLIT_X) {
				int leveldif = lj - level_j;
				int midsize = 1<<(leveldif - 1);
				int midj = j*(1<<leveldif) + midsize;
				int index = (nj < midj) ? 0 : 1;
				return children[index]->get(li, lj, ni, nj, allowEarly);
			} else {
				int leveldif = li - level_i;
				int midsize = 1<<(leveldif - 1);
				int midi = i*(1<<leveldif) + midsize;
				int index = (ni < midi) ? 0 : 1;
				return children[index]->get(li, lj, ni, nj, allowEarly);
			}
		}	
};


kdNode* root;
kdNode* oldRoot;


// debugging

kdNode* getLeaf(kdNode* node, double x, double y, int ml) {
	// TODO
	assert (false);
	/*
	int size = 1<<node->level;
	assert (x >= (node->j + 0.0)/size && x <= (node->j + 1.0)/size);
	assert (y >= (node->i + 0.0)/size && y <= (node->i + 1.0)/size);
	if (node->leaf || node->level == ml) {
		return node;
	}
	double midx = (node->j + 0.5)/size;
	double midy = (node->i + 0.5)/size;
	int newj = x < midx ? 0 : 1;
	int newi = y < midy ? 0 : 1;
	return getLeaf(node->children[2*newi + newj], x, y, ml);
	*/
}

/*kdNode*  getSemiLagrangianLookback(kdNode* r, double* x, double* y, int steps, int ml) {
	double newdt = dt / steps;
	kdNode* cur = getLeaf(r, *x, *y, ml);
	while (steps--) {
		std::pair<double, double> vel = cur->getVelocityAt(r, *x, *y);
		*x -= vel.first * newdt;
		*y -= vel.second * newdt;
		*x = fmin(1.0, fmax(0, *x));
		*y = fmin(1.0, fmax(0, *y));
		cur = getLeaf(r, *x, *y, ml);
	}
	return cur;
}*/

// not fast but who cares
void printValue(NodeValue v) {
	int size = 1<<(levels - 1);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			kdNode* n = root->get(levels, levels, i, j, true);
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

void drawMultilevel(kdNode* node, int ml) {
	assert (node->level_i + node->level_j <= ml);
	int size_i = 1 << node->level_i;
	int size_j = 1 << node->level_j;
	if (node->level_i + node->level_j == ml || node->leaf) {
		// draw it
		double x = -1.0 + (2.0 * node->j)/size_j;
		double y = -1.0 + (2.0 * node->i)/size_i;

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

		glRectf(x, y, x + 2.0/size_j, y + 2.0/size_i);

		if (drawVelocity) {
			if (maxMag >= eps) {
				glColor3f(0.0, 1.0, 0.0);
				double vx = node->vx;
				double vy = node->vy;
				double mag = sqrt(vx*vx + vy*vy);
				if (mag >= eps) {
					glBegin(GL_LINES);
					// max size is 1.0/side length, scaled by the max magnitude

					double scale = maxMag * std::min(size_j, size_i);

					// vx
					glVertex2f(x, y + 1.0/size_i);
					glVertex2f(x + vx / scale, y + 1.0/size_i);

					// vy
					glVertex2f(x + 1.0/size_j, y);
					glVertex2f(x + 1.0/size_j, y + vy / scale);

					glEnd();
				}
			}
		}
		// draw cell
		if (drawCells) {
			glColor3f(1.0, 1.0, 1.0);
			glBegin(GL_LINE_LOOP);
			glVertex2f(x, y);
			glVertex2f(x, y + 2.0/size_i);
			glVertex2f(x + 2.0/size_j, y + 2.0/size_i);
			glVertex2f(x + 2.0/size_j, y);
			glEnd();
		}
	} else {
		// draw children
		for (int k = 0; k < 2; k++) {
			drawMultilevel(node->children[k], ml);
		}
	}
}

void computeMax(kdNode* node, int ml) {
	if (node->leaf || node->level_i + node->level_j == ml) {
		if (!water) {
			minP = std::min(minP, node->p);
			maxP = std::max(maxP, node->p);
		}

		if (drawVelocity) {
			//maxMag = std::max(maxMag, sqrt(node->vx*node->vx + node->vy*node->vy));
			maxMag = std::max(maxMag, node->vx);
			maxMag = std::max(maxMag, node->vy);
		}
	} else {
		for (int k = 0; k < 2; k++) {
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
void computeNodalVelocity(kdNode* root) {
	//root->invalidateNodalValues();
	// loop to try to figure out if done
	int oldTotalLeft;
	int totalLeft = 10000000;
	while (totalLeft > 0) {
		oldTotalLeft = totalLeft;
		//totalLeft = root->attemptComputeNodalVals(root);
		assert (totalLeft < oldTotalLeft);
	}
}


/*void testNodalVelocity() {
	int oldLevels = levels;

	levels = 3;

	root = new kdNode(NULL, 0, 0, 0, 0);
	root->expand(false, SPLIT_X);
	root->children[0]->expand(false, SPLIT_Y);
	root->children[1]->expand(false, SPLIT_Y);
	root->children[0]->children[0]->expand(false, SPLIT_X);
	root->children[0]->children[0]->children[1]->expand(false, SPLIT_Y);
	root->children[0]->children[1]->expand(false, SPLIT_Y);
	root->children[0]->children[1]->children[0]->expand(false, SPLIT_Y);
	root->children[1]->children[1]->expand(false, SPLIT_X);

	root->children[0]->vx = 0.0;
	root->children[0]->vx2 = 26.0/3.0;
	root->children[0]->vy = 0.0;
	root->children[0]->vy2 = 0.0;

	root->children[1]->vx = 26.0/3.0;
	root->children[1]->vx2 = 0.0;
	root->children[1]->vy = 0.0;
	root->children[1]->vy2 = 0.0;

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
	assertf(24, root->children[0]->cvx[3], 44.0/3.0);
	assertf(25, root->children[1]->cvx[2], 44.0/3.0);
	assertf(26, root->children[2]->cvx[1], 44.0/3.0);
	assertf(27, root->children[3]->cvx[0], 44.0/3.0);
	assertf(28, root->children[0]->children[3]->cvx[3], 44.0/3.0);
	assertf(29, root->children[0]->children[3]->children[3]->cvx[3], 44.0/3.0);
	assertf(30, root->children[3]->children[0]->cvx[0], 44.0/3.0);
	assertf(31, root->children[0]->cvy[3], 46.0/3.0);
	assertf(32, root->children[1]->cvy[2], 46.0/3.0);
	assertf(33, root->children[2]->cvy[1], 46.0/3.0);
	assertf(34, root->children[3]->cvy[0], 46.0/3.0);
	assertf(35, root->children[0]->children[3]->cvy[3], 46.0/3.0);
	assertf(36, root->children[0]->children[3]->children[3]->cvy[3], 46.0/3.0);
	assertf(37, root->children[3]->children[0]->cvy[0], 46.0/3.0);

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
	assertf(49, root->children[0]->children[1]->cvy[3], 6.40625);
	assertf(50, root->children[0]->children[3]->cvy[1], 6.40625);
	assertf(51, root->children[0]->children[3]->children[1]->cvy[1], 6.40625);

	assertf(52, root->children[0]->children[3]->cvx[2], 9.5);
	assertf(53, root->children[0]->children[3]->children[2]->cvx[2], 9.5);
	assertf(54, root->children[0]->children[2]->cvx[3], 9.5);
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
	assertf(68, root->children[3]->children[0]->cvx[1], 145.0/12.0);
	assertf(69, root->children[3]->children[1]->cvx[0], 145.0/12.0);
	assertf(70, root->children[3]->children[0]->cvy[1], 16.5);
	assertf(71, root->children[3]->children[1]->cvy[0], 16.5);

	assertf(72, root->children[3]->children[0]->cvy[2], 325.0/24.0);
	assertf(73, root->children[3]->children[2]->cvy[0], 325.0/24.0);
	assertf(74, root->children[3]->children[0]->cvx[2], 20);
	assertf(75, root->children[3]->children[2]->cvx[0], 20);

	assertf(76, root->children[3]->children[1]->cvx[3], 0.0);
	assertf(77, root->children[3]->children[3]->cvx[1], 0.0);
	assertf(78, root->children[3]->children[1]->cvy[3], 21);
	assertf(79, root->children[3]->children[3]->cvy[1], 21);

	assertf(80, root->children[3]->children[3]->cvx[2], 23);
	assertf(81, root->children[3]->children[2]->cvx[3], 23);
	assertf(82, root->children[3]->children[3]->cvy[2], 0.0);
	assertf(83, root->children[3]->children[2]->cvy[3], 0.0);

	assertf(84, root->children[3]->children[0]->cvx[3], 21);
	assertf(85, root->children[3]->children[0]->cvy[3], 20.5);	
	assertf(86, root->children[3]->children[1]->cvx[2], 21);	
	assertf(87, root->children[3]->children[1]->cvy[2], 20.5);	
	assertf(88, root->children[3]->children[2]->cvx[1], 21);
	assertf(89, root->children[3]->children[2]->cvy[1], 20.5);	
	assertf(90, root->children[3]->children[3]->cvx[0], 21);	
	assertf(91, root->children[3]->children[3]->cvy[0], 20.5);	

	///
	assertf(92, root->children[0]->children[3]->children[0]->cvx[1], 31.0/6.0);
	assertf(93, root->children[0]->children[3]->children[1]->cvx[0], 31.0/6.0);
	assertf(94, root->children[0]->children[3]->children[0]->cvy[1], 4.5);
	assertf(95, root->children[0]->children[3]->children[1]->cvy[0], 4.5);

	assertf(96, root->children[0]->children[3]->children[0]->cvy[2], 25.0/3.0);
	assertf(97, root->children[0]->children[3]->children[2]->cvy[0], 25.0/3.0);
	assertf(98, root->children[0]->children[3]->children[0]->cvx[2], 8.5);
	assertf(99, root->children[0]->children[3]->children[2]->cvx[0], 8.5);

	assertf(100, root->children[0]->children[3]->children[1]->cvx[3], 10.5);
	assertf(101, root->children[0]->children[3]->children[3]->cvx[1], 10.5);
	assertf(102, root->children[0]->children[3]->children[1]->cvy[3], 10.859375);
	assertf(103, root->children[0]->children[3]->children[3]->cvy[1], 10.859375);

	assertf(104, root->children[0]->children[3]->children[3]->cvx[2], 12.0);
	assertf(105, root->children[0]->children[3]->children[2]->cvx[3], 12.0);
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
}*/



// tests neighbor pointers
void testNeighbors() {
	// TODO
	//kdNode* testRoot = new kdNode(NULL, 0, 0);
	//testRoot->expand(false);

}


void testAdvect() {
	// TODO finish later..
	int oldLevels = levels;
	int oldTotalLevels = totalLevels;
	double oldDt = dt;

	dt = 0.0;
	levels = 2;
	totalLevels = 4;

	root = new kdNode(NULL, 0, 0, 0, 0);
	oldRoot = new kdNode(NULL, 0, 0, 0, 0);

	delete root;
	delete oldRoot;
	levels = oldLevels;
	totalLevels = oldTotalLevels;
	dt = oldDt;
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
	int oldTotalLevels = totalLevels;
	bool oldPressureInterp = pressureInterp;
	levels = 2;
	totalLevels = 4;
	pressureInterp = true;

	kdNode* testRoot = new kdNode(NULL, 0, 0, 0, 0);
	testRoot->expand(false, SPLIT_X);
	testRoot->children[0]->expand(false, SPLIT_Y);
	testRoot->children[1]->expand(false, SPLIT_Y);
	testRoot->children[0]->children[0]->expand(false, SPLIT_Y);
	testRoot->children[1]->children[0]->expand(false, SPLIT_X);

	testRoot->children[0]->children[1]->p = 5;
	testRoot->children[1]->children[1]->p = 6;
	testRoot->children[0]->children[0]->children[0]->p = 1;
	testRoot->children[0]->children[0]->children[1]->p = 2;
	testRoot->children[1]->children[0]->children[0]->p = 3;
	testRoot->children[1]->children[0]->children[1]->p = 4;

	assertf(0, testRoot->children[0]->children[1]->getFaceGradient(3, 0, P), 0.0);
	assertf(1, testRoot->children[0]->children[1]->getFaceGradient(3, 1, P), -8.0);
	assertf(2, testRoot->children[0]->children[1]->getFaceGradient(3, 2, P), 2.0);
	assertf(3, testRoot->children[0]->children[1]->getFaceGradient(3, 3, P), 0.0);
	
	assertf(4, testRoot->children[1]->children[1]->getFaceGradient(3, 0, P), 0.0);
	assertf(5, testRoot->children[1]->children[1]->getFaceGradient(3, 1, P), -5.0);
	assertf(6, testRoot->children[1]->children[1]->getFaceGradient(3, 2, P), 0.0);
	assertf(7, testRoot->children[1]->children[1]->getFaceGradient(3, 3, P), -2.0);
	
	assertf(8, testRoot->children[0]->children[0]->children[0]->getFaceGradient(3, 0, P), 4.0);
	assertf(9, testRoot->children[0]->children[0]->children[0]->getFaceGradient(3, 1, P), 0.0);
	assertf(10, testRoot->children[0]->children[0]->children[0]->getFaceGradient(3, 2, P), 16.0/3.0);
	assertf(11, testRoot->children[0]->children[0]->children[0]->getFaceGradient(3, 3, P), 0.0);
		
	assertf(12, testRoot->children[0]->children[0]->children[1]->getFaceGradient(3, 0, P), 8.0);
	assertf(13, testRoot->children[0]->children[0]->children[1]->getFaceGradient(3, 1, P), -4.0);
	assertf(14, testRoot->children[0]->children[0]->children[1]->getFaceGradient(3, 2, P), 9.0/2.0);
	assertf(15, testRoot->children[0]->children[0]->children[1]->getFaceGradient(3, 3, P), 0.0);
	
	assertf(16, testRoot->children[1]->children[0]->children[0]->getFaceGradient(3, 0, P), 11.0/2.0);
	assertf(17, testRoot->children[1]->children[0]->children[0]->getFaceGradient(3, 1, P), 0.0);
	assertf(18, testRoot->children[1]->children[0]->children[0]->getFaceGradient(3, 2, P), 4.0);
	assertf(19, testRoot->children[1]->children[0]->children[0]->getFaceGradient(3, 3, P), -4.0);
		
	assertf(20, testRoot->children[1]->children[0]->children[1]->getFaceGradient(3, 0, P), 4.0);
	assertf(21, testRoot->children[1]->children[0]->children[1]->getFaceGradient(3, 1, P), 0.0);
	assertf(22, testRoot->children[1]->children[0]->children[1]->getFaceGradient(3, 2, P), 0.0);
	assertf(23, testRoot->children[1]->children[0]->children[1]->getFaceGradient(3, 3, P), -4.0);
	
	levels = oldLevels;
	totalLevels = oldTotalLevels;
	pressureInterp = oldPressureInterp;
	delete testRoot;
}



void testGradLap() {
	kdNode* testRoot = new kdNode(NULL, 0, 0, 0, 0);
	testRoot->expand(false, SPLIT_X);
	testRoot->children[0]->expand(false, SPLIT_X);
	testRoot->children[1]->expand(false, SPLIT_Y);
	testRoot->children[1]->children[0]->expand(false, SPLIT_Y);
	testRoot->children[1]->children[1]->expand(false, SPLIT_X);

	testRoot->children[0]->children[0]->p = 1.0;
	testRoot->children[0]->children[1]->p = 2.0;
	testRoot->children[1]->children[0]->children[0]->p = 3.0;
	testRoot->children[1]->children[0]->children[1]->p = 4.0;
	testRoot->children[1]->children[1]->children[0]->p = 5.0;
	testRoot->children[1]->children[1]->children[1]->p = 6.0;

	// gradient tests
	assertf (0, testRoot->children[0]->children[0]->getFaceGradient(3, 0, P), 0.0);
	assertf (1, testRoot->children[0]->children[0]->getFaceGradient(3, 1, P), 0.0);
	assertf (2, testRoot->children[0]->children[0]->getFaceGradient(3, 2, P), 4.0);
	assertf (3, testRoot->children[0]->children[0]->getFaceGradient(3, 3, P), 0.0);

	assertf (4, testRoot->children[0]->children[1]->getFaceGradient(3, 0, P), 0.0);
	assertf (5, testRoot->children[0]->children[1]->getFaceGradient(3, 1, P), 0.0);
	assertf (6, testRoot->children[0]->children[1]->getFaceGradient(3, 2, P), 8.0);//?
	assertf (7, testRoot->children[0]->children[1]->getFaceGradient(3, 3, P), -4.0);

	assertf (8, testRoot->children[1]->children[0]->children[0]->getFaceGradient(3, 0, P), 4.0);
	assertf (9, testRoot->children[1]->children[0]->children[0]->getFaceGradient(3, 1, P), 0.0);
	assertf (10, testRoot->children[1]->children[0]->children[0]->getFaceGradient(3, 2, P), 0.0);//?
	assertf (11, testRoot->children[1]->children[0]->children[0]->getFaceGradient(3, 3, P), -8.0/3.0);

	assertf (12, testRoot->children[1]->children[0]->children[1]->getFaceGradient(3, 0, P), 4.0);
	assertf (13, testRoot->children[1]->children[0]->children[1]->getFaceGradient(3, 1, P), -4.0);
	assertf (14, testRoot->children[1]->children[0]->children[1]->getFaceGradient(3, 2, P), 0.0);
	assertf (15, testRoot->children[1]->children[0]->children[1]->getFaceGradient(3, 3, P), -16.0/3.0);

	assertf (16, testRoot->children[1]->children[1]->children[0]->getFaceGradient(3, 0, P), 0.0);
	assertf (17, testRoot->children[1]->children[1]->children[0]->getFaceGradient(3, 1, P), -8.0/3.0);
	assertf (18, testRoot->children[1]->children[1]->children[0]->getFaceGradient(3, 2, P), 4.0);
	assertf (19, testRoot->children[1]->children[1]->children[0]->getFaceGradient(3, 3, P), -12.0);

	assertf (20, testRoot->children[1]->children[1]->children[1]->getFaceGradient(3, 0, P), 0.0);
	assertf (21, testRoot->children[1]->children[1]->children[1]->getFaceGradient(3, 1, P), -16.0/3.0);
	assertf (22, testRoot->children[1]->children[1]->children[1]->getFaceGradient(3, 2, P), 0.0);
	assertf (23, testRoot->children[1]->children[1]->children[1]->getFaceGradient(3, 3, P), -4.0);

	// laplacian tests

	double aSum, bSum, laplacian;
	kdNode* n;

	aSum = 0.0;
	bSum = 0.0;
	n = testRoot->children[0]->children[0];
	n->getLaplacian(3, &aSum, &bSum, P);
	laplacian = (aSum * n->p + bSum) * (1<<n->level_i) * (1<<n->level_j);
	assertf (24, laplacian, 16);

	aSum = 0.0;
	bSum = 0.0;
	n = testRoot->children[0]->children[1];
	n->getLaplacian(3, &aSum, &bSum, P);
	laplacian = (aSum * n->p + bSum) * (1<<n->level_i) * (1<<n->level_j);
	assertf (25, laplacian, 16);

	aSum = 0.0;
	bSum = 0.0;
	n = testRoot->children[1]->children[0]->children[0];
	n->getLaplacian(3, &aSum, &bSum, P);
	laplacian = (aSum * n->p + bSum) *  (1<<n->level_i) * (1<<n->level_j);
	assertf (26, laplacian, 32.0/3.0);

	aSum = 0.0;
	bSum = 0.0;
	n = testRoot->children[1]->children[0]->children[1];
	n->getLaplacian(3, &aSum, &bSum, P);
	laplacian = (aSum * n->p + bSum) *  (1<<n->level_i) * (1<<n->level_j);
	assertf (27, laplacian, -32.0/3.0);

	aSum = 0.0;
	bSum = 0.0;
	n = testRoot->children[1]->children[1]->children[0];
	n->getLaplacian(3, &aSum, &bSum, P);
	laplacian = (aSum * n->p + bSum) *  (1<<n->level_i) * (1<<n->level_j);
	assertf (28, laplacian, -14.0*8.0/3.0);

	aSum = 0.0;
	bSum = 0.0;
	n = testRoot->children[1]->children[1]->children[1];
	n->getLaplacian(3, &aSum, &bSum, P);
	laplacian = (aSum * n->p + bSum) *  (1<<n->level_i) * (1<<n->level_j);
	assertf (29, laplacian, -80.0/3.0);

	// cleanup
	delete testRoot;
}

// poisson solver functions

// returns max(abs(R))
double computeResidual(kdNode* node) {
	int size_i = 1<<node->level_i;
	int size_j = 1<<node->level_j;
	if (node->leaf) { // if leaf cell, compute residual
		// compute it: R = div(vx, vy) - 1/(ha)*sum of (s * grad) for each face
		double aSum = 0.0;
		double bSum = 0.0;
		node->getLaplacian(totalLevels, &aSum, &bSum, P);
		double laplacian = (aSum * node->p + bSum) * size_i * size_j;
		
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
		for (int k = 0; k < 2; k++) {
			maxR = std::max(maxR, computeResidual(node->children[k]));
			node->R += node->children[k]->R;
		}
		node->R /= 2.0;
		return maxR;
	}
}


// relax the node at the given multilevel
// puts the new value of dp in temp
bool relaxRecursive(kdNode* node, int ml) {
	assert (node->level_i + node->level_j <= ml);
	int size_i = 1<<node->level_i;
	int size_j = 1<<node->level_j;
	if (node->level_i + node->level_j == ml || node->leaf) {
		double aSum = 0.0;
        double bSum = 0.0;
        double dp;
        double h;
		double oldDif = node->dp;

		node->getLaplacian(ml, &aSum, &bSum, DP);
		
        // A*R = bSum - aSum*dp, dp = (bSum - A*R)/asum, A = h*h = 1/size^2
		node->temp = -(bSum - node->R/size_i/size_j)/aSum;

		double diff = oldDif - node->temp;

        //printf("relaxing[%d][%d][%d]: aSum: %f, bSumL %f, R: %f, result: %f, diff from old: %f\n", d, i, j, aSum, bSum, grid[d][i*size+j].R, grid[d][i*size+j].temp, diff);
		if (fabs(diff) > eps) {
			return false;
		}
		return true;
	} else {
		// relax children
		bool done = true;
		for (int k = 0; k < 2; k++) {
			done &= relaxRecursive(node->children[k], ml);
		}
		return done;
	}
}

void recursiveInject(kdNode* node, int ml) {
	if (node->level_i + node->level_j  == ml) {
		node->dp = node->parent->dp;
	}
	if (!node->leaf) {
		for (int k = 0; k < 2; k++) {
			recursiveInject(node->children[k], ml);
		}
	}
}

void recursiveUpdate(kdNode* node, int ml) {	
	node->dp = node->temp;
	if (!node->leaf) {
		for (int k = 0; k < 2; k++) {
			recursiveUpdate(node->children[k], ml);
		}
	}
}

void relax(int d, int r) {
	assert (d <= totalLevels);

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

double curlAdaptFunction(kdNode* node) {
	double curl = node->getCurl();
	return fabs(curl) > thresh;
}

double vnormAdaptFunction(kdNode* node) {
	int size_i = 1<<node->level_i;
	int size_j = 1<<node->level_j;
	double vnorm = fabs(node->vx2-node->vx)/size_j + fabs(node->vy2-node->vy)/size_i;
	return vnorm > thresh;
}

bool pGradAdaptFunction(kdNode* node) {
	assert (false);
	//std::pair<double, double> pgrad = node->getValueGradient(P);
    //return fabs(pgrad.first + pgrad.second) > pressureThresh/(1<<node->level);
}

bool adaptFunction(kdNode* node) {
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

SplitDir velAdaptDirFunction(kdNode* node) {
	// dimension with bigger velocity gradient
	double dvdx = fabs(node->vx2 - node->vx)*(1<<node->level_j);
	double dvdy = fabs(node->vy2 - node->vy)*(1<<node->level_i);
	if (dvdx > dvdy) {
		return SPLIT_X;
	} else {
		return SPLIT_Y;
	}
}

SplitDir pGradAdaptDirFunction(kdNode* node) {
	assert(false);
}	

SplitDir adaptDirFunction(kdNode* node) {
	assert (adaptScheme != ADAPTNONE);
	assert (node->level_i + node->level_j < totalLevels);
	if (node->level_i == levels || node->level_i - MAX_LEVEL_DIF == node->level_j) {
		return SPLIT_X;
	} else if (node->level_j == levels || node->level_j - MAX_LEVEL_DIF == node->level_i) {
		return SPLIT_Y;
	}
	if (adaptScheme == PGRAD) {
		return pGradAdaptDirFunction(node);
	} else if (adaptScheme == CURL || adaptScheme == VNORM) {
		return velAdaptDirFunction(node);
	}
}

// trees should have the same structure
// use old tree for adapt calculations
// copies new pressure over
// returns true if any nodes were changed, false otherwise
bool recursiveAdaptAndCopy(kdNode* node, kdNode* oldNode) {
	// TODO implement properly with nodal/face velocities
	if (node->level_i + node->level_j == totalLevels) return false;
	assert (node->leaf == oldNode->leaf);
	node->p = oldNode->p;
	if (node->leaf) {
		if (adaptFunction(oldNode)) {
			SplitDir dir = adaptDirFunction(oldNode);
			double oldP = oldNode->p;

			oldNode->expand(true, dir);
			node->expand(false, dir);
			for (int k = 0; k < 2; k++) {
				node->children[k]->p = oldNode->children[k]->p;
				node->children[k]->vx = oldNode->children[k]->vx;
				node->children[k]->vy =  oldNode->children[k]->vy;
				node->children[k]->vx2 = oldNode->children[k]->vx2;
				node->children[k]->vy2 = oldNode->children[k]->vy2;
				 // TODO corner values
			}
			// reset old node to old state so it is the same as before, so it can be used in other calculations
			oldNode->contract();
			oldNode->p = oldP;

			// now node is adapted, with values from oldnode's calculations
			return true;
		}
		return false;
	}
	bool allChildrenLeaves = true;
    for (int k = 0; k < 2; k++) {
        if (!node->children[k]->leaf) {
            allChildrenLeaves = false;
            break;
        }
    }
	if (allChildrenLeaves && !adaptFunction(oldNode)) {
        // this wouldn't be expanded if it was a leaf, so it shouldn't have leaf children
		node->contract();
		node->p = 0.0;
		for (int k = 0; k < 2; k++) {
			node->p += oldNode->children[k]->p;
		}
		node->p /= 2.0;
		return true;
    } else {
		bool anyChanged = false;
        for (int k = 0; k < 2; k++) {
            anyChanged |= recursiveAdaptAndCopy(node->children[k], oldNode->children[k]);
        }
		return anyChanged;
    }
	
}

void copy(kdNode* node, kdNode* oldNode) {
	// need to create/delete new nodes if it adapted
	if (node->leaf && !oldNode->leaf) {
		// old node expanded, expand new node
		node->expand(false, oldNode->splitDir);
	} else if (!node->leaf && oldNode->leaf) {
		// old node contracted
		node->contract(); // TODO this doesn't handle multiple contractions per step, but rn only do single contract/expand per step
	}
	
	node->p = oldNode->p;

	if (!node->leaf) {
		for (int k = 0; k < 2; k++) {
			copy(node->children[k], oldNode->children[k]);
		}
	}
}

void setNewAdvect(kdNode* node) {
	assert (false);
	/*if (node->leaf) {
		if (node->i == 0 || node->j == 0) {
			node->cvx = 0.0;
			node->cvy = 0.0;
			return;
		}
		int size = 1<<node->level;
		double x = (node->j + 0.5)/size;
		double y = (node->i + 0.5)/size;
		kdNode* last= getSemiLagrangianLookback(oldRoot, &x, &y, 1, node->level);
		// TODO implement
		std::pair<double, double> newvel = last->getVelocityAt(oldRoot, x, y);
		node->cvx = newvel.first;
		node->cvy = newvel.second;
	} else {
		for (int k = 0; k < 4; k++) {
			setNewAdvect(node->children[k]);
		}
		node->cvx = node->children[0]->cvx;
		node->cvy = node->children[0]->cvy;
	}*/
}


// TODO put on node class
/*void setNewFace(kdNode* node) {
	assert (false);
	// TODO this is innacurate at T-junctions on edge?
	std::pair<double, double> c00 = node->getNodalAt(root, 0, 0);
	std::pair<double, double> c01 = node->getNodalAt(root, 0, 1);
	std::pair<double, double> c10 = node->getNodalAt(root, 1, 0);
	std::pair<double, double> c11 = node->getNodalAt(root, 1, 1);
	node->vx = (c00.first + c10.first)/2.0;
	node->vy = (c00.second + c01.second)/2.0;
	node->vx2 = (c01.first + c11.first)/2.0;
	node->vy2 = (c10.second + c11.second)/2.0;

	if (!node->leaf) {
		for (int k = 0; k < 2; k++) {
			setNewFace(node->children[k]);
		}
	}
}*/

void advectAndCopy() {

	copy(root, oldRoot);

	// assume nodal values have been computed on old tree
	// set nodal values on new tree to advected nodal values on old tree
	//setNewAdvect(root);

	// set new face values from nodal values
	//setNewFace(root);

}

void correctPressure(kdNode* node) {
	//node->p += node->dp;
	if (node->leaf) {
		//printf("at node %d (%d, %d) correcting pressure %f by %f\n", node->level, node->i, node->j, node->p, node->dp);
		node->p += node->dp;
	} else {
		for (int k = 0; k < 2; k++) {
			correctPressure(node->children[k]);
		}
	}
}

void project(kdNode* node) {
	// correct velocity with updated pressure field to make non-divergent
	if (node->leaf) {
		//if (node->i == 0 || node->j == 0) return;
		double pgradX, pgradY;
		int size_i = 1<<node->level_i;
		int size_j = 1<<node->level_j;
		if (node->j > 0) {
			pgradX = -node->getFaceGradient(levels - 1, 3, P);
			node->vx -= pgradX;
		}
		if (node->i > 0) {
			pgradY = -node->getFaceGradient(levels - 1, 1, P);
			node->vy -= pgradY;
		}
		if (node->j < size_j - 1) {
			pgradX = node->getFaceGradient(levels - 1, 2, P);
			node->vx2 -= pgradX;
		}
		if (node->i < size_i - 1) {
			pgradY = node->getFaceGradient(levels - 1, 0, P);
			node->vy2 -= pgradY;
		}
		
	} else {
		for (int k = 0; k < 2; k++) {
			project(node->children[k]);
		}

		assert (node->splitDir != SPLIT_NONE);
		if (node->splitDir == SPLIT_X) {
			node->vx = node->children[0]->vx;
			node->vx2 = node->children[1]->vx2;
			node->vy = (node->children[0]->vy + node->children[1]->vy)/2.0;
			node->vy2 = (node->children[0]->vy2 + node->children[1]->vy2)/2.0;	
		} else {
			node->vy = node->children[0]->vy;
			node->vy2 = node->children[1]->vy;
			node->vx = (node->children[0]->vx + node->children[1]->vx)/2.0;
			node->vx2 = (node->children[0]->vx2 + node->children[1]->vx2)/2.0;
			
		}
	}
}

void clampVelocity(kdNode* node) {
	assert (false);
	// clamp velocity
	/*int size = 1<<node->level;
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
	}*/
}

void eulerAdvectParticle(Particle& p, std::pair<double, double> v) {
	p.x += v.first * dt;
	p.y += v.second * dt;
}

void advectParticles() {
	for (int i = 0; i < numParticles; i++) {
		if (particleAlgorithm == EULER) {
		}
	}
}

double runVCycle() {
	//relax(0, MAX_RELAX_COUNT);
	root->dp = 0;
	// do not relax level 0 since it makes no sense to do so...

	for (int d = 1; d <= totalLevels; d++) {
		relax(d, (d==1) ? MAX_RELAX_COUNT : OPTIMIZED_RELAX_COUNT);
	}
	/*relax(1, MAX_RELAX_COUNT);
	for (int d = 2; d <= totalLevels; d++) {
		if (d%2 == 0 || d == totalLevels)
			relax(d, OPTIMIZED_RELAX_COUNT);
		else
			recursiveInject(root, d);
	}*/

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

void poissonAverageR(kdNode* node, double* total) {
	if (node->leaf) {
		int size_i = 1<<node->level_i;
		int size_j = 1<<node->level_j;
		*total += node->R / size_i / size_j;
	} else {
		for (int k = 0; k < 2; k++) {
			poissonAverageR(node->children[k], total);
		}
	}
}

void poissonCorrectR(kdNode* node, double K) {
	node->R -= K;
	if (!node->leaf) {
		for (int k = 0; k < 2; k++) {
			poissonCorrectR(node->children[k], K);
		}
	}
}

double getMaxR(kdNode* node) {
	if (node->leaf)	return fabs(node->R);
	double max = 0.0;
	for (int k = 0; k < 2; k++)
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

void projectionCheck(kdNode* node, double* max, double* avg) {
	if (node->leaf) {
		//printf("divV after projection: %f\n", node->divV);
		*max = fmax(*max, fabs(node->divV));
		int size_i = 1<<node->level_i;
		int size_j = 1<<node->level_j;
		*avg += fabs(node->divV)/size_i/size_j;
	} else {
		for (int k = 0; k < 2; k++) {
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

	//clampVelocity(root);
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

void initNodeLeftUniform(kdNode* node) {
	int size_i = 1<<node->level_i;
	//node->p = 1.0/size/size;
	node->p = 0.0;
	node->vx = 0.0;
	node->vy = 0.0;
	node->phi = 0.0;
	if (node->i == size_i/2 || node->i == size_i/2-1) {
		node->vx = -1.0;
	}	
}

void initNodeSink(kdNode* node) {
	int size_i = 1<<node->level_i;
	int size_j = 1<<node->level_j;
	node->p = 1.0/size_j/size_i;
	node->phi = 0.0;

	double center = 0.5;
	node->vx = center - (node->j + 0.5)/size_j;
	node->vy = center - (node->i + 0.5)/size_i;
	double len = sqrt(node->vx * node->vx + node->vy * node->vy);
	node->vx /= len;
	node->vy /= len;
}

void initNodeSrcSink(kdNode* node) {
	int size_i = 1<<node->level_i;
	int size_j = 1<<node->level_j;
	node->p = 1.0/size_i/size_j;
	node->phi = 0.0;
	
	double src = .125;
	double sink = .875;

	double vxsrc = (node->j + 0.5)/size_j - src;
	double vysrc = (node->i + 0.5)/size_i - src;
	double lensrc = sqrt(vxsrc * vxsrc + vysrc * vysrc);

	double vxsink = sink - (node->j + 0.5)/size_j;
	double vysink = sink - (node->i + 0.5)/size_i;
	double lensink = sqrt(vxsink * vxsink + vysink * vysink);
	
	node->vx = vxsrc/lensrc + vxsink/lensink;
	node->vy = vysrc/lensrc + vysink/lensink;
}

void initPoissonTest(kdNode* node) {
	
	node->phi = 0;
	node->vx = 0;
	node->vy = 0;
	
	//int size = 1<<node->level;
	//double x = (node->j + 0.5)/size;
	//double y = (node->i + 0.5)/size;

	node->p = 0;
	
}

void initProjectTest(kdNode* node) {
	int size_i = 1<<node->level_i;
	int size_j = 1<<node->level_j;
	node->p = 0.0;
	node->dp = 0.0;

	assert (projectTestFunc != PROJECTFUNCNONE);

	double x = ((float)node->j)/size_j;
	double y = ((float)node->i)/size_i;

	if (projectTestFunc == PROJECTXY) {
		node->vx = (node->j+0.5)*(node->j+0.5)/size_j/size_j;
		node->vy = 1-(node->i+0.5)*(node->i+0.5)/size_i/size_i;
	
		node->vx2 = (node->j+1.5)*(node->j+1.5)/size_j/size_j;
		node->vy2 = 1-(node->i+1.5)*(node->i+1.5)/size_i/size_i;
	} else if (projectTestFunc == PROJECTSIN) {
		// vx
		double x2 = x;
		double y2 = y + 0.5/size_i;
		node->vx = 2*M_PI*sin(2*M_PI*x2) - 2*M_PI*sin(2*M_PI*x2)*cos(2*M_PI*y2);

		x2 = x+1.0/size_j;
		node->vx2 = 2*M_PI*sin(2*M_PI*x2) - 2*M_PI*sin(2*M_PI*x2)*cos(2*M_PI*y2);

		//vy
		x2 = x + 0.5/size_j;
		y2 = y;
		node->vy = 2*M_PI*sin(2*M_PI*y2) - 2*M_PI*cos(2*M_PI*x2)*sin(2*M_PI*y2);

		y2 = y + 1.0/size_i;
		node->vy2 = 2*M_PI*sin(2*M_PI*y2) - 2*M_PI*cos(2*M_PI*x2)*sin(2*M_PI*y2);	
	}
}

void initNodeFunction(kdNode* node) {
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
void initRecursive(kdNode* node, int d) {
	if (node->level_i + node->level_j == d) {
		initNodeFunction(node);
	} else {
		SplitDir dir = (node->level_i == node->level_j) ? SPLIT_X : SPLIT_Y;
		node->expand(false, dir);

		node->p = 0.0;
        node->vx = 0.0;
        node->vy = 0.0;
        //node->phi = 0.0;

        for (int k = 0; k < 2; k++) {
			initRecursive(node->children[k], d);
			node->p += node->children[k]->p;
			//node->vx += node->children[k]->vx;
			//node->vy += node->children[k]->vy;
			//node->phi += node->children[k]->phi;
			// TODO do velocity properly here
        }

        node->p /= 4.0;
        //node->vx /= 4.0;
        //node->vy /= 4.0;
        //node->phi /= 4.0;
	}
}

int poissonk = 3;
int poissonl = 3;

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

void poissonReset(kdNode* node) {
	node->p = 0;
	int size_i = 1<<node->level_i;
	int size_j = 1<<node->level_j;
	double x = (node->j + 0.0)/size_j;
	double y = (node->i + 0.0)/size_i;

	// vel
	double x2 = x;
	double y2 = y + 0.5/size_i;
	node->vx = getPoissonVX(x2, y2);
	x2 = x + 1.0/size_j;
	node->vx2 = getPoissonVX(x2, y2);

	x2 = x + 0.5/size_j;
	y2 = y;
	node->vy = getPoissonVY(x2, y2);
	y2 = y + 1.0/size_i;
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
		for (int k = 0; k < 2; k++) {
			poissonReset(node->children[k]);
		}
	}
}

void computePoissonError(kdNode* node, double* total) {
	if (node->leaf) {
		int size_i = 1<<node->level_i;
		int size_j = 1<<node->level_j;
		double x = (node->j + 0.5)/size_j;
		double y = (node->i + 0.5)/size_i;
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
		*total += fabs(node->p - correct)/size_i/size_j;
	} else {
		for (int k = 0; k < 2; k++) {
			computePoissonError(node->children[k], total);
		}
	}
}

void poissonAverage(kdNode* node, double* total) {
	if (node->leaf) {
		int size_i = 1<<node->level_i;
		int size_j = 1<<node->level_j;
		*total += node->p / size_i / size_j;
	} else {
		for (int k = 0; k < 2; k++) {
			poissonAverage(node->children[k], total);
		}
	}
}

void expandRadius(kdNode* node, double radius) {
	if (node->leaf) {
		int size_j = 1<<node->level_j;
		int size_i = 1<<node->level_i;
		double x = (node->j + 0.5)/size_j;
		double y = (node->i + 0.5)/size_i;
		double dist = sqrt((0.5-x)*(0.5-x) + (0.5-y)*(0.5-y));
		if (dist < radius) {
			SplitDir dir;
			if (node->level_i == node->level_j) {
				dir = SPLIT_X;
			} else {
				dir = SPLIT_Y;
			}
			node->expand(false, dir);
		}
	} else {
		for (int k = 0; k < 2; k++) {
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
		//printf("residual after %d vcycles: %f\n", i, newR);
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
		doneVCycle = getMaxR(root) < eps;
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
void checkPoissonErrorRecursive(kdNode* node, double* maxG, double* avgG, double* maxL, double* avgL) {
	if (node->leaf) {
		int size_i = 1<<node->level_i;
		int size_j = 1<<node->level_j;
		for (int k = 0; k < 4; k++) {
			double vel;
			if (k == 0) vel = node->vy2;
			else if (k == 1) vel = -node->vy;
			else if (k == 2) vel = node->vx2;
			else vel = -node->vx;
			double grad = node->getFaceGradient(totalLevels, k, P);
			double error = fabs(grad - vel);
			//printf("d: %d, velocity: %f, pressure gradient :%f, error: %f\n", k, vel, grad, error);
			*maxG = fmax(*maxG, error);
			*avgG += error/ size_i/size_j/4;
		}
		double aSum = 0.0;
		double bSum = 0.0;
		node->getLaplacian(totalLevels, &aSum, &bSum, P);
		double laplacian = (aSum * node->p + bSum) * size_i * size_j;
		double lap2 = 0.0;
		double error = fabs(laplacian - node->divV);
		*maxL = fmax(*maxL, error);
		*avgL += error / size_i/size_j;
		
	} else {
		for (int k = 0; k < 2; k++) {
			checkPoissonErrorRecursive(node->children[k], maxG, avgG, maxL, avgL);
		}
	}
}

void checkPoissonError() {
	// TODO pass values to sum
	double maxG = 0.0;
	double avgG = 0.0;
	double maxL = 0.0;
	double avgL = 0.0;
	checkPoissonErrorRecursive(root, &maxG, &avgG, &maxL, &avgL);	

	printf("gradient error. max: %f, avg: %f\n", maxG, avgG);
	printf("laplacian error. max: %f, avg: %f\n", maxL, avgL);
}

void setErrorPressure(kdNode* node) {
	int size_i = 1<<node->level_i;
	int size_j = 1<<node->level_j;
	double x = (node->j + 0.5)/size_j;
	double y = (node->i + 0.5)/size_i;

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
		for (int k = 0; k < 2; k++) {
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

double calculateError(kdNode* node, double* avgError) {
	if (node->leaf) {
		int size_i = 1<<node->level_i;
		int size_j = 1<<node->level_j;
		double max = 0.0;
		for (int k = 0; k < 4; k++) {
			// calculate gradient in direction k

			// 4/4: (4, 7) in dir 3
			if (node->level_i == 4 && node->level_j == 4 && node->i ==4 && node->j == 7 && k == 3) {
				printf("asdfasdfsafa\n");
			}
			int ml = totalLevels;
			double calc = node->getFaceGradient(ml, k, P);
			ml = node->level_i + node->level_j;
			
			double x = (node->j + 0.5 + 0.5*deltas[k][1])/size_j;
			double y = (node->i + 0.5 + 0.5*deltas[k][0])/size_i;
			double real = getRealDerivative(x, y, k);
			if (k % 2 == 1) {
				// switch directions because other way
				real = -real;
			}
			double error = fabs(real - calc);
			printf("Error for node %d/%d: (%d, %d) in dir %d, ", node->level_i, node->level_j, node->i, node->j, k);
			if (ml == node->level_i + node->level_j) {
				printf(" at the same level.");
			} else if (ml < node->level_i + node->level_j) {
				printf(" up %d levels.", node->level_i + node->level_j - ml);
			} else {
				printf(" down %d levels.", ml - node->level_i - node->level_j);
			}
			printf("real: %f, calc: %f, error: %f\n", real, calc, error);
			*avgError += error/size_i/size_j/4;
			max = fmax(max, error);
		}
		return max;
	} else {
		double maxError = 0.0;
		for (int k = 0; k < 2; k++) {
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

void adaptTestRecursive(kdNode* node) {
	node->p = 0.0;
	if (node->level_i + node->level_j == totalLevels) {
		return;
	}
	if (node->leaf) {
		int size_i = 1<<node->level_i;
		int size_j = 1<<node->level_j;
		double x1, y1, x2, y2;
		x1 = ((float)node->j)/size_j;
		y1 = ((float)node->i)/size_i;
		x2 = x1 + 1.0/size_j;
		y2 = y1 + 1.0/size_i;
		if (rectIntersectFunc(x1, y1, x2, y2)) {
			// determine direction of adaptation	

			SplitDir dir;
			if (node->level_i == levels || node->level_i - MAX_LEVEL_DIF == node->level_j) {
				dir =  SPLIT_X;
			} else if (node->level_j == levels || node->level_j - MAX_LEVEL_DIF == node->level_i) {
				dir = SPLIT_Y;
			} else {
				// x
				double midx = (x1+x2)/2.0;
				bool bothX = rectIntersectFunc(x1, y1, midx, y2) && rectIntersectFunc(midx, y1, x2, y2);
				double midy = (y1+y2)/2.0;
				bool bothY = rectIntersectFunc(x1, y1, x2, midy) && rectIntersectFunc(x1, midy, x2, y2); 

				if ((bothX && bothY) || (!bothX && !bothY)) {
					// both will end up being split next time, choice doesn't matter
					if (node->level_i > node->level_j) {
						dir = SPLIT_Y;
					} else {
						dir = SPLIT_X;
					}
				} else if (bothX) {
					dir = SPLIT_Y;
				} else if (bothY) {
					dir = SPLIT_X;
				} else {
					assert (false);
				}
				
			}
			node->expand(false, dir);
			for (int k = 0; k < 2; k++) {
				adaptTestRecursive(node->children[k]);
			}
		}
	} else {
		for (int k = 0; k < 2; k++) {
			adaptTestRecursive(node->children[k]);
		}
	}
}

void runAdaptTest() {	
	bool done = false;
	int numCycles = -1;

	//root = new kdNode(NULL, 0, 0);

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

	root = new kdNode(NULL, 0, 0, 0, 0);
	oldRoot = new kdNode(NULL, 0, 0, 0, 0);

	int startLevel = totalLevels - 2;
	if (startState == POISSONTEST || startState == ERRORTEST || startState == PROJECTTEST) {
		if (adaptScheme == ADAPTNONE)
			startLevel = totalLevels;
		else
			startLevel = 6;
	} else if (startState == ADAPTTEST) {
		startLevel = 2;
	}

	initRecursive(root, startLevel);

	//computeNodalValues(root);
	
	if (startState == POISSONTEST) {
		poissonReset(root);
		// adapt
		if (adaptScheme != ADAPTNONE) {
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

		root->computeVelocityDivergence();

		for (int i = 0; i < warmupRuns; i++) {
			runPoissonTest(false);
			printf("prerun %d\n", i+1);
		}

		if (numToRun == 0) numToRun++;
		double times[numToRun];
		for (int i = 0; i < numToRun; i++) {
			times[i] = runPoissonTest(true);
		}
		checkPoissonError();
		computeResultStats(numToRun, times);
		numToRun = 0;
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

	if (particleAlgorithm != PARTICLENONE) {

		particles = new Particle[numParticles];
		std::default_random_engine generator;
		double distmin = 0.0;
		double distmax = 1.0;

		/*if (startState == LEFT) {
    		int size = 1<<levelToDisplay;
			double h = 1.0 / size;
			distmin = 0.5 - 2 * h;
			distmax = 0.5 + 2 * h;
		}*/

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
			levelToDisplay = std::min(totalLevels, levelToDisplay + 1);
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

 
/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {
	glutInit(&argc, argv);          // Initialize GLUT
	headless = false;
	levels = 6;
	totalLevels = 12;
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
			levels = atoi(argv[++i]);
			totalLevels = 2 * levels;
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
				return 1;
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
		} else if (!strcmp("-maxleveldif", arg)) {
			MAX_LEVEL_DIF = atoi(argv[++i]);
			assert (MAX_LEVEL_DIF > 0);
		} else if (!strcmp("-adaptscale", arg)) {
			adaptScale = atof(argv[++i]);
		} else if (!strcmp("-wavePeriod",arg)) {
			wavePeriod = atof(argv[++i]);
		} else if (!strcmp("-offsetY", arg)) {
			offsetY = atof(argv[++i]);
		} else if (!strcmp("-thresh", arg)) {
			thresh = atof(argv[++i]);
		} else if (!strcmp("-eps", arg)) {
			eps = atof(argv[++i]);
		} else if (!strcmp("-k", arg)) {
			poissonk = atoi(argv[++i]);
		} else if (!strcmp("-l", arg)) {
			poissonl = atoi(argv[++i]);
		} else if (!strcmp("--pressureinterp", arg)) {
			pressureInterp = true;
		}
	}
	//levelToDisplay = levels/2;
    levelToDisplay = totalLevels;
	printf("headless: %d, levels: %d\n", headless, levels);

	// run tests
	//testNeighbors();
	//testMultilevelNeighbors();
	//testGradLap();
	testPressureInterp();

	initSim();
	
    printf("pre-running %d steps.\n", numToRun);
    for (int i = 0; i < numToRun; i++) {
        runStep();
    }

	// TODO don't do if headless
	glutInitWindowSize(windowWidth, windowHeight);   // Set the window's initial width & height
	glutInitWindowPosition(50, 50);
	windowid = glutCreateWindow("K-D Tree");  // Create window with the given title
	glutDisplayFunc(display);       // Register callback handler for window re-paint event
	glutKeyboardFunc(keyboard);
	initGL();                       // Our own OpenGL initialization
	glutMainLoop();                 // Enter the event-processing loop
	return 0;
}
