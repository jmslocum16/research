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
enum AdaptScheme { ADAPTNONE, CURL, PGRAD };
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
		//double cvx[4];
		//double cvy[4];
		double cvx, cvy;

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
				return (getVal(v) - original->getVal(v)) * d;
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
				*aSum -= d;
				*bSum += d*getVal(v);
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
					//qNode* n = this;
					//while (n->neighbors[2] == NULL) n = n->parent;
					//b = n->neighbors[2]->getOtherVelocityFace(true);
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
					//qNode* n = this;
					//while (n->neighbors[0] == NULL) n = n->parent;
					//b = n->neighbors[0]->getOtherVelocityFace(false);
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

		std::pair<double, double> getNodal(int nl, int ni, int nj) {
			int nsize = 1<<nl;
			double x, y;
			if (nj == 0 || nj == nsize) {
				x = 0.0;
			} else if (ni == 0 || ni ==nsize) {
				//x = 2*vx - get(nl, ni+1, nj)->getNodal(nl, ni+1, nj);
				x = vx;
			} else {
				x = cvx;
			}
			if (ni == 0 || ni == nsize) {
				y = 0.0;
			} else if (nj == 0 || nj == nsize) {
				y = vy;
			} else {
				y = cvy;
			}
			return std::make_pair(x, y);
		}

		std::pair<double, double> getNodalAt(qNode* r, int i, int j) {
			if (i == 0 && j == 0) {
				return getNodal(level, i, j);
			} else if (i == 0 && j == 0) {
				return r->get(level, i, j+1)->getNodal(level, i, j+1);
			} else if (i == 1 && j == 0) {
				return r->get(level, i+1, j)->getNodal(level, i+1, j);
			} else if (i == i && j == 1) {
				return r->get(level, i+1, j+1)->getNodal(level, i+1, j+1);
			} else {
				assert(false);
			}
		}

		std::pair<double, double> getVelocityAt(qNode* r, double x, double y) {
			//return std::make_pair(getValueAt(x, y, VX), getValueAt(x, y, VY));
			// TODO
			int size = 1<<level;
			double minX = ((float)j)/size;
			double minY = ((float)i)/size;
			assert (!(x < minX || y < minY || x > minX + 1.0/size || y > minY + 1.0/size));

			double dj = (x*size)-j;
			double di = (y*size)-i;
			
			std::pair<double, double> c00 = getNodalAt(r, 0, 0);
			std::pair<double, double> c01 = getNodalAt(r, 0, 1);
			std::pair<double, double> c10 = getNodalAt(r, 1, 0);
			std::pair<double, double> c11 = getNodalAt(r, 1, 1);

			// x
			double newvx = bilinearInterpolation(c00.first, c01.first, c10.first, c11.first, di, dj);
			double newvy = bilinearInterpolation(c00.second, c01.second, c10.second, c11.second, di, dj);	
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
			// TODO
    		//printf("expanding cell %d, (%d, %d)\n", level, i, j);
    		int size = 1<<level;
    		// std::pair<double, double> levelSetGrad;
    		for (int k = 0; k < 4; k++) {

				this->children[k] = new qNode(this, 2*i+(k/2), 2*j+(k%2));

				if (calculateNewVals) {
					//assert (false);
					/*
	        		int constX = k%2 == 0 ? -1 : 1;
	        		int constY = k/2 == 0 ? -1 : 1;
					double newx = (j + 0.5)/size + (constX*.25)/size;
					double newy = (i + 0.5)/size + (constY*.25)/size;
	
					//children[k]->p = getPressureAt(newx, newy);
					children[k]->p = getValueAt(newx, newy, P);
					// TODO incorporate full velocity gradients?
					std::pair<double, double> newVel = getVelocityAt(newx, newy);
        			children[k]->vx = newVel.first;
        			children[k]->vy = newVel.second;
        			// children[k]->phi = ...
					*/
				}
    		}
			setChildrenNeighbors();
			leaf = false;
		}

		// don't need to average values since adapt function takes care of that
		void contract() {
			assert(!leaf);
			for (int k = 0; k < 4; k++) {
				assert(children[k]->leaf);
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
			/*if (node->level == d) {
				assert(node->i == i);
				assert(node->j == j);
				return node;
			} else if (node->leaf) return NULL;
			*/
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

qNode*  getSemiLagrangianLookback(qNode* r, double* x, double* y, int steps, int ml) {
	double newdt = dt / steps;
	qNode* cur = getLeaf(r, *x, *y, ml);
	while (steps--) {
		std::pair<double, double> vel = cur->getVelocityAt(r, *x, *y);
		*x -= vel.first * newdt;
		*y -= vel.second * newdt;
		*x = fmin(1.0, fmax(0, *x));
		*y = fmin(1.0, fmax(0, *y));
		cur = getLeaf(r, *x, *y, ml);
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

double adaptScale = 0.5;
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
	return .5 + adaptScale*.5*cos(wavePeriod*2*M_PI*x);
}

bool rectIntersectFuncWave(double x1, double y1, double x2, double y2) {
	return false;
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
}

void drawWaveSubdivide(double x1, double x2) {
	double y1 = getWaveFunc(x1);
	double y2 = getWaveFunc(x2);
	double slope = (y2-y1)/(x2-x1);
	if (x2-x1 < .01) {
		glBegin(GL_LINES);
		glVertex2f(-1+2*x1, -1+2*y1 + 2*offsetY);
		glVertex2f(-1+2*x2, -1+2*y2 + 2*offsetY);
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


void testAdvect() {
	int oldLevels = levels;
	double oldDt = dt;

	dt = 0.0;
	levels = 2;

	root = new qNode(NULL, 0, 0);
	oldRoot = new qNode(NULL, 0, 0);
	root->expand(false);
	root->children[0]->expand(false);
	oldRoot->expand(false);
	oldRoot->children[0]->expand(false);

	oldRoot->children[0]->vx = 0.0;
	oldRoot->children[0]->vy = 0.0;
	oldRoot->children[0]->vx2 = 1.0;
	oldRoot->children[0]->vy2 = 2.0;

	oldRoot->children[1]->vx = 1.0;
	oldRoot->children[1]->vy = 0.0;
	oldRoot->children[1]->vx2 = 0.0;
	oldRoot->children[1]->vy2 = 4.0;

	oldRoot->children[2]->vx = 0.0;
	oldRoot->children[2]->vy = 2.0;
	oldRoot->children[2]->vx2 = 3.0;
	oldRoot->children[2]->vy2 = 0.0;

	oldRoot->children[3]->vx = 3.0;
	oldRoot->children[3]->vy = 4.0;
	oldRoot->children[3]->vx2 = 0.0;
	oldRoot->children[3]->vy2 = 0.0;

	// next level
	oldRoot->children[0]->children[0]->vx = 0.0;
	oldRoot->children[0]->children[0]->vy = 0.0;
	oldRoot->children[0]->children[0]->vx2 = -1;
	oldRoot->children[0]->children[0]->vy2 = 8;

	delete root;
	delete oldRoot;
	levels = oldLevels;
	dt = oldDt;
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

double curlThresh = 0.1;
double curlAdaptFunction(qNode* node) {
	double curl = node->getCurl();
	return fabs(curl) > curlThresh / (1<<node->level);
}

double pressureThresh = .01;
bool pGradAdaptFunction(qNode* node) {
    //std::pair<double, double> pgrad = node->getPressureGradient();
	std::pair<double, double> pgrad = node->getValueGradient(P);
    return fabs(pgrad.first + pgrad.second) > pressureThresh/(1<<node->level);
}

bool adaptFunction(qNode* node) {
	assert(adaptScheme != ADAPTNONE);
	if (adaptScheme == PGRAD) {
    	return pGradAdaptFunction(node);
	} else { // CURL
		return curlAdaptFunction(node);
	}
}

// trees should have the same structure
// use old tree for adapt calculations
// copies new pressure over
// returns true if any nodes were changed, false otherwise
bool recursiveAdaptAndCopy(qNode* node, qNode* oldNode) {
	// TODO implement properly with nodal/face velocities
	if (node->level == levels - 1) return false;
	assert (node->leaf == oldNode->leaf);
	node->p = oldNode->p;
	if (node->leaf) {
		if (adaptFunction(oldNode)) {
			double oldP = oldNode->p;
			double oldVx = oldNode->vx;
			double oldVy = oldNode->vy;

			oldNode->expand(true);
			node->expand(false);
			for (int k = 0; k < 4; k++) {
				node->children[k]->p = oldNode->children[k]->p;
				node->children[k]->vx = oldNode->children[k]->vx;
				node->children[k]->vy =  oldNode->children[k]->vy;
			}
			// reset old node to old state so it is the same as before, so it can be used in other calculations
			oldNode->contract();
			oldNode->p = oldP;
			oldNode->vx = oldVx;
			oldNode->vy = oldVy;

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
		node->contract();
		node->p = 0.0;
		node->vx = 0.0;
		node->vy = 0.0;
		for (int k = 0; k < 4; k++) {
			node->p += oldNode->children[k]->p;
			node->vx += oldNode->children[k]->vx;
			node->vy += oldNode->children[k]->vy;
		}
		node->p /= 4.0;
		node->vx /= 4.0;
		node->vy /= 4.0;
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
		node->contract(); // TODO this doesn't handle multiple contractions per step, but rn only do single contract/expand per step
	}
	
	node->p = oldNode->p;

	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			copy(node->children[k], oldNode->children[k]);
		}
	}
}

// assumes face values set
void computeNodalValues(qNode* node) {
	if (node->leaf) {
		if (node->i == 0 || node->j == 0) {
			node->cvx = 0.0;
			node->cvy = 0.0;
			return;
		}
		//x
		qNode* n = node;
		while (n->neighbors[3] == NULL) n = n->parent;
		node->cvx = (n->vx + n->neighbors[3]->vx)/2.0;

		// y
		n = node;
		while (n->neighbors[1] == NULL) n = n->parent;
		node->cvy = (n->vy + n->neighbors[1]->vy)/2.0;
	} else {
		for (int k = 0; k < 4; k++) {
			computeNodalValues(node->children[k]);
		}
		node->cvx = node->children[0]->cvx;
		node->cvy = node->children[0]->cvy;
	}
}

void setNewAdvect(qNode* node) {
	if (node->leaf) {
		if (node->i == 0 || node->j == 0) {
			node->cvx = 0.0;
			node->cvy = 0.0;
			return;
		}
		int size = 1<<node->level;
		double x = (node->j + 0.5)/size;
		double y = (node->i + 0.5)/size;
		qNode* last= getSemiLagrangianLookback(oldRoot, &x, &y, 1, node->level);
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
	}
}


// TODO put on node class
void setNewFace(qNode* node) {
	int l = node->level;
	int i = node->i;
	int j = node->j;
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
		for (int k = 0; k < 4; k++) {
			setNewFace(node->children[k]);
		}
	}
}

void advectAndCopy() {

	copy(root, oldRoot);

	// assume nodal values have been computed on old tree
	// set nodal values on new tree to advected nodal values on old tree
	setNewAdvect(root);

	// set new face values from nodal values
	setNewFace(root);

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
		std::pair<double, double> vGrad = getLeaf(root, particles[i].x, particles[i].y, levels-1)->getVelocityAt(root, particles[i].x, particles[i].y);
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

void poissonReset(qNode* node) {
	node->p = 0;
	int size = 1<<node->level;
	//node->p = 1.0/size/size/12;
	if (node->leaf) {
		// set pressure to integral over domain
		double x = (node->j + 0.5)/size;
		double y = (node->i + 0.5)/size;
		if (poissonTestFunc == POISSONXY) {
			node->p = .2667;
			node->divV = 96 * (2*x-1) * (y-1)*(y-1) * y*y  +  32 * (x-1)*(x-1) * (2*x+1) * (1 - 6*y + 6*y*y);
		} else if (poissonTestFunc == POISSONCOS) {
			node->p = 1.0;
			double cx = cos(M_PI*2*x);
			double cy = cos(M_PI*2*y);
			node->divV = 4*M_PI*M_PI * (cx + cy - 2*cx*cy);
		} else if (poissonTestFunc == POISSONCOSKL) {
			node->p = 0.0;
			int k = 2;
			int l = 2;
			node->divV = -M_PI*M_PI*(k*k + l*l)*cos(M_PI*k*x)*cos(M_PI*l*y);
		} else {
			assert(false);
		}
	} else {
		node->divV = 0;
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
			int k = 2;
			int l = 2;
			correct = cos(M_PI * k * x) * cos(M_PI * l * y);
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

void runPoissonTest(bool print) {
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
	if (print) {
		endTime("poisson test");	
		printf("poisson test took %d vcycles\n", i);
		printf("average error: %f\n", avgError);
	}

	//printPressure();
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



void initSim() {
	printf("init sim\n");

	root = new qNode(NULL, 0, 0);
	oldRoot = new qNode(NULL, 0, 0);

	int startLevel = levels - 2;
	if (startState == POISSONTEST || startState == ERRORTEST || startState == PROJECTTEST) {
		if (adaptScheme == ADAPTNONE)
			startLevel = levels - 1;
		else
			startLevel = 0;
	} else if (startState == ADAPTTEST) {
		//assert(levels > 3);
		//startLevel = std::max(levels - 4, levels/2 + 1);
		//startLevel = levels - 3;
		startLevel = 1;
	}


	initRecursive(root, startLevel);

	//computeNodalValues(root);
	
	if (startState == POISSONTEST) {
		// adapt
		if (adaptScheme != ADAPTNONE) {
			//expandRadius(root, .3);
			//expandRadius(root, .2);
			runAdaptTest();
		}

		for (int i = 0; i < warmupRuns; i++) {
			runPoissonTest(false);
		}

		if (numToRun == 0) numToRun++;
		for (int i = 0; i < numToRun; i++) {
			runPoissonTest(true);
		}
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
			char* adapt = argv[++i];
			if (!strcmp("xy", func)) {
				poissonTestFunc = POISSONXY;
			} else if (!strcmp("cos", func)) {
				poissonTestFunc = POISSONCOS;
			} else if (!strcmp("coskl", func)) {
				poissonTestFunc = POISSONCOSKL;
			} else {
				return 1;
			}
			if (!strcmp("none", adapt)) {
				adaptTestFunc = ADAPTFUNCNONE;
			} else if (!strcmp("circle", adapt)) {
				adaptTestFunc = ADAPTCIRCLE;
			} else if (!strcmp("wave", adapt)) {
				adaptTestFunc = ADAPTWAVE;
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
		}
	}
	//levelToDisplay = levels/2;
    levelToDisplay = levels - 1;
	printf("headless: %d, levels: %d\n", headless, levels);

	// run tests
	//testNeighbors();
	//testMultilevelNeighbors();
	testIntersect();

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
