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
enum StartState { LEFT, SINK, SRCSINK, POISSONTEST, ADAPTTEST , ERRORTEST };
StartState startState;
enum AdaptScheme { ADAPTNONE, CURL, PGRAD };
AdaptScheme adaptScheme = ADAPTNONE;
enum AdaptTestFunc { ADAPTFUNCNONE, ADAPTXY, ADAPTSIN, PARABOLA };
AdaptTestFunc adaptTestFunc = ADAPTFUNCNONE;
enum PoissonTestFunc { POISSONFUNCNONE, POISSONXY, POISSONCOS, POISSONCOSKL };
PoissonTestFunc poissonTestFunc = POISSONFUNCNONE;
enum ErrorTestFunc { ERRORFUNCNONE, ERRORXY, ERRORCOS, ERRORCOSKL };
ErrorTestFunc errorTestFunc = ERRORFUNCNONE;

double dt = .03;

// drawing stuff
double minP;
double maxP;
double maxMag;


// navigating 
double deltas[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
// for a given delta, the 2 corners to go to on finer levels when calculating face gradients/laplacian
double corner1coords[4][2] = {{0, 0}, {1, 0}, {0, 0}, {0, 1}};
double corner2coords[4][2] = {{0, 1}, {1, 1}, {1, 0}, {1, 1}};

// utility functions
std::pair<double, double> getGradient(double vals[], int sizes[], int size) {
	// grad = a - b / (dist) = (a - b) / ((1/sizea) / 2 + (1/sizeb) / 2 + 1/size) = (a - b) / (.5/sizea + .5/sizeb + size)
	return std::make_pair((vals[2] - vals[3])/(0.5/sizes[2] + 0.5/sizes[3] + 1.0/size), (vals[0]-vals[1])/(0.5/sizes[0] + 0.5/sizes[1] + 1.0/size));	
}

int nodeId = 0;
enum NodeValue { P, VX, VY, PHI, DP};
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
		double p, vx, vy, dp, R, phi, divV, temp; // pressure ,x velocity, y velocity, pressure correction, residual, level set value

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
			assert(false);
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

		double getValueInDir(int k, NodeValue v, int* ml) {
			std::pair<qNode*, qNode*> neighbor = getNeighborInDir(k, ml);
			double val;
			if (*ml >= level) {
				if (v == P) {
					if (neighbor.second == NULL) val = neighbor.first->p;
					else val = (neighbor.first->p + neighbor.second->p)/2.0;
				} else if (v == VX) {
					if (neighbor.second == NULL) val = neighbor.first->vx;
					else val = (neighbor.first->vx + neighbor.second->vx)/2.0;
				} else if (v == VY) {
					if (neighbor.second == NULL) val = neighbor.first->vy;
					else val = (neighbor.first->vy + neighbor.second->vy)/2.0;
				} else if (v == PHI) {
					// TODO
				} else if (v == DP) {
					if (neighbor.second == NULL) val = neighbor.first->dp;
					else val = (neighbor.first->dp + neighbor.second->dp)/2.0;
				}
			} else {
				assert(neighbor.second == NULL);
				int nsize = 1<<(neighbor.first->level);
				double targetX = (neighbor.first->j + 0.5)/nsize;
				double targetY = (neighbor.first->i + 0.5)/nsize;
				if (k < 2) { // in y dir
					targetX = (j + 0.5)/(1<<level);
				} else {
					targetY = (i + 0.5)/(1<<level);
				}
				val = neighbor.first->getValueAt(targetX, targetY, v);
			}
			return val;
		}
	

		// finds the  neighbor at the given multilevel
		// returns 2  cells if on the same (higher) level, otherwise returns one and null
		// d is level, (i, j), k is the direction (deltas index). Value initially in level is the target level, value returned in level is the level of the neighboring cell
		// guaranteed d <= target level, otherwise that wouldn't make any sense..
		std::pair<qNode*, qNode*> getNeighborInDir(int k, int* ml) {
			int size = 1<<level;
			assert (level <= *ml);
			int newi = i + deltas[k][0];
			int newj = j + deltas[k][1];
			if (newi < 0 || newi >= size || newj < 0 || newj >= size) {
				// not on grid, use boundary conditions
				*ml = level;
				return std::make_pair(this, static_cast<qNode*>(NULL));
			} else if (neighbors[k] == NULL) {
				// go up until find the cell with that neighbor
				qNode* node = parent;
				while (node->neighbors[k] == NULL) {
					node = node->parent;
				}
				*ml = node->level;
				return std::make_pair(node->neighbors[k], static_cast<qNode*>(NULL));
			} else if (level == *ml || neighbors[k]->leaf) {
				// simply just your neighbor
				*ml = level;
				return std::make_pair(neighbors[k], static_cast<qNode*>(NULL));
			} else {
				qNode* c1 = neighbors[k]->corner1(k);
				qNode* c2 = neighbors[k]->corner2(k);
				int d = level + 1;
				while (d < *ml && !c1->leaf && !c2->leaf) {
					c1 = c1->corner2(k);
					c2 = c2->corner1(k);
					d++;
				}
				if (*ml == d || (c1->leaf && c2->leaf)) {
					*ml = d;
					return std::make_pair(c1, c2);
				} else if (!c1->leaf) {
					// terminated because c2 not used anymore. Keep following c1 to the end then return it.
					while (d < *ml && !c1->leaf) {
						c1 = c1->corner2(k);
						d++;
					}
					*ml = d;
					return std::make_pair(c1, static_cast<qNode*>(NULL));
				} else if (!c2->leaf) {
					// terminated because c1 not used anymore. Keep following c2 to the end then return it.
					while (d < *ml && !c2->leaf) {
						c2 = c2->corner1(k);
						d++;
					}
					*ml = d;
					return std::make_pair(c2, static_cast<qNode*>(NULL));
				} else {
		            printf("THIS SHOULD NEVER HAPPEN\n");
				}
			}
		}

		qNode* corner1(int k) {
			if (leaf) return NULL;
			int i = corner1coords[k][0];
			int j = corner1coords[k][1];
			return this->children[2*i+j];
		}
		qNode* corner2(int k) {
			if (leaf) return NULL;
			int i = corner2coords[k][0];
			int j = corner2coords[k][1];
			return this->children[2*i+j];
		}

		// utility functions
		
		std::pair<double, double> getVelocityGradient() {
		    double vals[4];
			int sizes[4];
		
			for (int k = 0; k < 4; k++) {
				int ml = leaf ? levels - 1 : level;
				vals[k] = getValueInDir(k, (k < 2) ? VY : VX, &ml);
				sizes[k] = 1<<ml;
			}
			
			return getGradient(vals, sizes, 1<<level);
		}

		std::pair<double, double> getValueGradient(NodeValue v) {
			double vals[4];
			int sizes[4];

			for (int k = 0; k < 4; k++) {
				int ml = leaf? levels - 1: level;
				vals[k] = getValueInDir(k, v, &ml);
				sizes[k] = 1<<ml;
			}
			return getGradient(vals, sizes, 1<<level);
		}
		
		// curl(F) = d(Fy)/dx - d(Fx)/dy
		double getCurl() {
			if (adaptTestFunc == ADAPTFUNCNONE) {
				//std::pair<double, double> xgrad = getVelocitySingleGradient(true);
				//std::pair<double, double> ygrad = getVelocitySingleGradient(false);
				std::pair<double, double> xgrad = getValueGradient(VX);
				std::pair<double, double> ygrad = getValueGradient(VY);
				return ygrad.first - xgrad.second;
			}
			int size = 1<<level;
			double x = (j + 0.5)/size;
			double y = (i + 0.5)/size;
			if (adaptTestFunc == ADAPTXY) {
				return .5*x*y;
			} else if (adaptTestFunc == ADAPTSIN) {
				return .2*sin(M_PI*x)*sin(M_PI*y);
			} else if (adaptTestFunc == PARABOLA) {
				return (.5-x)*(.5-y);
			}
			assert(false);
		}
		
		void computeVelocityDivergence() {
			std::pair<double, double> grad = getVelocityGradient();
			divV = grad.first + grad.second;
			if (!leaf) {
				for (int k = 0; k < 4; k++) {
					children[k]->computeVelocityDivergence();
				}
			}
		}

		std::pair<double, double> getVelocityAt(double x, double y) {
			return std::make_pair(getValueAt(x, y, VX), getValueAt(x, y, VY));
		}
		
		double getValueAt(double x, double y, NodeValue v) {
			double val;
			if (v == P)
				val = p;
			else if (v == VX)
				val = vx;
			else if (v == VY)	
				val = vy;
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
    		//printf("expanding cell %d, (%d, %d)\n", level, i, j);
    		int size = 1<<level;
    		// std::pair<double, double> levelSetGrad;
    		for (int k = 0; k < 4; k++) {

				this->children[k] = new qNode(this, 2*i+(k/2), 2*j+(k%2));

				if (calculateNewVals) {
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
};


qNode* root;
qNode* oldRoot;


// debugging

// current node, target d/i/j
qNode* get(qNode* node, int d, int i, int j) {
	if (node->level == d) {
		assert(node->i == i);
		assert(node->j == j);
		return node;
	} else if (node->leaf) return NULL;
	int leveldif = d - node->level;
	int midsize = 1<<(leveldif-1);
	int midi = node->i*(1<<leveldif) + midsize;
	int midj = node->j*(1<<leveldif) + midsize;
	int newi = (i < midi) ? 0 : 1;
	int newj = (j < midj) ? 0 : 1;
	return get(node->children[newi*2+newj], d, i, j);
}

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
		std::pair<double, double> vel = cur->getVelocityAt(*x, *y);
		*x -= vel.first * newdt;
		*y -= vel.second * newdt;
		*x = fmin(1.0, fmax(0, *x));
		*y = fmin(1.0, fmax(0, *y));
		cur = getLeaf(r, *x, *y, ml);
	}
	return cur;
}

// not fast but who cares
void printPressure() {
	for (int d = 0; d < levels; d++) {
		int size = 1<<d;
		printf("level %d\n", d);
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				qNode* n = get(root, d, i, j);
				if (n != NULL) {
					printf(" %.4f", n->p);
				} else {
					printf("      ");
				}
			}
			printf("\n");
		}
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
					// max size is 1.0/side length, scaled by the max magnitude
					glVertex2f(x + 1.0/size, y + 1.0/size);
					double scale = maxMag * size;
					glVertex2f(x + 1.0/size + vx / scale, y + 1.0/size + vy /scale);
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

// test multilevel neighbor functions
// construct following quadtree with 4 levels:
/* - - - - - - - -
 *|       | | |   |
 *         - -     
 *|       | | |   |
 *         - - - -
 *|       | | |   |
 *         - -
 *|       | | |   |
 * - - - - - - - - 
 *|   |   |       |
 *
 *|   |   |       |
 * - - - -
 *|   |   |       |
 *
 *|   |   |       |
 * - - - - - - - -
 */
void testMultilevelNeighbors() {
	// construct tree
	qNode* testRoot = new qNode(NULL, 0, 0);
	testRoot->expand(false);
	testRoot->children[1]->expand(false);
	testRoot->children[2] ->expand(false);
	testRoot->children[1]->children[0]->expand(false);
	testRoot->children[1]->children[2]->expand(false);

	// actual tests

	// deltas: {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
	// level 1
	int d = 1;
	std::pair<qNode*, qNode*> n;
	// cell 0's neighbor to the right
	int level = d;
	n = testRoot->children[0]->getNeighborInDir(2, &level);
	assert (assertNeighbors(testRoot->children[1], NULL, n.first, n.second));
	assert (level == 1);
	// cell 1's neighbor to the left
	level = d;
	n = testRoot->children[1]->getNeighborInDir(3, &level);
	assert (assertNeighbors(testRoot->children[0],  NULL, n.first, n.second));
	assert (level == 1);

	// level 2
	d = 2;

	level = d;
	n = testRoot->children[0]->getNeighborInDir(2, &level);
	assert (assertNeighbors(testRoot->children[1]->children[0], testRoot->children[1]->children[2], n.first, n.second));
	assert (level == 2);

	level = d;
	n = testRoot->children[1]->children[2]->getNeighborInDir(3, &level);
	assert (assertNeighbors(testRoot->children[0], NULL, n.first, n.second));
	assert (level == 1);

	// level 3
	d = 3;

	level = d;
	n = testRoot->children[0]->getNeighborInDir(2, &level);
	assert (assertNeighbors(testRoot->children[1]->children[0]->children[2], testRoot->children[1]->children[2]->children[0], n.first, n.second));
	assert (level == 3);

	level = d;
	n = testRoot->children[0]->getNeighborInDir(0, &level);
	assert (assertNeighbors(testRoot->children[2]->children[0], testRoot->children[2]->children[1], n.first, n.second));
	assert (level == 2);

	level = d;
	n = testRoot->children[3]->getNeighborInDir(1, &level);
	assert (assertNeighbors(testRoot->children[1]->children[2]->children[3], NULL, n.first, n.second));
	assert (level == 3);

	// clean up
	printf("done testing multilevel neighbors\n");
	delete testRoot;
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

// poisson solver functions

// returns max(abs(R))
double computeResidual(qNode* node) {
	int size = 1<<node->level;
	if (node->leaf) { // if leaf cell, compute residual
		// compute it: R = div(vx, vy) - 1/(ha)*sum of (s * grad) for each face
		double faceGradSum = 0.0;
		for (int k = 0; k < 4; k++) {
			int level = levels - 1; // residual is computed at largest multilevel only.
			double neighborp = node->getValueInDir(k, P, &level);
			int neighborsize = 1 << level;
			// integral around the edge of flux, or side length * face gradient
			faceGradSum += (neighborp - node->p) / (0.5/size + 0.5/neighborsize);
		}
		
		// double flux = 1/ha * faceGradSum = 1/(1/size * 1) * faceGradSum = size * faceGradSum;
		double flux = size * faceGradSum;
		double R = node->divV - flux;
		
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
		double aSum2 = 0.0;
		double bSum2 = 0.0;
		for (int k = 0; k < 4; k++) {
			int level = ml;
			dp = node->getValueInDir(k, DP, &level);
            h = 0.5 / size + 0.5 /(1<<level);
            aSum2 -= 1/h;
            bSum2 += dp/h;
		}

		for (int k = 0; k < 4; k++) {
			int oppositeK = 1 - (k%2);
			qNode* n = node;
			while (n != NULL && n->neighbors[k] == NULL) n = n->parent;
			if (n != NULL) {
				n->neighbors[k]->addFaceToLaplacian(node, ml, oppositeK, &aSum, &bSum, DP);
			} else {
				aSum -= 1;
				bSum += node->dp;
			}
		}
		
		// node->temp = (node->R/size - bSum) / aSum;
        // A*R = bSum - aSum*dp, dp = (bSum - A*R)/asum, A = h*h = 1/size^2
		node->temp = -(bSum - node->R/size/size)/aSum;

		//printf("old dp: bSum: %f, aSum: %f, R: %f\n", bSum2, aSum2, node->R);
		//printf("new dp: bSum: %f, aSum: %f, R: %f\n", bSum, aSum, node->R);
		double oldDp = (node->R/size - bSum2)/aSum2;
		double newDp = node->temp;

		/*if (oldDp != newDp) {
			printf("[%d](%d,%d) old dp: %f, new dp: %f\n", node->level, node->i, node->j, (node->R/size - bSum2)/aSum2, node->temp);
		}*/
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
	if (adaptTestFunc == ADAPTFUNCNONE)
		return fabs(curl) > curlThresh / (1<<node->level);
	return fabs(curl) > curlThresh;
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

void advectAndCopy(qNode* node, qNode* oldNode) {
	// need to create/delete new nodes if it adapted
	if (node->leaf && !oldNode->leaf) {
		// old node expanded, expand new node
		node->expand(false);
	} else if (!node->leaf && oldNode->leaf) {
		// old node contracted
		node->contract(); // TODO this doesn't handle multiple contractions per step, but rn only do single contract/expand per step
	}
	
	// copy values
	node->p = oldNode->p;
	node->vx = oldNode->vx;
	node->vy = oldNode->vy;
	//node->phi = oldNode->phi;


	// semi lagrangian advection - copy value at t - dt
	int size = 1<<node->level;
	double x = (node->j + 0.5)/size;
	double y = (node->i + 0.5)/size;
	qNode* last= getSemiLagrangianLookback(oldRoot, &x, &y, 1, node->level);
	//printf("lookback from node %d (%d, %d) gives %d (%d, %d)", node->level, node->i, node->j, last->level, last->i, last->j);
	std::pair<double, double> velInterp = last->getVelocityAt(x, y);

	node->vx = velInterp.first;
	node->vy = velInterp.second;

	// add gravity TODO(make better)
	//grid[index].vy -= .98 * dt;
	// density of water is 1000kg/m^3, so 100kg/m^2?.. density = M/V,so mass = V * density = density/(size+2)
	//grid[index].vy -= dt * 9.8*100.0/size/size;

	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			advectAndCopy(node->children[k], oldNode->children[k]);
		}
	}
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
	// std::pair<double, double> grad = node->getPressureGradient();
	std::pair<double, double> grad = node->getValueGradient(P);
	//node->vx -= grad.first * dt;
	//node->vy -= grad.second * dt;
	node->vx -= grad.first;
	node->vy -= grad.second;

	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			project(node->children[k]);
		}
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
	advectAndCopy(root, oldRoot);
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

	int numCycles = 0;
	
	while (!doneVCycle) {
		numCycles++;
		//printf("start v cycle %d\n", numCycles);

		double residual = runVCycle();	
		printf("residual after %d cycles: %f\n", numCycles, residual);
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
	project(root);

	clampVelocity(root);
	totalTime += endTime("projecting and clamping velocity");

	printf("total run step time: %.3f\n", totalTime);
	
	resetProfiling();
	root->profile();
	printProfiling();

	if (particleAlgorithm != PARTICLENONE) {
		advectParticles();
	}
	//printf("advecting phi\n");
	// dphi/dt + v dot grad(phi) = 0
	// dphi/dt = -v dot grad(phi)
	/*size = 1 << (levels - 1);
	double dphi[size*size];
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			std::pair<double, double> grad = getLevelSetGradient(grid, levels - 1, i, j);
			//printf("level set gradient at [%d][%d] is (%.2f, %.2f), v is (%.2f, %.2f), dphi is %.2f\n", i, j, grad.first, grad.second, grid[i*size+j].vx, grid[i*size+j].vy, dphi[i*size+j]);
			//dphi[i*size+j] = -dt * (grad.first * grid[i*size+j].vx + grad.second * grid[i*size+j].vy);
			dphi[i*size+j] = dt * (-grad.first * grid[levels - 1][i*size+j].vx - grad.second * grid[levels-1][i*size+j].vy);
		}
	}*/
	/*for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			grid[i * size + j].phi += dphi[i*size+j];
			printf(" %.2f", grid[i*size+j].phi);
		}
		printf("\n");
	}*/
	
	
	/*printf("done step, took %d iterations\n", numCycles);
	for (int d = 0; d < levels; d++) {
		printf("level %d\n", d);
		int size = 1<<d;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (grid[d][i*size+j].used) {
					printf(" %.3f", grid[d][i*size+j].p);
				} else {
					printf("      ");
				}
			}
			printf("\n");
		}
	}*/
	//printPressure();

	// increase frame number here instead of in display so that you get 1 of each frame, not dependent on glut's redrawing, like on alt-tabbing or anything
	frameNumber++;
}

void initNodeLeftUniform(qNode* node) {
	int size = 1<<node->level;
	node->p = 1.0/size/size;
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
	
	//node->p = (.5-x)*(.5-x) + (.5-y)*(.5-y); // use pressure gradient split to make multilevel
	
}

void initNodeFunction(qNode* node) {
	if (startState == LEFT || startState == ADAPTTEST || startState == ERRORTEST)
		initNodeLeftUniform(node);
	else if (startState == SINK)
		initNodeSink(node);
	else if (startState == SRCSINK)
		initNodeSrcSink(node);
	else if (startState == POISSONTEST)
		initPoissonTest(node);
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
	//node->p = 0;
	int size = 1<<node->level;
	node->p = 1.0/size/size/12;
	if (node->leaf) {
		double x = (node->j + 0.5)/size;
		double y = (node->i + 0.5)/size;
		if (poissonTestFunc == POISSONXY) {
			node->divV = 96 * (2*x-1) * (y-1)*(y-1) * y*y  +  32 * (x-1)*(x-1) * (2*x+1) * (1 - 6*y + 6*y*y);
		} else if (poissonTestFunc == POISSONCOS) {
			double cx = cos(M_PI*2*x);
			double cy = cos(M_PI*2*y);
			node->divV = 4*M_PI*M_PI * (cx + cy - 2*cx*cy);
		} else if (poissonTestFunc == POISSONCOSKL) {
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

void computePoissonError(qNode* node, double K, double* total) {
	if (node->leaf) {
		int size = 1<<node->level;
		double x = (node->j + 0.5)/size;
		double y = (node->i + 0.5)/size;
		double correct = K;
		if (poissonTestFunc == POISSONXY) {
			double newy = 2*y+1;
			correct += (2*x*x*x - 3*x*x + 1)  *  (newy*newy*newy*newy - 2*newy*newy + 1);
		} else if (poissonTestFunc == POISSONCOS) {
			correct += (1-cos(M_PI * 2 * x)) * (1-cos(M_PI * 2 * y));
		} else if (poissonTestFunc == POISSONCOSKL) {
			int k = 2;
			int l = 2;
			correct += cos(M_PI * k * x) * cos(M_PI * l * y);
		} else {
			assert(false);
		}
		*total += fabs(node->p - correct)/size/size;
	} else {
		for (int k = 0; k < 4; k++) {
			computePoissonError(node->children[k], K, total);
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

void runPoissonTest() {
	// adapt
	if (adaptScheme != ADAPTNONE) {
		expandRadius(root, .3);
		expandRadius(root, .2);
	}

	// set all pressure to 0 again
	assert(poissonTestFunc != POISSONFUNCNONE);
	poissonReset(root);
	
	double initialResidual = computeResidual(root);
	printf("initial residual: %f\n", initialResidual);

	// try manually fixing residual
	double avgR = 0.0;
	poissonAverageR(root, &avgR);
	printf("new average residual: %f\n", avgR);
	poissonCorrectR(root, avgR);
	
	// run V cycles, compute stuff
	int i = 0;

	double newR, oldR;
	newR = initialResidual;

	startTime();
	while (!doneVCycle) {
		i++;
		oldR = newR;
		newR = runVCycle();
		printf("residual after %d vcycles: %f\n", i, newR);
		double total = 0.0;
		double K = 0.0;
		poissonAverage(root, &K);
		double avgError = 0.0;
		computePoissonError(root, K, &avgError);
		printf("average error after %d vcycles: %f\n", i, avgError);

		// try manually fixing residual
		avgR = 0.0;
		poissonAverageR(root, &avgR);
		printf("new average residual: %f\n", avgR);
		poissonCorrectR(root, avgR);
		avgR = 0.0;
		poissonAverageR(root, &avgR);
		printf("average R after correction: %f\n", avgR);
		//doneVCycle |= fabs(newR-oldR) < eps/100;
		doneVCycle = getMaxR(root) < eps;
	}
	endTime("poisson test");
	
	printf("poisson test took %d vcycles\n", i);

	printPressure();
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
		int k = 3;
		int l = 3;
		if (k < 2) // y partial
			return -M_PI*l*cos(M_PI*k*x)*sin(M_PI*l*y);
		else //x partial
			return -M_PI*k*sin(M_PI*k*x)*cos(M_PI*l*y);
	}
	assert(false);
}

double calculateError(qNode* node, double* avgError) {
	if (node->leaf) {
		int size = 1<<node->level;
		for (int k = 0; k < 4; k++) {
			// calculate gradient in direction k
			int ml = levels - 1;
			double otherp = node->getValueInDir(k, P, &ml);
			double h = 0.5/size + 0.5/(1<<ml);
			double calc = (otherp - node->p) / h;
			
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
				printf(" at the same level = %f\n", error);
			} else if (ml < node->level) {
				printf(" up %d levels = %f\n", node->level - ml, error);
			} else {
				printf(" down %d levels = %f\n", ml - node->level, error);
			}
			*avgError += error/size/size;
			return error;
		}
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

void setAdaptTestValues(qNode* node) {
	// sets velocity for rendering n shit
	int size = 1<<node->level;
	double x = (node->j + 0.5)/size;
	double y = (node->i + 0.5)/size;
	if (adaptTestFunc == ADAPTXY) {
		node->vx = -.125*x*y*y;
		node->vy = .125*x*x*y;
	} else if (adaptTestFunc == ADAPTSIN) {
		node->vx = -.8/M_PI*cos(M_PI*x)*sin(M_PI*y);
		node->vy = -1/M_PI*sin(M_PI*x)*cos(M_PI*y);
	} else if (adaptTestFunc == PARABOLA) {
		node->vx = .5*x*y + .25*y*y;
		node->vy = .25*y + .5*x*x*y;
	}
	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			setAdaptTestValues(node->children[k]);
		}
	}
}

void runAdaptTest() {	
	setAdaptTestValues(root);
	assert(adaptScheme == CURL);
	std::swap(root, oldRoot);
	advectAndCopy(root, oldRoot);
	bool done = false;
	int numCycles = -1;
	while (!done) {
		std::swap(root, oldRoot);
		done = recursiveAdaptAndCopy(root, oldRoot);
		setAdaptTestValues(root);
		numCycles++;
	}
	printf("number of adapts done: %d\n", numCycles);

	resetProfiling();
	root->profile();
	printProfiling();	

}

// TODO add different possible starting states
void initSim() {
	printf("init sim\n");

	root = new qNode(NULL, 0, 0);
	oldRoot = new qNode(NULL, 0, 0);

	int startLevel = levels - 2;
	if (startState == POISSONTEST || startState == ERRORTEST) {
		if (adaptScheme == ADAPTNONE)
			startLevel = levels - 1;
		else
			startLevel = levels - 3;
	} else if (startState == ADAPTTEST) {
		assert(levels > 3);
		startLevel = std::max(levels - 4, levels/2 + 1);
	}

	initRecursive(root, startLevel);
	if (startState == POISSONTEST) {
		runPoissonTest();
		return;
	} else if (startState == ADAPTTEST) {
		runAdaptTest();
		return;
	} else if (startState == ERRORTEST) {
		runErrorTest();
		return;
	}

	clampVelocity(root);

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
			if (!strcmp("xy", func)) {
				adaptTestFunc = ADAPTXY;
			} else if (!strcmp("sin", func)) {
				adaptTestFunc = ADAPTSIN;
			} else if (!strcmp("parabola", func)) {
				adaptTestFunc = PARABOLA;
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
			
		}
	}
	//levelToDisplay = levels/2;
    levelToDisplay = levels - 1;
	printf("headless: %d, levels: %d\n", headless, levels);

	// run tests
	//testNeighbors();
	//testMultilevelNeighbors();

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
