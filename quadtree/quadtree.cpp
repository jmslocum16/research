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
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

// profiling
int numLeaves, numNodes;
int* leafLevels;
long frameMS = 0;

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

double dt = .02;

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
class qNode {
	public:
		// tree stuffs
		int id;
		bool leaf;
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
	

		// finds the highest level used neighbor at the given multilevel
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
		std::pair<double, double> getPressureGradient(int targetLevel) {
			int size = 1<<level;
			double vals[4];	
			int sizes[4];
			for (int k = 0; k < 4; k++) {
				int ml = targetLevel;
				std::pair<qNode*, qNode*> neighbor = getNeighborInDir(k, &ml);
				
				if (neighbor.second == NULL) {
					// only 1
					vals[k] = neighbor.first->p;
				} else {
					vals[k] = (neighbor.first->p + neighbor.second->p)/2.0;
				}
				sizes[k] = 1<<ml;
			}
			return getGradient(vals, sizes, size);
		}
		
		// TODO
		/*std::pair<double, double> getLevelSetGradient(int targetLevel) {}*/
		
		std::pair<double, double> getVelocityGradient(int targetLevel) {
		    double vals[4];
			int sizes[4];
		
			for (int k = 0; k < 4; k++) {
				int ml = targetLevel;
				std::pair<qNode*, qNode*> neighbor = getNeighborInDir(k, &ml);
				if (neighbor.second == NULL) {
				vals[k] = (k < 2) ? neighbor.first->vy : neighbor.first->vx;
				} else {
					if (k < 2) {
						vals[k] = (neighbor.first->vy + neighbor.second->vy) / 2.0;
					} else {
						vals[k] = (neighbor.first->vx + neighbor.second->vx) / 2.0;
					}
				}	
				sizes[k] = 1<<ml;
			}
			
			return getGradient(vals, sizes, 1<<level);
		}
		
		std::pair<double, double> getVelocitySingleGradient(int targetLevel, bool x) {
		    double vals[4];	
			int sizes[4];
		
			for (int k = 0; k < 4; k++) {
				int ml = targetLevel;
				std::pair<qNode*, qNode*> neighbor = getNeighborInDir(k, &ml);
				if (neighbor.second == NULL) {
					vals[k] = x ? neighbor.first->vx : neighbor.first->vy;
				} else {
					if (x) {
						vals[k] = (neighbor.first->vx + neighbor.second->vx) / 2.0;
					} else {
						vals[k] = (neighbor.first->vy + neighbor.second->vy) / 2.0;
					}
				}	
				sizes[k] = 1<<ml;
			}
			
			return getGradient(vals, sizes, 1<<level);
		}

		std::pair<double, double> getVelocityProductGradient(int targetLevel) {
			 double vals[4];	
			int sizes[4];
		
			for (int k = 0; k < 4; k++) {
				int ml = targetLevel;
				std::pair<qNode*, qNode*> neighbor = getNeighborInDir(k, &ml);
				if (neighbor.second == NULL) {
					vals[k] = neighbor.first->vx * neighbor.first->vy;
				} else {
					double x = (neighbor.first->vx + neighbor.second->vx) / 2.0;
					double y = (neighbor.first->vy + neighbor.second->vy) / 2.0;
					vals[k] = x * y;
				}	
				sizes[k] = 1<<ml;
			}
			
			return getGradient(vals, sizes, 1<<level);
		}
		
		// curl(F) = d(Fy)/dx - d(Fx)/dy
		double getCurl() {
			std::pair<double, double> xgrad = getVelocitySingleGradient(levels - 1, true);
			std::pair<double, double> ygrad = getVelocitySingleGradient(levels - 1, false);
			return ygrad.first - xgrad.second;
		}
		
		void computeVelocityDivergence() {
			std::pair<double, double> grad = getVelocityGradient(levels - 1);
			divV = grad.first + grad.second;
			if (!leaf) {
				for (int k = 0; k < 4; k++) {
					children[k]->computeVelocityDivergence();
				}
			}
		}

		// adaptivity
		void expand() {
			assert(leaf);
    		printf("expanding cell %d, (%d, %d)\n", level, i, j);
    		std::pair<double, double> pressureGrad = getPressureGradient(levels - 1); // TODO should these be this level?.. maybe if it's a leaf do multilevel otherwise its level?
    		std::pair<double, double> velocityGrad = getVelocityGradient(levels - 1);
    		int size = 1<<level;
    		// std::pair<double, double> levelSetGrad;
    		for (int k = 0; k < 4; k++) {
        		int constX = k%2 == 0 ? -1 : 1;
        		int constY = k/2 == 0 ? -1 : 1;
				this->children[k] = new qNode(this, 2*i+(k/2), 2*j+(k%2));

        		children[k]->p = p + (constX * .25/size) * pressureGrad.first + (constY * .25/size) * pressureGrad.second;
				// TODO incorporate full velocity gradients?
        		children[k]->vx = vx + (constX * .25/size) * velocityGrad.first;
        		children[k]->vy = vy + (constY * .25/size) * velocityGrad.second;
        		// children[k]->phi = ...
    		}
			setChildrenNeighbors();
			leaf = false;
		}

		void contract() {
			assert(!leaf);
			for (int k = 0; k < 4; k++) {
				assert(children[k]->leaf);
			}
		    printf("contracting cell %d (%d, %d)\n", level, i, j);
			// average child values and set used/leaf
		    int size = 1<<level;
			p = 0.0;
			vx = 0.0;
			vy = 0.0;
			//phi = 0.0;
		    
			for (int k = 0; k < 4; k++) {
				p += children[k]->p;
				vx += children[k]->vx;
				vy += children[k]->vy;
				//phi += children[k]->phi;
				delete children[k];
			}
			
			p /= 4.0;
			vx /= 4.0;
			vy /= 4.0;
			//phi /= 4.0;
		    leaf = true;
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

qNode* getLeaf(qNode* node, double x, double y) {
	if (node->leaf) {
		return node;
	}
	int size = 1<<node->level;
	double midx = ((double)node->j)/size + 0.5/size;
	double midy = ((double)node->i)/size + 0.5/size;
	int newj = x < midx ? 0 : 1;
	int newi = y < midy ? 0 : 1;
	return getLeaf(node->children[2*newi + newj], x, y);
}

qNode* getSemiLagrangianLookback(qNode* r, double x, double y, int steps) {
	double newdt = dt / steps;
	qNode* cur = getLeaf(r, x, y);
	while (steps--) {
		x -= cur->vx * newdt;
		y -= cur->vy * newdt;
		if (x < 0 || x > 1)
		cur = getLeaf(r, x, y);
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
			if (node->p > 0) {
				redPercent = node->p/maxP;
			}
			if (node->p < 0) {
				bluePercent = node->p/-minP;
			}
			
			glColor3f(redPercent, 0, bluePercent);
		}

		glRectf(x, y, x + 2.0/size, y + 2.0/size);

		if (drawVelocity && node->level <= 5) { // if size is > 40 this is pretty useless
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
		if (drawCells && node->level <= 5) {
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

/* Handler for window-repaint event. Call back when the window first appears and
   whenever the window needs to be re-painted. */
void display() {
	glClear(GL_COLOR_BUFFER_BIT);   // Clear the color buffer with current clearing color

	printf("display\n");

	minP = 100000000.0;
	maxP = -100000000.0;
	maxMag = 0.0;

	computeMax(root, levelToDisplay);

	drawMultilevel(root, levelToDisplay);

	
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
	testRoot->expand();
	testRoot->children[1]->expand();
	testRoot->children[2] ->expand();
	testRoot->children[1]->children[0]->expand();
	testRoot->children[1]->children[2]->expand();

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
	testRoot->expand();

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			qNode* c = testRoot->children[i * 2 + j];
			c->expand();
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

void computeResidual(qNode* node) {
	int size = 1<<node->level;
	if (node->leaf) { // if leaf cell, compute residual
		// compute it: R = div(vx, vy) - 1/(ha)*sum of (s * grad) for each face
		double faceGradSum = 0.0;
		for (int k = 0; k < 4; k++) {
			int level = levels - 1; // residual is computed at largest multilevel only.
			std::pair<qNode*, qNode*> neighbor = node->getNeighborInDir(k, &level);
			double neighborp = (neighbor.second == NULL) ? neighbor.first->p : (neighbor.first->p + neighbor.second->p) / 2.0;
			int neighborsize = 1 << level;
			// integral around the edge of flux, or side length * face gradient
			//faceGradSum += 1 * (neighbor.p - this.p)/ (0.5/size + 0.5/neighborsize);
			faceGradSum += (neighborp - node->p) / (0.5/size + 0.5/neighborsize);
		}
		// h = length of cell = 1.0/size
		// a = "fluid volume fraction of the cell". Since no boundaries cutting through cell, always 1
		
		// double flux = 1/ha * faceGradSum = 1/(1/size * 1) * faceGradSum = size * faceGradSum;
		double flux = size * faceGradSum;
		double R = node->divV - flux;
		
		node->R = R;
		//printf("at [%d][%d], divV: %f, flux: %f, R: %f\n", i, j, grid[d][i*size+j].divV, flux, R);
		// if a*R > e, not done
		if (fabs(R) > eps) {
			doneVCycle = false;
			//printf("more work to do at [%d][%d], %f is bigger than epsilon of %f\n", i, j, grid[d][i*size+j].R, eps);
		} else {
			//printf("done with this cell already, %f is smaller than epsilon of %f\n", fabs(grid[d][i*size+j].R), eps);
		}
	} else {
        node->R = 0.0;
		for (int k = 0; k < 4; k++) {
			computeResidual(node->children[k]);
			node->R += node->children[k]->R;
		}
		node->R /= 4.0;
	}
}


// relax the node at the given multilevel
// puts the new value of dp in temp
bool relaxRecursive(qNode* node, int ml) {
	assert (node->level <= ml);
	int size = 1<<node->level;
	if (node->level == ml || node->leaf) {
        // dp = (h*a*divV - bsum)/asum
		double aSum = 0.0;
        double bSum = 0.0;
        double dp;
        double h;
		double oldDif = node->dp;
		for (int k = 0; k < 4; k++) {
			int level = ml;
			std::pair<qNode*, qNode*> neighbor = node->getNeighborInDir(k, &level);
			if (neighbor.second == NULL) {
				dp = neighbor.first->dp;
			} else {
				dp = (neighbor.first->dp + neighbor.second->dp)/2.0;
			}
            h = 0.5 / size + 0.5 /(1<<level);
            aSum -= 1/h;
            bSum += dp/h;
		}
		
        // dp = (h*R - bsum)/asum
        node->temp = (node->R/size - bSum) / aSum;
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

	printf("relaxing level %d\n", d);
	// get initial gress from previous level, if possible
	recursiveInject(root, d);

	bool done = false;
    int totalCycles = 0;
	while (r-- > 0 && !done) {
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

double curlThresh = .01;
double curlAdaptFunction(qNode* node) {
	double curl = node->getCurl();
	return fabs(curl) > curlThresh / (1<<node->level);
}

double pressureThresh = .01;
bool pGradAdaptFunction(qNode* node) {
    std::pair<double, double> pgrad = node->getPressureGradient(levels - 1);
    return fabs(pgrad.first + pgrad.second) > pressureThresh/(1<<node->level);
}

bool adaptFunction(qNode* node) {
    return pGradAdaptFunction(node);
	// return curlAdaptFunction(node);
}

void recursiveAdaptivity(qNode* node) {
    if (node->leaf) {
        // see if should node cell
        if (adaptFunction(node)) {
            node->expand();
        }
        return;
    }
    bool allChildrenLeaves = true;
    for (int k = 0; k < 4; k++) {
        if (!node->children[k]->leaf) {
            allChildrenLeaves = false;
            break;
        }
    }
    if (allChildrenLeaves) {
        if (!adaptFunction(node)) {
            // this wouldn't be expanded if it was a leaf, so it shouldn't have leaf children
            node->contract();
        }
    } else {
        for (int k = 0; k < 4; k++) {
            recursiveAdaptivity(node->children[k]);
        }
    }
}

void advectAndCopy(qNode* node, qNode* oldNode) {
	// need to create/delete new nodes if it adapted
	if (node->leaf && !oldNode->leaf) {
		// old node expanded, expand new node
		node->expand();
	} else if (!node->leaf && oldNode->leaf) {
		// old node contracted
		node->contract(); // TODO this doesn't handle multiple contractions per step, but rn only do single contract/expand per step
	}
	
	// copy values
	node->p = oldNode->p;
	node->vx = oldNode->vx;
	node->vy = oldNode->vy;
	//node->phi = oldNode->phi;


	// advect: U** = Un - An*dt
	// An = U*gradU, U*grad = du/dx + du/dy, so U* that is du*U/dx + du*U/dy
	// An = <dVx^2/dx + dVyVy/dy, dXyVy/dx + dVy^2/dy>
	int size = 1<<node->level;
	double x = (node->j + 0.5)/size;
	double y = (node->i + 0.5)/size;
	qNode* last = getSemiLagrangianLookback(oldNode, x, y, 1);
	printf("lookback from node %d (%d, %d) gives %d (%d, %d)", node->level, node->i, node->j, last->level, last->i, last->j);
	std::pair<double, double> vGrad = last->getVelocityGradient(levels - 1);
	std::pair<double, double> vProductGrad = last->getVelocityProductGradient(levels - 1);

	double advectX = 2 * node->vx * vGrad.first + vProductGrad.second;
	double advectY = vProductGrad.first + 2 * node->vy * vGrad.second;
	
	//printf("advecting grid[%d][%d][%d] by (%f, %f)\n", k, i, j, advectX, advectY);	
	//printf("pressure at [%d][%d][%d]: %f\n", k, i, j, grid[k][i*size+j].p);

	node->vx -= advectX * dt;
	node->vy -= advectY * dt;

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
	node->p += node->dp;
	if (!node->leaf) {
		for (int k = 0; k < 4; k++) {
			correctPressure(node->children[k]);
		}
	}
}

void project(qNode* node) {
	// correct velocity with updated pressure field to make non-divergent
	std::pair<double, double> grad = node->getPressureGradient(levels - 1);
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

void runStep() {
	printf("running step\n");
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


	// compute velocity divergence for advection calculation
	// TODO fix?
	//oldRoot->computeVelocityDivergence();
	
	// advect velocity, copy old stuff
	//printf("advecting velocity\n");
	advectAndCopy(root, oldRoot);

	// grid's vx and vy are now provisional velocity
	// compute velocity divergence of new grid to use in residual calcuation
	root->computeVelocityDivergence();

	// compute all residuals, so that the whole bottom multilevel is computed
	doneVCycle = true;
	computeResidual(root);

	int numCycles = 0;
	
	while (!doneVCycle) {
		numCycles++;
		//printf("start v cycle %d\n", numCycles);

		
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
		computeResidual(root);

		// TODO REMOVE
		// doneVCycle = true;
		//if (numCycles == 1) break;

		//printf("end v cycle\n");
	}
	
    // given new state, do adaptivity
    //recursiveAdaptivity(root);

	project(root);

	clampVelocity(root);

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
	printPressure();

	// increase frame number here instead of in display so that you get 1 of each frame, not dependent on glut's redrawing, like on alt-tabbing or anything
	frameNumber++;
}

void initNodeLeftUniform(qNode* node) {
	int size = 1<<node->level;
	node->p = 1.0 / size / size;
	node->vx = 0.0;
	node->vy = 0.0;
	node->phi = 0.0;
	if (node->i == size/2 || node->i == size/2-1) {
		node->vx = -1;
	}
	
}

void initNodeFunction(qNode* node) {
	initNodeLeftUniform(node);
	// TODO
}

// inits leaves according to test and averages non-leaves
void initRecursive(qNode* node, int d) {
	if (node->level == d) {
		initNodeFunction(node);
	} else {
		node->expand();

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

// TODO add different possible starting states
void initSim() {
	printf("init sim\n");
	root = new qNode(NULL, 0, 0);
	oldRoot = new qNode(NULL, 0, 0);

	initRecursive(root, levels - 2);

	clampVelocity(root);

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
	for (int i = 1; i < argc; i++) {
		char* arg = argv[i];
		if (!strcmp("--headless", arg)) {
			headless = true;
		} else if (!strcmp("-levels", arg)) {
			levels= atoi(argv[++i]);
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
        }
	}
	//levelToDisplay = levels/2;
    levelToDisplay = levels - 1;
	printf("headless: %d, levels: %d\n", headless, levels);

	// seed random
	srandom(time(NULL));

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
	windowid = glutCreateWindow("MultiGrid");  // Create window with the given title
	glutDisplayFunc(display);       // Register callback handler for window re-paint event
	glutKeyboardFunc(keyboard);
	initGL();                       // Our own OpenGL initialization
	glutMainLoop();                 // Enter the event-processing loop
	return 0;
}
