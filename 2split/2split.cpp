#include <GL/glut.h>  // GLUT, include glu.h and gl.h

#include <algorithm>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

double epsilon = .001;
int changed;
double maxP, minP;
int maxheight;

// profiling
int numLeaves, numNodes;
int* leafLevels;
long frameMS = 0;

int nodeId = 0;
class s2Node {
	public:
		int id;
		bool leaf;
		bool wasExpanded;
		//double P, Pold, Vx, Vxold, Vy, Vyold;
		double values[3];
		double oldValues[3];
		s2Node* parent;
		int split = -2; // -2: none, -1 all ways, 0 xsplit, 1 ysplit
		s2Node* children[4]; // 2D array in 1D array
		s2Node* splitChildren[2];
		s2Node* neighbors[2][2];


		int xlevel, ylevel;
		double x, y, x2, y2;
		bool hasFaces;
		double faces[2][2][3]; // x/y, direction, P, Vx, Vy

		s2Node(s2Node *p, double x1, double y1, double x2, double y2, int xlevel, int ylevel): parent(p), x(x1), y(y1), x2(x2), y2(y2), xlevel(xlevel), ylevel(ylevel) {
			id = nodeId++;
			leaf = true;
			wasExpanded = false;
			for (int i = 0; i < 4; i++) {
				children[i] = NULL;
				neighbors[i/2][i%2] = NULL;
			}
			splitChildren[0] = NULL;
			splitChildren[1] = NULL;
		}
		~s2Node() {
			// set your neighbors' neighbor pointer at you to be null
			for (int i = 0; i < 2; i++) {
				if (neighbors[0][i] != NULL) {
					neighbors[0][i]->neighbors[0][1-i] = NULL;
				}
				if (neighbors[1][i] != NULL) {
					neighbors[1][i]->neighbors[1][1-i] = NULL;
				}
			}
			if (parent != NULL) parent->leaf = true;

			// TODO don't need to delete parent/child/neighbor arrays right?...
		}

	
		void contractLeaves() {
			if (leaf) {
				return;
			}
			bool allLeaves = true;
			double max = -1000000000.0;
			double min = -max;
			double sums[2][2][3];
			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
					sums[i][j] = 0.0;
			if (split == -1) {
				for (int i = 0; i < 4; i++) {
					this->children[i]->contractLeaves();
					allLeaves &= children[i]->leaf;
					max = std::max(max, children[i]->values[0]);
					min = std::min(min, children[i]->values[0]);
					for (int k = 0; k < 3; k++) {
						sums[0][i/2][k] += children[i]->values[k];
						sums[1][i%2][k] += children[i]->values[k];
					}
				}
				for (int i = 0; i < 4; i++)
					for (int k = 0; k < 3; k++)
						sums[i/2][i%2][k] /= 2.0;
				
			}
			double mult = (x2-x) * (y2-y) / 4.0;
			if (allLeaves && ((max - min) * mult) < epsilon) {
				contract();
			} else {
				double xdif = abs(sums[0][0][0] - sums[0][1][0]);
				double ydif = abs(sums[1][0][0] - sums[1][1][0]);
				if (split == -1 && xdif > ydif) {
					contractToXSplit(sums[0][0], sums[0][1]);
				} else if (split == -1) {
					contractToYSplit(sums[1][0], sums[1][1]);
				}
				if (wasExpanded) {
					changed++;
				}
			}
			wasExpanded = false;
		}

		// should only be called 1 level up from leaves
		void contract() {
			if (!leaf) {
				for (int k = 0; k < 3; k++) this->values[k] = 0.0;
				for (int i = 0; i < 4; i++) {
					for (int k = 0; k < 3; k++) {
						this->values[k] += this->children[i]->values[k];
					}	
				}
				for (int k = 0; k < 3; k++) {
					this->values[k] /= 4;
				}
				
				// delete children
				for (int i = 0; i < 4; i++) {
					delete(children[i]);
				}
			}
			leaf = true;
			//basically delete children and average values
		}

		void contractToXSplit(double* x0vals, double* x1vals) {
			if (!leaf) {
				contract();
			}
			double xMid = (x+x2)/2;
			this->splitChildren[0] = new s2Node(this, x, y, xmid, y2, xlevel+1, ylevel);
			this->splitChildren[1] = new s2Node(this, xmid, y, x2, y2, xlevel+1, ylevel);
			for (int k = 0; k < 3; k++) {
				this->splitChildren[0]->values[k] = x0vals[k];
				this->splitChildren[1]->values[k] = x1vals[k];
			}
			setXSplitNeighbors();
		}

		void contractToYSplit(double* y0vals, double* y1vals) {
			if (!leaf) {
				contract();
			}
			double yMid = (y+y2)/2;
			this->splitChildren[0] = new s2Node(this, x, y, x2, ymid, xlevel, ylevel+1);
			this->splitChildren[1] = new s2Node(this, x, ymid, x2, y2, xlevel, ylevel+1);
			for (int k = 0; k < 3; k++) {
				this->splitChildren[0]->values[k] = y0vals[k];
				this->splitChildren[1]->values[k] = y1vals[k];
			}
			setYSplitNeighbors();
		}


		void setXSplitNeighbors() {
			// x neighbors
			splitChildren[0]->neighbors[0][0] = this->neighbors[0][0];
			splitChildren[1]->neighbors[0][1] = this->neighbors[0][1];
			splitChildren[0]->neighbors[0][1] = splitChildren[1];
			splitChildren[1]->neighbors[0][0] = splitChildren[0];

			// y neighbors
			s2Node* bottomParent = NULL;
			if (this->neighbors[1][0] != NULL)
				bottomParent = this->neighbors[1][0]->getSplitFace(1, 1);

			if (bottomParent != NULL) {
				s2Node* other = (bottomParent->split == -1) ? bottomParent->children[1] : bottomParent->splitChildren[0];
				splitChildren[0]->neighbors[1][0] = other;
				//if (other->neighbors[1][1] == NULL)
					other->neighbors[1][1] = splitChildren[0];
				
				
				other = (bottomParent->split == -1) ? bottomParent->children[3] : bottomParent->splitChildren[1];
				splitChildren[1]->heighbors[1][0] = other;
				//if (other->neighbors[1][1] == NULL)
					other->neighbors[1][1] = splitChildren[1];
			}
			
			s2Node* topParent = NULL;
			if (this->neighbors[1][1] != NULL)
				topParent = this->neighbors[1][1]->getSplitFace(1, 0);

			if (this->neighbors[1][1] != NULL) {
				s2Node* other = (topParent->split == -1) ? topParent->children[0] : topParent->splitChildren[0];
				splitChildren[0]->neighbors[1][1] = other;
				//if(other->neighbors[1][0] == NULL)
					other->neighbors[1][0] = splitChildren[0];

				other = (topParent->split == -1) ? topParent->children[2] : topParent->splitChildren[1];
				// if (other->neighbors[1][0] == NULL)
					other->neighbors[1][0] = splitChildren[1];
			}
		}

		void setYSplitNeighbors() {
			// y neighbors
			splitChildren[0]->neighbors[1][0] = this->neighbors[1][0];
			splitChildren[1]->neighbors[1][1] = this->neighbors[1][1];
			splitChildren[0]->neighbors[1][1] = splitChildren[1];
			splitChildren[1]->neighbors[1][0] = splitChildren[0];

			// y neighbors
			s2Node* bottomParent = NULL;
			if (this->neighbors[1][0] != NULL)
				bottomParent = this->neighbors[1][0]->getSplitFace(1, 1);

			if (bottomParent != NULL) {
				s2Node* other = (bottomParent->split == -1) ? bottomParent->children[1] : bottomParent->splitChildren[0];
				splitChildren[0]->neighbors[1][0] = other;
				//if (other->neighbors[1][1] == NULL)
					other->neighbors[1][1] = splitChildren[0];
				
				
				other = (bottomParent->split == -1) ? bottomParent->children[3] : bottomParent->splitChildren[1];
				splitChildren[1]->heighbors[1][0] = other;
				//if (other->neighbors[1][1] == NULL)
					other->neighbors[1][1] = splitChildren[1];
			}
			
			s2Node* topParent = NULL;
			if (this->neighbors[1][1] != NULL)
				topParent = this->neighbors[1][1]->getSplitFace(1, 0);

			if (this->neighbors[1][1] != NULL) {
				s2Node* other = (topParent->split == -1) ? topParent->children[0] : topParent->splitChildren[0];
				splitChildren[0]->neighbors[1][1] = other;
				//if(other->neighbors[1][0] == NULL)
					other->neighbors[1][0] = splitChildren[0];

				other = (topParent->split == -1) ? topParent->children[2] : topParent->splitChildren[1];
				// if (other->neighbors[1][0] == NULL)
					other->neighbors[1][0] = splitChildren[1];
			}
		}



		s2Node* getSplitFace(int f, int fi) {
			if (leaf) return NULL;
			if (split == -1 || split == f) {
				return this;
			} else {
				splitChildren[fi]->getSplitFace(f, fi);
			}
		}

		// functions for expanding

		void expandLeaves() {
			if (leaf) {
				if (this->level < maxheight) {
					expand(true);
				}
			} else {
				for (int i = 0; i < 4; i++) {
					this->children[i]->expandLeaves();
				}
			}
		}

		void expand(bool set) {
			if (!leaf) return;
			leaf = false;
			double sideLenX = (x2-x)/2;
			double sideLenY = (y2-y)/2;
			for (int i = 0; i < 4; i++) {
				double newx = this->x + (i/2) * sideLenX;
				double newy = this->y + (i%2) * sideLenY;
				children[i] = new s2Node(this, newx, newy, newx + sideLenX, newy + sideLenY);
				for (int k = 0; k < 3; k++) {
					children[i]->oldValues[k] = this->oldValues[k];
				}
				// TODO REMOVE
				double midx = newx + sideLenX/2;
				double midy = newy + sideLenY/2;
				double theta = atan(midy/midx);
				
				if (set) {
					if (sqrt(midx*midx+midy*midy) < 0.5 + /*.25*sin(theta * 8.0)*/0.0)  {
						children[i]->oldValues[0] = 10.0;
					} else {
						children[i]->oldValues[0] = 0.0;
					}
					children[i]->oldValues[1] = 0.0;
					children[i]->oldValues[2] = 0.0;
				}
			}
			setChildrenNeighbors();
			wasExpanded = true;
		}

		void setChildrenNeighbors() {
			// set both ways because the first time a node is expanded the other one may not have been
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					s2Node* c = children[i * 2 + j];
					c->neighbors[0][i] = NULL;
					if (this->neighbors[0][i] != NULL
							&& !this->neighbors[0][i]->leaf) {
						c->neighbors[0][i] = this->neighbors[0][i]->children[(1 - i) * 2 + j]; 
						this->neighbors[0][i]->children[(1 - i) * 2 + j]->neighbors[0][1-i] = c; 
					}
					c->neighbors[0][1-i] = this->children[(1-i) * 2 + j];
					this->children[(1-i) * 2 + j]->neighbors[0][i] = c;

					c->neighbors[1][j] = NULL;
					if (this->neighbors[1][j] != NULL && !this->neighbors[1][j]->leaf) {
						c->neighbors[1][j] = this->neighbors[1][j]->children[i * 2 + (1 - j)]; 
						this->neighbors[1][j]->children[i * 2 + (1 - j)]->neighbors[1][1-j] = c; 
					}
					c->neighbors[1][1-j] = this->children[i * 2 + (1 - j)];
					this->children[i * 2 + (1 - j)]->neighbors[1][j] = c;
				}
			}
		}

		void solveFaces() {
			if (leaf) {
				computeFaces();
			} else {
				for (int i = 0; i < 4; i++) {
					this->children[i]->solveFaces();
				}
			}
		}

		void computeFaces() {
			//printf("computing faces for %d\n", this->id);
			double xp = (x+x2)/2;
			double yp = (y+y2)/2;
			
			computeFace(0, 0, xp, x2-x);
			computeFace(0, 1, xp, x2-x);
			computeFace(1, 0, yp, y2-y);
			computeFace(1, 1, yp, y2-y);
		}
		
		// TODO change to use old vals
		void computeFace(int f, int fi, double p, double h) {
			
			if (this->neighbors[f][fi] != NULL && this->neighbors[f][fi]->leaf) {
				// simple interpolation
				//printf("simple interp\n");
				for (int k = 0; k < 3; k++) {
					this->faces[f][fi][k] = (this->neighbors[f][fi]->oldValues[k] - this->oldValues[k])/h;
				}
				return;
			}
			double nh; //height of neighbors
			double v[3]; // values
			if (neighbors[f][fi] != NULL) {
				neighbors[f][fi]->getInterpVals(this->level, f, 1-fi,  &nh, v);
				double dist = h/2 + nh/2;
				for (int k = 0; k < 3; k++) {
					this->faces[f][fi][k] = (v[k] - this->oldValues[k]) / dist; // using distance between pts
				}
			} else {
				s2Node* par = this->parent;
				while (par != NULL && par->neighbors[f][fi] == NULL) {
					par = par->parent;
				}
				if (par == NULL) {
					// no neighbor at that point, use edge value
					//printf("edge\n");
					for (int k = 0; k < 3; k++) {
						this->faces[f][fi][k] = 0.0;
					}
				} else {
					//printf("interp at higher node\n");
					// TODO make this not suck
					// use cell center
					double ph = (f==0) ? par->x2-par->x : par->y2 - par->y;
					double dist = h/2 + ph/2;
					for (int k = 0; k < 3; k++) {
						this->faces[f][fi][k] = (par->neighbors[f][fi]->oldValues[k] - this->oldValues[k]) / dist;
					}
				}
			}
		}

		void getInterpVals(double level, int f, int fi, double* h, double* nv) {
			if (leaf) {
				*h = (f == 0) ? x2-x : y2-y;
				for (int k = 0; k < 3; k++) {
					nv[k] = this->oldValues[k];
				}
			} else {
				s2Node* best = NULL;
				s2Node* a = NULL;
				s2Node* b = NULL;
				if (f == 0) {
					a = this->children[2*fi]->getCorner(fi, 1);
					b = this->children[2*fi + 1]->getCorner(fi, 0);
				} else {
					a = this->children[fi]->getCorner(1, fi);
					b = this->children[2+fi]->getCorner(0,fi);
				}
				if (a->level > b->level) {
					best = a;
				} else if (b->level > a->level) {
					best = b;
				} else {
					// same level
					*h = (f==0 ? a->x2 - a->x : a->y2 - a->y);
					for (int k = 0; k < 3; k++) {
						nv[k] = (a->oldValues[k] + b->oldValues[k])/2.0;
					}
					return;
				}
				*h = (f==0) ? best->x2 - best->x : best->y2 - best->y;
				for (int k = 0; k < 3; k++) {
					nv[k] = best->oldValues[k];
				}
			}
		}

		s2Node* getCorner(int i, int j) {
			if (leaf) return this;
			return this->children[2*i+j]->getCorner(i, j);
		}


		/// reel functions
		void advect(double dt) {
			if (leaf) {
				double dot = this->faces[0][1][1] - this->faces[0][0][1] + this->faces[1][1][2] - this->faces[1][0][2];
				double change = 1 - dt * dot;
				this->values[1] = this->oldValues[1] * change;
				this->values[2] = this->oldValues[2] * change;
			} else {
				for (int i = 0; i < 4; i++) children[i]->advect(dt);
			}
		}

		void diffuse(double dt) {
			if (leaf) {
				/*double div = this->faces[0][1][0] + this->faces[0][0][0] + this->faces[1][1][0] + this->faces[1][0][0];
				this->values[0] = this->oldValues[0] + div * dt;*/
				this->values[0] = this->oldValues[0];
				minP = std::min(minP, this->values[0]);
				maxP = std::max(maxP, this->values[0]);
			} else {
				for (int i = 0; i < 4; i++) children[i]->diffuse(dt);
			}
		}

		void project(double dt) {
			if (leaf) {
				double xgrad = this->faces[0][1][0]  - this->faces[0][0][0];
				double ygrad = this->faces[1][1][0] - this->faces[1][0][0];
				this->values[1] -= xgrad * dt;
				this->values[2] -= ygrad * dt;
			}
		}

		void reset() {
			if (!leaf) {
				for (int i = 0; i < 4; i++) this->children[i]->reset();
			}
			std::swap(oldValues, values);
			this->hasFaces = false;
			this->wasExpanded = false;
		}

		void resetExpanded() {
			if (!leaf) {
				for (int i = 0; i < 4; i++) children[i]->resetExpanded();
			}
			this->wasExpanded = false;
		}
		
		void drawCell() {
			if (leaf) {
				double p = this->values[0];
				double percent = (p - minP) / (maxP - minP);
				if (abs(minP - maxP) < .001) {
					percent = 0.5;
				}
				glColor3f(percent, 0, 1.0 - percent);
				glRectf(x, y, x2, y2);
			} else {
				for (int i = 0; i < 4; i++) {
					this->children[i]->drawCell();
				}
			}
		}
		void drawOutline() {
			if (!leaf) {
				for (int i = 0; i < 4; i++) {
					this->children[i]->drawOutline();
				}
			}
			glColor3f(.1, .1, .1);
			glBegin(GL_LINE_LOOP);
				glVertex2f(x, y);
				glVertex2f(x2, y);
				glVertex2f(x2, y2);
				glVertex2f(x, y2);
			glEnd();
		}

		// profiling
		void profileNode() {
			if (!leaf) {
				for (int i = 0; i < 4; i++) children[i]->profileNode();
			} else {
				numLeaves++;	
				leafLevels[this->level]++;
			}
			numNodes++;
		}
};


 
/* Initialize OpenGL Graphics */
void initGL() {
   // Set "clearing" or background color
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
}

bool headless; // run without ui for profiling

s2Node* root;

double visc, kS, aS;

int sourceX, sourceY;
double sourceVal;
 
/* Handler for window-repaint event. Call back when the window first appears and
   whenever the window needs to be re-painted. */
void display() {
	glClear(GL_COLOR_BUFFER_BIT);   // Clear the color buffer with current clearing color
	
	/*double sideLen = 2.0 / size;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double p = P[i * size + j];
			double percent = (p - minP) / (maxP - minP);
			if (abs(minP - maxP) < .001) {
				percent = 0.5;
			}
			glColor3f(percent, 0, 1.0 - percent);
			double x = -1.0  + j * sideLen;
			double y = 1.0 - i * sideLen; 
			glRectf(x, y, x + sideLen, y - sideLen);
		}
	}*/

	// draw tree
	root->drawCell();
	root->drawOutline();

	glFlush();  // Render now
}

/*double linInterp(double x, double y, double* arr) {
	int i = (int)y;
	int j = (int)x;
	double di = y - i;
	double dj = x - j;
	double result = 0;
	for (int ii = 0; ii <= 1; ii++) {
		for(int jj = 0; jj <= 1; jj++) {
			double disti = (ii == 0) ? 1-di : di;
			double distj = (jj == 0) ? 1-dj : dj;
			result += disti * distj * get(i + ii, j + jj, arr);
		}
	}
}*/

/*void traceParticle(double x, double y, double* vx, double* vy, double delta, double* resultX, double* resultY) {
	// TODO runge-kutta 2

	// linear euler
	*resultX = linInterp(x, y, vx) * delta;
	*resultY = linInterp(x, y, vy) * delta;
}*/

void velocityStep(double visc, double dt) {
	// add forces
	// TODO

	// advect
	/*transport(Vx, Vxold, Vxold, Vyold, dt);
	transport(Vy, Vyold, Vxold, Vyold, dt);

	// diffuse
	diffuse(Vx, Vxold, visc, dt);
	diffuse(Vy, Vyold, visc, dt);

	// project
	project(dt);*/
}

void scalarStep(double k, double a, double dt) {
	// add forces (source)
	/*Pold[sourceX * size + sourceY] += sourceVal * dt;	

	// transport by velocity
	transport(P, Pold, Vx, Vy, dt);

	// diffuse
	diffuse(P, Pold, k, dt);

	// dissipate
	dissipate(a, dt);*/
}

long getMS() {
	struct timeval time;
	gettimeofday(&time, NULL);
	return 1000L*time.tv_sec + time.tv_usec / 1000L;
}

int curstep = -1;

void initFrameProfile() {
	memset(leafLevels, 0, maxheight + 1);
	numLeaves = 0;
	numNodes = 0;
	frameMS = getMS();
}

void endFrameProfile() {
	frameMS = getMS() - frameMS;
	
	printf("Profiling frame--------------\n");

	printf("frame ms: %ld\n", frameMS);
	printf("total leaves: %d\n", numLeaves);
	printf("total non-leaves: %d\n", numNodes - numLeaves);
	printf("total nodes: %d\n", numNodes);
	
	printf("height statistics\n");
	printf("height totals:\n");

	int sum = 0;
	int num = 0;
	for (int i = 0; i <= maxheight; i++) {
		printf("%d: %d\n", i, leafLevels[i]);
		sum += i * leafLevels[i];
		num += leafLevels[i];
	}

	double avg = 1.0 * sum / num;
	printf("average height: %lf\n", avg);

	printf("Done profiling frame----------\n\n");
}



void step() {
	// do simulation step in grid
	/*long t = time(NULL);
	double dt = .05;
	std::swap(P, Pold);
	std::swap(Vx, Vxold);
	std::swap(Vy, Vyold);

	velocityStep(visc, dt);
	scalarStep(kS, aS, dt);
	
	t = time(NULL) - t;
	printf("# milliseconds for frame: %ld\n", t);*/

	initFrameProfile();
	
	double dt = .05;
	
	changed = 1;
	root->reset();
	int i  = 0;
	while (changed && i++ < 20) {
		maxP = -1000000000.0;
		minP = -maxP;
		changed = 0;
		root->expandLeaves();
		root->solveFaces();
		root->advect(dt);
		root->diffuse(dt);
		root->project(dt);
		root->contractLeaves();
	}

	root->profileNode();
	endFrameProfile();

	printf("converged in %d iterations\n", i);
	/*if (curstep == -1) {
		root->reset();
		printf("reset\n");
	} else if (curstep % 6 == 0) {
		maxP = -1000000000.0;
		minP = -maxP;
		changed = 0;
		root->resetExpanded();
		root->expandLeaves();
		printf("expand\n");
	} else if (curstep % 6 == 1) {
		root->solveFaces();
		printf("solvefaces\n");
	} else if (curstep % 6 == 2) {
		//root->advect(dt);
		printf("advect\n");
	} else if (curstep % 6 == 3) {
		root->diffuse(dt);
		printf("diffuse\n");
	} else if (curstep % 6 == 4) {
		//root->project(dt);
		printf("project\n");
	} else {
		root->contractLeaves();
		printf("contract: %d\n", changed);
	}*/

	curstep++;	
}

double frand(double l, double h) {
	double f = (double)rand() / RAND_MAX;
	return l + f * (h - l);
}

//tests
void testNeighbors() {
	root = new s2Node(NULL, -1, -1, 1, 1);
	root->expand(false);

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			s2Node* c = root->children[i * 2 + j];
			c->expand(false);
			if (c->neighbors[0][1- i] != root->children[(1-i) * 2 + j]) {
				printf("WRONG INNER X NEIGHBOR (%d, %d)\n", i, j);
			}
			if (c->neighbors[1][1-j] != root->children[i * 2 + (1-j)]) {
				printf("WRONG INNER Y NEIGHBOR (%d,%d)\n", i, j);
			}
		}
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			s2Node* c = root->children[i * 2 + j];
			for (int ii = 0; ii < 2; ii++) {
				for (int jj = 0; jj < 2; jj++) {
					s2Node* cc = c->children[ii* 2 + jj];
					if (cc->neighbors[0][1- ii] != c->children[(1-ii) * 2 + jj]) {
						printf("WRONG INNER Y NEIGHBOR of %d    - (%d, %d): (%d, %d)\n", cc->id, i, j,  ii, jj);
						printf("expected: %d, actual: %d\n", cc->neighbors[0][1-ii]->id, c->children[(1-ii) * 2 + jj]->id);
					}			
					if (cc->neighbors[1][1-jj] != c->children[ii * 2 + (1-jj)]) {
						printf("WRONG INNER X NEIGHBOR of %d    - (%d, %d): (%d, %d)\n", cc->id, i, j,  ii, jj);
						printf("expected: %d, actual: %d\n", cc->neighbors[1][1-jj]->id, c->children[ii * 2 + (1-jj)]->id);
					}

					if (i == ii) {
						if (cc->neighbors[0][i] != NULL) {
							printf("null x neighbor incorrect: (%d,%d)(%d,%d)\n", i, j, ii, jj);
						}
					} else {
						s2Node* shouldNeighbor = root->children[(1-i) * 2 + j]->children[(1-ii) * 2 + jj];
						if (cc->neighbors[0][ii] != shouldNeighbor) {
							printf("cross-node neighbor of %d    incorrect: (%d,%d)(%d,%d)\n", cc->id, i, j, ii, jj);
							printf("expected: %d, actual: %d\n", cc->neighbors[0][ii]->id, shouldNeighbor->id);
						}
					}
					if (j == jj) {
						if (cc->neighbors[1][j] != NULL) {
							printf("null y neighbor incorrect: (%d,%d)(%d,%d)", i, j, ii, jj);
						}
					} else {
						s2Node* shouldNeighbor = root->children[i * 2 + (1-j)]->children[ii * 2 + (1-jj)];
						if (cc->neighbors[1][jj] != shouldNeighbor) {
							printf("cross-node neighbor of %d    incorrect: (%d,%d)(%d,%d)\n", cc->id, i, j, ii, jj);
							printf("expected: %d, actual: %d\n", cc->neighbors[1][jj]->id, shouldNeighbor->id);
						}

					}
				}
			}
		}
	}
}

void testInterp() {
	root = new s2Node(NULL, -1, -1, 1, 1);
	root->expand(false);

	root->children[0]->oldValues[0] = 1.0;
	root->children[1]->oldValues[0] = 2.0;
	root->children[2]->oldValues[0] = 4.0;
	root->children[3]->oldValues[0] = 8.0;

	/*root->solveFaces();


	// top level
	if (root->children[0]->faces[0][0][0] != 0.0 || root->children[0]->faces[1][0][0] != 0.0) {
		printf("children[0] boundary wrong\n");
	}
	if (root->children[0]->faces[0][1][0] != 3.0) {
		printf("children[0] inner x face wrong: %lf\n", root->children[0]->faces[0][1][0]);
	}
	if (root->children[0]->faces[1][1][0] != 1.0) {
		printf("children[0] inner y face wrong: %lf\n", root->children[0]->faces[1][1][0]);
	}

	if (root->children[1]->faces[0][0][0] != 0.0 || root->children[1]->faces[1][1][0] != 0.0) {
		printf("children[1] boundary wrong\n");
	}
	if (root ->children[1]->faces[0][1][0] != 6.0) {
		printf("children[1] inner x face wrong: %lf\n", root->children[1]->faces[0][1][0]);
	}
	if (root ->children[1]->faces[1][0][0] != -1.0) {
		printf("children[1] inner y face wrong: %lf\n", root->children[1]->faces[1][0][0]);
	}

	if (root->children[2]->faces[0][1][0] != 0.0 || root->children[2]->faces[1][0][0] != 0.0) {
		printf("children[2] boundary wrong\n");
	
	if (root ->children[2]->faces[0][0][0] != -3.0) {
		printf("children[2] inner x face wrong: %lf\n", root->children[2]->faces[0][0][0]);
	}
	if (root ->children[2]->faces[1][1][0] != 4.0) {
		printf("children[2] inner y face wrong: %lf\n", root->children[2]->faces[1][1][0]);
	}

	if (root->children[3]->faces[0][1][0] != 0.0 || root->children[3]->faces[1][1][0] != 0.0) {
		printf("children[3] boundary wrong\n");
	}
	if (root ->children[3]->faces[0][0][0] != -6.0) {
		printf("children[3] inner x face wrong: %lf\n", root->children[3]->faces[0][0][0]);
	}
	if (root ->children[3]->faces[1][0][0] != -4.0) {
		printf("children[3] inner y face wrong: %lf\n", root->children[3]->faces[1][0][0]);
	}*/

	// next level
	root->children[0]->expand(false);
	root->children[0]->children[0]->oldValues[0] = 1.0;
	root->children[0]->children[1]->oldValues[0] = 2.0;
	root->children[0]->children[2]->oldValues[0] = 4.0;
	root->children[0]->children[3]->oldValues[0] = 8.0;

	root->solveFaces();
	
	if (root->children[0]->children[1]->faces[1][1][0] != 0.0) {
		printf("node 0-1's 1-1 face is incorrect: %lf\n", root->children[0]->children[1]->faces[1][1][0]);
	}
	if (root->children[0]->children[2]->faces[0][1][0] != 0.0) {
		printf("node 0-2's 0-1 face is incorrect: %lf\n", root->children[0]->children[2]->faces[0][1][0]);
	}

	if (root->children[0]->children[0]->faces[1][1][0] != 2.0) {
		printf("inner face of 2nd level wrong: %lf\n", root->children[0]->children[0]->faces[1][1][0]);
	}

	if (root->children[0]->children[3]->faces[0][1][0] != -4/.75) {
		printf("node 0-3's 0-1 face is incorrect: %lf\n", root->children[0]->children[3]->faces[0][1][0]);
	}
	if (root->children[0]->children[3]->faces[1][1][0] != -6/.75) {
		printf("node 0-3's 1-1 face is incorrect: %lf\n", root->children[0]->children[3]->faces[1][1][0]);
	}

	if (root->children[1]->faces[1][0][0] != 3/.75) {
		printf("node 1's face did down interp wrong: %lf\n", root->children[1]->faces[1][0][0]);
	}
	if (root->children[2]->faces[0][0][0] != 2/.75) {
		printf("node 2's face did down interp wrong: %lf\n", root->children[2]->faces[0][0][0]);
	}
}

void initSim() {
	/*P = new double[size * size];
	Pold = new double[size * size];
	Vx = new double[size * size];
	Vxold = new double[size * size];
	Vy = new double[size * size];
	Vyold = new double[size * size];*/

	// TESTS
	// testNeighbors();
	testInterp();
	
	printf("done tests\n");

	root = new s2Node(NULL, -1, -1, 1, 1);
	root->expand(false);
	/*root->children[0]->values[0] = 1.0;
	root->children[2]->values[0] = 4.0;
	root->children[3]->values[0] = 8.0;
	root->children[1]->expand();
	root->children[1]->children[0]->values[0] = 1.0;
	root->children[1]->children[1]->values[0] = 2.0;
	root->children[1]->children[2]->values[0] = 4.0;
	root->children[1]->children[3]->values[0] =	8.0;*/


	// TODO REMOVE
	//root->reset();
	
	minP = 1.0;
	maxP = 8.0;

	visc = 1.0;
	kS = 1.0;
	aS = 1.0; // TODO configurable

	//sourceX = size / 2;
	//sourceY = size / 2;

	sourceVal = 1.0;
}

void initProfile() {
	leafLevels = new int[maxheight + 1];
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
		case 32:     // space key
			step();
			if (!headless) {
				glutPostRedisplay();
			}
		break;
	}
}

 
/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {
	glutInit(&argc, argv);          // Initialize GLUT
	headless = false;
	maxheight = 5;
	for (int i = 1; i < argc; i++) {
		char* arg = argv[i];
		if (!strcmp("--headless", arg)) {
			headless = true;
		} else if (!strcmp("-maxheight", arg)) {
			maxheight = atoi(argv[++i]);
		}
	}
	printf("headless: %d, maxheight: %d\n", headless, maxheight);

	// seed random
	srandom(time(NULL));

	initSim();
	initProfile();

	glutInitWindowSize(800, 800);   // Set the window's initial width & height
	glutCreateWindow("Quadtree");  // Create window with the given title
	glutInitWindowPosition(50, 50); // Position the window's initial top-left corner
	glutDisplayFunc(display);       // Register callback handler for window re-paint event
	glutKeyboardFunc(keyboard);
	initGL();                       // Our own OpenGL initialization
	glutMainLoop();                 // Enter the event-processing loop
	return 0;
}
