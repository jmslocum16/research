#include <GL/glut.h>  // GLUT, include glu.h and gl.h

#include <algorithm>
#include <ctime>
#include <stdio.h>
#include <string.h>

int nodeId = 0;
class qNode {
	public:
		int id;
		bool leaf;
		double P, Pold, Vx, Vxold, Vy, Vyold;
		qNode* parent;
		qNode* children[4]; // 2D array in 1D array
		//qNode* xNeighbors[2]; // 2 nodes
		//qNode* yNeighbors[2]; // 2 nodes
		qNode* neighbors[2][2]; //TODO  rewrite neighbors
		bool hasFaces;
		double Faces[2][2][3]; // x/y, direction, P, Vx, Vy

		qNode(qNode *p): parent(p) {
			id = nodeId++;
			leaf = true;
			for (int i = 0; i < 4; i++) {
				children[i] = NULL;
			}
			xNeighbors[0] = xNeighbors[1] = NULL;
			yNeighbors[0] = yNeighbors[1] = NULL;
		}
		~qNode() {
			// set your neighbors' neighbor pointer at you to be null
			for (int i = 0; i < 2; i++) {
				if (xNeighbors[i] != NULL) {
					xNeighbors[i]->xNeighbors[1-i] = NULL;
				}
				if (yNeighbors[i] != NULL) {
					yNeighbors[i]->yNeighbors[1-i] = NULL;
				}
			}
			if (parent != NULL) parent->leaf = true;

			// TODO don't need to delete parent/child/neighbor arrays right?...
		}

		double getAvgPressure() {
			if (leaf) {
				return P;
			} else {
				return (children[0]->getAvgPressure() + children[1]->getAvgPressure() + children[2]->getAvgPressure() + children[3]->getAvgPressure())/4.0;
			}
		}
		double getAvgVx() {
			if (leaf) {
				return Vx;
			} else {
				return (children[0]->getAvgVx() + children[1]->getAvgVx() + children[2]->getAvgVx() + children[3]->getAvgVx())/4.0;
			}
		}

		double getAvgVy() {
				if (leaf) {
				return Vy;
			} else {
				return (children[0]->getAvgVy() + children[1]->getAvgVy() + children[2]->getAvgVy() + children[3]->getAvgVy())/4.0;
			}

		}

		// should only be called 1 level up from leaves
		void contract() {
			if (!leaf) {
				//average values
				P = getAvgPressure();
				Vx = getAvgVx();
				Vy = getAvgVy();
				
				// delete children
				for (int i = 0; i < 4; i++) {
					delete(children[i]);
				}
			}
			leaf = true;
			//basically delete children and average values
		}

		// functions for expanding
		void expand() {
			if (!leaf) return;
			leaf = false;
			for (int i = 0; i < 4; i++) {
				children[i] = new qNode(this);
			}
			setChildrenNeighbors();
		}
		void setChildrenNeighbors() {
			// set both ways because the first time a node is expanded the other one may not have been
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					qNode* c = children[i * 2 + j];
					c->xNeighbors[i] = NULL;
					if (this->xNeighbors[i] != NULL
							&& !this->xNeighbors[i]->leaf) {
						c->xNeighbors[i] = this->xNeighbors[i]->children[(1 - i) * 2 + j]; 
						this->xNeighbors[i]->children[(1 - i) * 2 + j]->xNeighbors[1-i] = c; 
					}
					c->xNeighbors[1-i] = this->children[(1-i) * 2 + j];
					this->children[(1-i) * 2 + j]->xNeighbors[i] = c;

					c->yNeighbors[j] = NULL;
					if (this->yNeighbors[j] != NULL && !this->yNeighbors[j]->leaf) {
						c->yNeighbors[j] = this->yNeighbors[j]->children[i * 2 + (1 - j)]; 
						this->yNeighbors[j]->children[i * 2 + (1 - j)]->yNeighbors[1-j] = c; 
					}
					c->yNeighbors[1-j] = this->children[i * 2 + (1 - j)];
					this->children[i * 2 + (1 - j)]->yNeighbors[j] = c;
				}
			}
		}

		//
		void clearFaces() {
			if (!leaf) {
				for (it i = 0; i < 4; i++) {
					children[i].clearFaces();
				}
				this.hasFaces = false;
			}
		}
};


 
/* Initialize OpenGL Graphics */
void initGL() {
   // Set "clearing" or background color
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
}

bool headless; // run without ui for profiling
int size; // size of initial simulation grid

double maxP, minP;
qNode* root;

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

	glFlush();  // Render now
}

//returns 0 if out of bounds, otherwise returns value at arr[i][j]
double get(int i, int j, double* arr) {
	if (i < 0 || i >= size || j < 0 || j >= size) return 0;
	return arr[i * size + j];
}


double linInterp(double x, double y, double* arr) {
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
}

void traceParticle(double x, double y, double* vx, double* vy, double delta, double* resultX, double* resultY) {
	// TODO runge-kutta 2

	// linear euler
	*resultX = linInterp(x, y, vx) * delta;
	*resultY = linInterp(x, y, vy) * delta;
}

void linearSolver() {
	// TODO ????? Use actual library?
}

void transport(double* newArr, double* oldArr, double* vx, double*vy, double dt) {
	double x, y, oldX, oldY;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)  {
			x = j + 0.5;	
			y = i + 0.5;
			traceParticle(x, y, vx, vy, -dt, &oldX, &oldY);
			newArr[i * size + j] = linInterp(oldX, oldY, oldArr);
		}
	}	
}

void diffuse(double * newArr, double * oldArr, double constant, double dt) {
	// TODO linear solver
}

void project(double dt) {
	// Vx, Vxold, Vy, Vyold

	// TODO linear solver
}

void dissipate(double constant, double dt) {
	// P, Pold
	double denom = 1 + constant * dt;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			// TODO wat
			// P[i * size + j] += Pold[i * size + j] / denom;
			// P[i * size + j] * = 1 / denom;
		}
	}
}

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
}

double frand(double l, double h) {
	double f = (double)rand() / RAND_MAX;
	return l + f * (h - l);
}


void initSim() {
	/*P = new double[size * size];
	Pold = new double[size * size];
	Vx = new double[size * size];
	Vxold = new double[size * size];
	Vy = new double[size * size];
	Vyold = new double[size * size];*/
	root = new qNode(NULL);
	root->expand();
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			qNode* c = root->children[i * 2 + j];
			c->expand();
			if (c->xNeighbors[1- i] != root->children[(1-i) * 2 + j]) {
				printf("WRONG INNER X NEIGHBOR (%d, %d)\n", i, j);
			}
			if (c->yNeighbors[1-j] != root->children[i * 2 + (1-j)]) {
				printf("WRONG INNER Y NEIGHBOR (%d,%d)\n", i, j);
			}
		}
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			qNode* c = root->children[i * 2 + j];
			for (int ii = 0; ii < 2; ii++) {
				for (int jj = 0; jj < 2; jj++) {
					qNode* cc = c->children[ii* 2 + jj];
					if (cc->xNeighbors[1- ii] != c->children[(1-ii) * 2 + jj]) {
						printf("WRONG INNER Y NEIGHBOR of %d    - (%d, %d): (%d, %d)\n", cc->id, i, j,  ii, jj);
						printf("expected: %d, actual: %d\n", cc->xNeighbors[1-ii]->id, c->children[(1-ii) * 2 + jj]->id);
					}			
					if (cc->yNeighbors[1-jj] != c->children[ii * 2 + (1-jj)]) {
						printf("WRONG INNER X NEIGHBOR of %d    - (%d, %d): (%d, %d)\n", cc->id, i, j,  ii, jj);
						printf("expected: %d, actual: %d\n", cc->yNeighbors[1-jj]->id, c->children[ii * 2 + (1-jj)]->id);
					}

					if (i == ii) {
						if (cc->xNeighbors[i] != NULL) {
							printf("null x neighbor incorrect: (%d,%d)(%d,%d)\n", i, j, ii, jj);
						}
					} else {
						qNode* shouldNeighbor = root->children[(1-i) * 2 + j]->children[(1-ii) * 2 + jj];
						if (cc->xNeighbors[ii] != shouldNeighbor) {
							printf("cross-node neighbor of %d    incorrect: (%d,%d)(%d,%d)\n", cc->id, i, j, ii, jj);
							printf("expected: %d, actual: %d\n", cc->xNeighbors[ii]->id, shouldNeighbor->id);
						}
					}
					if (j == jj) {
						if (cc->yNeighbors[j] != NULL) {
							printf("null y neighbor incorrect: (%d,%d)(%d,%d)", i, j, ii, jj);
						}
					} else {
						qNode* shouldNeighbor = root->children[i * 2 + (1-j)]->children[ii * 2 + (1-jj)];
						if (cc->yNeighbors[jj] != shouldNeighbor) {
							printf("cross-node neighbor of %d    incorrect: (%d,%d)(%d,%d)\n", cc->id, i, j, ii, jj);
							printf("expected: %d, actual: %d\n", cc->yNeighbors[jj]->id, shouldNeighbor->id);
						}

					}
				}
			}
		}
	}
	
	printf("done\n");

	visc = 1.0;
	kS = 1.0;
	aS = 1.0; // TODO configurable

	sourceX = size / 2;
	sourceY = size / 2;

	sourceVal = 1.0;
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
	size = 50;
	for (int i = 1; i < argc; i++) {
		char* arg = argv[i];
		if (!strcmp("--headless", arg)) {
			headless = true;
		} else if (!strcmp("-size", arg)) {
			size = atoi(argv[++i]);
		}
	}
	printf("headless: %d, size: %d\n", headless, size);

	// seed random
	srandom(time(NULL));

	initSim();

	glutCreateWindow("Vertex, Primitive & Color");  // Create window with the given title
	glutInitWindowSize(320, 320);   // Set the window's initial width & height
	glutInitWindowPosition(50, 50); // Position the window's initial top-left corner
	glutDisplayFunc(display);       // Register callback handler for window re-paint event
	glutKeyboardFunc(keyboard);
	initGL();                       // Our own OpenGL initialization
	glutMainLoop();                 // Enter the event-processing loop
	return 0;
}
