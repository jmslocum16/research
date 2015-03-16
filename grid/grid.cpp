#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#include <bitmap.h>

#include <algorithm>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <string.h>
 
/* Initialize OpenGL Graphics */
void initGL() {
	// Set "clearing" or background color
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
}

// options
int windowWidth = 640;
int windowHeight = 640;

bool headless; // run without ui for profiling
int size; // number of simulation grids
bool water;
bool drawVelocity;
bool screenshot;
int frameNumber = 0;

double dt = .05;


// GRID[R][C] maps to first quadrant by R->Y, C->X!!

struct Cell {
	double p, vx, vy, dp, R, phi; // pressure ,x velocity, y velocity, pressure correction, residual, level set value
};

struct Cell* grid;
struct Cell* oldGrid;

bool doneVCycle = false;

int MAX_RELAX_COUNT = 20;
int OPTIMIZED_RELAX_COUNT = 4;


double eps = .0001;

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

void statScreenshotDir() {
	 // meh
}


/* Handler for window-repaint event. Call back when the window first appears and
   whenever the window needs to be re-painted. */
void display() {
	glClear(GL_COLOR_BUFFER_BIT);   // Clear the color buffer with current clearing color

	printf("display\n");

	if (water) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (grid[i*size+j].phi < 0) continue;
				double val = fmin(1.0, grid[i*size+j].phi);
				glColor3f(0.0, 0.0, val);
				
				double x = -1.0 + (2.0 * j) / size;
				double y = -1.0 + (2.0 * i) / size;
				glRectf(x, y, x + 2.0/size, y + 2.0/size);
			}
		}
	} else {
		double minP = 100000000.0;
		double maxP = -100000000.0;
	
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				minP = std::min(minP, grid[i*size+j].p);
				maxP = std::max(maxP, grid[i*size+j].p);
			}
		}
		
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				double redPercent = 0.0;
				double bluePercent = 0.0;
				if (grid[i*size+j].p > 0) {
					redPercent = grid[i*size+j].p/maxP;
				}
				if (grid[i*size+j].p < 0) {
					bluePercent = grid[i*size+j].p/-minP;
				}
				
				//printf("percent: %lf\n", percent);
				glColor3f(redPercent, 0, bluePercent);
				double x = -1.0 + (2.0 * j) / size;
				double y = -1.0 + (2.0 * i) / size;
				//double y = 1.0 - (2.0 * i) / size; 
				glRectf(x, y, x + 2.0/size, y + 2.0/size);
			}
		}
	}

	if (drawVelocity && size <= 40) { // if size is > 40 this is pretty useless
		double maxMag = 0.0;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				double x = grid[i*size+j].vx;
				double y = grid[i*size+j].vy;
				maxMag = std::max(maxMag, sqrt(x*x + y*y));
			}
		}

		if (maxMag >= eps) {
			glColor3f(0.0, 1.0, 0.0);
			glBegin(GL_LINES);
			for (int i = 0; i < size; i++) {
				for (int j = 0; j < size; j++) {
					double vx = grid[i*size+j].vx;
					double vy = grid[i*size+j].vy;
					double mag = sqrt(vx*vx + vy*vy);
					if (mag < eps) continue;
					// max size is 1.0/side length, scaled by the max magnitude
					double x = -1.0 + (2.0*j) / size + 1.0/size;
					double y = -1.0 + (2.0*i) / size + 1.0/size;
					glVertex2f(x, y);
					double scale = maxMag * size;
					glVertex2f(x + vx / scale, y + vy /scale);
				}
			}
			glEnd();
		}
	}

	glFlush();  // Render now

	// capture screen if necessary
	if (screenshot) {
		saveScreen();
	}

	frameNumber++;
}

double deltas[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

std::pair<double, double> getPressureGradient(Cell* g, int i, int j) {
	double vals[4];	
	for (int k = 0; k < 4; k++) {
		int newi = i + deltas[k][0];
		int newj = j + deltas[k][1];
		if (newi < 0 || newi >= size || newj <0 || newj >= size) {
			vals[k] = g[i*size+j].p;
		} else {
			vals[k] = g[newi*size+newj].p;
		}
	}
	return std::make_pair((vals[0] - vals[1])/(2*size), (vals[2]-vals[3])/(2*size));
}

std::pair<double, double> getLevelSetGradient(Cell* g, int i, int j) {
	double vals[4];	
	for (int k = 0; k < 4; k++) {
		int newi = i + deltas[k][0];
		int newj = j + deltas[k][1];
		if (newi < 0 || newi >= size || newj <0 || newj >= size) {
			vals[k] = g[i*size+j].phi;
		} else {
			vals[k] = g[newi*size+newj].phi;
		}
	}
	return std::make_pair((vals[0] - vals[1])/(2*size), (vals[2]-vals[3])/(2*size));
}

double determinant(double col1[], double col2[], double col3[]) {
	//printf("determinant of:\n");
	//printf("%f %f %f\n", col1[0], col2[0], col3[0]);
	//printf("%f %f %f\n", col1[1], col2[1], col3[1]);
	//printf("%f %f %f\n", col1[2], col2[2], col3[2]);
	double first = col1[0] * (col2[1] * col3[2] - col2[2] * col3[1]);
	double second = -col1[1] * (col2[0] * col3[2] - col2[2] * col3[0]);
	double third = col1[2] * (col2[0] * col3[1] - col2[1] * col3[0]);
	//printf("first: %f\n", first);
	//printf("second: %f\n", second);
	//printf("third: %f\n", third);
	//printf("is: %f\n", first + second + third);
	return first + second + third;
}

// gets derivative at middle point as if they are a parabola
// solves system of equations for a, b, c
/*
 * y1 = a * x1^2 + b * x1 + c
 * y2 = a * x2^2 + b * x2 + c
 * y3 = a * x3^2 + b * x3 + c
 *
 * and derivative at x2 is 2 * a * x2 + b
 */
double quadraticDerivative(double x1, double y1, double x2, double y2, double x3, double y3) {
	double col1[3] = {x1*x1, x2*x2, x3*x3};
	double col2[3] = {x1, x2, x3};
	double col3[3] = {1, 1, 1};
	double col4[3] = {y1, y2, y3};
	
	// cramer's rule
	double D = determinant(col1, col2, col3);
	double a = determinant(col4, col2, col3) / D;
	double b = determinant(col1, col4, col3) / D;
	// double c = determinant(col1, col2, col4) / D;
	return 2 * a * x2 + b;
}

void testDeterminant() {
	double col1[3] = {2, 1, 1};
	double col2[3] = {1, -1, 2};
	double col3[3] = {1, -1, 1};
	double col4[3] = {3, 0, 0};
	
	double det0 = determinant(col1, col2, col3);
	if (det0 != 3.0) {
		printf("det0 is %f\n", det0);
	}
	double det1 = determinant(col4, col2, col3);
	if (det1 != 3.0) {
		printf("det1 is %f\n", det1);
	}
	double det2 = determinant(col1, col4, col3);
	if (det2 != -6) {
		printf("det2 is %f\n", det2);
	}
	double det3 = determinant(col1, col2, col4);
	if (det3 != 9) {
		printf("det3 is \n", det3);
	}
	printf("done testing determinant\n");
}

// returns better but slower approximation
double getQuadraticVelocityDivergence(double vals[], Cell& center) {
	double dx = quadraticDerivative(-1, vals[2], 0, center.vx, 1, vals[3]);
	double dy = quadraticDerivative(-1, vals[0], 0, center.vy, 1, vals[1]);
	return dx + dy;
}

// returns shitty approximation
double getLinearVelocityDivergence(double vals[], Cell& center) {
	//double dx = (vals[2] - vals[3]) / (size + sizes[2]/2.0 + sizes[3]/2.0);
	//double dy = (vals[0] - vals[1]) / (size + sizes[0]/2.0 + sizes[1]/2.0);
	double dx = (vals[2] - vals[3]) / (2.0/size);
	double dy = (vals[0] - vals[1]) / (2.0/size);
	return dx + dy;
}

double getVelocityDivergence(Cell* g, int i, int j) {
	double vals[4];	
	for (int k = 0; k < 4; k++) {
		int newi = i + deltas[k][0];
		int newj = j + deltas[k][1];
		if (newi < 0 || newi >= size || newj <0 || newj >= size) {
			vals[k] = (k < 2) ? g[i*size+j].vy : g[i*size+j].vx;
			continue;
		}
		Cell c = g[newi*size+newj];
		vals[k] = (k < 2) ? c.vy : c.vx;
	}

	double lin = getLinearVelocityDivergence(vals, g[i*size+j]);
	double quad = getQuadraticVelocityDivergence(vals, g[i*size+j]);
	//printf("velocity divergence of grid[%d][%d]: linear: %f, quadratic: %f\n", i,j, lin, quad);
	return lin;
	//return quad;
}

void computeResidual() {
	//printf("computing residual\n");
	
	doneVCycle = true;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			// compute it: R = div(vx, vy) - 1/(ha)*sum of (s * grad) for each face
			double faceGradSum = 0.0;
			for (int k = 0; k < 4; k++) {
				int newi = i + deltas[k][0];
				int newj = j + deltas[k][1];
				if (newi < 0 || newi >= size || newj <0 || newj >= size) continue; // 0 face gradient at boundary condition
				// integral around the edge of flux, or side length * face gradient
				//faceGradSum += 1 * (grid[newi*size+newj].p - grid[i*size+j].p)/ (1.0/size); // s is always 1
				faceGradSum += size * (grid[newi*size+newj].p - grid[i*size+j].p);
			}
			// h = length of cell = 1.0/size
			// a = "fluid volume fraction of the cell". Since no boundaries cutting through cell, always 1
			
			double divV = getVelocityDivergence(grid, i, j);
			// double flux = 1/ha * faceGradSum = 1/(1/size * 1) * faceGradSum = size * faceGradSum;
			double flux = size * faceGradSum;
			double R = divV - flux;
			
			grid[i*size+j].R = R;

			// if a*R > e, not done
			if (fabs(R) > eps) {
				doneVCycle = false;
				//printf("more work to do here, %f is bigger than epsilon of %f\n", fabs(grid[d][i*size+j].R)/(size*size), eps);
				//printf("for cell[%d][%d], divV: %f, flux: %f, R: %f\n", i, j, divV, flux, R);
			} else {
				//printf("done with this cell already, %f is smaller than epsilon of %f\n", fabs(grid[i*size+j].R), eps);
			}
		}
	}
	//printf("done residual\n");
}

void relaxJacobi(int r) {
	//printf("relaxing %d times. residual matrix:\n", r);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			grid[i*size+j].dp = 0; // TODO what is better initial guess?
			//printf(" %.3f", grid[i*size+j].R);
		}
		//printf("\n");
	}
	
	bool done = false;
	int cycles = 0;
	double newdp[size*size];
	while (/*r-- > 0 && !done */!done) {
		cycles++;
		done = true;
		// Relax(dp, R): dp <-- (sum of adjacent dp - h^2*R)/4
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				double pSum = 0.0;
				double oldDif = grid[i*size+j].dp;
				for (int k = 0; k < 4; k++) {
					int newi = i + deltas[k][0];
					int newj = j + deltas[k][1];
					if (newi < 0 || newi >= size || newj <0 || newj >= size) {
						pSum += grid[i*size+j].dp; // if off, use this pressure
						continue;
					}
					pSum += grid[newi*size+newj].dp;
				}
				
				newdp[i*size+j] = (pSum - ( (grid[i*size+j].R)/(size*size) ) )/4.0;
				double diff = oldDif - newdp[i*size+j];
				if (fabs(diff) > eps) {
					done = false;
					//printf("relaxing[%d][%d]: pSum: %f, R: %f, result: %f, diff from old: %f\n", i, j, pSum, grid[i*size+j].R, newdp[i*size+j], diff);
					//printf("get out\n");
				}
			}
		}
		//printf("dp matrix with %d cycles left: \n", r);
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				grid[i*size+j].dp = newdp[i*size+j];
				//printf(" %.3f", newdp[i*size+j]);
			}
			//printf("\n");
		}
	}
	//printf("relaxing took %d cycles\n", cycles);
}

void relaxGaussSiedel(int r) {
	//printf("relaxing %d times. residual matrix:\n", r);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			grid[i*size+j].dp = 0; // TODO what is better initial guess?
			//printf(" %.3f", grid[i*size+j].R);
		}
		//printf("\n");
	}
	
	bool done = false;
	int cycles = 0;
	while (/*r-- > 0 && !done */!done) {
		cycles++;
		done = true;
		// Relax(dp, R): dp <-- (sum of adjacent dp - h^2*R)/4
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				double pSum = 0.0;
				double oldDif = grid[i*size+j].dp;
				for (int k = 0; k < 4; k++) {
					int newi = i + deltas[k][0];
					int newj = j + deltas[k][1];
					if (newi < 0 || newi >= size || newj <0 || newj >= size) {
						pSum += grid[i*size+j].dp; // if off, use this pressure
						continue;
					}
					pSum += grid[newi*size+newj].dp;
				}
				
				grid[i*size+j].dp = (pSum - ( (grid[i*size+j].R)/(size*size) ) )/4.0;
				double diff = oldDif - grid[i*size+j].dp;
				if (fabs(diff) > eps) {
					done = false;
					//printf("relaxing[%d][%d]: pSum: %f, R: %f, result: %f, diff from old: %f\n", i, j, pSum, grid[i*size+j].R, newdp[i*size+j], diff);
					//printf("get out\n");
				}
			}
		}
	}
	//printf("relaxing took %d cycles\n", cycles);
}

void relax(int r ){
	//relaxJacobi(r);
	relaxGaussSiedel(r);
}

void runStep() {
	printf("running step\n");
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			printf(" %.3f", grid[i*size+j].p);
		}
		printf("\n");
	}
	std::swap(grid, oldGrid);
	
	// advect velocity, copy old stuff
	//printf("advecting velocity\n");
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			// copy pressure, used
			int index = i*size + j;
			grid[index].p = oldGrid[index].p;
			grid[index].vx = oldGrid[index].vx;
			grid[index].vy = oldGrid[index].vy;
			grid[index].phi = oldGrid[index].phi;
			
			// U** = Un - An*dt, An = div(V) * V
			double divV = getVelocityDivergence(oldGrid, i, j);
			double advectX = divV * grid[index].vx;
			double advectY = divV * grid[index].vy;
			//printf("advecting grid[%d][%d][%d] by (%f, %f)\n", k, i, j, advectX, advectY);
			//printf("pressure at [%d][%d][%d]: %f\n", k, i, j, grid[k][i*size+j].p);

			advectX *= dt;
			advectY *= dt;
			grid[index].vx -= advectX;
			grid[index].vy -= advectY;

			// add gravity TODO(make better)
			//grid[index].vy -= .98 * dt;
			// density of water is 1000kg/m^3, so 100kg/m^2?.. density = M/V,so mass = V * density = density/(size+2)
			//grid[index].vy -= dt * 9.8*100.0/size/size;
		}
	}
	// grid's vx and vy are now provisional velocity

	// TODO clamp velocity along border positions?

	// project velocity to get divergent free velocity and pressure
	computeResidual();

	int numCycles = 0;
	
	while (!doneVCycle) {
		numCycles++;
		//printf("start v cycle %d\n", numCycles);

		relax(MAX_RELAX_COUNT);

		//printf("correcting pressure\n");
		// use corrections to improve pressure
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				grid[i*size+j].p += grid[i*size+j].dp;
			}
		}

		// TODO REMOVE
		/*printf("printing pressure after corrections\n");
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				printf("pressure at [%d][%d][%d]: %f\n", k, i, j, grid[k][i*size+j].p);
			}
		}
		*/

		
		// re-compute residual
		computeResidual();

		// TODO REMOVE
		//doneVCycle = true;
		//if (numCycles == 10) break;

		//printf("end v cycle\n");
	}

	// correct velocity with updated pressure field to make non-divergent
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			std::pair<double, double> grad = getPressureGradient(grid, i, j);
			//grid[i*size+j].vx -= grad.first * dt;
			//grid[i*size+j].vy -= grad.second * dt;



			// clamp velocity on boundary condition
			if (i == 0 || i == size-1) {
				grid[i*size+j].vy = 0.0;
			}
			if (j == 0 || j == size-1) {
				grid[i*size+j].vx = 0.0;
			}

		}
	}


	printf("advecting phi\n");
	// dphi/dt + v dot grad(phi) = 0
	// dphi/dt = -v dot grad(phi)
	double dphi[size*size];
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			std::pair<double, double> grad = getLevelSetGradient(grid, i, j);
			printf("level set gradient at [%d][%d] is (%.2f, %.2f), v is (%.2f, %.2f), dphi is %.2f\n", i, j, grad.first, grad.second, grid[i*size+j].vx, grid[i*size+j].vy, dphi[i*size+j]);
			//dphi[i*size+j] = -dt * (grad.first * grid[i*size+j].vx + grad.second * grid[i*size+j].vy);
			dphi[i*size+j] = dt * (-grad.first * grid[i*size+j].vx - grad.second * grid[i*size+j].vy);
		}
	}
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			grid[i * size + j].phi += dphi[i*size+j];
			printf(" %.2f", grid[i*size+j].phi);
		}
		printf("\n");
	}
	
	

	
	printf("done step, took %d iterations\n", numCycles);
}

void initSim() {
	printf("init sim\n");
	grid = new Cell[size*size];
	oldGrid = new Cell[size*size];
	//int start = 2;
	int start = size/2;

	int sinkR = 2;
	int sinkC = 2;

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			//grid[i*size+j].p = i*size+j;
			grid[i*size+j].p = 0;
			if (i == start && j == start) {
				grid[i*size+j].p = 1.0;
			}
			//grid[i*size+j].vx = (sinkC - j)/(1.0*size);
			//grid[i*size+j].vy = (sinkR - i)/(1.0*size);
			if (i == start)
				grid[i*size+j].vx = -1;
			// normalize
			double len = sqrt(grid[i*size+j].vx*grid[i*size+j].vx + grid[i*size+j].vy*grid[i*size+j].vy);
			if (len > 0) {
				grid[i*size+j].vx /= len;
				grid[i*size+j].vy /= len;
			}
			
			//printf(" %.2f", grid[i*size+j].p);
			// phi is signed distance function
			grid[i*size+j].phi = 1 - (abs(i-start) + abs(j-start));
			grid[i*size+j].phi /= size;
			printf (" %.2f", grid[i*size+j].phi);
		}
		printf("\n");
	}
}



void keyboard(unsigned char key, int x, int y) {
	switch (key) {
		case 32:     // space key
			runStep();
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
	size = 5;
	water = false;
	screenshot = false;
	for (int i = 1; i < argc; i++) {
		char* arg = argv[i];
		if (!strcmp("--headless", arg)) {
			headless = true;
		} else if (!strcmp("-size", arg)) {
			size = atoi(argv[++i]);
		} else if (!strcmp("--water", arg)) {
			water = true;
		} else if (!strcmp("--vel", arg)) {
			drawVelocity = true;
		} else if (!strcmp("--screenshot", arg)) {
			screenshot = true;
			statScreenshotDir();
		}
	}
	printf("headless: %d, size: %d\n", headless, size);

	// seed random
	srandom(time(NULL));

	initSim();

	// test determinant
	//testDeterminant();

	// TODO don't do if headless
	glutInitWindowSize(windowWidth, windowHeight);   // Set the window's initial width & height
	glutInitWindowPosition(50, 50);
	glutCreateWindow("Single Grid");  // Create window with the given title
	glutDisplayFunc(display);       // Register callback handler for window re-paint event
	glutKeyboardFunc(keyboard);
	initGL();                       // Our own OpenGL initialization
	glutMainLoop();                 // Enter the event-processing loop
	return 0;
}
