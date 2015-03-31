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
int windowid;
int windowWidth = 640;
int windowHeight = 640;

bool headless; // run without ui for profiling
int levels;
int levelToDisplay; // level to draw
bool water;
bool drawVelocity;
bool drawCells;
bool screenshot;
int frameNumber = 0;
int numToRun = 0;

double dt = .02;

// GRID[R][C] maps to first quadrant by R->Y, C->X!!

struct Cell {
	double p, vx, vy, dp, R, phi, divV, temp; // pressure ,x velocity, y velocity, pressure correction, residual, level set value
	bool used, leaf;
};

struct Cell** grid;
struct Cell** oldGrid;

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

double minP;
double maxP;
double maxMag;


void drawMultilevel(int d, int i, int j, int ml) {
	assert (d <= ml);
	int size = 1<<d;
	if (!grid[d][i*size+j].used) {
		return;
	} else if (d == ml || grid[d][i*size+j].leaf) {
		// draw it
		double x = -1.0 + (2.0 * j) / size;
		double y = -1.0 + (2.0 * i) / size;

		if (water) {
			if (grid[d][i*size+j].phi >= 0) {
				double val = fmin(1.0, grid[d][i*size+j].phi);
				glColor3f(0.0, 0.0, val);
			}
		} else {
			double redPercent = 0.0;
			double bluePercent = 0.0;
			if (grid[levelToDisplay][i*size+j].p > 0) {
				redPercent = grid[levelToDisplay][i*size+j].p/maxP;
			}
			if (grid[levelToDisplay][i*size+j].p < 0) {
				bluePercent = grid[levelToDisplay][i*size+j].p/-minP;
			}
			
			//printf("percent: %lf\n", percent);
			glColor3f(redPercent, 0, bluePercent);
		}

		glRectf(x, y, x + 2.0/size, y + 2.0/size);

		if (drawVelocity && d <= 5) { // if size is > 40 this is pretty useless
			if (maxMag >= eps) {
				glColor3f(0.0, 1.0, 0.0);
				double vx = grid[d][i*size+j].vx;
				double vy = grid[d][i*size+j].vy;
				double mag = sqrt(vx*vx + vy*vy);
				if (mag >= eps) {
					glBegin(GL_LINES);
					// max size is 1.0/side length, scaled by the max magnitude
					double x = -1.0 + (2.0*j) / size + 1.0/size;
					double y = -1.0 + (2.0*i) / size + 1.0/size;
					glVertex2f(x, y);
					double scale = maxMag * size;
					glVertex2f(x + vx / scale, y + vy /scale);
					glEnd();
				}
			}
		}
		// draw cell
		if (drawCells && d <= 5) {
			glColor3f(1.0, 1.0, 1.0);
			glBegin(GL_LINE_LOOP);
			glVertex2f(x, y);
			glVertex2f(x, y + 2.0/size);
			glVertex2f(x + 2.0/size, y+2.0/size);
			glVertex2f(x + 2.0/size, y);
			glEnd();
		}
	} else {
		// draw children
		for (int k = 0; k < 4; k++) {
			drawMultilevel(d + 1, i*2 + (k/2), j*2 + (k%2), ml);
		}
	}
}

/* Handler for window-repaint event. Call back when the window first appears and
   whenever the window needs to be re-painted. */
void display() {
	glClear(GL_COLOR_BUFFER_BIT);   // Clear the color buffer with current clearing color

	printf("display\n");

	int size = 1<<levelToDisplay;

	if (!water) {
		minP = 100000000.0;
		maxP = -100000000.0;

		for (int d = 0; d <= levelToDisplay; d++) {
			int size = 1<<d;
			for (int i = 0; i < size; i++) {
				for (int j = 0; j < size; j++) {
					if (grid[d][i*size+j].leaf || d == levelToDisplay) {
						minP = std::min(minP, grid[d][i*size+j].p);
						maxP = std::max(maxP, grid[d][i*size+j].p);
					}
				}
			}
		}
	}
	if (drawVelocity && levelToDisplay <= 5) {
		maxMag = 0.0;
		for (int d = 0; d <= levelToDisplay; d++) {
			int size = 1<<d;
			for (int i = 0; i < size; i++) {
				for (int j = 0; j < size; j++) {
					if (grid[d][i*size+j].leaf || d == levelToDisplay) {
						double x = grid[d][i*size+j].vx;
						double y = grid[d][i*size+j].vy;
						maxMag = std::max(maxMag, sqrt(x*x + y*y));
					}
				}
			}
		}
	}


	drawMultilevel(0, 0, 0, levelToDisplay);

	
	glFlush();  // Render now

	// capture screen if necessary
	if (screenshot) {
		saveScreen();
	}
}



double deltas[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
// for a given delta, the 2 corners to go to on finer levels when calculating face gradients/laplacian
double corner1[4][2] = {{0, 0}, {1, 0}, {0, 0}, {0, 1}};
double corner2[4][2] = {{0, 1}, {1, 1}, {1, 0}, {1, 1}};

Cell* nullCell = NULL;

// finds the highest level used neighbor at the given multilevel
// returns 2  cells if on the same (higher) level, otherwise returns one and null
// d is level, (i, j), k is the direction (deltas index). Value initially in level is the target level, value returned in level is the level of the neighboring cell
// guaranteed d <= target level, otherwise that wouldn't make any sense..
std::pair<Cell*, Cell*> getNeighborInDir(Cell** g, int d, int i, int j, int k, int* level) {
	int size = 1<<d;
	assert (d <= *level);
	int newi = i + deltas[k][0];
	int newj = j + deltas[k][1];
	if (newi < 0 || newi >= size || newj <0 || newj >= size) {
		// not on grid, use boundary conditions
		*level = d;
		return std::make_pair(&g[d][i*size+j], nullCell);
	} else if (!g[d][newi*size+newj].used) {
		// go up until find the used cell
		int newd = d;
		while (!g[newd][newi * (1<<newd) + newj].used) {
			newi >>= 1;
			newj >>= 1;
			newd -= 1;
		}
		*level = newd;
		return std::make_pair(&g[newd][newi * (1<<newd) + newj], nullCell);
	} else if (d == *level || g[d][newi*size+newj].leaf) {
		// simply just your neighbor
		*level = d;
		return std::make_pair(&g[d][newi * size + newj], nullCell);
	} else {
		int c1i = 2*newi + corner1[k][0];
		int c1j = 2*newj + corner1[k][1];
		
		int c2i = 2*newi + corner2[k][0];
		int c2j = 2*newj + corner2[k][1];
		int cd = d + 1;
		int csize = 1<<cd;
		while (cd < *level && !g[cd][c1i * csize + c1j].leaf && !g[cd][c2i * csize + c2j].leaf) {
			cd++;
			csize <<= 1;
			c1i = 2*c1i + corner2[k][0];
			c1j = 2*c1j + corner2[k][1];
			c2i = 2*c2i + corner1[k][0];
			c2j = 2*c2j + corner1[k][1];
		}
		if (*level == cd || g[cd][c1i * csize + c1j].leaf && g[cd][c2i * csize + c2j].leaf) {
			*level = cd;
			return std::make_pair(&g[cd][c1i * csize + c1j], &g[cd][c2i * csize + c2j]);
		} else if (!g[cd][c1i * csize + c1j].leaf) {
			// terminated because c2 not used anymore. Keep following c1 to the end then return it.
			while (cd < *level && !g[cd][c1i * csize + c1j].leaf) {
				cd++;
				csize <<= 1;
				c1i = 2*c1i + corner2[k][0];
				c1j = 2*c1j + corner2[k][1];
			}
			*level = cd;
			return std::make_pair(&g[cd][c1i * csize + c1j], nullCell);
		} else if (!g[cd][c2i * csize + c2j].leaf) {
			// terminated because c1 not used anymore. Keep following c2 to the end then return it.
			while (cd < *level && !g[cd][c2i * csize + c2j].leaf) {
				cd++;
				csize <<= 1;
				c2i = 2*c2i + corner1[k][0];
				c2j = 2*c2j + corner1[k][1];
			}
			*level = cd;
			return std::make_pair(&g[cd][c2i * csize + c2j], nullCell);
		} else {
			// neither c1 nor c2 used, but both were on previous level, so back both up 1 level and return both
            printf("THIS SHOULD NEVER HAPPEN\n");
		}
	}
}

bool assertNeighbors(Cell* n1, Cell* n2, Cell* realN1, Cell* realN2) {
	if (n2 == NULL) {
		return realN2 == NULL && n1 == realN1;
	} else {
		return (n1 == realN1 && n2 == realN2) || (n1 == realN2 && n2 == realN1);
	}
}

// test multilevel neighbor functions
// construct following multigrid with 4 levels:
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
	// construct grid
	Cell** testGrid = new Cell*[4];
	for (int i = 0; i < 4; i++) {
		int size = 1<<i;
		testGrid[i] = new Cell[size*size]; 
		for (int j = 0; j < size*size; j++) {
			testGrid[i][j].used = false;
            testGrid[i][j].leaf = false;
		}
	}

	testGrid[0][0].used = true;
	// level 1
	testGrid[1][0].used = true;
    testGrid[1][0].leaf = true;
	testGrid[1][1].used = true;
	testGrid[1][2].used = true;
	testGrid[1][3].used = true;
    testGrid[1][3].leaf = true;
	// level 2
	testGrid[2][2].used = true;
	testGrid[2][3].used = true;
    testGrid[2][3].leaf = true;
	testGrid[2][6].used = true;
	testGrid[2][7].used = true;
    testGrid[2][7].leaf = true;
	testGrid[2][8].used = true;
    testGrid[2][8].leaf = true;
	testGrid[2][9].used = true;
    testGrid[2][9].leaf = true;
	testGrid[2][12].used = true;
    testGrid[2][12].leaf = true;
	testGrid[2][13].used = true;
    testGrid[2][14].leaf = true;
	// level 3	
	testGrid[3][4].used = true;
    testGrid[3][4].leaf = true;
	testGrid[3][5].used = true;
    testGrid[3][5].leaf = true;
	testGrid[3][12].used = true;
    testGrid[3][12].leaf = true;
	testGrid[3][13].used = true;
    testGrid[3][13].leaf = true;
	testGrid[3][20].used = true;
    testGrid[3][20].leaf = true;
	testGrid[3][21].used = true;
    testGrid[3][21].leaf = true;
	testGrid[3][28].used = true;
    testGrid[3][28].leaf = true;
	testGrid[3][29].used = true;
    testGrid[3][29].leaf = true;

	// actual tests

	// deltas: {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
	// level 1
	int d = 1;
	std::pair<Cell*, Cell*> n;
	// cell 0's neighbor to the right
	int level = d;
	n = getNeighborInDir(testGrid, d, 0, 0, 2, &level);
	assert (assertNeighbors(&testGrid[1][1], NULL, n.first, n.second));
	assert (level == 1);
	// cell 1's neighbor to the left
	level = d;
	n = getNeighborInDir(testGrid, d, 0, 1, 3, &level);
	assert (assertNeighbors(&testGrid[1][0], NULL, n.first, n.second));
	assert (level == 1);

	// level 2
	d = 2;

	level = d;
	n = getNeighborInDir(testGrid, 1, 0, 0, 2, &level);
	assert (assertNeighbors(&testGrid[2][2], &testGrid[2][6], n.first, n.second));
	assert (level == 2);

	level = d;
	n = getNeighborInDir(testGrid, 2, 0, 2, 3, &level);
	assert (assertNeighbors(&testGrid[1][0], NULL, n.first, n.second));
	assert (level == 1);

	// level 3
	d = 3;

	level = d;
	n = getNeighborInDir(testGrid, 1, 0, 0, 2, &level);
	assert (assertNeighbors(&testGrid[3][12], &testGrid[3][20], n.first, n.second));
	assert (level == 3);

	level = d;
	n = getNeighborInDir(testGrid, 1, 0, 0, 0, &level);
	assert (assertNeighbors(&testGrid[2][8], &testGrid[2][9], n.first, n.second));
	assert (level == 2);

	level = d;
	n = getNeighborInDir(testGrid, 1, 1, 1, 1, &level);
	assert (assertNeighbors(&testGrid[3][29], NULL, n.first, n.second));
	assert (level == 3);

	// clean up
	printf("done testing multilevel neighbors\n");
	for (int i = 0; i < 4; i++) {
		delete testGrid[i];
	}
	delete testGrid;
}

std::pair<double, double> getGradient(double vals[], int sizes[], int size) {
	// grad = a - b / (dist) = (a - b) / ((1/sizea) / 2 + (1/sizeb) / 2 + 1/size) = (a - b) / (.5/sizea + .5/sizeb + size)
	return std::make_pair((vals[2] - vals[3])/(0.5/sizes[2] + 0.5/sizes[3] + 1.0/size), (vals[0]-vals[1])/(0.5/sizes[0] + 0.5/sizes[1] + 1.0/size));	
}

std::pair<double, double> getPressureGradient(Cell** g, int d, int i, int j, int targetLevel) {
	int size = 1<<d;
	double vals[4];	
	int sizes[4];
	double h;
	for (int k = 0; k < 4; k++) {
		int level = targetLevel;
		std::pair<Cell*, Cell*> neighbor = getNeighborInDir(grid, d, i, j, k, &level);
		
		if (neighbor.second == NULL) {
			// only 1
			vals[k] = neighbor.first->p;
		} else {
			vals[k] = (neighbor.first->p + neighbor.second->p)/2.0;
		}
		sizes[k] = 1<<level;
	}
	return getGradient(vals, sizes, size);
}

// TODO
/*std::pair<double, double> getLevelSetGradient(Cell** g, int d, int i, int j) {
	int size = 1<<d;
	double vals[4];
	for (int k = 0; k < 4; k++) {
		int newi = i + deltas[k][0];
		int newj = j + deltas[k][1];
		if (newi < 0 || newi >= size || newj <0 || newj >= size) {
			vals[k] = g[d][i*size+j].phi;
		} else {
			vals[k] = g[d][newi*size+newj].phi;
		}
	}
	return std::make_pair(size * (vals[2] - vals[3])/2, size * (vals[0]-vals[1])/2);
}*/

std::pair<double, double> getVelocityGradient(Cell** g, int d, int i, int j) {
    double vals[4];	
	int sizes[4];

	for (int k = 0; k < 4; k++) {
		int level = levels - 1;
		std::pair<Cell*, Cell*> neighbor = getNeighborInDir(g, d, i, j, k, &level);
		if (neighbor.second == NULL) {
		vals[k] = (k < 2) ? neighbor.first->vy : neighbor.first->vx;
		} else {
			if (k < 2) {
				vals[k] = (neighbor.first->vy + neighbor.second->vy) / 2.0;
			} else {
				vals[k] = (neighbor.first->vx + neighbor.second->vx) / 2.0;
			}
		}	
		sizes[k] = 1<<level;
	}
	
	return getGradient(vals, sizes, 1<<d);
}

std::pair<double, double> getVelocitySingleGradient(Cell** g, int d, int i, int j, bool x) {
    double vals[4];	
	int sizes[4];

	for (int k = 0; k < 4; k++) {
		int level = levels - 1;
		std::pair<Cell*, Cell*> neighbor = getNeighborInDir(g, d, i, j, k, &level);
		if (neighbor.second == NULL) {
		vals[k] = x ? neighbor.first->vx : neighbor.first->vy;
		} else {
			if (x) {
				vals[k] = (neighbor.first->vx + neighbor.second->vx) / 2.0;
			} else {
				vals[k] = (neighbor.first->vy + neighbor.second->vy) / 2.0;
			}
		}	
		sizes[k] = 1<<level;
	}
	
	return getGradient(vals, sizes, 1<<d);
}

// curl(F) = d(Fy)/dx - d(Fx)/dy
double getCurl(Cell** g, int d, int i, int j) {
	std::pair<double, double> xgrad = getVelocitySingleGradient(g, d, i, j, true);
	std::pair<double, double> ygrad = getVelocitySingleGradient(g, d, i, j, false);
	return ygrad.first - xgrad.second;
	
}

void computeVelocityDivergence(Cell** g) {
    for (int d = 0; d < levels; d++) {
		int size = 1<<d;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				for (int k = 0; k < 4; k++) {
				    std::pair<double, double> grad = getVelocityGradient(g, d, i, j);
			    	g[d][i*size+j].divV = grad.first + grad.second;
			   }   
		    }
    	}
    }
}

void computeResidual(int d) {
	printf("computing residual for level %d\n", d);

	int size = 1<<d;
	
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (!grid[d][i*size+j].used) {
				//printf("       ");
				continue;
			} else if (d == levels - 1 || grid[d][i*size + j].leaf) { // if leaf cell, compute residual
				// compute it: R = div(vx, vy) - 1/(ha)*sum of (s * grad) for each face
				double faceGradSum = 0.0;
				for (int k = 0; k < 4; k++) {
					int level = levels - 1; // residual is computed at largest multilevel only.
					std::pair<Cell*, Cell*> neighbor = getNeighborInDir(grid, d, i, j, k, &level);
					double neighborp = (neighbor.second == NULL) ? neighbor.first->p : (neighbor.first->p + neighbor.second->p) / 2.0;
					int neighborsize = 1 << level;
					// integral around the edge of flux, or side length * face gradient
					//faceGradSum += 1 * (neighbor.p - this.p)/ (0.5/size + 0.5/neighborsize);
					faceGradSum += (neighborp - grid[d][i*size+j].p) / (0.5/size + 0.5/neighborsize);
				}
				// h = length of cell = 1.0/size
				// a = "fluid volume fraction of the cell". Since no boundaries cutting through cell, always 1
				
				// double flux = 1/ha * faceGradSum = 1/(1/size * 1) * faceGradSum = size * faceGradSum;
				double flux = size * faceGradSum;
				double R = grid[d][i*size+j].divV - flux;
				
				grid[d][i*size+j].R = R;
				//printf("at [%d][%d], divV: %f, flux: %f, R: %f\n", i, j, grid[d][i*size+j].divV, flux, R);
				// if a*R > e, not done
				if (fabs(R) > eps) {
					doneVCycle = false;
					//printf("more work to do at [%d][%d], %f is bigger than epsilon of %f\n", i, j, grid[d][i*size+j].R, eps);
				} else {
					//printf("done with this cell already, %f is smaller than epsilon of %f\n", fabs(grid[d][i*size+j].R), eps);
				}
			} else {
                grid[d][i*size+j].R = 0.0;
				for (int k = 0; k < 4; k++) {
					grid[d][i*size+j].R += grid[d+1][(2*i+(k/2))*2*size+2*j+(k%2)].R;
				}
				grid[d][i*size+j].R /= 4.0;
			}

			//printf(" %.4f", grid[d][i*size+j].R);
		}
		//printf("\n");
	}
	
	printf("done residual\n");
}


// relax the cell at [d][i][j] at the given multilevel
// puts the new value of dp in temp
bool relaxRecursive(int d, int i, int j, int ml) {
	// d should never be > ml
	assert (d <= ml);
	int size = 1<<d;
	if (!grid[d][i*size+j].used) {
		return true;
	} else if (d == ml || grid[d][i*size+j].leaf) {
        // dp = (h*a*divV - bsum)/asum
		double aSum = 0.0;
        double bSum = 0.0;
        double dp;
        double h;
		double oldDif = grid[d][i*size+j].dp;
		for (int k = 0; k < 4; k++) {
			int level = ml;
			std::pair<Cell*, Cell*> neighbor = getNeighborInDir(grid, d, i, j, k, &level);
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
        grid[d][i*size+j].temp = (grid[d][i*size+j].R/size - bSum) / aSum;
		double diff = oldDif - grid[d][i*size+j].temp;

        //printf("relaxing[%d][%d][%d]: aSum: %f, bSumL %f, R: %f, result: %f, diff from old: %f\n", d, i, j, aSum, bSum, grid[d][i*size+j].R, grid[d][i*size+j].temp, diff);
		if (fabs(diff) > eps) {
			return false;
			//printf("get out\n");
		}
		return true;
	} else {
		// relax children
		bool done = true;
		for (int k = 0; k < 4; k++) {
			done &= relaxRecursive(d + 1, 2*i + (k/2), 2*j + (k%2), ml);
		}
		return done;
	}
}

void relax(int d, int r) {
	int size = 1<<d;
	assert (d < levels);

	printf("relaxing level %d\n", d);
	// get initial gress from previous level, if possible
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (!grid[d][i*size+j].used) {
				continue; // not part of its multilevel, higher cell already has guess
			} else if (d == 0) {
				grid[d][i*size+j].dp = 0; // TODO what is better initial guess?
			} else {
				// inject from higher level
				grid[d][i*size+j].dp = grid[d-1][i/2 * size/2 + j/2].dp;
				//printf("initial guess for [%d][%d].dp is %f\n", i, j, grid[d][i*size+j].dp);
			}
			//printf(" %.3f", grid[d][i*size+j].dp);
			//printf(" %.3f", grid[i*size+j].R);
		}
		//printf("\n");
	}
	bool done = false;
	int maxlevel = d;
    int totalCycles = 0;
	while (r-- > 0 && !done) {
        totalCycles++;
		done = relaxRecursive(0, 0, 0, maxlevel);
		for (d = 0; d <= maxlevel; d++) {
			size = 1<<d;
			for (int i = 0; i < size; i++) {
				for (int j = 0; j < size; j++) {
					if (grid[d][i*size+j].used) {
						grid[d][i*size+j].dp = grid[d][i*size+j].temp;
                        //printf(" %.4f", grid[d][i*size+j].dp);
					}
				}
                //printf("\n");
			}
		}
	}
	//printf("dp matrix with %d cycles left: \n", r);
    printf("relaxing took %d cycles\n", totalCycles);
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
    }
}


void expand(int d, int i, int j) {
    printf("expanding cell [%d][%d][%d]\n", d, i, j);
    std::pair<double, double> pressureGrad = getPressureGradient(grid, d, i, j, levels - 1);
    std::pair<double, double> velocityGrad = getVelocityGradient(grid, d, i, j);
    int size = 1<<d;
    Cell cur = grid[d][i*size+j];
    // std::pair<double, double> levelSetGrad;
    for (int k = 0; k < 4; k++) {
        Cell& child = grid[d+1][(2*i+(k/2))*2*size + 2*j + (k%2)];
        child.used = true;
        child.leaf = true;
        int constX = k%2 == 0 ? -1 : 1;
        int constY = k/2 == 0 ? -1 : 1;
        child.p = cur.p + (constX * .25/size) * pressureGrad.first + (constY * .25/size) * pressureGrad.second;
        child.vx = cur.vx + (constX * .25/size) * velocityGrad.first;
        child.vy = cur.vy + (constY * .25/size) * velocityGrad.second;
        // child.phi = ...
    }
    grid[d][i*size+j].used = true;
	grid[d][i*size+j].leaf = false;
}

void contract(int d, int i, int j) {
    printf("contracting cell [%d][%d][%d]\n", d, i, j);
    int size = 1<<d;
    for (int k = 0; k < 4; k++) {
        Cell& child = grid[d+1][(2*i+(k/2))*2*size + 2*j + (k%2)];
        child.used = false;
        child.leaf = false;
    }
    grid[d][i*size+j].leaf = true;
}

// returns true if leaf should be expanded, false if it should not

double curlThresh = .01;
double curlAdaptFunction(int d, int i, int j) {
	double curl = getCurl(grid, d, i, j);
	return fabs(curl) > curlThresh / (1<<d);
}

double pressureThresh = .01;
bool pGradAdaptFunction(int d, int i, int j) {
    std::pair<double, double> pgrad = getPressureGradient(grid, d, i, j, levels - 1);
    return fabs(pgrad.first + pgrad.second) > pressureThresh/(1<<d);
}

bool adaptFunction(int d, int i, int j) {
    return pGradAdaptFunction(d, i, j);
}

void recursiveAdaptivity(int d, int i, int j) {
    int size = 1<<d;
    if (d == levels - 1 || !grid[d][i*size+j].used) return;
    if (grid[d][i*size+j].leaf) {
        // see if should adapt cell
        if (adaptFunction(d, i, j)) {
            expand(d, i, j);
        }
        return;
    }
    bool allChildrenLeaves = true;
    for (int k = 0; k < 4; k++) {
        if (!grid[d+1][(2*i+(k/2))*2*size + 2*j + (k%2)].leaf) {
            allChildrenLeaves = false;
            break;
        }
    }
    if (allChildrenLeaves) {
        if (!adaptFunction(d, i, j)) {
            // this wouldn't be expanded if it was a leaf, so it shouldn't have leaf children
            contract(d, i, j);
        }
    } else {
        for (int k = 0; k < 4; k++) {
            recursiveAdaptivity(d+1, i*2 + (k/2), j*2 + (k%2));
        }
    }
}

void runStep() {
	printf("running step\n");
	for (int d = 0; d < levels; d++) {
		printf("level %d\n", d);
		int size = 1<<d;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				printf(" %.3f", grid[d][i*size+j].p);
			}
			printf("\n");
		}
	}
	std::swap(grid, oldGrid);


	// compute velocity divergence for advection calculation
	computeVelocityDivergence(oldGrid);
	
	// advect velocity, copy old stuff
	//printf("advecting velocity\n");
	for (int d = 0; d < levels; d++) {
		int size = 1<<d;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// copy pressure, used
				int index = i*size + j;
				grid[d][index].p = oldGrid[d][index].p;
				grid[d][index].vx = oldGrid[d][index].vx;
				grid[d][index].vy = oldGrid[d][index].vy;
				grid[d][index].phi = oldGrid[d][index].phi;
				grid[d][index].used = oldGrid[d][index].used;
                grid[d][index].leaf = oldGrid[d][index].leaf;

				if (!grid[d][index].used) continue;
				
				// U** = Un - An*dt, An = div(V) * V
				
				double advectX = oldGrid[d][index].divV * grid[d][index].vx;
				double advectY = oldGrid[d][index].divV * grid[d][index].vy;
				//printf("advecting grid[%d][%d][%d] by (%f, %f)\n", k, i, j, advectX, advectY);
				//printf("pressure at [%d][%d][%d]: %f\n", k, i, j, grid[k][i*size+j].p);
	
				advectX *= dt;
				advectY *= dt;
				grid[d][index].vx -= advectX;
				grid[d][index].vy -= advectY;
	
				// add gravity TODO(make better)
				//grid[index].vy -= .98 * dt;
				// density of water is 1000kg/m^3, so 100kg/m^2?.. density = M/V,so mass = V * density = density/(size+2)
				//grid[index].vy -= dt * 9.8*100.0/size/size;
			}
		}
	}
	// grid's vx and vy are now provisional velocity
	// compute velocity divergence of new grid to use in residual calcuation
	computeVelocityDivergence(grid);

	// compute all residuals, so that the whole bottom multilevel is computed
	doneVCycle = true;
	computeResidual(levels - 1);
	for (int d = levels - 2; d >= 0; d--) {
			computeResidual(d);
	}

	int numCycles = 0;
	
	while (!doneVCycle) {
		numCycles++;
		//printf("start v cycle %d\n", numCycles);

		
		//relax(0, MAX_RELAX_COUNT);
		grid[0][0].dp = 0;
		// Trying to not relax level 0 since it makes no sense to do so...

		for (int d = 1; d < levels; d++) {
			relax(d, (d==1) ? MAX_RELAX_COUNT : OPTIMIZED_RELAX_COUNT);
		}

		//printf("correcting pressure\n");
		// use corrections to improve pressure
		for (int d = 0; d < levels; d++) {
			int size = 1<<d;
			for (int i = 0; i < size; i++) {
				for (int j = 0; j < size; j++) {
					if (grid[d][i*size+j].used) {
						grid[d][i*size+j].p += grid[d][i*size+j].dp;
					}
				}
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
		doneVCycle = true;
		computeResidual(levels - 1);
		for (int d = levels - 2; d >= 0; d--) {
			computeResidual(d);
		}

		// TODO REMOVE
		// doneVCycle = true;
		//if (numCycles == 1) break;

		//printf("end v cycle\n");
	}
	
    // given new state, do adaptivity
    //recursiveAdaptivity(0, 0, 0);

	// correct velocity with updated pressure field to make non-divergent
	for (int d = 0; d < levels; d++) {
		int size = 1<<d;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (!grid[d][i*size+j].used) continue;
				std::pair<double, double> grad = getPressureGradient(grid, d, i, j, levels - 1);
				grid[d][i*size+j].vx -= grad.first * dt;
				grid[d][i*size+j].vy -= grad.second * dt;
				//grid[d][i*size+j].vx -= grad.first;
				//grid[d][i*size+j].vy -= grad.second;
	
	
	
				// clamp velocity on boundary condition
				if (i == 0 || i == size-1) {
					grid[d][i*size+j].vy = 0.0;
				}
				if (j == 0 || j == size-1) {
					grid[d][i*size+j].vx = 0.0;
				}	
			}
		}
	}


	printf("advecting phi\n");
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
	
	
	printf("done step, took %d iterations\n", numCycles);
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
	}

	// increase frame number here instead of in display so that you get 1 of each frame, not dependent on glut's redrawing, like on alt-tabbing or anything
	frameNumber++;
}

void initSim() {
	printf("init sim\n");
	grid = new Cell*[levels];
	oldGrid = new Cell*[levels];
	for (int d = 0; d < levels; d++) {
		int size = 1<<d;
		grid[d] = new Cell[size*size];
		oldGrid[d] = new Cell[size*size];
	}
	//int level = levels/2;
	int level = levels - 2;

    for (int d = levels - 1; d > level; d--) {
        int size = 1<<d;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                grid[d][i*size+j].used = false;
                grid[d][i*size+j].leaf = false;
                grid[d][i*size+j].p = 0;
                grid[d][i*size+j].vx = 0;
                grid[d][i*size+j].vy = 0;
                grid[d][i*size+j].phi = 0.0;
            }
        }
    }

	int size = 1<<level;
	int start = size/2;
	printf("deepest level: %d\n", level);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			//grid[i*size+j].p = i*size+j;
			grid[level][i*size+j].p = 0;
			if ((i == start-1 || i == start) && (j == start || j == start-1)) {
				grid[level][i*size+j].p = 1.0;
			}
            grid[level][i*size+j].vx = 0.0;
            grid[level][i*size+j].vy = 0.0;
			//grid[i*size+j].vx = (sinkC - j)/(1.0*size);
			//grid[i*size+j].vy = (sinkR - i)/(1.0*size);
			if (i == start || i == start-1)
				grid[level][i*size+j].vx = -1;
			// normalize
			/*double len = sqrt(grid[level][i*size+j].vx*grid[level][i*size+j].vx + grid[level][i*size+j].vy*grid[level][i*size+j].vy);
			if (len > 0) {
				grid[level][i*size+j].vx /= len;
				grid[level][i*size+j].vy /= len;
			}*/

			//printf(" %.2f", grid[i*size+j].p);
			// phi is signed distance function
			grid[level][i*size+j].phi = 1 - (abs(i-start) + abs(j-start));
			grid[level][i*size+j].phi /= size;
			grid[level][i*size+j].used = true;
            grid[level][i*size+j].leaf = true;
			printf (" %.2f", grid[level][i*size+j].p);
		}
		printf("\n");
	}

	for (int d = level - 1; d >= 0; d--) {
		printf("level %d\n", d);
		int size = 1<<d;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
	            grid[d][i*size+j].p = 0.0;
				grid[d][i*size+j].vx = 0.0;
				grid[d][i*size+j].vy = 0.0;
				grid[d][i*size+j].phi = 0.0;
                
				for (int k = 0; k < 4; k++) {
					grid[d][i*size+j].p += grid[d+1][(2*i+(k%2))*2*size + 2*j+k/2].p;
					grid[d][i*size+j].vx += grid[d+1][(2*i+(k%2))*2*size + 2*j+k/2].vx;
					grid[d][i*size+j].vy += grid[d+1][(2*i+(k%2))*2*size + 2*j+k/2].vy;
					grid[d][i*size+j].phi += grid[d+1][(2*i+(k%2))*2*size + 2*j+k/2].phi;
				}
				
				grid[d][i*size+j].p /= 4.0;
				grid[d][i*size+j].vx /= 4.0;
				grid[d][i*size+j].vy /= 4.0;
				grid[d][i*size+j].phi /= 4.0;
				grid[d][i*size+j].used = true;
                grid[d][i*size+j].leaf = false;
				printf(" %.2f", grid[d][i*size+j].p);
			
			}
			printf("\n");
		}
	}

	// MAKE MULTILEVEL
    contract(level - 1, 0, 1<< (level-2));
	//expand(level, (1<<level)-1, 1<<(level-1));

	// clamp starting velocity
	for (int d = 0; d < levels; d++) {
		int size = 1<<d;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (i == 0 || i == size-1) {
					grid[d][i*size+j].vy = 0.0;
				}
				if (j == 0 || j == size-1) {
					grid[d][i*size+j].vx = 0.0;
				}
			}
		}
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

	initSim();

	// test neighbors
	//testMultilevelNeighbors();



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
