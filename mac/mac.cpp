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
int size;
bool drawVelocity;
bool drawCells;
bool screenshot;
int numToRun = 0;

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
enum PoissonTestFunc { POISSONFUNCNONE, POISSONXY, POISSONCOS, POISSONCOSKL };
PoissonTestFunc poissonTestFunc = POISSONFUNCNONE;

double dt = .03;

// drawing stuff
double minP;
double maxP;
double maxMag;


// navigating 
double deltas[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

struct Cell {
	// cell-centered
	double p, dp, R, divV, temp; // pressure ,x velocity, y velocity, pressure correction, residual, level set value
	// velocity (vx, vy) on faces, cvx, cvy on corners
	double vx, vy, cvx, cvy;
};

Cell* grid;
Cell* oldGrid;

/* Initialize OpenGL Graphics */
void initGL() {
   // Set "clearing" or background color
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
}

void draw() {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {

			double x = -1.0 + 2.0*j/size;
			double y = -1.0 + 2.0*i/size;
			double redPercent = 0.0;
			double bluePercent = 0.0;
			if (grid[i*size+j].p > 0 && maxP > 0) {
				redPercent = grid[i*size+j].p/maxP;
			}
			if (grid[i*size+j].p < 0 && minP < 0) {
				bluePercent = grid[i*size+j].p/minP;
			}
			glColor3f(redPercent, 0, bluePercent);
			glRectf(x, y, x + 2.0/size, y + 2.0/size);

			if (drawVelocity) {
				if (maxMag >= eps) {
					glColor3f(0.0, 1.0, 0.0);
					double vx = grid[i*size+j].vx;
					double vy = grid[i*size+j].vy;
					if (vx >= eps || vy > eps) {
						glBegin(GL_LINES);
						// max size is 1.0/side length, scaled by the max magnitude
						double scale = maxMag * size;
						glVertex2f(x, y + 1.0/size);
						glVertex2f(x + vx / scale, y + 1.0/size);
						glVertex2f(x + 1.0/size, y);
						glVertex2f(x + 1.0/size, y + vy/scale);
						glEnd();
					}
				}
			}
			// draw cell
			if (drawCells) {
				glColor3f(1.0, 1.0, 1.0);
				glBegin(GL_LINE_LOOP);
				glVertex2f(x, y);
				glVertex2f(x, y + 2.0/size);
				glVertex2f(x + 2.0/size, y + 2.0/size);
				glVertex2f(x + 2.0/size, y);
				glEnd();
			}
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

		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				minP = fmin(minP, grid[i*size+j].p);
				maxP = fmax(maxP, grid[i*size+j].p);
				maxMag = fmax(maxMag, fabs(grid[i*size+j].vx));
				maxMag = fmax(maxMag, fabs(grid[i*size+j].vy));
			}
		}

		draw();
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
// TODO velocity interp tests


double computeResidual() {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double faceGradSum = 0.0;
			for (int k = 0; k < 4; k++) {

				int newi = i + deltas[k][0];
				int newj = j + deltas[k][1];
				double neighborp = grid[i*size+j].p;
				if (newi >= 0 && newi < size && newj >= 0 && newj < size) {
					neighborp = grid[newi*size+newj].p;
				}
				
				// integral around the edge of flux, or side length * face gradient
				faceGradSum += size * (neighborp - grid[i*size+j].p);
			}
			
			// double flux = 1/ha * faceGradSum = 1/(1/size * 1) * faceGradSum = size * faceGradSum;
			double flux = size * faceGradSum;
			double R = grid[i*size+j].divV - flux;
			
			grid[i*size+j].R = R;
			//printf("at [%d][%d], divV: %f, flux: %f, R: %f\n", node->i, node->j, node->divV, flux, R);
			// if a*R > e, not done
			if (fabs(R) > eps) {
				doneVCycle = false;
				//printf("more work to do at [%d][%d], %f is bigger than epsilon of %f\n", i, j, grid[d][i*size+j].R, eps);
			} else {
				//printf("done with this cell already, %f is smaller than epsilon of %f\n", fabs(grid[d][i*size+j].R), eps);
			}
		}
	}
}

// puts the new value of dp in temp
bool relaxIter() {
	bool done = true;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {	
			double oldDif = grid[i*size+j].dp;
			double aSum2 = 0.0;
			double bSum2 = 0.0;
			for (int k = 0; k < 4; k++) {
				double dp = grid[i*size+j].dp;
				int newi = i + deltas[k][0];
				int newj = j + deltas[k][1];
				if (newi >= 0 && newi < size && newj >= 0 && newj < size) {
					dp = grid[newi*size+newj].dp;
				}
           		double h = 1.0 / size;
       		    aSum2 -= 1/h;
           		bSum2 += dp/h;
			}

			// TODO make sure scaling residual right
			grid[i*size+j].temp = (grid[i*size+j].R/size/size - bSum2)/aSum2;
	
			double diff = oldDif - grid[i*size+j].temp;
	
	        //printf("relaxing[%d][%d][%d]: aSum: %f, bSumL %f, R: %f, result: %f, diff from old: %f\n", d, i, j, aSum, bSum, grid[d][i*size+j].R, grid[d][i*size+j].temp, diff);
			if (fabs(diff) > eps) {
				done = false;
			}
		}
	}
	return done;
}

void relax(int r) {

	bool done = false;
    int totalCycles = 0;
	while (r-- > 0/* && !done*/) {
        totalCycles++;
		done = relaxIter();
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				grid[i*size+j].dp = grid[i*size+j].temp;
			}
		}
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

std::pair<double, double> getVelocityAt(Cell* g, double x, double y) {
	if (x < 0) x = 0;
	if (x > 1) x = .99999;
	if (y < 0) y = 0;
	if (y > 1) y = .99999;
	// get cell these values are in

	int i = (int)(y * size);
	int j = (int)(x * size);

	double di = (y*size)-i;
	double dj = (x*size)-j;
	
	// do bilinear interpolation
	
	// x
	double a, b, c, d;
	if (j == 0) {
		a = 0.0;
		c = 0.0;
	} else {
		a = g[i*size+j].cvx;
		if (i == size-1) // extrapolate
			c = 2*g[i*size+j].vx-a;
		else
			c = g[(i+1)*size+j].cvx;
	}
	if (j == size-1) {
		b = 0.0;
		d = 0.0;
	} else {
		if (i == 0) {
			b = -g[(i+1)*size+j+1].cvx + 2*g[i*size+j+1].vx;
		} else {
			b = g[i*size+j+1].cvx;	
		}
		if (i == size-1) {
			d = 2*g[i*size+j+1].vx - b;
		} else {
			d = g[(i+1)*size+j+1].cvx;
		}
	}
	double newvx = bilinearInterpolation(a, b, c, d, di, dj);
	

	// y
	if (i == 0) {
		a = 0.0;
		b = 0.0;
	} else {
		a = g[i*size+j].cvy;
		if (j == size - 1)
			b = 2*g[i*size+j].vy - a;
		else
			b = g[i*size+j+1].cvy;
	}
	if (i == size - 1) {
		c = 0.0;
		d = 0.0;
	} else {
		if (j == 0) {
			c = -g[i*size+j+1].cvy + 2*g[(i+1)*size+j].vy;
		} else {
			c = g[(i+1)*size+j].cvy;
		}
		if (j == size-1) {
			d = 2*g[(i+1)*size+j].vy - c;
		} else {
			d = g[(i+1)*size+j+1].cvy;
		}
	}
	double newvy = bilinearInterpolation(a, b, c, d, di, dj);

	return std::make_pair(newvx, newvy);
}

void getSemiLagrangianLookback(Cell* g, double* x, double*y) {
	int steps = 1;
	double newdt = dt/steps;
	while (steps -- > 0) {
		std::pair<double, double> newvel = getVelocityAt(g, *x, *y);
		*x -= newvel.first * newdt;
		*y -= newvel.second * newdt;
	}
}

/*
 * Copies old quantities and sets new velocity quantities.
 * First, computes centered nodal velocities on old grid by averaging face velocities
 * Then, does semi-lagrangian advection, reading from old grid and writing to new grid
 * Finally, averages nodal velocities back to cells
 */
void advectAndCopy() {
	// compute nodal velocities except on edges
	for (int i = 1; i < size; i++) {
		for (int j = 1; j < size; j++) {
			oldGrid[i*size+j].cvx = (oldGrid[i*size+j].vx  + oldGrid[(i-1)*size+j].vx)/2;
			oldGrid[i*size+j].cvy = (oldGrid[i*size+j].vy  + oldGrid[i*size+j-1].vy)/2;
		}
	}
	
	// computes SL nodal velocities, writes them to new grid nodal velocities
	// also does face stuff for nodes on edges with no nodal velocity
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double x = float(j)/size;
			double y = float(i)/size;
			getSemiLagrangianLookback(oldGrid, &x, &y);
			std::pair<double, double> newvel = getVelocityAt(oldGrid, x, y);
			grid[i*size+j].cvx = newvel.first;
			grid[i*size+j].cvy = newvel.second;

		}
	}

	// set new face velocities
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double vx2, vy2;
			if (i == size-1) {
				// extrapolate and set x value
				double x = float(j)/size;
				double y = float(i+1)/size;
				getSemiLagrangianLookback(oldGrid, &x, &y);
				std::pair<double, double> newvel = getVelocityAt(oldGrid, x, y);
				vx2 = newvel.first;
			} else {
				vx2 = oldGrid[(i+1)*size+j].cvx;
			}
			if (j == size-1) {
				double x = float(j+1)/size;
				double y = float(i)/size;
				getSemiLagrangianLookback(oldGrid, &x, &y);
				std::pair<double, double> newvel = getVelocityAt(oldGrid, x, y);
				vy2 = newvel.second;
			} else {
				vy2 = oldGrid[i*size+j+1].cvy;
			}
			grid[i*size+j].vx = (grid[i*size+j].cvx + vx2)/2.0;
			grid[i*size+j].vy = (grid[i*size+j].cvy + vy2)/2.0;
		}
	}
	
}

void correctPressure() {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			grid[i*size+j].p += grid[i*size+j].dp;
		}
	}
}

void project() {
	// correct velocity with updated pressure field to make non-divergent
	for (int i = 1; i < size; i++) {
		for (int j = 1; j < size; j++) {
			double pxGrad = (grid[i*size+j].p - grid[i*size+j-1].p)*size;
			double pyGrad = (grid[i*size+j].p - grid[(i-1)*size+j].p)*size;
			grid[i*size+j].vx -= pxGrad;
			grid[i*size+j].vy -= pyGrad;
		}
	}
}
void eulerAdvectParticle(Particle& p, std::pair<double, double> v) {
	p.x += v.first * dt;
	p.y += v.second * dt;
}

void advectParticles() {
	// TODO implement
	return;
	/*for (int i = 0; i < numParticles; i++) {
		std::pair<double, double> vGrad = getLeaf(root, particles[i].x, particles[i].y, levels-1)->getVelocityAt(particles[i].x, particles[i].y);
		if (particleAlgorithm == EULER) {
			eulerAdvectParticle(particles[i], vGrad);
		}
	}*/
}

double runVCycle() {
	//relax(0, MAX_RELAX_COUNT);
	// do not relax level 0 since it makes no sense to do so...

	relax(MAX_RELAX_COUNT);

	// printf("correcting pressure\n");
	// use corrections to improve pressure
	correctPressure();

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
	return computeResidual();

	// TODO REMOVE
	// doneVCycle = true;
	//if (numCycles == 1) break;

	//printf("end v cycle\n");

}

void poissonAverageR(double* total) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			*total += grid[i*size+j].R/size/size;
		}
	}
}

void poissonCorrectR(double K) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			grid[i*size+j].R -= K;
		}
	}
}

double getMaxR() {
	double max = 0.0;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			max = fmax(max, fabs(grid[i*size+j].R));
		}
	}
	return max;
}

double fixAndCheckResidual() {
	double avgR = 0.0;
	poissonAverageR(&avgR);
	poissonCorrectR(avgR);

	double newMaxR = getMaxR();
	if (newMaxR < eps) {
		doneVCycle = true;
	}
	return newMaxR;
}

void computeVelocityDivergence() {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double a,b;

			// x
			if (j == 0)
				a = 0;
			else
				a = grid[i*size+j].vx;
			if (j == size-1)
				b = 0;
			else
				b = grid[i*size+j+1].vx;
			grid[i*size+j].divV = b-a;


			// y
			if (i == 0)
				a = 0;
			else
				a = grid[i*size+j].vy;
			if (i == size - 1)
				b = 0;
			else
				b = grid[(i+1)*size+j].vy;				

			grid[i*size+j].divV += b-a;
		}
	}
}

void projectionCheck(double* max, double* avg) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			//printf("divV after projection: %f\n", node->divV);
			*max = fmax(*max, fabs(grid[i*size+j].divV));
			*avg += fabs(grid[i*size+j].divV)/size/size;
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
	//std::swap(root, oldRoot);

	double totalTime  = 0.0;

	// compute velocity divergence for advection calculation
	// TODO fix?
	
	// advect velocity, copy old stuff
	printf("advecting velocity\n");

	advectAndCopy();
	
	printf("computing divV\n");

	// grid's vx and vy are now provisional velocity
	// compute velocity divergence of new grid to use in residual calcuation
	computeVelocityDivergence();

	printf("starting poisson solver\n");


	// compute all residuals, so that the whole bottom multilevel is computed
	doneVCycle = true;
	computeResidual();

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
	

	double max = 0.0;
	double avg = 0.0;

	computeVelocityDivergence();
	projectionCheck(&max, &avg);
	printf("max divV before project: %f, avg divV before project: %f\n", max, avg);

	project();

	max = 0.0;
	avg = 0.0;

	computeVelocityDivergence();
	projectionCheck(&max, &avg);
	printf("max divV afer project: %f, avg divV after project: %f\n", max, avg);


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

void initLeftUniform(int i, int j) {
	grid[i*size+j].p = 0.0;
	grid[i*size+j].vx = 0.0;
	grid[i*size+j].vy = 0.0;
	if (i == size/2 || i == size/2-1) {
		grid[i*size+j].vx = -1.0;
	}	
}

void initSink(int i, int j) {
	grid[i*size+j].p = 1.0/size/size;

	double center = 0.5;
	grid[i*size+j].vx = center - (j + 0.5)/size;
	grid[i*size+j].vy = center - (i + 0.5)/size;
	double len = sqrt(grid[i*size+j].vx * grid[i*size+j].vx + grid[i*size+j].vy * grid[i*size+j].vy);
	grid[i*size+j].vx /= len;
	grid[i*size+j].vy /= len;
}

void initSrcSink(int i, int j) {
	grid[i*size+j].p = 1.0/size/size;
	
	double src = .125;
	double sink = .875;

	double vxsrc = (j + 0.5)/size - src;
	double vysrc = (i + 0.5)/size - src;
	double lensrc = sqrt(vxsrc * vxsrc + vysrc * vysrc);

	double vxsink = sink - (j + 0.5)/size;
	double vysink = sink - (i + 0.5)/size;
	double lensink = sqrt(vxsink * vxsink + vysink * vysink);
	
	grid[i*size+j].vx = vxsrc/lensrc + vxsink/lensink;
	grid[i*size+j].vy = vysrc/lensrc + vysink/lensink;
}

void initPoissonTest(int i, int j) {
	
	grid[i*size+j].vx = 0;
	grid[i*size+j].vy = 0;
	
	grid[i*size+j].p = 0;
	
	//node->p = (.5-x)*(.5-x) + (.5-y)*(.5-y); // use pressure gradient split to make multilevel
	
}

void initGrid() {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (startState == LEFT || startState == ADAPTTEST || startState == ERRORTEST)
				initLeftUniform(i, j);
			else if (startState == SINK)
				initSink(i, j);
			else if (startState == SRCSINK)
				initSrcSink(i, j);
			else if (startState == POISSONTEST)
				initPoissonTest(i, j);
			else
				printf("invalid start state\n");
		}
	}
}

void poissonReset() {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			// set pressure to integral over domain
			double x = (j + 0.5)/size;
			double y = (i + 0.5)/size;
			if (poissonTestFunc == POISSONXY) {
				grid[i*size+j].p = .2667;
				grid[i*size+j].divV = 96 * (2*x-1) * (y-1)*(y-1) * y*y  +  32 * (x-1)*(x-1) * (2*x+1) * (1 - 6*y + 6*y*y);
			} else if (poissonTestFunc == POISSONCOS) {
				grid[i*size+j].p = 1.0;
				double cx = cos(M_PI*2*x);
				double cy = cos(M_PI*2*y);
				grid[i*size+j].divV = 4*M_PI*M_PI * (cx + cy - 2*cx*cy);
			} else if (poissonTestFunc == POISSONCOSKL) {
				grid[i*size+j].p = 0.0;
				int k = 2;
				int l = 2;
				grid[i*size+j].divV = -M_PI*M_PI*(k*k + l*l)*cos(M_PI*k*x)*cos(M_PI*l*y);
			} else {
				assert(false);
			}
		}
	}
}

void computePoissonError(double* total) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double x = (j + 0.5)/size;
			double y = (i + 0.5)/size;
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
			*total += fabs(grid[i*size+j].p - correct)/size/size;
		}
	}
}

void runPoissonTest() {
	// adapt
	assert(poissonTestFunc != POISSONFUNCNONE);
	poissonReset();
	
	double initialResidual = computeResidual();
	printf("initial residual: %f\n", initialResidual);

	// try manually fixing residual
	double avgR = 0.0;
	poissonAverageR(&avgR);
	printf("new average residual: %f\n", avgR);
	poissonCorrectR(avgR);
	
	// run V cycles, compute stuff
	int i = 0;

	double newR, oldR;
	newR = initialResidual;

	while (!doneVCycle) {
		i++;
		oldR = newR;
		newR = runVCycle();
		double total = 0.0;
		//double K = 0.0;
		//poissonAverage(root, &K);
		double avgError = 0.0;
		computePoissonError(&avgError);
		printf("average error after %d vcycles: %f\n", i, avgError);

		// try manually fixing residual
		avgR = 0.0;
		poissonAverageR(&avgR);
		poissonCorrectR(avgR);
		avgR = 0.0;
		poissonAverageR(&avgR);
		//doneVCycle |= fabs(newR-oldR) < eps/100;
		double newR = getMaxR();

		printf("residual after %d vcycles: %f\n", i, newR);
		doneVCycle = newR < eps;
	}
	
	printf("poisson test took %d vcycles\n", i);

	//printPressure();

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			printf(" %.4f", grid[i*size+j].p);
		}
		printf("\n");
	}
}




void initSim() {
	printf("init sim\n");

	grid = new Cell[size*size];
	oldGrid = new Cell[size*size];

	initGrid();

	if (startState == POISSONTEST) {
		runPoissonTest();
	}

	if (particleAlgorithm != PARTICLENONE) {

		particles = new Particle[numParticles];
		std::default_random_engine generator;
		double distmin = 0.0;
		double distmax = 1.0;

		if (startState == LEFT) {
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
	size = 32;
	drawCells = false;
	screenshot = false;
	particleAlgorithm = PARTICLENONE;
	poissonTestFunc = POISSONFUNCNONE;
	startState = LEFT;
	for (int i = 1; i < argc; i++) {
		char* arg = argv[i];
		if (!strcmp("--headless", arg)) {
			headless = true;
		} else if (!strcmp("-size", arg)) {
			size = atoi(argv[++i]);
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
		}
	}
	printf("headless: %d, size: %d\n", headless, size);

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
	windowid = glutCreateWindow("MAC Grid");  // Create window with the given title
	glutDisplayFunc(display);       // Register callback handler for window re-paint event
	glutKeyboardFunc(keyboard);
	initGL();                       // Our own OpenGL initialization
	glutMainLoop();                 // Enter the event-processing loop
	return 0;
}
