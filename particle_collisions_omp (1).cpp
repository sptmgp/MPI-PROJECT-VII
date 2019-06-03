#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include <math.h>       /* asin */
#include <iostream>
#include <fstream>
#include <ctime>

#include <omp.h>


#define NBALLS 1000

using namespace std;
////////////////////////////////////////////////////////////////////////////

  //  Initialize sphere location, velocity, and
  //  coefficient of restitution.
float  e = 0.9;
float  sphereRadius = 5.0;
float  timeIncrement = 0.07;

float tmp;
float vp1, vp2, vp1_prime,vp2_prime, vn1, vn2, PI = 3.14159265359f;;

float thetha;

int count=0,i,j, xmax=500, ymax=500, cycles=400, method=-1;

double wtime;
////////////////////////////////////////////////////////////////////////////

class point
{
    public:
    float x,y;
    float vx,vy;
    float mass;
    float angle;
    bool isColliding(point point_to_compare)
	{
        return sqrt(pow(x - point_to_compare.x,2) + pow(y - point_to_compare.y,2)) <= 2*sphereRadius;
	}
};

point positions[NBALLS];

void init() 
{
    for(i=0;i < NBALLS; i++ )
    {
        positions[i].x = (rand()%xmax*2-xmax);
        positions[i].y = (rand()%ymax-ymax/2);
        positions[i].vx = (rand()%20-10);
        positions[i].vy = (rand()%20-10);
        positions[i].mass = rand() % 10 + 1;
    }
}
void calculate_new_velocities()
{
    
    #pragma omp parallel for shared(positions,i) private (j, thetha, vp1, vp2, vn1, vn2, tmp,vp1_prime,vp2_prime ) ordered schedule(dynamic)
    for(i=0;i< NBALLS; i++ )
    {
        for(j=0;j < NBALLS; j++ )
        {      
            if ( (i!=j) && positions[i].isColliding(positions[j]))
            {
                
                thetha = asin(  (positions[j].y-positions[i].y)/ (sqrt(pow(positions[j].x-positions[i].x,2) + pow(positions[j].y-positions[i].y,2)) ) ) *180/PI;

                vp1 = positions[i].vx*cos(thetha) + positions[i].vy*sin(thetha);
                vn1 = -positions[i].vx*sin(thetha) + positions[i].vy*cos(thetha);

                vp2 = positions[j].vx*cos(thetha) +positions[j].vy*sin(thetha);
                vn2 = -positions[j].vx*sin(thetha) +positions[j].vy*cos(thetha);            

                tmp = 1.0/(positions[i].mass + positions[j].mass);

                vp1_prime = (positions[i].mass - e*positions[j].mass)*vp1*tmp + (1.0 + e)*positions[j].mass*vp2*tmp;
                vp2_prime = (1.0 + e)*positions[i].mass*vp1*tmp + (positions[j].mass - e*positions[i].mass)*vp2*tmp;
                
                positions[i].vx = vp1_prime*cos(thetha) - vn1*sin(thetha);
                positions[i].vy = vp1_prime*sin(thetha) + vn1*cos(thetha);

                positions[j].vx = vp2_prime*cos(thetha) - vn2*sin(thetha);
                positions[j].vy = vp2_prime*sin(thetha) + vn2*cos(thetha);

            }
        }
    }
}

void update_positions_Verlet()
{
    # pragma omp parallel shared ( positions, timeIncrement ) private ( i )
    // # pragma omp simd
    # pragma omp for
    for(i=0;i < NBALLS; i++ )
    {
        positions[i].x += timeIncrement*positions[i].vx;
        positions[i].y += timeIncrement*positions[i].vy;     
    }
}

void update_positions_EulerCromerMethod()
// Euler-Cromer method
// Vn+1 = Vn + g(tn,Xn)delta(t)
// Xn+1 = Xn + f(tn,Vn+1)delta(t)
{
    # pragma omp parallel shared ( positions, timeIncrement ) private ( i)
    # pragma omp for
    for(i=0;i < NBALLS; i++ )
    {
        positions[i].vx += timeIncrement*positions[i].x*timeIncrement;
        positions[i].vy += timeIncrement*positions[i].y*timeIncrement;

        positions[i].x += timeIncrement*positions[i].vx*timeIncrement;
        positions[i].y += timeIncrement*positions[i].vy*timeIncrement;     
    }
}

void update_positions_EulerMethod()
// Euler method
// Vn+1 = Vn + g(tn,Xn)delta(t)
// Xn+1 = Xn + f(tn,Vn)delta(t)
{
    float aux_x, aux_y;
    # pragma omp parallel shared ( positions, timeIncrement ) private ( i, aux_x, aux_y )
    # pragma omp for
    for(i=0;i < NBALLS; i++ )
    {
        aux_x=positions[i].x;
        aux_y=positions[i].y;

        positions[i].x += timeIncrement*positions[i].vx*timeIncrement;
        positions[i].y += timeIncrement*positions[i].vy*timeIncrement;     

        positions[i].vx += timeIncrement*aux_x*timeIncrement;
        positions[i].vy += timeIncrement*aux_y*timeIncrement;

    }
}

int main(int argc, char** argv)
{
    if(argc>1) method= atoi(argv[1]);
    srand(42);
    
    wtime = omp_get_wtime ( );    
    
    init();
    while(true)
    {
        calculate_new_velocities();
        if(method==1) update_positions_Verlet();
        else if (method==2) update_positions_EulerMethod();
        else if (method==3)  update_positions_EulerCromerMethod();
        else update_positions_Verlet();     
        count+=1;            
        if(cycles==count) break;
    }
    
    wtime = omp_get_wtime ( ) - wtime;

    cout << "Elapsed cpu time for main computation: " <<wtime << " seconds.\n";
    cout << positions[NBALLS-1].x << " "<<positions[NBALLS-1].y <<endl;
    return 0;
}
