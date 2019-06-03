#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include <math.h>       /* asin */
#include <iostream>
#include <fstream>
#include <ctime>


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
bool perfTiming = true;	// performance timing?


int count=0,i,j, xmax=500, ymax=500, cycles=400;

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
    // #pragma omp parallel for
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
    // cout<<"IN"<<endl;
    // std::ofstream outf("collisions_serial.dat", std::ios::app);
    // if (!outf)
    // {
    //     std::cerr << "Uh oh, collisions_serial.dat could not be opened for writing!" << std::endl;
    //     exit(1);
    // }

    // #pragma omp parallel for collapse(2)
    for(i=0;i< NBALLS; i++ )
    {
        for(j=0;j < NBALLS; j++ )
        {       
            if ( i!=j )
            {
                if (positions[i].isColliding(positions[j]))
                {
                    // printf("x1 %f, y1 %f \n x2 %f, y2 %f\n",positions[i].x, positions[i].y, positions[j].x, positions[j].y);
                    thetha = asin(  (positions[j].y-positions[i].y)/ (sqrt(pow(positions[j].x-positions[i].x,2) + pow(positions[j].y-positions[i].y,2)) ) ) *180/PI;
                    // printf("thethaa %f\n", thetha);
                    vp1 = positions[i].vx*cos(thetha) + positions[i].vy*sin(thetha);
                    vn1 = -positions[i].vx*sin(thetha) + positions[i].vy*cos(thetha);

                    vp2 = positions[j].vx*cos(thetha) +positions[j].vy*sin(thetha);
                    vn2 = -positions[j].vx*sin(thetha) +positions[j].vy*cos(thetha);            
                    // printf("\nCollision occurred\n\n");

                    tmp = 1.0/(positions[i].mass + positions[j].mass);

                    vp1_prime = (positions[i].mass - e*positions[j].mass)*vp1*tmp + (1.0 + e)*positions[j].mass*vp2*tmp;
                    vp2_prime = (1.0 + e)*positions[i].mass*vp1*tmp + (positions[j].mass - e*positions[i].mass)*vp2*tmp;
                    
                    positions[i].vx = vp1_prime*cos(thetha) - vn1*sin(thetha);
                    positions[i].vy = vp1_prime*sin(thetha) + vn1*cos(thetha);

                    positions[j].vx = vp2_prime*cos(thetha) - vn2*sin(thetha);
                    positions[j].vy = vp2_prime*sin(thetha) + vn2*cos(thetha);

                    // printf("vp1 %f, vp2 %f\n vn1 %f, vn2 %f\n v1prim %f, v2prim%f\n positions[i].vx %f, positions[i].vy %f\n positions[j].vx %f,positions[j].vy %f\n", vp1,vp2, vn1,vn2, 
                    // vp1_prime, vp2_prime,  positions[i].vx, positions[i].vy, positions[j].vx,positions[j].vy);
                } 
            }
        }
    }
    // cout<<"out"<<endl;
}

void update_positions()
{
    for(i=0;i < NBALLS; i++ )
    {
        positions[i].x += timeIncrement*positions[i].vx;
        positions[i].y += timeIncrement*positions[i].vy;     
    }
}

int main(int argc, char** argv)
{
    clock_t begin = clock();
    srand(42);
    
    init();
    while(true)
    {
        calculate_new_velocities();
        update_positions();
        count+=1;            
        if(cycles==count) break;
    }
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    cout << "  Elapsed cpu time for main computation: "<< elapsed_secs << " seconds.\n";
    return 0;
}
