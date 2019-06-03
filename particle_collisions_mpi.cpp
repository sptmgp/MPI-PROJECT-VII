#include <stdio.h>
#include <stdlib.h>
#include <math.h>       /* asin */
#include <iostream>
#include <mpi.h>
// #include <omp.h>
#include <chrono> // for std::chrono functions


#define NBALLS 1000

using namespace std;
////////////////////////////////////////////////////////////////////////////

float  e = 0.9; //  coefficient of restitution.
float  sphereRadius = 5.0;
float  timeIncrement = 0.07;

float PI = 3.14159265359f;
float vp1, vp2, vp1_prime,vp2_prime, vn1, vn2, tmp;
float thetha;

int count=0, xmax=500, ymax=500, cycles=400;
int i,j,k,from,to, id, nproc;
bool loop=true;
////////////////////////////////////////////////////////////////////////////
class Timer
{
    private:
	// Type aliases to make accessing nested type easier
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1> >;
 
	std::chrono::time_point<clock_t> m_beg;
 
    public:
	Timer() : m_beg(clock_t::now())
	{
	}

	void reset()
	{
		m_beg = clock_t::now();
	}
 
	double elapsed() const
	{
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	}
};
struct point
{
    public:
    float x,y;      //positions
    float vx,vy;    //velocities
    float mass;     //mass
    bool isColliding(point point_to_compare)
	{
        return sqrt(pow(x - point_to_compare.x,2) + pow(y - point_to_compare.y,2)) <= 2*sphereRadius;
	}
};

point *positions = new point[NBALLS];

void init() 
//Initializes the particles mass, velocities, positions
{
    for(i=0;i<NBALLS;i++)
    {
        positions[i].x = (rand()%xmax*2-xmax);
        positions[i].y = (rand()%ymax-ymax/2);
        positions[i].vx = (rand()%20-10);
        positions[i].vy = (rand()%20-10);
        positions[i].mass = rand() % 10 + 1;
    }
}
point* calculate_new_velocities(MPI_Datatype point_type)
{
    point *aux = new point[NBALLS];

    MPI_Scatter(positions, NBALLS/nproc, point_type, positions, NBALLS/nproc, point_type, 0, MPI_COMM_WORLD);
    for(i=from;i< to; i++ )
    {
        aux[i].x = positions[i].x;
        aux[i].y = positions[i].y;        
        aux[i].vx =positions[i].vx;
        aux[i].vy =positions[i].vy;     
        aux[i].mass =positions[i].mass;     
    }
    MPI_Gather(&aux[from], NBALLS/nproc, point_type, aux, NBALLS/nproc, point_type, 0, MPI_COMM_WORLD);

    MPI_Scatter(positions, NBALLS/nproc, point_type, positions, NBALLS/nproc, point_type, 0, MPI_COMM_WORLD);
    for(i=from;i< to; i++ )
    {
        for(j=from;j < to; j++ )
        {       
            if ( (i!=j) && positions[i].isColliding(positions[j]) )
            {        
                thetha = asin(  (positions[j].y-positions[i].y)/ (sqrt(pow(positions[j].x-positions[i].x,2) + pow(positions[j].y-positions[i].y,2)) ) ) *180/PI;
        
                vp1 = positions[i].vx*cos(thetha) + positions[i].vy*sin(thetha);
                vn1 = -positions[i].vx*sin(thetha) + positions[i].vy*cos(thetha);

                vp2 = positions[j].vx*cos(thetha) +positions[j].vy*sin(thetha);
                vn2 = -positions[j].vx*sin(thetha) +positions[j].vy*cos(thetha);            
             
                tmp = 1.0/(positions[i].mass + positions[j].mass);

                vp1_prime = (positions[i].mass - e*positions[j].mass)*vp1*tmp + (1.0 + e)*positions[j].mass*vp2*tmp;
                vp2_prime = (1.0 + e)*positions[i].mass*vp1*tmp + (positions[j].mass - e*positions[i].mass)*vp2*tmp;
                
                aux[i].vx = vp1_prime*cos(thetha) - vn1*sin(thetha);
                aux[i].vy = vp1_prime*sin(thetha) + vn1*cos(thetha);

                aux[j].vx = vp2_prime*cos(thetha) - vn2*sin(thetha);
                aux[j].vy = vp2_prime*sin(thetha) + vn2*cos(thetha);
            }
        }
    }

    MPI_Gather(&aux[from], NBALLS/nproc, point_type, aux, NBALLS/nproc, point_type, 0, MPI_COMM_WORLD);
    return aux;
    delete aux;
}

point* update_positions(MPI_Datatype point_type)
{
    point *aux = new point[NBALLS];
    
    MPI_Bcast(&timeIncrement, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(positions, NBALLS/nproc, point_type, positions, NBALLS/nproc, point_type, 0, MPI_COMM_WORLD);
    for(i=from;i< to; i++ )
    {
        aux[i].x = positions[i].x + timeIncrement*positions[i].vx;
        aux[i].y = positions[i].y + timeIncrement*positions[i].vy;     
        
        aux[i].vx =positions[i].vx;
        aux[i].vy =positions[i].vy;     
        aux[i].mass =positions[i].mass;     
    }
    MPI_Gather(&aux[from], NBALLS/nproc, point_type, aux, NBALLS/nproc, point_type, 0, MPI_COMM_WORLD);

    return aux;
    delete aux;
}

int main(int argc, char** argv)
{
    srand(42);
    Timer t;    
    init();

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc); //get number of total nodes
    MPI_Comm_rank(MPI_COMM_WORLD, &id);    // get id of my node

    /* create a type for struct point */
    const int nitems=5;
    int          blocklengths[5] = {1,1,1,1,1};
    MPI_Datatype types[5] = {MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT};
    MPI_Datatype point_type;
    MPI_Aint     offsets[5];

    offsets[0] = offsetof(point, x);
    offsets[1] = offsetof(point, y);
    offsets[2] = offsetof(point, vy);
    offsets[3] = offsetof(point, vx);
    offsets[4] = offsetof(point, mass);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &point_type);
    MPI_Type_commit(&point_type);
    
    from = id * NBALLS/nproc;
    to = (id+1) * NBALLS/nproc;  

    while (loop)
    {
        positions = calculate_new_velocities(point_type);
        positions = update_positions(point_type);
        count+=1;
        if(count ==cycles){ loop=false;}
    }    

    MPI_Finalize();
    if(id==0) cout << "Elapsed cpu time for main computation: " << t.elapsed() << " seconds.\n";
    return 0;
}
