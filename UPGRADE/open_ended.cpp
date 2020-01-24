//
// Created by HYH on 2020/1/4.
//

//
// Created by HYH on 2019/12/5.
//

#pragma GCC optimize(2)
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
using namespace std;

/*Parameters*/
#define pop_robot 50
#define generation 10000

#define pi  3.14159265358979323846;  /* pi */
#define NUM_THREADS 6


/*  Simulation Parameters  */
double delta_t = 0.0001;          // Simulation time-step
double g[3] = {0, 0, -9.81};      // Gravity acceleration
double k_c = 1000000;              // Ground Constraint
double t = 0.1;                     // Global time variable
double u_s = 1;                   // Static friction coefficient
double u_k = 0.6;                 // Kinetic friction coefficient
double damping = 0.999;           // Velocity dampening
double init_high =0;             // initial height
double spring_f[3] = {0,0,0};
/*  Robot Parameters  */
int cycle = 5 *10000;  // Cycles for every generation
double w = 2 * pi;                // Frequency parameter w
double k_hard = 1000000;          // k value of hard spring fixing the shape
double bdy_l = 0.3;               // the length of robot
double bdy_w = 0.08;              // the with of robot
double bdy_h = 0.05;              // the hight of robot

/*  Evolutionary Parameters  */

double crossover_rate = 0.2;      // The probability of crossover for each individual in generation
double selection_rate = 0.5;      // Selection probability
double mutation_rate = 0.2;       // Mutation probability
double great_mutation_rate = 0.1; // Each generation, there is few individuals will be mutated significantly.

int protect_elite = 1;            // keep the best one and avoid mutation
double gen_pop[pop_robot][8];     // sorted by toe1,knee1,toe2,knee2...
double spr_change[pop_robot][8];  // sorted by toe1,knee1,toe2,knee2...
double fitness[pop_robot];        // Evolution fitness
double fitness_cut[10][pop_robot];        // Evolution fitness
double best_fitness;              // the best fitness of all

/*  Other Parameters  */
double knee_shift = 0;
double step_v =0.05;
double step_h = 1;
double step = 0;
int type =0;
double step_length = 0.4;
double pose[pop_robot];
int ni =-1;
int capture =0;
double data[generation][pop_robot];

struct Mass{
    double m;
    double p[3];
    double v[3];
    double a[3];
    double color[3];
};
struct Spring{
    double k;
    double lgh;
    double const_lgh;
    int cnn1;
    int cnn2;
};

template <class T>
int getArrayLen(T &array){
    return sizeof(array) / sizeof(array[0]);
}




double robot_masses[40][3] = {
        {0,-0.02,0.2},{0, 0, 0.2}, {0, bdy_l/3, 0.2}, {0, bdy_l*2/3, 0.2}, {0, bdy_l, 0.2},
        {0,-0.02,0.2+bdy_h},{0, 0, 0.2+bdy_h}, {0, bdy_l/3, 0.2+bdy_h}, {0, bdy_l*2/3, 0.2+bdy_h}, {0, bdy_l, 0.2+bdy_h},
        {bdy_w,-0.02,0.2}, {bdy_w, 0, 0.2}, {bdy_w, bdy_l/3, 0.2}, {bdy_w, bdy_l*2/3, 0.2}, {bdy_w, bdy_l, 0.2},
        {bdy_w,-0.02,0.2+bdy_h},{bdy_w, 0, 0.2+bdy_h}, {bdy_w, bdy_l/3, 0.2+bdy_h}, {bdy_w, bdy_l*2/3, 0.2+bdy_h}, {bdy_w, bdy_l, 0.2+bdy_h}
// 0-19 total 20 masses for body
};

double shape_gens[pop_robot][20][3];




int robot_hard_springs[96][2] = {
        {0,1},{1,2},{2,3},{3,4},{5,6},{6,7},{7,8},{8,9},{0,5},{1,6},{2,7},{3,8},{4,9},{1,5},{1,7},{2,6},{2,8},{3,7},{3,9},{4,8},
        {10,11},{11,12},{12,13},{13,14},{15,16},{16,17},{17,18},{18,19},{10,15},{11,16},{12,17},{13,18},{14,19},{11,15},{11,17},{12,16},{12,18},{13,17},{13,19},{14,18},
        {0,10},{1,11},{2,12},{3,13},{4,14},{5,15},{6,16},{7,17},{8,18},{9,19},{0,15},{10,5},{4,19},{9,14},{0,19},{9,10}
        //0-55 total 56 fixed body springs
        //48 ; maximum springs are 104 when there are 4 points in the region. (55~103)
};

int legs_gens[pop_robot][40][2];

int robot_soft_springs[8][2] = {
        {2, 20}, {5, 20},
        {12,22}, {15,22},
        {4, 21}, {8, 21},
        {14,23}, {18,23}
};


int mass_num = getArrayLen(robot_masses);
int hard_spring = getArrayLen(robot_hard_springs);
int ttl_spring = getArrayLen(robot_soft_springs) + hard_spring;

Mass masses_struct[pop_robot][40];
Spring springs_struct[pop_robot][104];




double linalg_norm(double a,double b,double c);
void BubbleSort(double  *p, int length, int * ind_diff);

void result();
double * f_spring_get(int nf, int i);
double * friction_get(double z_udr0, double a_x, double a_y, double a_z);



int init_robot(){
    // initialize the masses
    for(int num = 0; num < pop_robot; num++){

        //embed gens
        for (int trans = 0; trans<20; trans++){
            for (int po = 0; po < 3; po++) {
                robot_masses[20 + trans][po] = shape_gens[num][trans][po];
            }
        }
        for (int trans = 0; trans<40; trans++){
            for (int po = 0; po < 2; po++) {
                robot_hard_springs[56 + trans][po] = legs_gens[num][trans][po];
            }
        }


        for(int i = 0; i < mass_num; i++){
            masses_struct[num][i] = {1,robot_masses[i][0]+0.2* num,robot_masses[i][1],robot_masses[i][2]+init_high};
        }


        for(int i = 0; i < ttl_spring; i++){

            // initialize the hard springs
            if (i<hard_spring) {
                double l = linalg_norm((masses_struct[num][robot_hard_springs[i][0]].p[0] - masses_struct[num][robot_hard_springs[i][1]].p[0]),
                                       (masses_struct[num][robot_hard_springs[i][0]].p[1] - masses_struct[num][robot_hard_springs[i][1]].p[1]),
                                       (masses_struct[num][robot_hard_springs[i][0]].p[2] - masses_struct[num][robot_hard_springs[i][1]].p[2]));
                springs_struct[num][i] = {k_hard, l,l, robot_hard_springs[i][0], robot_hard_springs[i][1]};}

                // initialize the soft springs
            else {
                double l = linalg_norm((masses_struct[num][robot_soft_springs[i - hard_spring][0]].p[0] -
                                        masses_struct[num][robot_soft_springs[i - hard_spring][1]].p[0]),
                                       (masses_struct[num][robot_soft_springs[i - hard_spring][0]].p[1] -
                                        masses_struct[num][robot_soft_springs[i - hard_spring][1]].p[1]),
                                       (masses_struct[num][robot_soft_springs[i - hard_spring][0]].p[2] -
                                        masses_struct[num][robot_soft_springs[i - hard_spring][1]].p[2]));

                if (i == hard_spring or i == hard_spring + 6) {
                    springs_struct[num][i] = {gen_pop[num][0], l + gen_pop[num][2] * sin(gen_pop[num][4]),
                                              l,
                                              robot_soft_springs[i - hard_spring][0],
                                              robot_soft_springs[i - hard_spring][1]};
                }
                else if (i == hard_spring+1 or i == hard_spring + 7) {
                    springs_struct[num][i] = {gen_pop[num][1], l + gen_pop[num][3] * sin(gen_pop[num][5])+ knee_shift,
                                              l,
                                              robot_soft_springs[i - hard_spring][0],
                                              robot_soft_springs[i - hard_spring][1]};
                }
                else if (i == hard_spring+2 or i == hard_spring + 4){
                    springs_struct[num][i] = {gen_pop[num][0], l+gen_pop[num][2]* sin(gen_pop[num][6]),
                                              l,
                                              robot_soft_springs[i - hard_spring][0],
                                              robot_soft_springs[i - hard_spring][1]};
                }
                else if (i == hard_spring+3 or i == hard_spring + 5){
                    springs_struct[num][i] = {gen_pop[num][1], l + gen_pop[num][3] * sin(gen_pop[num][7])+ knee_shift,
                                              l,
                                              robot_soft_springs[i - hard_spring][0],
                                              robot_soft_springs[i - hard_spring][1]};
                }
            }

        }
    }
//    cout<<"Initialize Successfully"<<endl;
    return 0;
}
void create(int start, int end){
//    cout<<"Creating the New Individual from:"<<start<<"to"<<end<<endl;
    //  ||||||k_toe, k_knee, b_toe, b_knee, c1, c2, c3, c4 |||||| number of dots
    //    int toe = 55; //because leg length = 1.414 * 2
    //    int knee = 35;
    //    double knee_shift = -0.029;
    for (int i = start; i<end; i++) {
        double c[4];
        for (int j = 0; j < 4; j++) {
            c[j] = rand()%62800 *0.01;// generate a random value from 0-628

            gen_pop[i][4 + j] = c[j];

        }
        gen_pop[i][0] = (rand()%(1000000-100000)+100000)*0.1 ; //k_toe generate a random value from 10000-100000
        gen_pop[i][1] = (rand()%(1000000-100000)+100000)*0.1 ; //k_knee
        gen_pop[i][2] = (rand()%550)*0.0001;  //b_toe get a random number(0~0.055)
        gen_pop[i][3] = (rand()%350)*0.0001;  //b_knee get a random number(0~0.035)


//        cout<<"Legs assemble"<<endl;
        /*                   CREATE LEGS                    */
        //clear gens
        for (int k = 0; k<20; k++) {
            for (int f = 0; f<3; f++) {
                shape_gens[i][k][f] = 0;
            }
            for (int f=0;f<2;f++){
                legs_gens[i][k][f] = 0;
                legs_gens[i][20+k][f] = 0;
            }
        }
//        cout<<"1";
        //ground points
        double grd_f = (rand()%(4000-1500)+1500) * 0.0001; //0.15-0.4
        double grd_b = rand()%2500 * 0.0001 - 0.1;        //-0.1- 0.15

        for (int k = 0; k<4; k++) {
            shape_gens[i][k][0] = int(k/2) * bdy_w;
            shape_gens[i][k][1] = grd_f*(k%2) + grd_b*(1-k%2);
            shape_gens[i][k][2] = 0;
        }

//        cout<<"2"<<endl;
        int cre_num_dots_b = rand() % 3 + 2;  // get 2,3,4
//        cout<<cre_num_dots_b<<endl;
        int cre_num_dots_f = rand() % 3 + 2;  // get less or equal to back
//        cout<<cre_num_dots_f<<endl;

//        cout<<"Create Legs"<<endl;
//        cout<< "Back:"<< cre_num_dots_b<<';'<< "Front:"<<cre_num_dots_f<<endl;
        for (int k = 0; k<cre_num_dots_b; k++){
            double iy = rand()%2500 * 0.0001 - 0.1; // -1000~3000 -> -0.1~0.15
            double iz = rand()%4000 * 0.0001; // 0~4000 -> 0~0.4
//            cout<<"back_y"<<iy<<";back_z"<<iz<<endl;
            shape_gens[i][4+k][0] = 0;
            shape_gens[i][12+k][0] = bdy_w;
            shape_gens[i][4+k][1] = iy;
            shape_gens[i][12+k][1] = iy;
            shape_gens[i][4+k][2] = iz;
            shape_gens[i][12+k][2] = iz;
        }

        for (int k = 0; k<cre_num_dots_f; k++){
            double i2y = rand()%2500 * 0.0001 - 0.1; // 1500~7000 -> 0.15~0.4
            double i2z = rand()%2000 * 0.0001; // 0~2000 -> 0~0.2

            shape_gens[i][8+k][0] = 0;
            shape_gens[i][16+k][0] = bdy_w;
            shape_gens[i][8+k][1] = i2y;
            shape_gens[i][16+k][1] = i2y;
            shape_gens[i][8+k][2] = i2z;
            shape_gens[i][16+k][2] = i2z;
        }
        mass_num = getArrayLen(robot_masses);
        hard_spring = getArrayLen(robot_hard_springs);
        ttl_spring = getArrayLen(robot_soft_springs)+hard_spring;

        // connect leg masses
        //region 1
//        cout<<"Legs assemble"<<endl;
//        cout<<cre_num_dots_f<<endl;
//        cout<<cre_num_dots_b<<endl;

        for( int k = 0; k <( 2 * cre_num_dots_b-4); k++){
            if (k ==3){
                legs_gens[i][k][0] = 4 + k + 20;
                legs_gens[i][k][1] = 4 + 20;
                legs_gens[i][k+20][0] = 12 + k + 20;
                legs_gens[i][k+20][1] = 12 + 20;
            }
            else {
                legs_gens[i][k][0] = 4 + k + 20;
                legs_gens[i][k][1] = 4 + k +1 + 20;
                legs_gens[i][k+20][0] = 12 + k + 20;
                legs_gens[i][k+20][1] = 12 + k +1 + 20;
            }
        }


        //region 2

        for (int k = 0; k < (2 * cre_num_dots_f - 4); k++) {
            if (k == 3) {
                legs_gens[i][k + 10][0] = 8 + k + 20;
                legs_gens[i][k + 10][1] = 8 + 20;
                legs_gens[i][k + 30][0] = 16 + k + 20;
                legs_gens[i][k + 30][1] = 16 + 20;
            }
            else {
                legs_gens[i][k + 10][0] = 8 + k + 20;
                legs_gens[i][k + 10][1] = 8 + k + 1 + 20;
                legs_gens[i][k + 30][0] = 16 + k + 20;
                legs_gens[i][k + 30][1] = 16 + k + 1 + 20;
            }
        }

        //connect to body
//        cout<<"Legs to body"<<endl;
//        cout<<cre_num_dots_d<<endl;
        for (int k =0; k<4; k++){
            legs_gens[i][10 * k + 5][0] = 1 + k/2*10 + k%2*2;
            legs_gens[i][10 * k + 5][1] = 4 + 4 * k + 20;
            legs_gens[i][10 * k + 6][0] = 1 + k/2*10 + k%2*2;
            legs_gens[i][10 * k + 6][1] = 4 + 4 * k + 1 + 20;

            legs_gens[i][10 * k + 7][0] = k + 20;
            legs_gens[i][10 * k + 7][1] = 4 + 4 * k + k%2 * cre_num_dots_f + (1-k%2) * cre_num_dots_b - 2 + 20;
            legs_gens[i][10 * k + 8][0] = k + 20;
            legs_gens[i][10 * k + 8][1] = 4 + 4 * k + k%2 * cre_num_dots_f + (1-k%2) * cre_num_dots_b - 1 + 20;;
        }

//        cout<<"Add actuators"<<endl;
        // Add actuators
        //back
        for (int k =0; k<2; k++){
            robot_soft_springs[2*k][1] = 4 + k*8 +20;
//            robot_soft_springs[2*k+1][1] = 4 + cre_num_dots_b-1 + k*8 +20;
        }

        //front
        for (int k =0; k<2; k++){
            robot_soft_springs[2*k+4][1] = 8 + k*8 +20;
//            robot_soft_springs[2*k+5][1] = 8+ cre_num_dots_f-1 + k*8 +20;
        }

//        cout<<"Over assemble"<<endl;

    }

}
void moveit(){
#pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i<pop_robot; i++){

        for (int j = 0; j<8; j+=2) {
            springs_struct[i][j+hard_spring].k = gen_pop[i][0];
            springs_struct[i][j+1+hard_spring].k = gen_pop[i][1];
            double b_toe = gen_pop[i][2];
            double b_knee = gen_pop[i][3];
            double c1 = gen_pop[i][4];
            double c2 = gen_pop[i][5];
            double c3 = gen_pop[i][6];
            double c4 = gen_pop[i][7];
            if (j == 0 or j == 6){
                spr_change[i][j] = (b_toe * sin(w * t + c1));
                spr_change[i][j+1] = (b_knee * sin(w * t + c2) + knee_shift);

                springs_struct[i][j + hard_spring].lgh = springs_struct[i][j + hard_spring].const_lgh + spr_change[i][j];
                springs_struct[i][j + 1 + hard_spring].lgh = springs_struct[i][j + 1 + hard_spring].const_lgh + spr_change[i][j + 1];

            }
            else{
                spr_change[i][j] = (b_toe * sin(w * t + c3));
                spr_change[i][j+1] = (b_knee * sin(w * t + c4) + knee_shift);

                springs_struct[i][j + hard_spring].lgh = springs_struct[i][j + hard_spring].const_lgh + spr_change[i][j];
                springs_struct[i][j + 1 + hard_spring].lgh = springs_struct[i][j + 1 + hard_spring].const_lgh + spr_change[i][j + 1];

            }
        }

    }
}
void physic(){
#pragma omp parallel for num_threads(NUM_THREADS)
    for (int num = 0; num < pop_robot; num++) {
        //each masses will
        double f_spring[mass_num][3];
        double *f_spring_get_p;
        // get the distance change of every mass
        for (int i = 0; i < mass_num; i++) {

            f_spring_get_p = f_spring_get(num, i);
            for (int j = 0; j < 3; j++) {
                f_spring[i][j] = *(f_spring_get_p + j);
//                cout<<i<<'i'<<f_spring[i][j]<<endl;
            }

        }

        for (int i = 0; i < mass_num; i++) {

            if (masses_struct[num][i].p[2] >= 0) {
                double col[3] = {0.8, 0.8, 0.8};
                for (int j = 0; j < 3; j++) {
                    masses_struct[num][i].a[j] = g[j] + f_spring[i][j]/masses_struct[num][i].m;//
                    masses_struct[num][i].color[j] = *(col + j);
                }
            }
//            when the cube fall on the floor
            if (masses_struct[num][i].p[2] < 0) {
                if (masses_struct[num][i].a[0]< 0.00001)
                    masses_struct[num][i].a[0] = 0;
                if (masses_struct[num][i].a[1]< 0.00001)
                    masses_struct[num][i].a[1] = 0;
                double col[3] = {0.9, 0, 0};
                double f_f[3];
                double a_h = sqrt(masses_struct[num][i].a[0] * masses_struct[num][i].a[0] + masses_struct[num][i].a[1] * masses_struct[num][i].a[1]);

                if (a_h <= abs(masses_struct[num][i].a[2]) * u_s){
                    f_f[2] = -k_c * masses_struct[num][i].a[2];
                    masses_struct[num][i].v[0] = 0;
                    masses_struct[num][i].v[1] = 0;

                }


                for (int h = 0; h < 3; h++){
                    f_f[h] = *(friction_get(masses_struct[num][i].p[2],
                                            masses_struct[num][i].a[0],
                                            masses_struct[num][i].a[1],
                                            masses_struct[num][i].a[2])+h);}


                for (int j = 0; j < 3; j++) {
                    masses_struct[num][i].a[j] = g[j]+f_f[j]+ f_spring[i][j]/masses_struct[num][i].m;//

                    masses_struct[num][i].color[j] = *(col + j);}
            }


            double dir[8][3];
            for (int ad = 0; ad < 8; ad++){
                for (int j = 0; j<3;j++){
                    dir[ad][j] = masses_struct[num][springs_struct[num][ad+hard_spring].cnn1].p[j] - masses_struct[num][springs_struct[num][ad+hard_spring].cnn2].p[j];
                }
                double length_c = linalg_norm(dir[ad][0],dir[ad][1],dir[ad][2]);
                dir[ad][0] /= length_c;
                dir[ad][1] /= length_c;
                dir[ad][2] /= length_c;
//                for (int j = 1; j<3; j++) {
//                    masses_struct[num][springs_struct[num][ad+24].cnn1].p[j] += dir[ad][j] * spr_change[num][ad]/2;
//                    masses_struct[num][springs_struct[num][ad+24].cnn2].p[j] -= dir[ad][j] * spr_change[num][ad]/2;
//                }
            }
            //move the cube based on the new acceleration.
            for (int j = 1; j < 3; j++) {
                masses_struct[num][i].v[j] = (masses_struct[num][i].v[j] + masses_struct[num][i].a[j] * delta_t) * damping;
                masses_struct[num][i].p[j] = masses_struct[num][i].p[j] + masses_struct[num][i].v[j] * delta_t;
            }


//                lock(num);



        }

        pose[num] = - abs(  (masses_struct[num][2].p[2]+masses_struct[num][6].p[2])/2 - 0.2)*0.0001;
    }

}
void get_fitness(){

    for (int i= 0; i<pop_robot; i++){

        double move_dis = masses_struct[i][0].p[1];
        fitness[i] = move_dis;

        //PUNISH
        if (masses_struct[i][5].p[2]<masses_struct[i][0].p[2])   //if the robot is turn over --> die
            fitness[i] = 0;
        if (move_dis < 0)   //if the robot retrogress --> die
            fitness[i] = 0;
        if(masses_struct[i][0].p[2]>0.3) //if the robot jump too high --> die
            fitness[i] = 0;
        if(masses_struct[i][4].p[2]<0.1) //if the robot lie down --> die
            fitness[i] = 0;
        if(masses_struct[i][4].p[2]>0.3) //if the robot jump too high --> die
            fitness[i] = 0;

        //during the time:




    }
}

void crossover(int choose){

    // parents are chosen from good individuals.
    int n = pop_robot*selection_rate;
    int parent1 = 0;
    int parent2 = 0;
    if (rand()%2 == 0) {
        parent1 = rand() % (n);
        parent2 = rand() % (pop_robot);
    }
    else{
        parent2 = rand() % (n);
        parent1 = rand() % (pop_robot);
    }
//     cout<<"MATE:"<<parent1<<"and"<<parent2<<endl;
    // Here is the actuator gens inherit.

    if (rand()%2 == 0) {
        for (int i = 0; i < 8; i++) {
            gen_pop[choose][i] = gen_pop[parent1][i];
//            cout<<"actuator gens inherit: "<<parent1<<" to "<<choose<<endl;
        }
    }
    else{
        for (int i = 0; i<8; i++){
            gen_pop[choose][i] = gen_pop[parent2][i];
//            cout<<"actuator gens inherit: "<<parent2<<" to "<<choose<<endl;
        }
    }


//  the hint leg comes from parent1, the front leg comes from parent2.
    for(int i = 4; i < 8; i++) {
        for (int j = 0; j<3; j++) {
//            hint leg from one parent
            shape_gens[choose][i][j] = shape_gens[parent1][i][j];
            shape_gens[choose][i+8][j] = shape_gens[parent1][i+8][j];
//            front leg from the other parent
            shape_gens[choose][i+4][j] = shape_gens[parent2][i+4][j];
            shape_gens[choose][i+12][j] = shape_gens[parent2][i+12][j];
        }
    }

    for (int j = 0; j<3; j++) {
//        hint leg ground points
        shape_gens[choose][0][j] = shape_gens[parent1][0][j];
        shape_gens[choose][2][j] = shape_gens[parent1][2][j];
//        front leg
        shape_gens[choose][1][j] = shape_gens[parent2][1][j];
        shape_gens[choose][3][j] = shape_gens[parent2][3][j];
    }


    //leg_genres: the connection of leg masses
    for (int j = 0; j < 10; j++){
        for (int k = 0; k < 2; k++){
            legs_gens[choose][j+20+0][k] = legs_gens[parent1][j+20+0][k];
            legs_gens[choose][j+20+10][k] = legs_gens[parent2][j+20+10][k];
        }
        for (int k = 0; k < 2; k++){
            legs_gens[choose][j+0][k] = legs_gens[parent1][j+0][k];
            legs_gens[choose][j+10][k] = legs_gens[parent2][j+10][k];
        }
    }

//    cout<<"the hint leg comes from "<<parent1<<", the front leg comes from "<<parent2<<endl;
}

void mutation(int choose){

    if (rand()%2 == 0){
//        cout<< "muscle actuator will be changed."<<endl;
        for (int i = 0; i<8;i++){
            if (i ==2 or i ==3) {
                if (rand()% 3 == 0) {
                    gen_pop[choose][i] = gen_pop[choose][i] + (rand()%10)*0.0001;
//                    cout<<"1+++"<<endl;
                }
                if (rand() % 3 == 1) {
                    //            cout << "mutation_gen:" << i <<"divide"<< endl;
                    gen_pop[choose][i] = gen_pop[choose][i] - (rand()%10)*0.0001;
//                    cout<<"1----"<<endl;
                }
            }
            else{
                if (rand()% 3 == 0) {
                    gen_pop[choose][i] = gen_pop[choose][i]+(rand()%10);
                }
                if (rand()%3 == 1){

                    gen_pop[choose][i] = gen_pop[choose][i]-(rand()%10);
                }
            }
        }
    }

    else{
//        cout<< "the shape of leg will be changed a little bit randomly."<<endl;
        for (int i = 4; i<11;i++)
            if (shape_gens[choose][i][1] > -0.02 and shape_gens[choose][i][2] > 0){
                double direction1 = ((rand()%11)-5 ) * 0.005;
                double direction2 = ((rand()%11)-5 ) * 0.005;
                if(shape_gens[choose][i][1] <0.4 and shape_gens[choose][i][2] <0.4) {
                    shape_gens[choose][i][1] += direction1;
                    shape_gens[choose][i + 8][1] += direction1;
                    shape_gens[choose][i][2] += direction2;
                    shape_gens[choose][i + 8][2] += direction2;
                }
            }

    }



}

void evolution(){
    int indx[pop_robot]={0};
    double selected[pop_robot][8]={0};
    double selected_shape[pop_robot][20][3];
    int selected_legs[pop_robot][40][2];
//    selection
    BubbleSort(fitness, pop_robot, indx);
    best_fitness = fitness[0];
    cout<<fitness[0]<<endl;
    for (int i = 0; i < pop_robot; i++){
        for (int j = 0; j<8; j++)
            selected[i][j] = gen_pop[indx[i]][j];
        for (int j = 0; j < 20; j++){
            for (int k = 0; k < 3; k++){
                selected_shape[i][j][k] = shape_gens[indx[i]][j][k];
            }
        }
        for (int j = 0; j < 40; j++){
            for (int k = 0; k < 2; k++){
                selected_legs[i][j][k] = legs_gens[indx[i]][j][k];
            }
        }
    }

    for (int i = 0; i < pop_robot; i++) {
        for (int j = 0; j < 8; j++)
            gen_pop[i][j] = selected[i][j];
        for (int j = 0; j < 20; j++){
            for (int k = 0; k < 3; k++){
                shape_gens[i][j][k] = selected_shape[i][j][k];
            }
        }
        for (int j = 0; j < 40; j++){
            for (int k = 0; k < 2; k++){
                legs_gens[i][j][k] = selected_legs[i][j][k];
            }
        }
    }

    // the bad individuals will be recreated, which means great mutation.
    create(int((1-great_mutation_rate) * pop_robot), pop_robot);


//    cout<<"crossover"<<endl;
    for (int i = selection_rate * pop_robot; i < pop_robot * (1 - great_mutation_rate); i++) {
//            cout<<"Number "<<i<<" becomes child"<<endl;
        crossover(i);
    }

//    cout<<"mutation"<<endl;
//    mutation
    for (int i = protect_elite; i < pop_robot*(1-great_mutation_rate); i++) {
        if  (rand()%10 <=  mutation_rate*10){
//            cout<<"Number "<<i<<"robot is ready for mutation!"<<endl;
            mutation(i);
        }
    }
}


int main(){

    srand(int(time(0)));
    float threshold = 0.4;
    create(0, pop_robot);
    for (int i= 0; i< generation; i++){

        int catcher = 0;
        init_robot();
        t=0;
        // loop five cycle to get fitness
//        cout<<"GENERATION:"<<i<<endl;
        for(int j = 0; j< cycle; j++){

            physic();
            moveit();
            t+=delta_t;

            if (j%7000 ==0){
                get_fitness();
                for (int indi=0; indi<pop_robot; indi ++)
                    fitness_cut[catcher][indi] = fitness[indi];
                catcher+=1;
//                cout<<"t:"<<t<<"; catcher"<<endl;
            }
        }
        get_fitness();
        for (auto & j : fitness_cut) {
            if (j[i] == 0){
                fitness[i] = 0;
            }
        }
        evolution();
        for (int k = 0; k< pop_robot; k++){
            data[i][k] = fitness[k];
        }
        if (data[i][0] > threshold){
            threshold += 0.1;
            cout<<"GENERATION "<<i<<endl;
            result();
        }

    }
    result();

    std::ofstream outfile;
    outfile.open("C:\\Users\\yh3187\\Desktop\\EA_Project\\all_individuals.txt", std::ios_base::app);
    for (int tg = 0; tg< generation; tg++){
        for (int al = 0; al< pop_robot; al++) {
            outfile <<data[tg][al]<<endl;
        }
    }
    return 0;
}


double linalg_norm(double a,double b,double c){
    double ln = sqrt(a*a+b*b+c*c);
    return ln;
}
double *friction_get(double z_udr0, double a_x, double a_y, double a_z){
    static double f_f[3];
    if (a_z<0){
        f_f[0] = -a_x - a_z * u_k;
        f_f[1] = -a_y - a_z * u_k;
    }
    else{
        f_f[0] = -a_x + a_z * u_k;
        f_f[1] = -a_y + a_z * u_k;}

    f_f[2] = -k_c * z_udr0;
    return f_f;
}
double * f_spring_get(int nf, int i)
{
    spring_f[0] = 0;
    spring_f[1] = 0;
    spring_f[2] = 0;
    for(int j = 0; j <ttl_spring; j ++)
    {
        if (springs_struct[nf][j].cnn1 == i)
        {
            double x1 = masses_struct[nf][springs_struct[nf][j].cnn1].p[0] - masses_struct[nf][springs_struct[nf][j].cnn2].p[0];
            double y1 = masses_struct[nf][springs_struct[nf][j].cnn1].p[1] - masses_struct[nf][springs_struct[nf][j].cnn2].p[1];
            double z1 = masses_struct[nf][springs_struct[nf][j].cnn1].p[2] - masses_struct[nf][springs_struct[nf][j].cnn2].p[2];
            double l_c = linalg_norm(x1, y1, z1);
            double f0 = springs_struct[nf][j].k * (l_c - springs_struct[nf][j].lgh);
            if (l_c!=0){
                spring_f[0] -= (x1 / l_c) *f0;
                spring_f[1] -= (y1 / l_c) *f0;
                spring_f[2] -= (z1 / l_c )*f0;
            }
        }

        if (springs_struct[nf][j].cnn2 == i)
        {
            double x1 = masses_struct[nf][springs_struct[nf][j].cnn1].p[0] - masses_struct[nf][springs_struct[nf][j].cnn2].p[0];
            double y1 = masses_struct[nf][springs_struct[nf][j].cnn1].p[1] - masses_struct[nf][springs_struct[nf][j].cnn2].p[1];
            double z1 = masses_struct[nf][springs_struct[nf][j].cnn1].p[2] - masses_struct[nf][springs_struct[nf][j].cnn2].p[2];
            double l_c = linalg_norm(x1, y1, z1);
            double f0 = springs_struct[nf][j].k * (l_c - springs_struct[nf][j].lgh);
            if (l_c!=0){
                spring_f[0] += (x1 / l_c) *f0;
                spring_f[1] += (y1 / l_c) *f0;
                spring_f[2] += (z1 / l_c )*f0;
            }
        }
    }
    return spring_f;
}
void BubbleSort(double  *p, int length, int * ind_diff) {
    for (int m = 0; m < length; m++) {
        ind_diff[m] = m;
    }

    for (int i = 0; i < length; i++) {
        for (int j = 0; j < length - i - 1; j++) {
            if (p[j] < p[j + 1]) {
                double temp = p[j];
                p[j] = p[j + 1];
                p[j + 1] = temp;

                int ind_temp = ind_diff[j];
                ind_diff[j] = ind_diff[j + 1];
                ind_diff[j + 1] = ind_temp;
            }
        }
    }
}
void result(){

    cout<<"The Result Page"<<endl;
    cout<<"Best Fitness:"<<fitness[0]<<endl;
    cout<<"Control Genres:"<<endl;
    for(int i = 0; i <8; i ++){
        cout<<gen_pop[0][i]<<",";
    }
    cout<<endl;
    cout<<"Shape Genres:"<<endl;
    for(int i = 0; i <20; i ++){
        cout<<"{"<<shape_gens[0][i][0]<<","<<shape_gens[0][i][1]<<","<<shape_gens[0][i][2]<<"},";
        if ((i+1)%5==0){
            cout<<endl;
        }
    }
    cout<<endl;

    //k_toe, k_knee, b_toe, b_knee, c1, c2, c3, c4
}


