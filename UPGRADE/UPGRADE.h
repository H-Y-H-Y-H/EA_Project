//
// Created by HYH on 2020/1/4.
//

#ifndef OPENGL_UPGRADE_ROBOT_H
#define OPENGL_UPGRADE_ROBOT_H

#define mode 1
#define pi  3.14159265358979323846;  /* pi */
#define pop_robot 13              // the number of robot
#define generation  5000          // the number of generation

/*  Simulation Parameters  */
double delta_t = 0.0001;          // Simulation time-step
double g[3] = {0, 0, -9.81};      // Gravity acceleration
double k_c = 1000000;              // Ground Constraint
double t = 0;                     // Global time variable
double u_s = 1;                   // Static friction coefficient
double u_k = 0.6;                 // Kinetic friction coefficient
double damping = 0.999;           // Velocity dampening
double init_high =0;             // initial height
int ti = 0;
/*  Robot Parameters  */
#define cycle 5*10000             // Cycles for every generation
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
double fitness_cut[10][pop_robot];        // Evolution fitness during the period
double best_fitness;              // the best fitness of all

/*  Other Parameters  */
double knee_shift = 0;
double step_v =0.05;
double step_h = 0.9;
double step = 0;
int type =0;
double step_length = 0.5;
double pose[pop_robot];
int ni =-1;
int capture = 0;
float ice_start = 2;
float ice_over = 3;
double data[generation][pop_robot];

//data collection

double fric_data_front[pop_robot][cycle];
double fric_data_hint[pop_robot][cycle];


double fric_data_front2[pop_robot][cycle];
double fric_data_hint2[pop_robot][cycle];

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
        {bdy_w,-0.02,0.2+bdy_h},{bdy_w, 0, 0.2+bdy_h}, {bdy_w, bdy_l/3, 0.2+bdy_h}, {bdy_w, bdy_l*2/3, 0.2+bdy_h}, {bdy_w, bdy_l, 0.2+bdy_h},
// 0-19 total 20 masses for body
};

double shape_gens[pop_robot][20][3];
// 0~3 for the ground points
// 4~7 for the first region
// 8~11 for the second region
// 12~15 for the 2first region
// 16~19 for the 2second region


int robot_hard_springs[96][2] = {
        {0,1},{1,2},{2,3},{3,4},{5,6},{6,7},{7,8},{8,9},{0,5},{1,6},{2,7},{3,8},{4,9},{1,5},{1,7},{2,6},{2,8},{3,7},{3,9},{4,8},
        {10,11},{11,12},{12,13},{13,14},{15,16},{16,17},{17,18},{18,19},{10,15},{11,16},{12,17},{13,18},{14,19},{11,15},{11,17},{12,16},{12,18},{13,17},{13,19},{14,18},
        {0,10},{1,11},{2,12},{3,13},{4,14},{5,15},{6,16},{7,17},{8,18},{9,19},{0,15},{10,5},{4,19},{9,14},{0,19},{9,10}
        //0-55 total 56 fixed body springs
        //48 ; maximum springs are 104 when there are 4 points in the region. (56~104)
};

//how leg masses collected
int legs_gens[pop_robot][40][2];// 5 + 4
// 0~9 first region
// 10~19 second region
// 20~29 2first region
// 30~39 2second region

//initialize the soft_springs
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



#endif //OPENGL_UPGRADE_ROBOT_H


