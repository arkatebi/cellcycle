

void estimate_stable_expression(int NUM_NODES,int NUM_EDGES,\
                                int ITER_FOR_ODE,\
                                int TOTAL_EULER_STEPS,\
                                double EULER_STEP_SIZE,\
                                long double CONVERGENCE_PROXIMITY,\
                                double TRANS_RATE_FACTOR,\
                                const double *MPR_arrv,\
                                const double *DNR_arrv,\
                                const int *node_type_arrv,\
                    		        const int *edge_source_arrv,\
                    		        const int *edge_target_arrv,\
                    		        const int *edge_type_arrv,\
	                    	        const double *TSH_arrv,\
                    		        const int *HCO_arrv,\
                    		        const double *FCH_arrv,\
				                    const double *IC_arr,double *exp_arr,\
                                double *fX_arr);

void cluster_solutions(char *WORK_DIR, \
                       char *FNAME_STATES, char *FNAME_LIMITCYCLES, \
                       char *FNAME_SUMMARY,\
                       const int MODEL_NO,\
                       int NUM_NODES,int NUM_EDGES,\
                       int NUM_RANDOM_ICS,int ITER_FOR_ODE,\
                       double EULER_SIM_TIME,\
                       double EULER_STEP_SIZE,\
                       long double CONVERGENCE_PROXIMITY,\
                       double TRANS_RATE_FACTOR,\
                       const double *MPR_arrv,\
                       const double *DNR_arrv,\
                       const int *node_type_arrv,\
                       const int *edge_source_arrv,\
                       const int *edge_target_arrv,\
                       const int *edge_type_arrv,\
                       const double *TSH_arrv,\
                       const int *HCO_arrv,\
                       const double *FCH_arrv,\
                       double *fX_arr_norm,\
                       double EXP_arr[][NUM_NODES]);

int find_states(char *WORK_DIR, char *FNAME_STATES, const int MODEL_NO,\
                int NUM_NODES, int NUM_RANDOM_ICS,\
                long double CONVERGENCE_PROXIMITY,\
                double *fX_arr_norm,\
                double EXP_arr[][NUM_NODES],\
                bool fX_arr_bool[NUM_RANDOM_ICS]);

int find_limitcycles(char *WORK_DIR, char *FNAME_LIMITCYCLES, \
                     const int MODEL_NO,\
                     int NUM_NODES,int NUM_EDGES,\
                     int NUM_RANDOM_ICS,int ITER_FOR_ODE,\
                     double EULER_SIM_TIME,\
                     double EULER_STEP_SIZE,\
                     long double CONVERGENCE_PROXIMITY,\
                     double TRANS_RATE_FACTOR,\
                     const double *MPR_arrv,\
                     const double *DNR_arrv,\
                     const int *node_type_arrv,\
                     const int *edge_source_arrv,\
                     const int *edge_target_arrv,\
                     const int *edge_type_arrv,\
                     const double *TSH_arrv,\
                     const int *HCO_arrv,\
                     const double *FCH_arrv,\
                     double *fX_arr_norm,\
                     double EXP_arr[][NUM_NODES],\
                     bool fX_arr_bool[NUM_RANDOM_ICS]);

void write_limitcycle(FILE *fh_LCs, int MODEL_NO, int NUM_NODES,\
                      int count_LC, int SIZE_OF_PERIOD, int SIZE_OF_LIMIT_CYCLE,\
                      double **LC_exp_arr);

void cal_euler_approximation(int NUM_NODES,int NUM_EDGES,\
			                    int TOTAL_EULER_STEPS,double EULER_STEP_SIZE,\
                             double TRANS_RATE_FACTOR,\
                             const double *MPR_arrv,const double *DNR_arrv,\
                             const int *node_type_arrv,\
                             const int *edge_source_arrv,\
                             const int *edge_target_arrv,\
                             const int *edge_type_arrv,\
                             const double *TSH_arrv,\
                             const int *HCO_arrv,\
                             const double *FCH_arrv,\
                             double const *curr_exp_arrv,\
                             double *next_exp_arr,
                             double *fX_arr);

int detect_limitcycle(int NUM_NODES,int NUM_EDGES,\
                      int SIM_STEPS,double SIM_STEP_SIZE,\
                      long double CONVERGENCE_PROXIMITY,\
                      double TRANS_RATE_FACTOR,\
                      const double *MPR_arrv,const double *DNR_arrv,\
                      const int *node_type_arrv,\
                      const int *edge_source_arrv,\
                      const int *edge_target_arrv,\
                      const int *edge_type_arrv,\
                      const double *TSH_arrv,\
                      const int *HCO_arrv,\
                      const double *FCH_arrv,\
                      double *start_exp_arr,\
                      double *LC_start_exp_arr);

int cal_period(int NUM_NODES,\
               int MAX_PERIODS,\
               int count_min_exp,\
               int *min_idx_arr,\
               double min_exp_arr[][NUM_NODES],\
               double *LC_start_exp_arr);

double cal_limitcycle(int NUM_NODES,int NUM_EDGES,\
                      int SIM_STEPS,double SIM_STEP_SIZE,\
                      long double CONVERGENCE_PROXIMITY,\
                      double TRANS_RATE_FACTOR,\
                      const double *MPR_arrv,const double *DNR_arrv,\
                      const int *node_type_arrv,\
                      const int *edge_source_arrv,\
                      const int *edge_target_arrv,\
                      const int *edge_type_arrv,\
                      const double *TSH_arrv,\
                      const int *HCO_arrv,\
                      const double *FCH_arrv,\
                      double *start_exp_arr,\
                      double **EXP_arr);


void cal_fX(int NUM_NODES,int NUM_EDGES,\
            double EULER_STEP_SIZE,\
            double TRANS_RATE_FACTOR,\
            const double *MPR_arrv,
            const double *DNR_arrv,\
            const int *node_type_arrv,\
            const int *edge_source_arrv,\
            const int *edge_target_arrv,\
            const int *edge_type_arrv,\
            const double *TSH_arrv,\
            const int *HCO_arrv,\
            const double *FCH_arrv,\
            double const *curr_exp_arr,\
            double *fX_arr);

void cal_fX2(int NUM_NODES,int NUM_EDGES,\
            double EULER_STEP_SIZE,\
            double TRANS_RATE_FACTOR,\
            const double *MPR_arrv,
            const double *DNR_arrv,\
            const int *node_type_arrv,\
            const int *edge_source_arrv,\
            const int *edge_target_arrv,\
            const int *edge_type_arrv,\
            const double *TSH_arrv,\
            const int *HCO_arrv,\
            const double *FCH_arrv,\
            double const *curr_exp_arr,\
            double *fX_arr);

void cal_ode(int NUM_NODES,int NUM_EDGES,\
			    int TOTAL_EULER_STEPS,\
             double EULER_STEP_SIZE,\
             double TRANS_RATE_FACTOR,\
             const double *MPR_arrv,const double *DNR_arrv,\
             const int *node_type_arrv,\
             const int *edge_source_arrv,\
             const int *edge_target_arrv,\
             const int *edge_type_arrv,\
             const double *TSH_arrv,\
             const int *HCO_arrv,\
             const double *FCH_arrv,\
             double const *curr_exp_arrv);

void cal_EXPrange(int NUM_NODES,int NUM_EDGES,const double *MPR_arrv,\
                  const double *DNR_arrv,const int *edge_source_arrv,\
                  const int *edge_target_arrv,const int *edge_type_arrv,\
                  const double *FCH_arrv,double *minEXP_arr, double *maxEXP_arr);

void set_ICs(int NUM_NODES, const double *minEXP_arr, \
             const double *maxEXP_arr, double *IC_arr);

void set_ICs_by_std_st(int NUM_NODES, const double *std_st_exp_arr, \
                       double *IC_arr);


double eval_shiftedHill_fn (double x, double x0, double nx, double lamda);
double randu(double minV, double maxV);
double cal_norm (double *curr_exp, int NUM_NODES);
double sum_delta (double *curr_exp, double *next_exp, int NUM_NODES);
bool is_element(int key, const int *arr, int size); 
double **allocate_LCmemory(int LC_SIZE, int NUM_NODES);
void free_LCmemory(double ** EXP_arr, int LC_SIZE, int NUM_NODES);
void read_config(char *WORK_DIR, char* fname_cfg);
void print_config();

/********-------------- common stochastic simulation  -------- ***********/
/*------------------------------- START --------------------------------*/

void write_state_expressions(FILE *fh_states, int MODEL_NO, \
                             int NUM_NODES, double *exp_arr);

void cal_fX_gnw(int NUM_NODES,int NUM_EDGES,\
            double EULER_STEP_SIZE,\
            double TRANS_RATE_FACTOR,\
            const double *MPR_arrv,
            const double *DNR_arrv,\
            const int *node_type_arrv,\
            const int *edge_source_arrv,\
            const int *edge_target_arrv,\
            const int *edge_type_arrv,\
            const double *TSH_arrv,\
            const int *HCO_arrv,\
            const double *FCH_arrv,\
            double const *curr_exp_arr,\
            double *fX_arr,\
            double *fX2_arr);

/********-------------- common stochastic simulation  -------- ***********/
/*------------------------------- END  --------------------------------*/



/*--- probing models using stochastic simulation with GNW constant noise --*/
/*------------------------------- START -----------------------------------*/

void simulate_network_cnoise_gnw(char *WORK_DIR,\
                             char *FNAME_STATES, \
                             char *FNAME_LIMITCYCLES, \
                             char *FNAME_SUMMARY, \
                             const int MODEL_NO,\
                             int NUM_NODES,\
                             int NUM_EDGES,\
                             int NUM_RANDOM_ICS,\
                             int ITER_FOR_ODE,\
                             double EULER_SIM_TIME,\
                             double EULER_STEP_SIZE,\
                             long double CONVERGENCE_PROXIMITY,\
                             double TRANS_RATE_FACTOR,\
                             const double *MPR_arrv,\
                             const double *DNR_arrv,\
                             const int *node_type_arrv,\
                             const int *edge_source_arrv,\
                             const int *edge_target_arrv,\
                             const int *edge_type_arrv,\
                             const double *TSH_arrv,\
                             const int *HCO_arrv,\
                             const double *FCH_arrv,\
                             double *EXP_dict_arrv);

void estimate_stable_expression_cnoise_gnw(int NUM_NODES,int NUM_EDGES,\
                                           int ITER_FOR_ODE,int TOTAL_EULER_STEPS,\
                                           double EULER_STEP_SIZE,\
                                           long double CONVERGENCE_PROXIMITY,\
                                           double TRANS_RATE_FACTOR,\
                                           const double *MPR_arrv,\
                                           const double *DNR_arrv,\
                                           const int *node_type_arrv,
                                           const int *edge_source_arrv,\
                                           const int *edge_target_arrv,\
                                           const int *edge_type_arrv,\
                                           const double *TSH_arrv,\
                                           const int *HCO_arrv,\
                                           const double *FCH_arrv,\
                                           const double *IC_arr, double *exp_arr,\
                                           double *fX_arr,\
                                           double NOISE);

void cal_euler_approximation_cnoise_gnw(int NUM_NODES,int NUM_EDGES,\
                             int TOTAL_EULER_STEPS,double EULER_STEP_SIZE,\
                             double TRANS_RATE_FACTOR,\
                             const double *MPR_arrv,\
                            const double *DNR_arrv,\
                             const int *node_type_arrv,\
                             const int *edge_source_arrv,\
                             const int *edge_target_arrv,\
                             const int *edge_type_arrv,\
                             const double *TSH_arrv,\
                             const int *HCO_arrv,\
                             const double *FCH_arrv,\
                             double const *curr_exp_arrv,\
                             double *next_exp_arr,\
                            double *fX_arr,\
                            double NOISE);

/*--- probing models using stochastic simulation with GNW constant noise --*/
/*------------------------------- END  ----------------------------------*/



/********----- stochastic simulation with constant noise ---- ***********/
/*------------------------------- START --------------------------------*/
void find_solutions_stochastic(char *WORK_DIR, const int MODEL_NO,\
                               int NUM_NODES,int NUM_EDGES,\
                               int NUM_RANDOM_ICS,int ITER_FOR_ODE,\
                               double EULER_SIM_TIME,\
                               double EULER_STEP_SIZE,\
                               long double CONVERGENCE_PROXIMITY,\
                               double TRANS_RATE_FACTOR,\
                               const double *MPR_arrv,\
                               const double *DNR_arrv,\
                               const int *node_type_arrv,\
                               const int *edge_source_arrv,\
                               const int *edge_target_arrv,\
                               const int *edge_type_arrv,\
                               const double *TSH_arrv,\
                               const int *HCO_arrv,\
                               const double *FCH_arrv,\
                               double *EXP_dict_arrv,\
                               double *NOISE_strength,\
                               double *NOISE_strength_shot,\
                               double NOISE, double NOISE_SHOT);


void estimate_stable_expression_stochastic(int NUM_NODES,int NUM_EDGES,\
                                           int ITER_FOR_ODE,int TOTAL_EULER_STEPS,\
                                           double EULER_STEP_SIZE,\
                                           long double CONVERGENCE_PROXIMITY,\
                                           double TRANS_RATE_FACTOR,\
                                           const double *MPR_arrv,\
                                           const double *DNR_arrv,\
                                           const int *node_type_arrv,
                                           const int *edge_source_arrv,\
                                           const int *edge_target_arrv,\
                                           const int *edge_type_arrv,\
                                           const double *TSH_arrv,\
                                           const int *HCO_arrv,\
                                           const double *FCH_arrv,\
                                           const double *IC_arr, \
                                           double *exp_arr,\
                                           double *fX_arr,\
                                           double *NOISE_strength,\
                                           double *NOISE_strength_shot,\
                                           double NOISE,\
                                           double NOISE_SHOT);


void cal_euler_approximation_stochastic (int NUM_NODES,\
                                        int NUM_EDGES,\
                                        int TOTAL_EULER_STEPS,\
                                        double EULER_STEP_SIZE,\
                                        double TRANS_RATE_FACTOR,\
                                        
                                        const double *MPR_arrv,\
                                        const double *DNR_arrv,\
                                        const int *node_type_arrv,\
                                        const int *edge_source_arrv,\
                                        const int *edge_target_arrv,\
                                        const int *edge_type_arrv,\
                                        const double *TSH_arrv,\
                                        const int *HCO_arrv,\
                                        const double *FCH_arrv,\
                                        double const *curr_exp_arrv,\
                                        double *next_exp_arr,\
                                        double *fX_arr,\
                                        double *NOISE_strength,\
                                        double *NOISE_strength_shot,\
                                        double NOISE, double NOISE_SHOT );

/********----- stochastic simulation with constant noise ---- ***********/
/*------------------------------- END  --------------------------------*/




/********----- stochastic simulation with annealing ---- *****************/
/*************------------START-------------****************************/
void estimate_stable_expression_stochastic_anneal_network(\
                                           char *FNAME_STATES, \
                                           const int MODEL_NO,\
                                           int NUM_NODES,int NUM_EDGES,\
                                           int ITER_FOR_ODE,int TOTAL_EULER_STEPS,\
                                           double EULER_STEP_SIZE,\
                                           long double CONVERGENCE_PROXIMITY,\
                                           double TRANS_RATE_FACTOR,\
                                           const double *MPR_arrv,\
                                           const double *DNR_arrv,\
                                           const int *node_type_arrv,
                                           const int *edge_source_arrv,\
                                           const int *edge_target_arrv,\
                                           const int *edge_type_arrv,\
                                           const double *TSH_arrv,\
                                           const int *HCO_arrv,\
                                           const double *FCH_arrv,\
                                           const double *IC_arr, double *exp_arr,\
                                           double *fX_arr,\
                                           double *NOISE_strength,\
                                           double *NOISE_strength_shot,\
                                           double NOISE, double NOISE_SHOT,\
                                           double SCALING);




void find_solutions_stochastic_annealing(char *WORK_DIR, const int MODEL_NO,\
                                         int NUM_NODES,int NUM_EDGES,\
                                         int NUM_RANDOM_ICS,int ITER_FOR_ODE,\
                                         double EULER_SIM_TIME,\
                                         double EULER_STEP_SIZE,\
                                         long double CONVERGENCE_PROXIMITY,\
                                         double TRANS_RATE_FACTOR,\
                                         const double *MPR_arrv,\
                                         const double *DNR_arrv,\
                                         const int *node_type_arrv,\
                                         const int *edge_source_arrv,\
                                         const int *edge_target_arrv,\
                                         const int *edge_type_arrv,\
                                         const double *TSH_arrv,\
                                         const int *HCO_arrv,\
                                         const double *FCH_arrv,\
                                         double *EXP_dict_arrv,\
                                         double *NOISE_strength,\
                                         double *NOISE_strength_shot,\
                                         double NOISE, double NOISE_SHOT,\
                                         double SCALING);


void estimate_stable_expression_stochastic_annealing (const int MODEL_NO,\
                                                     int NUM_NODES,int NUM_EDGES,\
                                                     int ITER_FOR_ODE,\
                                                     int TOTAL_EULER_STEPS,\
                                                     double EULER_STEP_SIZE,\
                                                     long double CONVERGENCE_PROXIMITY,\
                                                     double TRANS_RATE_FACTOR,\
                                                     const double *MPR_arrv,\
                                                     const double *DNR_arrv,\
                                                     const int *node_type_arrv,
                                                     const int *edge_source_arrv,\
                                                     const int *edge_target_arrv,\
                                                     const int *edge_type_arrv,\
                                                     const double *TSH_arrv,\
                                                     const int *HCO_arrv,\
                                                     const double *FCH_arrv,\
                                                     const double *IC_arr, \
                                                     double *exp_arr,\
                                                     double *fX_arr,\
                                                     double *NOISE_strength,\
                                                     double *NOISE_strength_shot,\
                                                     double NOISE, double NOISE_SHOT,\
                                                     double SCALING);


void cal_rk_approximation(int NUM_NODES,int NUM_EDGES,\
                          double RK_TOLERANCE,\
                          double *ptr_TOTAL_TIME,\
                          double *ptr_RK_STEP_SIZE,\
                          double TRANS_RATE_FACTOR,\
                          const double *MPR_arrv,const double *DNR_arrv,\
                          const int *node_type_arrv,\
                          const int *edge_source_arrv,\
                          const int *edge_target_arrv,\
                          const int *edge_type_arrv,\
                          const double *TSH_arrv,\
                          const int *HCO_arrv,\
                          const double *FCH_arrv,\
                          double const *curr_exp_arrv,\
                          double *next_exp_arr,double *fX_arr);

/********----- stochastic simulation with annealing ---- *****************/
/*************------------END-------------*****************************/

/*** functions for stochastic annealing with GeneNetworkWeaver(GNW) ****/
/*************------------START-------------****************************/
void stochastic_anneal_network_gnw(char *WORK_DIR, 
                               char *FNAME_STATES, \
                               char *FNAME_LIMITCYCLES, \
                               char *FNAME_SUMMARY, \
                               const int MODEL_NO,\
                               int NUM_NODES,int NUM_EDGES,\
                               int NUM_RANDOM_ICS,int ITER_FOR_ODE,\
                               double EULER_SIM_TIME,\
                               double EULER_STEP_SIZE,\
                               long double CONVERGENCE_PROXIMITY,\
                               double TRANS_RATE_FACTOR,\
                               const double *steady_state_exp_arrv,\
                               const double *MPR_arrv,\
                               const double *DNR_arrv,\
                               const int *node_type_arrv,\
                               const int *edge_source_arrv,\
                               const int *edge_target_arrv,\
                               const int *edge_type_arrv,\
                               const double *TSH_arrv,\
                               const int *HCO_arrv,\
                               const double *FCH_arrv,\
                               double *EXP_dict_arrv,\
                               double NOISE, double SCALING);

void estimate_stable_expression_stochastic_anneal_network_gnw(\
                                           char *FNAME_STATES, \
                                           const int MODEL_NO,\
                                           int NUM_NODES,int NUM_EDGES,\
                                           int ITER_FOR_ODE,int TOTAL_EULER_STEPS,\
                                           double EULER_STEP_SIZE,\
                                           long double CONVERGENCE_PROXIMITY,\
                                           double TRANS_RATE_FACTOR,\
                                           const double *MPR_arrv,\
                                           const double *DNR_arrv,\
                                           const int *node_type_arrv,
                                           const int *edge_source_arrv,\
                                           const int *edge_target_arrv,\
                                           const int *edge_type_arrv,\
                                           const double *TSH_arrv,\
                                           const int *HCO_arrv,\
                                           const double *FCH_arrv,\
                                           const double *IC_arr, double *exp_arr,\
                                           double *fX_arr,\
                                           double NOISE, \
                                           double SCALING);

void cal_euler_approximation_stochastic_gnw (int NUM_NODES,int NUM_EDGES,\
                             int TOTAL_EULER_STEPS,double EULER_STEP_SIZE,\
                             double TRANS_RATE_FACTOR,\
                             const double *MPR_arrv,\
                            const double *DNR_arrv,\
                             const int *node_type_arrv,\
                             const int *edge_source_arrv,\
                             const int *edge_target_arrv,\
                             const int *edge_type_arrv,\
                             const double *TSH_arrv,\
                             const int *HCO_arrv,\
                             const double *FCH_arrv,\
                             double const *curr_exp_arrv,\
                             double *next_exp_arr,\
                            double *fX_arr,\
                            double NOISE);

/*****----- functions for stochastic annealing with gnew ----- ********/
/*************------------END -------------*****************************/

