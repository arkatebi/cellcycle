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
                      int SIMSTEPS_SIG_INCREASING,\
                      double SIM_STEP_SIZE,\
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
                      double *start_exp_arr,\
                      double **EXP_arr);

void cal_traj(int NUM_NODES,int NUM_EDGES,\
                    int SIM_STEPS,\
                    double SIM_STEP_SIZE,\
                    long double CONVERGENCE_PROXIMITY,\
                    double TRANS_RATE_FACTOR,\
                    double *MPR_arrv,\
                    const double *DNR_arrv,\
                    const int *node_type_arrv,\
                    const int *edge_source_arrv,\
                    const int *edge_target_arrv,\
                    const int *edge_type_arrv,\
                    const double *TSH_arrv,\
                    const int *HCO_arrv,\
                    const double *FCH_arrv,\
                    const int SIGNALING_NODE_ID,\
                    double *MPR_along_traj,\
                    double *start_exp_arr,\
                    double **TRAJ_exp_arr);

void cal_long_traj(FILE *fh_states,\
                   int MODEL_NO,\
                   int NUM_NODES,int NUM_EDGES,\
                    int SIM_STEPS,\
                    double SIM_STEP_SIZE,\
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
                    const double *start_exp_arr);


void cal_long_traj_with_updated_MPR(FILE *fh_states,\
                   int MODEL_NO,\
                   int NUM_NODES,int NUM_EDGES,\
                    int SIM_STEPS,\
                    double SIM_STEP_SIZE,\
                    long double CONVERGENCE_PROXIMITY,\
                    double TRANS_RATE_FACTOR,\
                    double *MPR_arrv,\
                    const double *DNR_arrv,\
                    const int *node_type_arrv,\
                    const int *edge_source_arrv,\
                    const int *edge_target_arrv,\
                    const int *edge_type_arrv,\
                    const double *TSH_arrv,\
                    const int *HCO_arrv,\
                    const double *FCH_arrv,\
                    const double *start_exp_arr);



void cal_traj_given_MPR_along(FILE *fh_states,\
                              int MODEL_NO,\
                              int NUM_NODES,int NUM_EDGES,\
                              int SIM_STEPS,\
                              double SIM_STEP_SIZE,\
                              long double CONVERGENCE_PROXIMITY,\
                              double TRANS_RATE_FACTOR,\
                              double *MPR_arrv,\
                              const double *DNR_arrv,\
                              const int *node_type_arrv,\
                              const int *edge_source_arrv,\
                              const int *edge_target_arrv,\
                              const int *edge_type_arrv,\
                              const double *TSH_arrv,\
                              const int *HCO_arrv,\
                              const double *FCH_arrv,\
                              const int SIGNALING_NODE_ID,\
                              double *MPR_along_traj,\
                              double *start_exp_arr,\
                              double *end_exp_arr);


void cal_trajectory(int NUM_NODES,int NUM_EDGES,\
                    int SIM_STEPS,\
                    double SIM_STEP_SIZE,\
                    long double CONVERGENCE_PROXIMITY,\
                    double TRANS_RATE_FACTOR,\
                    double *MPR_arrv,\
                    const double *DNR_arrv,\
                    const int *node_type_arrv,\
                    const int *edge_source_arrv,\
                    const int *edge_target_arrv,\
                    const int *edge_type_arrv,\
                    const double *TSH_arrv,\
                    const int *HCO_arrv,\
                    const double *FCH_arrv,\
                    const int SIGNALING_NODE_ID,\
                    double *MPR_along_traj,\
                    double *start_exp_arr,\
                    double **TRAJ_exp_arr);

void save_mpr_traj_sstate(FILE *fh_MPR, int MODEL_NO, double *MPR_arr, \
                    int SIGNALING_NODE);

//void save_mpr_ld_traj_sstate(FILE *fh_MPR, int MODEL_NO, \
//                             long double *MPR_arr, 
//                             int SIGNALING_NODE);

void save_mpr_at_traj_end(FILE *fh_MPR, int MODEL_NO, int SIGNALING_NODE,\
                          double *MPR_arr);

void save_expr_traj_sstate(FILE *fh_TRAJ, int MODEL_NO, int NUM_NODES,\
                          int TRAJ_SIZE, double **TRAJ_exp_arr);

//void save_expr_traj_spoints_sstate(FILE *fh_TRAJ, int MODEL_NO, int NUM_NODES,\
//                                   int TRAJ_SIZE, double **TRAJ_exp_arr);

void save_expr_traj_spoints_sstate(FILE *fh_TRAJ, int MODEL_NO, int NUM_NODES,\
                                   int TRAJ_SIZE, double **TRAJ_exp_arr, \
                                   int SAMPLING_DIST);

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

void set_ICs_by_sstate(int NUM_NODES, const double *std_st_exp_arr, \
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
