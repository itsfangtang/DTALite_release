// TAPLite.cpp : This file contains the 'main' function. Program execution begins and ends here.

// This code is built based on the implementation available at:
// http://www.bgu.ac.il/~bargera/tntp/FW.zip
// tntp from Dr. Hillel Bar-Gera. It includes one of the most efficient Deque implementations of 
// the modified label correcting (MLC) algorithm in C. For more details, see mincostroutes.cpp. 
// The enhanced C++ counterpart has served as the path engine for Path4GMNS and TransOMS.
// https://github.com/jdlph/TAP101

#define FLOAT_ACCURACY 1.0E-15
#define NO_COSTPARAMETERS 4
#define IVTT 0
#define OVTT 1
#define MONETARY 2
#define DIST 3

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define SQR(x) ((x) * (x))
#define FABS(x) ((x) >= 0 ? (x) : (-x))

#define INVALID -1

#define VALID(x) ((x) != -1)

#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <chrono>  // for high_resolution_clock
#define MAX_MODE_TYPES  10

//#ifndef _win32
//void fopen_s(FILE** file, const char* fileName, const char* mode)
//{
//    *file = fopen(fileName, mode);
//}
//#endif

typedef double cost_vector[NO_COSTPARAMETERS];

struct link_record {
	int internal_from_node_id;
	int internal_to_node_id;
	int link_id;
	int link_type;
	int external_from_node_id;
	int external_to_node_id;

	double Lane_Capacity;
	double Link_Capacity;
	double lanes;
	double FreeTravelTime;
	double free_speed;
	double Cutoff_Speed;


	double VDF_Alpha;
	double BoverC;
	double VDF_Beta;
	double VDF_plf;
	double Q_cd, Q_n, Q_cp, Q_s;

	double length;
	double Speed;
	std::string allowed_uses;
	int mode_allowed_use[MAX_MODE_TYPES];
	double mode_MainVolume[MAX_MODE_TYPES];
	double mode_Base_Volume[MAX_MODE_TYPES];
	double mode_SubVolume[MAX_MODE_TYPES];
	double mode_SDVolume[MAX_MODE_TYPES];

	double mode_Toll[MAX_MODE_TYPES];
	double mode_AdditionalCost[MAX_MODE_TYPES];

	double Travel_time;  // final travel time used in assignment 
	double BPR_TT;  // BPR_based travel time  for the entire assignment period  
	double QVDF_TT;  // QVDF_based travel time for the entire assignment period  


	double GenCost;
	double GenCostDer;
	double Ref_volume;
	double Base_volume;
	double Obs_volume;
	std::string geometry;
	link_record()
	{
		link_type = 1;

		VDF_Alpha = 0.15;
		VDF_Beta = 4;
		VDF_plf = 1;

		Q_cd = 1.0;
		Q_n = 1.0;
		Q_cp = 0.28125 /*0.15*15/8*/;
		Q_s = 4;

		Travel_time = 0;
		BPR_TT = 0;
		QVDF_TT = 0;
		Ref_volume = 0;
		Base_volume = 0;
		Obs_volume = -1;
	}
	void setup(int num_of_modes)
	{
		link_type = 1;
		VDF_Alpha = 0.15;
		VDF_Beta = 4;
		VDF_plf = 1;

		Q_cd = 1.0;
		Q_n = 1.0;
		Q_cp = 0.28125 /*0.15*15/8*/;
		Q_s = 4;

		Travel_time = 0;
		BPR_TT = 0;
		QVDF_TT = 0;
		Ref_volume = 0;
		Base_volume = 0;
		Obs_volume = -1;
		for (int m = 1; m <= num_of_modes; m++)
		{
			mode_Base_Volume[m] = 0;
		}
	}
};

struct mode_type {
	std::string mode_type;
	float vot;
	float pce;
	float occ;
	int dedicated_shortest_path;
	std::string demand_file;
};


mode_type g_mode_type_vector[MAX_MODE_TYPES];
int g_metric_system_flag = 0; 
void StatusMessage(const char* group, const char* format, ...);
void ExitMessage(const char* format, ...);
#include "TAPlite.h"

struct link_record* Link;
int* FirstLinkFrom;
int* LastLinkFrom;
sorted_list* LinksTo;

double Link_GenCost(int k, double* Volume);
double LinkCost_Integral(int k, double* Volume);
double Link_GenCostDer(int k, double* Volume);

void* Alloc_1D(int dim1, size_t size);
void** Alloc_2D(int dim1, int dim2, size_t size);
void Free_2D(void** Array, int dim1, int dim2);
void*** Alloc_3D(int dim1, int dim2, int dim3, size_t size);
void Free_3D(void*** Array, int dim1, int dim2, int dim3);

/* LINK VECTOR FUNCTIONS */
void ClearVolume(double* VolumeArray);
void VolumeDifference(double* Volume1, double* Volume2, double* Difference);
void UpdateVolume(double* MainVolume, double* SubVolume, double Lambda);

/* LINK COST FUNCTIONS */
void UpdateLinkAdditionalCost(void);
double UpdateLinkCost(double* Volume);
void UpdateLinkCostDer(double* Volume);
void GetLinkTravelTimes(double* Volume, double* TravelTime);
double TotalLinkCost(double* Volume);

/*  LINK OBJECTIVE FUNCTION */
double OF_Links(double* MainVolume);
double OF_LinksDirectionalDerivative(double* MainVolume, double* SDVolume, double Lambda);

void InitLinks(void);
void CloseLinks(void);

int Read_ODflow(double* TotalODflow, int* number_of_modes, int* no_zones);

int Minpath(int m, int Orig, int* PredLink, double* Cost_to);
double FindMinCostRoutes(int** MinPathPredLink);
void All_or_Nothing_Assign(double** ODflow, int** MinPathPredLink, double* Volume);


double LinksSDLineSearch(double* MainVolume, double* SDVolume);


/* Gloabal variables */

int no_zones, number_of_modes, no_nodes, number_of_links, FirstThruNode;
int TotalAssignIterations = 20;
int demand_period_starting_hours = 7;
int	demand_period_ending_hours = 8;
int g_tap_log_file = 0;
int g_base_demand_mode = 1;
int g_ODME_mode = 0;
double g_ODME_obs_VMT = -1;
double g_System_VMT = 0;

double g_ODME_link_volume_penalty = 0.01;  // relative weight on volume , convert the deviation of link volume to the travel time 
double g_ODME_VMT_penalty = 0.01;  // relative weight on VMT , convert the deviation of VMT to the travel time 
FILE* summary_log_file;


double** ODflow, TotalODflow;
double* TotalOFlow;
int* zone_outbound_link_size;

double*** MDODflow;
double*** MDDiffODflow;  // D^c - D^b

double*** MDRouteCost;
/* Local Declarations */
/* void FW(void); Should there be a function for fw, or should it be included in main? */
static void Init(int input_number_of_modes, int input_no_zones);
static void Close();
static void InitODflow(int input_number_of_modes, int input_no_zones);
static void CloseODflow(void);
/* End of local declarations. */

FILE* logfile;
int shortest_path_log_flag = 0;
int baseODDemand_loaded_flag = 0;
int baselinkvolume_loaded_flag = 0;

FILE* link_performance_file;

#define INVALID -1      /* Represents an invalid value. */
#define BIGM 9999999      /* Represents an invalid value. */
#define WAS_IN_QUEUE -7 /* Shows that the node was in the queue before. (7 is for luck.) */

double Link_QueueVDF(int k, double Volume, double& IncomingDemand, double& DOC, double& P, double& t0, double& t2, double& t3, double& vt2, double& Q_mu, double& Q_gamma, double& congestion_ref_speed,
	double& avg_queue_speed, double& avg_QVDF_period_speed, double& Severe_Congestion_P, double model_speed[300]);

int Minpath(int mode, int Orig, int* PredLink, double* CostTo)
{
	int node, now, NewNode, k, Return2Q_Count = 0;
	// Orig is the zone number 
	// now is the internal node id (Seq. no)
	double NewCost;
	int* QueueNext;
	int QueueFirst, QueueLast;

	QueueNext = (int*)Alloc_1D(no_nodes, sizeof(int));

	for (node = 1; node <= no_nodes; node++)
	{
		QueueNext[node] = INVALID;
		CostTo[node] = BIGM;
		PredLink[node] = INVALID;
	}

	now = g_map_external_node_id_2_node_seq_no[Orig];  // mapping from external zone id of Orig (which is defined in demand.csv_ to the corresponding node id (== zone_id) and then to the node internal number 
	int internal_node_id_for_origin_zone = now;
	QueueNext[now] = WAS_IN_QUEUE;
	PredLink[now] = INVALID;
	CostTo[now] = 0.0;

	QueueFirst = QueueLast = INVALID;

	while ((now != INVALID) && (now != WAS_IN_QUEUE))
	{
		if (now >= FirstThruNode || now == internal_node_id_for_origin_zone)  // this is the key implementation for FirstThruNode on connector
		{
			for (k = FirstLinkFrom[now]; k <= LastLinkFrom[now]; k++)
			{

				if (Link[k].mode_allowed_use[mode] == 0)  // implementation for allowed uses 
					continue;
				/* For every link that terminate at "now": */

				NewNode = Link[k].internal_to_node_id;
				NewCost = CostTo[now] + Link[k].Travel_time + Link[k].mode_AdditionalCost[mode];

				if (CostTo[NewNode] > NewCost)
				{
					/* If the new label is better than the old one, correct it, and make sure that
					 * the new node to the queue. */

					CostTo[NewNode] = NewCost;
					PredLink[NewNode] = k;  // PredLink is coded in terms of internal node id 

					/* If the new node was in the queue before, add it as the first in the queue. */
					if (QueueNext[NewNode] == WAS_IN_QUEUE)
					{
						QueueNext[NewNode] = QueueFirst;
						QueueFirst = NewNode;
						if (QueueLast == INVALID)
							QueueLast = NewNode;
						Return2Q_Count++;
					}

					/* If the new node is not in the queue, and wasn't there before, add it at the
					 * end of the queue. */
					else if (QueueNext[NewNode] == INVALID && NewNode != QueueLast)
					{
						if (QueueLast != INVALID)
						{ /*Usually*/
							QueueNext[QueueLast] = NewNode;
							QueueLast = NewNode;
						}
						else
						{ /* If the queue is empty, initialize it. */
							QueueFirst = QueueLast = NewNode;
							QueueNext[QueueLast] = INVALID;
						}
					}

					/* If the new node is in the queue, just leave it there. (Do nothing) */
				}
			}
		}

		/* Get the first node out of the queue, and use it as the current node. */
		now = QueueFirst;
		if ((now == INVALID) || (now == WAS_IN_QUEUE))
			break;

		QueueFirst = QueueNext[now];
		QueueNext[now] = WAS_IN_QUEUE;
		if (QueueLast == now)
			QueueLast = INVALID;
	}

	free(QueueNext);

	return (Return2Q_Count);
}

/* Find minimum cost routes .
Input: 	None
Output:	RouteCost - route generalized cost, by origin and destination
		MinPathSuccLink - trees of minimum cost routes, by destination and node. */
std::vector<int> Processor_origin_zones[50];

int g_number_of_processors = 4;

double FindMinCostRoutes(int*** MinPathPredLink)
{
	double** CostTo;

	CostTo = (double**)Alloc_2D(no_zones, no_nodes, sizeof(double));
	StatusMessage("Minpath", "Starting the minpath calculations.");
	double* system_least_travel_time_org_zone = (double*)Alloc_1D(no_zones, sizeof(double));

#pragma omp parallel for
	for (int p = 0; p < g_number_of_processors; p++)
	{

		for (int i = 0; i < Processor_origin_zones[p].size(); i++)
		{

		int Orig = Processor_origin_zones[p][i];  // get origin zone id


		if(TotalOFlow[Orig] < 0.00001)  // only work on positive zone flow 
			continue;

		system_least_travel_time_org_zone[Orig] = 0;  // reset it before mode based computing 

		for (int m = 1; m <= number_of_modes; m++)
		{
			if (g_mode_type_vector[m].dedicated_shortest_path == 0)  // skip the shortest path computing
				continue;

			Minpath(m, Orig, MinPathPredLink[m][Orig], CostTo[Orig]);
			if (MDRouteCost != NULL)
			{
				for (int Dest = 1; Dest <= no_zones; Dest++)
				{
					MDRouteCost[m][Orig][Dest] = BIGM; // initialization 

					if (MDODflow[m][Orig][Dest] > 0.000001)
					{
						if (CostTo[Orig][Dest] < 0.0)
							ExitMessage("Negative cost %lg from Origin %d to Destination %d.",
								(double)CostTo[Orig][Dest], Orig, Dest);

						// CostTo is coded as internal node id 
						int  internal_node_id_for_destination_zone  = g_map_external_node_id_2_node_seq_no[Dest]; 

						if (CostTo[Orig][internal_node_id_for_destination_zone] <= BIGM - 1)  // feasible cost 
						{
							MDRouteCost[m][Orig][Dest] = CostTo[Orig][internal_node_id_for_destination_zone];

							system_least_travel_time_org_zone[Orig] += MDRouteCost[m][Orig][Dest] * MDODflow[m][Orig][Dest] * g_mode_type_vector[m].pce;

						}
						else
						{
							int debug_flag = 1;
						}
					}
				}
			}
		}
	}

	}

	double system_least_travel_time = 0;
	for (int Orig = 1; Orig <= no_zones; Orig++)
	{
		system_least_travel_time += system_least_travel_time_org_zone[Orig];
	}


	// free(CostTo);
	Free_2D((void**)CostTo, no_zones, no_nodes);
	free(system_least_travel_time_org_zone);
	StatusMessage("Minpath", "Found all minpath.");
	return system_least_travel_time;
}

/* Assign OD flows to links according to the routes in MinPathPredLink. */

// Define a global 3D vector to store link indices for each OD pair

// Global 5D vector for storing link sequences
std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> linkIndices;

void AddLinkSequence(int m, int Orig, int Dest, int route_id, const std::vector<int>& linkIDs)
{
	if (linkIndices.size() == 0)
		return; 

	// Ensure we are within bounds before adding the link sequence
	if (Orig > 0 && Orig < linkIndices[m].size() &&
		Dest > 0 && Dest < linkIndices[m][Orig].size() &&
		route_id >= 0 && route_id < linkIndices[m][Orig][Dest].size())
	{
		linkIndices[m][Orig][Dest][route_id] = linkIDs;  // Add link sequence to the 5D vector
	}
	else
	{
		std::cerr << "Error: Invalid indices for adding link sequence." << std::endl;
	}
}





void All_or_Nothing_Assign(int Assignment_iteration_no, double*** ODflow, int*** MinPathPredLink, double* Volume)
{
	printf("All or nothing assignment\n");

	auto start0 = std::chrono::high_resolution_clock::now();
	double** ProcessorVolume;
	double*** ProcessorModeVolume;


	ProcessorVolume = (double**)Alloc_2D(number_of_links, g_number_of_processors, sizeof(double));

	if (ProcessorVolume == nullptr) {
		std::cerr << "Error: Memory allocation for ProcessorVolume failed." << std::endl;
		// Handle the error, e.g., exit the program, clean up resources, etc.
		exit(EXIT_FAILURE);
	}

	ProcessorModeVolume = (double***)Alloc_3D(number_of_links, number_of_modes, g_number_of_processors, sizeof(double));
	if (ProcessorModeVolume == nullptr) {
		std::cerr << "Error: Memory allocation for ProcessorModeVolume failed." << std::endl;
		// Handle the error, e.g., exit the program, clean up resources, etc.
		exit(EXIT_FAILURE);
	}

#pragma omp parallel for
	for (int k = 1; k <= number_of_links; k++)
	{


		for (int p = 0; p < g_number_of_processors; p++)
		{
			ProcessorVolume[k][p] = 0;

			for (int m = 1; m <= number_of_modes; m++)
				ProcessorModeVolume[k][m][p] = 0;
		}
	}
	//StatusMessage("Assign", "Starting assign.");

	if(Assignment_iteration_no == 0)
	{
		printf("The list of zero-volume zones:");
		for (int Orig = 1; Orig <= no_zones; Orig++)
		{

			if (TotalOFlow[Orig] < 0.00001)  // only work on positive zone flow 
			{ 
				printf("%d,", Orig);
			}

		}

		printf("\n");

		printf("The list of zones without outbound connecting links:");
		for (int Orig = 1; Orig <= no_zones; Orig++)
		{

			if (zone_outbound_link_size[Orig] == 0)  // there is no outbound link from the origin 
			{
				printf("%d,", Orig);
			}
		}
		printf("\n");
	}
#pragma omp parallel for
	for (int p = 0; p < g_number_of_processors; p++)
	{

		for (int i = 0; i < Processor_origin_zones[p].size(); i++)
		{

			int Orig = Processor_origin_zones[p][i];  // get origin zone id
			
			if (TotalOFlow[Orig] < 0.00001)  // only work on positive zone flow 
				continue;

			if (zone_outbound_link_size[Orig] == 0)  // there is no outbound link from the origin 
				continue; 


			int Dest, k;
			int CurrentNode;
			double RouteFlow;
			std::vector<int> currentLinkSequence; // Temporary vector to store link indices


			for (int m = 1; m <= number_of_modes; m++)
			{



				//	printf("Assign", "Assigning origin %6d.", Orig);
				for (Dest = 1; Dest <= no_zones; Dest++)  // Dest zone 
				{
					if (Dest == Orig)
						continue;

					if (zone_outbound_link_size[Dest] == 0)  // there is no outbound or inbound link from the origin zone
						continue;

					if (shortest_path_log_flag || Assignment_iteration_no == 0)
						currentLinkSequence.clear();

					RouteFlow = ODflow[m][Orig][Dest];  //test
					if (RouteFlow == 0)
						continue;

					if (g_mode_type_vector[m].dedicated_shortest_path == 0)  // skip the shortest path computing
					{
						if (MDRouteCost[1][Orig][Dest] >= BIGM - 1)  // if feasible 
							continue;
					}
					else
					{
						if (MDRouteCost[m][Orig][Dest] >= BIGM - 1)  // if feasible 
							continue;
					}

					CurrentNode = Dest;
					CurrentNode = g_map_external_node_id_2_node_seq_no[Dest];  // mapping from external  zone id of Dest (which is defined in demand.csv_ to the corresponding node id (== zone_id) and then to the node internal number 
					int internal_node_for_origin_node = g_map_external_node_id_2_node_seq_no[Orig];  // mapping from external  zone id of Orig (which is defined in demand.csv_ to the corresponding node id (== zone_id) and then to the node internal number 
					// MinPathPredLink is coded as internal node id 
					// 
					//double total_travel_time = 0;
					//double total_length = 0;
					//double total_FFTT = 0;

					while (CurrentNode != internal_node_for_origin_node)
					{
						if (g_mode_type_vector[m].dedicated_shortest_path == 0)  // skip the shortest path computing
							k = MinPathPredLink[1][Orig][CurrentNode];  // default to mode 1
						else
							k = MinPathPredLink[m][Orig][CurrentNode];

						if (k <= 0 || k > number_of_links || k == INVALID)
						{
							printf("A problem in All_or_Nothing_Assign() Invalid pred for node seq no %d Orig zone= %d \n\n", CurrentNode, Orig);
							break;
						}
						ProcessorVolume[k][p] += RouteFlow * g_mode_type_vector[m].pce;
						ProcessorModeVolume[k][m][p] += RouteFlow;  // pure volume 

						CurrentNode = Link[k].internal_from_node_id;

						if (CurrentNode <= 0 || CurrentNode > no_nodes )
						{
							printf("A problem in All_or_Nothing_Assign() Invalid node seq no %d Orig zone = %d \n\n", CurrentNode, Orig);
							break;
						}

						if(linkIndices.size() >0)
						{
						if (shortest_path_log_flag || Assignment_iteration_no == 0)
						{
#pragma omp critical
							{
								currentLinkSequence.push_back(k); // Store the link index
							}

						}
					}

							if (linkIndices.size() > 0)
							{
								if (shortest_path_log_flag || Assignment_iteration_no == 0)
								{
									AddLinkSequence(m, Orig, Dest, Assignment_iteration_no, currentLinkSequence);
									// Store the link sequence for this OD pair

								}
						}
						}

				}
			}

		}


	}

	if (Assignment_iteration_no == 0)  // with base demand
	{
		for (int k = 1; k <= number_of_links; k++)
		{
			Volume[k] = Link[k].Base_volume;

			for (int p = 0; p < g_number_of_processors; p++)
			{
				Volume[k] += ProcessorVolume[k][p];
			}

			for (int m = 1; m <= number_of_modes; m++)
			{
				Link[k].mode_MainVolume[m] = Link[k].mode_Base_Volume[m];

				for (int p = 0; p < g_number_of_processors; p++)
				{
					Link[k].mode_MainVolume[m] += ProcessorModeVolume[k][m][p];
				}
			}

		}

	}
	else
	{  // all or nothing 
#pragma omp parallel for
		for (int k = 1; k <= number_of_links; k++)
		{
			Volume[k] = 0;

			for (int p = 0; p < g_number_of_processors; p++)
			{
				Volume[k] += ProcessorVolume[k][p];
			}

			for (int m = 1; m <= number_of_modes; m++)
			{
				Link[k].mode_SubVolume[m] = 0.0;

				for (int p = 0; p < g_number_of_processors; p++)
				{
					Link[k].mode_SubVolume[m] += ProcessorModeVolume[k][m][p];
				}
			}

		}
	}


	Free_2D((void**)ProcessorVolume, number_of_links, g_number_of_processors);
	Free_3D((void***)ProcessorModeVolume, number_of_links, number_of_modes, g_number_of_processors);

	auto end0 = std::chrono::high_resolution_clock::now();

	// Calculate the duration in milliseconds

	// Calculate the duration in seconds
	auto duration = end0 - start0;

	// Convert to hours, minutes, seconds
	auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1));
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration % std::chrono::seconds(1));

	printf("All or nothing assignment: %lld hours %lld minutes %lld seconds %lld ms\n", hours.count(), minutes.count(), seconds.count(), milliseconds.count());

}




#include <fstream>  // for file output

#include <unordered_map>  // For hash table (unordered_map)
#include <fstream>        // For file output
#include <string>
#include <vector>
#include <iostream>

// Hash table to store unique (node sum, link sum) combinations


void OutputRouteDetails(const std::string& filename)
{
	std::ofstream outputFile(filename);  // Open the file for writing

	if (linkIndices.size() == 0)
		return; 
	// Write the CSV header in lowercase
	outputFile << "mode,route_id,o_zone_id,d_zone_id,unique_route_id,node_ids,link_ids,total_distance_mile,total_distance_km,total_free_flow_travel_time,total_travel_time,route_key,volume,\n";
	for (int m = 1; m < linkIndices.size(); ++m)
	{
		for (int Orig = 1; Orig < linkIndices[m].size(); ++Orig)
		{
			for (int Dest = 1; Dest < linkIndices[m][Orig].size(); ++Dest)
			{
				std::unordered_map<std::string, bool> uniqueRoutes;
				int unique_route_id = 1;
				for (int route_id = 0; route_id < linkIndices[m][Orig][Dest].size(); ++route_id)
				{
					if (!linkIndices[m][Orig][Dest][route_id].empty())
					{

						//for (int route_id_2 = 0; route_id_2 < route_id; ++route_id_2)
						//{  // mimic the route swiching machanisum, from route_id_2 to route_id, using the step size 
						//}
					
						double totalDistance = 0.0;
						double totalFreeFlowTravelTime = 0.0;
						double totalTravelTime = 0.0;
						std::string nodeIDsStr;
						std::string linkIDsStr;

						int nodeSum = 0;  // Sum of node IDs
						int linkSum = 0;  // Sum of link IDs

						// Collect node IDs, link indices, compute total distance, travel times, and calculate sums
						for (int i = linkIndices[m][Orig][Dest][route_id].size() - 1; i >= 0; --i)
						{
							int k = linkIndices[m][Orig][Dest][route_id][i];

							// Append the from_node_id for each link and calculate the node sum
							int fromNodeID = Link[k].external_from_node_id;
							nodeIDsStr += std::to_string(fromNodeID) + ";";
							nodeSum += fromNodeID;

							// Append the link index (link ID) to the string and calculate the link sum
							linkIDsStr += std::to_string(k) + ";";
							linkSum += k;

							// Sum up the total distance and travel times
							totalDistance += Link[k].length;
							totalFreeFlowTravelTime += Link[k].FreeTravelTime;
							totalTravelTime += Link[k].Travel_time;

							// For the last link, also add the to_node_id
							if (i == 0)
							{
								int toNodeID = Link[k].external_to_node_id;
								nodeIDsStr += std::to_string(toNodeID);
								nodeSum += toNodeID;
							}
						}

						// Create a unique key based on the node sum and link sum
						std::string routeKey = std::to_string(nodeSum) + "_" + std::to_string(linkSum);

						// Check if this route (based on node and link sums) is already output
						if (uniqueRoutes.find(routeKey) == uniqueRoutes.end())
						{
							// This is a unique route, store it in the hash table
							uniqueRoutes[routeKey] = true;


							unique_route_id++;
						}
						//else
						//{
						//    // Duplicate path found, skipping output
						//    //std::cout << "Duplicate route skipped for Origin: " << Orig << ", Destination: " << Dest << "\n";
						//}
					}
				}



				///////////////////////////////////// print out 
				uniqueRoutes.clear();
				int unique_route_id_size = unique_route_id; 
				unique_route_id = 1; 

				for (int route_id = 0; route_id < linkIndices[m][Orig][Dest].size(); ++route_id)
				{
					if (!linkIndices[m][Orig][Dest][route_id].empty())
					{
						double totalDistance = 0.0;
						double totalFreeFlowTravelTime = 0.0;
						double totalTravelTime = 0.0;
						std::string nodeIDsStr;
						std::string linkIDsStr;

						int nodeSum = 0;  // Sum of node IDs
						int linkSum = 0;  // Sum of link IDs

						// Collect node IDs, link indices, compute total distance, travel times, and calculate sums
						for (int i = linkIndices[m][Orig][Dest][route_id].size() - 1; i >= 0; --i)
						{
							int k = linkIndices[m][Orig][Dest][route_id][i];

							// Append the from_node_id for each link and calculate the node sum
							int fromNodeID = Link[k].external_from_node_id;
							nodeIDsStr += std::to_string(fromNodeID) + ";";
							nodeSum += fromNodeID;

							// Append the link index (link ID) to the string and calculate the link sum
							linkIDsStr += std::to_string(k) + ";";
							linkSum += k;

							// Sum up the total distance and travel times
							totalDistance += Link[k].length;
							totalFreeFlowTravelTime += Link[k].FreeTravelTime;
							totalTravelTime += Link[k].Travel_time;

							// For the last link, also add the to_node_id
							if (i == 0)
							{
								int toNodeID = Link[k].external_to_node_id;
								nodeIDsStr += std::to_string(toNodeID);
								nodeSum += toNodeID;
							}
						}

						// Create a unique key based on the node sum and link sum
						std::string routeKey = std::to_string(nodeSum) + "_" + std::to_string(linkSum);

						// Check if this route (based on node and link sums) is already output
						if (uniqueRoutes.find(routeKey) == uniqueRoutes.end())
						{
							// This is a unique route, store it in the hash table
							uniqueRoutes[routeKey] = true;

							// Remove trailing space from the link IDs string
							if (!linkIDsStr.empty())
								linkIDsStr.pop_back();

							float od_volume = MDODflow[m][Orig][Dest];
							float route_volume = od_volume / unique_route_id_size;




							// Write the data for this OD pair and route to the CSV file
							outputFile << g_mode_type_vector[m].mode_type.c_str() << ","
								<< route_id << "," << Orig << "," << Dest << "," << unique_route_id << ","
								<< nodeIDsStr << "," << linkIDsStr << ","
								<< totalDistance << "," << totalDistance*1.609 << "," << totalFreeFlowTravelTime << ","
								<< totalTravelTime << "," << routeKey.c_str() << "," << route_volume << "\n";

							unique_route_id++;
						}
						//else
						//{
						//    // Duplicate path found, skipping output
						//    //std::cout << "Duplicate route skipped for Origin: " << Orig << ", Destination: " << Dest << "\n";
						//}
					}
				}

			}
		}
	}

	// Close the file after writing
	outputFile.close();


	std::cout << "Output written to " << filename << std::endl;

}


void OutputODPerformance(const std::string& filename)
{
	std::ofstream outputFile(filename);  // Open the file for writing

	// Write the CSV header in lowercase
	outputFile << "mode,o_zone_id,d_zone_id,total_distance_mile,total_distance_km,total_free_flow_travel_time,total_congestion_travel_time,volume,\n";
	double grand_totalDistance = 0.0;
	double grand_totalFreeFlowTravelTime = 0.0;
	double grand_totalTravelTime = 0.0;
	double grand_total_count = 0;


	for (int m = 1; m < linkIndices.size(); ++m)
	{
		for (int Orig = 1; Orig < linkIndices[m].size(); ++Orig)
		{
			for (int Dest = 1; Dest < linkIndices[m][Orig].size(); ++Dest)
			{
				std::unordered_map<std::string, bool> uniqueRoutes;
				int unique_route_id = 1;
				for (int route_id = 0; route_id < linkIndices[m][Orig][Dest].size(); ++route_id)
				{
					if (!linkIndices[m][Orig][Dest][route_id].empty())
					{
						double totalDistance = 0.0;
						double totalFreeFlowTravelTime = 0.0;
						double totalTravelTime = 0.0;
						std::string nodeIDsStr;
						std::string linkIDsStr;

						long nodeSum = 0;  // Sum of node IDs
						long linkSum = 0;  // Sum of link IDs

						// Collect node IDs, link indices, compute total distance, travel times, and calculate sums
						for (int i = linkIndices[m][Orig][Dest][route_id].size() - 1; i >= 0; --i)
						{
							long k = linkIndices[m][Orig][Dest][route_id][i];

							// Append the from_node_id for each link and calculate the node sum
							int fromNodeID = Link[k].external_from_node_id;
							nodeIDsStr += std::to_string(fromNodeID) + ";";
							nodeSum += fromNodeID;

							// Append the link index (link ID) to the string and calculate the link sum
							linkIDsStr += std::to_string(k) + ";";
							linkSum += k;

							// Sum up the total distance and travel times
							totalDistance += Link[k].length;
							totalFreeFlowTravelTime += Link[k].FreeTravelTime;
							totalTravelTime += Link[k].Travel_time;

							// For the last link, also add the to_node_id
							if (i == 0)
							{
								int toNodeID = Link[k].external_to_node_id;
								nodeIDsStr += std::to_string(toNodeID);
								nodeSum += toNodeID;
							}
						}

						// Create a unique key based on the node sum and link sum
						std::string routeKey = std::to_string(nodeSum) + "_" + std::to_string(linkSum);

						// Check if this route (based on node and link sums) is already output
						if (uniqueRoutes.find(routeKey) == uniqueRoutes.end())
						{
							// This is a unique route, store it in the hash table
							uniqueRoutes[routeKey] = true;

							// Remove trailing space from the link IDs string
							if (!linkIDsStr.empty())
								linkIDsStr.pop_back();

							float volume = MDODflow[m][Orig][Dest]; 


							grand_totalDistance += totalDistance * volume;
							grand_totalFreeFlowTravelTime += totalFreeFlowTravelTime * volume;
							grand_totalTravelTime += totalTravelTime * volume;
							grand_total_count += volume;
							// Write the data for this OD pair and route to the CSV file
							outputFile << g_mode_type_vector[m].mode_type.c_str() << "," << Orig << "," << Dest << ","
								<< totalDistance << "," << totalDistance*1.609 << "," << totalFreeFlowTravelTime << ","
								<< totalTravelTime << "," << volume << "\n";

							unique_route_id++;
						}
						//else
						//{
						//    // Duplicate path found, skipping output
						//    //std::cout << "Duplicate route skipped for Origin: " << Orig << ", Destination: " << Dest << "\n";
						//}
					}
					break; // only compute travel time for the first route 
				}
			}
		}
	}
	if (grand_total_count < 0.001)
		grand_total_count = 0.001;

	std::cout << "OD performance summary: avg distance = "
		<< grand_totalDistance / grand_total_count << " miles, "
		<< ", avg free-flow travel time = "
		<< grand_totalFreeFlowTravelTime / grand_total_count << " min, "
		<< ", avg total travel time = "
		<< grand_totalTravelTime / grand_total_count << " min"
		<< ", avg travel time index = "
		<< grand_totalTravelTime / grand_totalFreeFlowTravelTime << " "
		<< std::endl;
	// Close the file after writing
	outputFile.close();
	std::cout << "Output written to " << filename << std::endl;
}


int get_number_of_nodes_from_node_file(int& number_of_zones, int& l_FirstThruNode)
{
	number_of_zones = 0;
	CDTACSVParser parser_node;
	l_FirstThruNode = 1;
	int number_of_nodes = 0;

	if (parser_node.OpenCSVFile("node.csv", true))
	{
		while (parser_node.ReadRecord())  // if this line contains [] mark, then we will also read
			// field headers.
		{
			// Read node id
			int node_id = 0;
			int zone_id = 0;
			parser_node.GetValueByFieldName("node_id", node_id);
			parser_node.GetValueByFieldName("zone_id", zone_id);

			if (zone_id >= 1 && zone_id != node_id)
			{
				printf("Error: zone_id should be the same as node_id but zone_id  = %d, node_id = %d\n", zone_id, node_id);
			}

			g_map_node_seq_no_2_external_node_id[number_of_nodes + 1] = node_id;
			g_map_external_node_id_2_node_seq_no[node_id] =
				number_of_nodes + 1;  // this code node sequential number starts from 1

			if (zone_id >= 1 && zone_id > number_of_zones)
				number_of_zones = zone_id;

			if (zone_id == 0 && l_FirstThruNode == 1 /* not initialized*/)
				l_FirstThruNode = number_of_nodes + 1;  //use sequential node id


			if (g_tap_log_file == 1)
			{
				fprintf(logfile, "node_id = %d, node_seq_no = %d\n", node_id, g_map_external_node_id_2_node_seq_no[node_id]);

			}

			number_of_nodes++;
		}


		parser_node.CloseCSVFile();
	}



	return number_of_nodes;
}

int get_number_of_links_from_link_file()
{
	CDTACSVParser parser_link;

	int number_of_links = 0;

	if (parser_link.OpenCSVFile("link.csv", true))
	{
		while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read
			// field headers.
		{
			int link_id = 0;
			int internal_from_node_id = 0;
			// Read link_id
			parser_link.GetValueByFieldName("link_id", link_id);
			parser_link.GetValueByFieldName("from_node_id", internal_from_node_id);

			number_of_links++;
		}

		parser_link.CloseCSVFile();
	}

	return number_of_links;
}

void createSettingsFile(const std::string& fileName) {
	std::ofstream file(fileName);
	if (!file.is_open()) {
		std::cerr << "Could not create the file: " << fileName << std::endl;
		return;
	}

	// Writing the headers of the CSV
	file << "metric_system,number_of_iterations,number_of_processors,demand_period_starting_hours,demand_period_ending_hours,base_demand_mode,route_output,log_file,odme_mode,odme_vmt\n";

	// Writing the sample data (from your provided file)
	file << "0,10,8,14,18,0,0,0,0,0\n";

	file.close();
	std::cout << "sample_settings.csv file created successfully!" << std::endl;
}

void read_settings_file()
{

	createSettingsFile("sample_settings.csv");

	CDTACSVParser parser_settings;

	if (parser_settings.OpenCSVFile("settings.csv", true))
	{
		while (parser_settings.ReadRecord())  // if this line contains [] mark, then we will also read
			// field headers.
		{
			g_number_of_processors = 4;
			parser_settings.GetValueByFieldName("metric_system", g_metric_system_flag);
			parser_settings.GetValueByFieldName("number_of_iterations", TotalAssignIterations);
			parser_settings.GetValueByFieldName("number_of_processors", g_number_of_processors);
			parser_settings.GetValueByFieldName("demand_period_starting_hours", demand_period_starting_hours);
			parser_settings.GetValueByFieldName("demand_period_ending_hours", demand_period_ending_hours);
			parser_settings.GetValueByFieldName("log_file", g_tap_log_file);
			parser_settings.GetValueByFieldName("base_demand_mode", g_base_demand_mode);
			parser_settings.GetValueByFieldName("odme_mode", g_ODME_mode);
			parser_settings.GetValueByFieldName("odme_vmt", g_ODME_obs_VMT);
			parser_settings.GetValueByFieldName("route_output", shortest_path_log_flag);



			

			break;
		}

		parser_settings.CloseCSVFile();
	}

	return;
}

void createModeTypeFile(const std::string& fileName) {
	std::ofstream file(fileName);
	if (!file.is_open()) {
		std::cerr << "Could not create the file: " << fileName << std::endl;
		return;
	}

	// Writing the headers of the CSV
	file << "mode_type,name,vot,pce,occ,demand_file,dedicated_shortest_path,\n";

	// Writing the sample data (from your provided file)
	file << "sov,DRIVE, 10, 1, 1,demand_1400_1800_sov.csv, 1\n";
	file << "hov,HOV, 10, 1, 2,demand_1400_1800_hov.csv, 1\n";
	file << "trk,truck, 10, 2, 1,demand_1400_1800_trk.csv, 0\n";

	std::cout << "sample_mode_type.csv file created successfully!" << std::endl;
}


void read_mode_type_file()
{
	createModeTypeFile("sample_mode_type.csv");
	CDTACSVParser parser_mode_type;

	if (parser_mode_type.OpenCSVFile("mode_type.csv", true))
	{

		number_of_modes = 0;
		while (parser_mode_type.ReadRecord())  // if this line contains [] mark, then we will also read
			// field headers.
		{
			number_of_modes += 1;
			if (number_of_modes < MAX_MODE_TYPES)
			{
				parser_mode_type.GetValueByFieldName("mode_type", g_mode_type_vector[number_of_modes].mode_type);
				parser_mode_type.GetValueByFieldName("vot", g_mode_type_vector[number_of_modes].vot);
				parser_mode_type.GetValueByFieldName("pce", g_mode_type_vector[number_of_modes].pce);
				parser_mode_type.GetValueByFieldName("occ", g_mode_type_vector[number_of_modes].occ);
				parser_mode_type.GetValueByFieldName("demand_file", g_mode_type_vector[number_of_modes].demand_file);

				g_mode_type_vector[number_of_modes].dedicated_shortest_path = 1;
				parser_mode_type.GetValueByFieldName("dedicated_shortest_path", g_mode_type_vector[number_of_modes].dedicated_shortest_path);

				if (number_of_modes == 1)
					g_mode_type_vector[1].dedicated_shortest_path = 1;  // reset; 

			}

		}

		parser_mode_type.CloseCSVFile();
	}

	if (number_of_modes == 0)  //default 
	{
		g_mode_type_vector[1].demand_file = "demand.csv";
		g_mode_type_vector[1].mode_type = "auto";
		g_mode_type_vector[1].vot = 10;
		g_mode_type_vector[1].pce = 1;
		g_mode_type_vector[1].occ = 1;
		g_mode_type_vector[1].dedicated_shortest_path = 1;  // reset; 
		number_of_modes = 1;
	}

	printf("number_of_modes = %d\n", number_of_modes);
	return;
}
// Initialize the 5D vector
void InitializeLinkIndices(int num_modes, int no_zones, int max_routes)
{


	auto start = std::chrono::high_resolution_clock::now();  // Start timing



	linkIndices = std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>(
		num_modes + 1,
		std::vector<std::vector<std::vector<std::vector<int>>>>(
			no_zones + 1,
			std::vector<std::vector<std::vector<int>>>(
				no_zones + 1,
				std::vector<std::vector<int>>(
					max_routes + 1
				)
			)
		)
	);

	// Record the end time
	auto end = std::chrono::high_resolution_clock::now();

	// Calculate the duration in milliseconds

	// Calculate the duration in seconds
	auto duration = end - start;

	// Convert to hours, minutes, seconds
	auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1));
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration % std::chrono::seconds(1));

	printf("Memmory creation time for 5D link path matrix: %lld hours %lld minutes %lld seconds %lld ms\n", hours.count(), minutes.count(), seconds.count(), milliseconds.count());
}

int main(int argc, char** argv)
{

	fopen_s(&summary_log_file, "summary_log_file.txt", "w");

	double* MainVolume, * SubVolume, * SDVolume, Lambda;
	int*** MDMinPathPredLink;

	read_settings_file();
	read_mode_type_file();
	fopen_s(&logfile, "TAP_log.csv", "w");  // Open the log file for writing.
	no_nodes = get_number_of_nodes_from_node_file(no_zones, FirstThruNode);
	number_of_links = get_number_of_links_from_link_file();

	printf("no_nodes= %d, no_zones = %d, FirstThruNode (seq No) = %d, number_of_links = %d\n", no_nodes, no_zones,
		FirstThruNode, number_of_links);

	fprintf(summary_log_file, "no_nodes= %d, no_zones = %d, FirstThruNode (seq No) = %d, number_of_links = %d\n", no_nodes, no_zones,
		FirstThruNode, number_of_links);

	fopen_s(&link_performance_file, "link_performance.csv", "w");
	if (link_performance_file == NULL)
	{
		printf("Error opening file!\n");
		return 1;
	}
	fclose(link_performance_file);

	fopen_s(&link_performance_file, "link_performance.csv", "a+");


	double system_wide_travel_time = 0;
	double system_least_travel_time = 0;

	//fprintf(logfile,
	//    "iteration_no,link_id,internal_from_node_id,internal_to_node_id,volume,capacity,voc,"
	//    "fftt,travel_time,delay\n");



	Init(number_of_modes, no_zones);


    InitializeLinkIndices(number_of_modes, no_zones, TotalAssignIterations);


		for (int Orig = 1; Orig <= no_zones; Orig++)  // initialization 
		{
			int p = Orig % g_number_of_processors;

			Processor_origin_zones[p].push_back(Orig);
		}


	int iteration_no = 0;
	MainVolume = (double*)Alloc_1D(number_of_links, sizeof(double));
	SDVolume = (double*)Alloc_1D(
		number_of_links, sizeof(double)); /* Compute search direction and sub-volume in the same place. */
	SubVolume = (double*)Alloc_1D(
		number_of_links, sizeof(double)); /* Compute search direction and sub-volume in the same place. */
	MDMinPathPredLink = (int***)Alloc_3D(number_of_modes, no_zones, no_nodes, sizeof(int));

	// Record the start time
	auto start = std::chrono::high_resolution_clock::now();

	for (int k = 1; k <= number_of_links; k++)
	{
		MainVolume[k] = Link[k].Base_volume;  // assign the base volume  to main volume 
	}

	system_wide_travel_time = UpdateLinkCost(MainVolume);  // set up the cost first using FFTT

	fprintf(link_performance_file,
		"iteration_no,link_id,from_node_id,to_node_id,volume,ref_volume,base_volume,obs_volume,"
		"capacity,D,doc,fftt,travel_time,VDF_alpha,VDF_beta,VDF_plf,speed_mph,speed_kmph,VMT,VHT,PMT,PHT,VHT_QVDF,PHT_QVDF,geometry,");

	fprintf(logfile, "iteration_no,link_id,from_node_id,to_node_id,volume,ref_volume,obs_volume,capacity,doc,fftt,travel_time,delay,");

	for (int m = 1; m <= number_of_modes; m++)
		fprintf(logfile, "mod_vol_%s,", g_mode_type_vector[m].mode_type.c_str());

	fprintf(logfile, "sub_main_vol,");

	for (int m = 1; m <= number_of_modes; m++)
		fprintf(logfile, "sub_mod_vol_%s,", g_mode_type_vector[m].mode_type.c_str());

	fprintf(logfile, "\n");


	for (int m = 1; m <= number_of_modes; m++)
		fprintf(link_performance_file, "mod_vol_%s,", g_mode_type_vector[m].mode_type.c_str());

	fprintf(link_performance_file, "P,t0,t2,t3,vt2_mph,vt2_kmph,mu,Q_gamma,free_speed_mph,cutoff_speed_mph,free_speed_kmph,cutoff_speed_kmph,congestion_ref_speed_mph,avg_queue_speed_mph,avg_QVDF_period_speed_mph,congestion_ref_speed_kmph,avg_queue_speed_kmph,avg_QVDF_period_speed_kmph,avg_QVDF_period_travel_time,Severe_Congestion_P,");

	for (int t = demand_period_starting_hours * 60; t < demand_period_ending_hours * 60; t += 5)
	{
		int hour = t / 60;
		int minute = t - hour * 60;

		fprintf(link_performance_file, "spd_mph_%02d:%02d,", hour, minute);
	}

	fprintf(link_performance_file, "\n");

	system_least_travel_time = FindMinCostRoutes(MDMinPathPredLink);
	All_or_Nothing_Assign(iteration_no, MDDiffODflow, MDMinPathPredLink, MainVolume);  // here we use MDDiffODflow as our OD search direction of D^c - D^b,



	system_wide_travel_time = UpdateLinkCost(MainVolume);
	double gap = (system_wide_travel_time - system_least_travel_time) /
		(fmax(0.1, system_least_travel_time)) * 100;

	if (gap < 0)
	{
		int ii = 0;
	}

	printf("iter No = %d, sys. TT =  %lf, least TT =  %lf, gap = %f %%\n", iteration_no,
		system_wide_travel_time, system_least_travel_time, gap);
	fprintf(summary_log_file, "iter No = %d, sys. TT =  %lf, least TT =  %lf, gap = %f %%\n", iteration_no,
		system_wide_travel_time, system_least_travel_time, gap);

	if (g_tap_log_file == 1)
	{
		for (int k = 1; k <= number_of_links; k++)
		{
			fprintf(logfile, "%d,%d,%d,%d,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,",
				iteration_no, k, Link[k].external_from_node_id, Link[k].external_to_node_id,
				MainVolume[k], Link[k].Ref_volume, Link[k].Obs_volume, Link[k].Link_Capacity,
				MainVolume[k] / fmax(0.01, Link[k].Link_Capacity), Link[k].FreeTravelTime,
				Link[k].Travel_time, Link[k].Travel_time - Link[k].FreeTravelTime);

			for (int m = 1; m <= number_of_modes; m++)
				fprintf(logfile, "%2lf,", Link[k].mode_MainVolume[m]);

			fprintf(logfile, "%2lf,", MainVolume[k]);

			for (int m = 1; m <= number_of_modes; m++)
				fprintf(logfile, "%2lf,", Link[k].mode_SubVolume[m]);

			fprintf(logfile, "\n");

		}
	}
	std::chrono::duration<double, std::milli> duration_FindMinCostRoutes;
	std::chrono::duration<double, std::milli> duration_All_or_Nothing_Assign;
	std::chrono::duration<double, std::milli> duration_LinksSDLineSearch;
	auto start0 = std::chrono::high_resolution_clock::now();  // Start timing


	for (iteration_no = 1; iteration_no < TotalAssignIterations; iteration_no++)
	{
		system_least_travel_time = FindMinCostRoutes(MDMinPathPredLink);  // the one right before the assignment iteration 
		

		g_System_VMT = 0;

		for (int k = 1; k <= number_of_links; k++)
		{
			if (Link[k].link_type >= 1)  // only include physical links
			{
				g_System_VMT += MainVolume[k] * Link[k].length;
			}
		}

		
		All_or_Nothing_Assign(iteration_no, MDODflow, MDMinPathPredLink, SubVolume); // this uees D^c, subvolume is Y, 
		auto end1 = std::chrono::high_resolution_clock::now();    // End timing
		auto start2 = std::chrono::high_resolution_clock::now();  // Start timing

		VolumeDifference(SubVolume, MainVolume, SDVolume); /* Which yields the search direction. SDVolume = y-X */


		Lambda = LinksSDLineSearch(MainVolume, SDVolume);
	
		// MSA options
	 //   Lambda = 1.0 / (iteration_no + 1);



		UpdateVolume(MainVolume, SDVolume, Lambda);  // x(k+1) = x(k) +lambda * (y-x) MainVolume is MainVolume
		system_wide_travel_time = UpdateLinkCost(MainVolume);

		//system_least_travel_time = FindMinCostRoutes(MinPathPredLink);  // the one right after the updated link cost 

		gap = (system_wide_travel_time - system_least_travel_time) /
			(fmax(0.1, system_least_travel_time)) * 100;
		auto end2 = std::chrono::high_resolution_clock::now();    // End timing

		duration_LinksSDLineSearch += end2 - start2;  // Compute duration in ms


		if (gap < 0)
		{
			int ii = 0;
		}

		printf("iter No = %d, Lambda = %f, g_System_VMT = %.1f, sys. TT =  %.1f, least TT =  %.1f, gap = %f %%\n",
			iteration_no, Lambda, g_System_VMT, system_wide_travel_time, system_least_travel_time, gap);
		fprintf(summary_log_file, "iter No = %d, Lambda = %f, g_System_VMT = %f, sys. TT =  %lf, least TT =  %lf, gap = %f %% \n",
			iteration_no, Lambda, g_System_VMT, system_wide_travel_time, system_least_travel_time, gap);

		if (g_tap_log_file == 1)
		{
			for (int k = 1; k <= number_of_links; k++)
			{
				fprintf(logfile, "%d,%d,%d,%d,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,",
					iteration_no, k, Link[k].external_from_node_id, Link[k].external_to_node_id,
					MainVolume[k], Link[k].Ref_volume, Link[k].Link_Capacity,
					MainVolume[k] / fmax(0.01, Link[k].Link_Capacity), Link[k].FreeTravelTime,
					Link[k].Travel_time, Link[k].Travel_time - Link[k].FreeTravelTime);

				for (int m = 1; m <= number_of_modes; m++)
					fprintf(logfile, "%2lf,", Link[k].mode_MainVolume[m]);

				fprintf(logfile, "%2lf,", SDVolume[k]);

				for (int m = 1; m <= number_of_modes; m++)
					fprintf(logfile, "%2lf,", Link[k].mode_SDVolume[m]);

				fprintf(logfile, "\n");

			}
		}


	}


	// Record the end time
	auto end = std::chrono::high_resolution_clock::now();

	// Calculate the duration in milliseconds

	// Calculate the duration in seconds
	auto duration = end - start;

	// Convert to hours, minutes, seconds
	auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1));
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration % std::chrono::seconds(1));

	printf("CPU running time: %lld hours %lld minutes %lld seconds %lld ms\n",
		hours.count(), minutes.count(), seconds.count(), milliseconds.count());

	// Log the result with hours, minutes, seconds, and milliseconds
	fprintf(summary_log_file, "CPU running time: %lld hours %lld minutes %lld seconds %lld ms\n",
		hours.count(), minutes.count(), seconds.count(), milliseconds.count());

	//for (int k = 1; k <= number_of_links; k++)
	//{
	//    fprintf(logfile, "%d,%d,%d,%d,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n", iteration_no, k,
	//        Link[k].external_from_node_id, Link[k].external_to_node_id, MainVolume[k],
	//        Link[k].Capacity, MainVolume[k] / fmax(0.01, Link[k].Capacity),
	//        Link[k].FreeTravelTime, Link[k].Travel_time,
	//        Link[k].Travel_time - Link[k].FreeTravelTime);
	//}

	for (int k = 1; k <= number_of_links; k++)
	{




		double P = 0;
		double vt2 = Link[k].Cutoff_Speed;
		double mu = Link[k].Lane_Capacity;
		double Severe_Congestion_P;
		double model_speed[300];
		double t0 = 0;
		double t2 = 0;
		double t3 = 0;
		double Q_gamma = 0;
		double congestion_ref_speed = 0;
		double avg_queue_speed = 0;
		double avg_QVDF_period_speed = 0;
		double IncomingDemand = 0;
		double DOC = 0;
		Link_QueueVDF(k, MainVolume[k], IncomingDemand, DOC, P, t0, t2, t3, vt2,  mu, Q_gamma, congestion_ref_speed, avg_queue_speed, avg_QVDF_period_speed, Severe_Congestion_P, model_speed);


		double VMT, VHT, PMT, PHT, VHT_QVDF, PHT_QVDF;
		VMT = 0; VHT = 0;  PMT = 0; PHT = 0; VHT_QVDF = 0; PHT_QVDF = 0;
		for (int m = 1; m <= number_of_modes; m++)
		{
			VMT += MainVolume[k] * Link[k].length;
			VHT += MainVolume[k] * Link[k].Travel_time / 60.0;

			PMT += Link[k].mode_MainVolume[m] * g_mode_type_vector[m].occ * Link[k].length;
			PHT += Link[k].mode_MainVolume[m] * g_mode_type_vector[m].occ * Link[k].Travel_time / 60.0;

			VHT_QVDF += MainVolume[k] * Link[k].QVDF_TT / 60.0;
			PHT_QVDF += Link[k].mode_MainVolume[m] * g_mode_type_vector[m].occ * Link[k].QVDF_TT / 60.0;


		}

        fprintf(link_performance_file,
            "%d,%d,%d,%d,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n",
            iteration_no, Link[k].link_id, Link[k].external_from_node_id, Link[k].external_to_node_id,
            MainVolume[k], Link[k].Ref_volume, Link[k].Base_volume, Link[k].Obs_volume, Link[k].Link_Capacity, IncomingDemand, DOC, Link[k].FreeTravelTime,
            Link[k].Travel_time, Link[k].VDF_Alpha, Link[k].VDF_Beta, Link[k].VDF_plf,
            Link[k].length / fmax(Link[k].Travel_time / 60.0, 0.001),
            Link[k].length / fmax(Link[k].Travel_time / 60.0, 0.001) * 1.609,
            Link[k].Travel_time - Link[k].FreeTravelTime);

		fprintf(link_performance_file, "%2lf,%2lf,%2lf,%2lf,%2lf,%2lf,", VMT, VHT, PMT, PHT, VHT_QVDF, PHT_QVDF);

		fprintf(link_performance_file, "\"%s\",",
			Link[k].geometry.c_str());

		for (int m = 1; m <= number_of_modes; m++)
			fprintf(link_performance_file, "%2lf,", Link[k].mode_MainVolume[m]);

		fprintf(link_performance_file, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,", P, t0, t2, t3, vt2, vt2 * 1.609, mu, Q_gamma,
			Link[k].free_speed, Link[k].Cutoff_Speed, Link[k].free_speed*1.609, Link[k].Cutoff_Speed * 1.609,
			congestion_ref_speed, avg_queue_speed, avg_QVDF_period_speed, 
			congestion_ref_speed * 1.609, avg_queue_speed * 1.609, avg_QVDF_period_speed * 1.609,
			Link[k].QVDF_TT, Severe_Congestion_P);
		for (int t = demand_period_starting_hours * 60; t < demand_period_ending_hours * 60; t += 5)
		{
			int t_interval = t / 5;
			double speed = model_speed[t_interval];
			fprintf(link_performance_file, "%.3f,", speed);
		}



		fprintf(link_performance_file, "\n");
	}

	OutputODPerformance("od_performance.csv");

	if (shortest_path_log_flag)
	{
		OutputRouteDetails("route_assignment.csv");
	}
	free(MainVolume);
	free(SubVolume);
	free(SDVolume);

	Free_3D((void***)MDMinPathPredLink, number_of_modes, no_zones, no_nodes);

	Close();

	fclose(link_performance_file);
	fclose(logfile);  // Close the log file when you're done with it.
	fclose(summary_log_file);
	return 0;
}

static void Init(int int_number_of_modes, int input_no_zones)
{
	// tuiInit(tuiFileName);
	InitLinks();
	baseODDemand_loaded_flag = Read_ODflow(&TotalODflow, &int_number_of_modes, &input_no_zones);

	if (baseODDemand_loaded_flag == 0)
	{
		// reset 
		for (int k = 1; k <= number_of_links; k++)
		{

			Link[k].Base_volume = 0;

			for (int m = 1; m <= int_number_of_modes; m++)
			{
				Link[k].mode_Base_Volume[m] = 0;
			}
		}
	}
}

static void Close()
{
	StatusMessage("General", "Closing all modules");
	// tuiClose(tuiFileName);
	CloseLinks();
	CloseODflow();

}


static void CloseODflow(void)
{

	free(TotalOFlow);

	Free_3D((void***)MDODflow, number_of_modes, no_zones, no_zones);
	Free_3D((void***)MDDiffODflow, number_of_modes, no_zones, no_zones);
	Free_3D((void***)MDRouteCost, number_of_modes, no_zones, no_zones);

}






void ReadLinks()
{
	CDTACSVParser parser_link;
	std::vector<CLink> links;

	std::string line;

	if (parser_link.OpenCSVFile("link.csv", true))
	{
		float total_base_link_volume = 0;
		int line_no = 0;

		int k = 1;  // link start from 1

		while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read
			// field headers.
		{
			// CLink link;
			int lanes = 1;
			float capacity = 0;
			float free_speed = 10;

			parser_link.GetValueByFieldName("lanes", lanes);
			parser_link.GetValueByFieldName("capacity", capacity);
			parser_link.GetValueByFieldName("free_speed", free_speed);

			if (g_metric_system_flag == 1)
				free_speed = free_speed / 1.609;

			if (lanes <= 0 || capacity < 0.0001 || free_speed < 0.0001)
				continue; 


			Link[k].setup(number_of_modes);
			std::string value;



			// Read link_id
			parser_link.GetValueByFieldName("from_node_id", Link[k].external_from_node_id);
			parser_link.GetValueByFieldName("to_node_id", Link[k].external_to_node_id);
			parser_link.GetValueByFieldName("link_id", Link[k].link_id);
			parser_link.GetValueByFieldName("link_type", Link[k].link_type);



			Link[k].internal_from_node_id = Link[k].external_from_node_id;
			Link[k].internal_to_node_id = Link[k].external_to_node_id;

			if (g_map_external_node_id_2_node_seq_no.find(Link[k].external_from_node_id) !=
				g_map_external_node_id_2_node_seq_no.end())
			{
				Link[k].internal_from_node_id =
				g_map_external_node_id_2_node_seq_no[Link[k].external_from_node_id];
			}
			else
			{
				printf("Error in from_node_id =%d for link_id = %d\n", Link[k].external_from_node_id, Link[k].link_id);

				continue; 
			}

			if (Link[k].internal_from_node_id == 0)
			{
				printf("Error in Link[k].internal_from_node_id\n"); 
			}
			if (g_map_external_node_id_2_node_seq_no.find(Link[k].external_to_node_id) !=
				g_map_external_node_id_2_node_seq_no.end())
			{
				Link[k].internal_to_node_id =
				g_map_external_node_id_2_node_seq_no[Link[k].external_to_node_id];
			}
			else
			{
				printf("Error in to_node_id =%d for link_id = %d\n", Link[k].external_to_node_id, Link[k].link_id);
				continue;
			}


			parser_link.GetValueByFieldName("length", Link[k].length);

			if (g_metric_system_flag == 1)
				Link[k].length = Link[k].length / 1609;

			parser_link.GetValueByFieldName("ref_volume", Link[k].Ref_volume);

			if (g_ODME_mode == 1)
			{
				parser_link.GetValueByFieldName("obs_volume", Link[k].Obs_volume);

			}
			if (g_base_demand_mode == 1)
			{

				parser_link.GetValueByFieldName("base_volume", Link[k].Base_volume);
				total_base_link_volume += Link[k].Base_volume;

				if (number_of_modes == 1)  // single mode 
					Link[k].mode_Base_Volume[1] = Link[k].Base_volume;

				for (int m = 1; m <= number_of_modes; m++)
				{
					std::string field_name = "base_vol_" + std::string(g_mode_type_vector[m].mode_type);
					parser_link.GetValueByFieldName(field_name.c_str(), Link[k].mode_Base_Volume[m]);

				}
			}



			Link[k].lanes = lanes;
			Link[k].Lane_Capacity = capacity;
			Link[k].Link_Capacity = lanes * capacity;


			parser_link.GetValueByFieldName("allowed_use", Link[k].allowed_uses);

			if (g_tap_log_file == 1)
			{
				fprintf(logfile, "link %d->%d, node_seq_no %d->%d\n", Link[k].external_from_node_id, Link[k].external_to_node_id, Link[k].internal_from_node_id, Link[k].internal_to_node_id);

			}


			for (int m = 1; m <= number_of_modes; m++)
			{
				Link[k].mode_allowed_use[m] = 1;
				Link[k].mode_MainVolume[m] = 0;
				Link[k].mode_SubVolume[m] = 0;
				Link[k].mode_Toll[m] = 0;
				Link[k].mode_AdditionalCost[m] = 0;

			}

			if (Link[k].allowed_uses.size() > 0 && Link[k].allowed_uses != "all")
			{
				for (int m = 1; m <= number_of_modes; m++)
				{
					if (Link[k].allowed_uses.find(g_mode_type_vector[m].mode_type) != std::string::npos)  // otherwise, only an agent type is listed in this "allowed_uses", then this agent type is allowed to travel on this link
					{
						Link[k].mode_allowed_use[m] = 1;  // found 
					}
					else
					{
						Link[k].mode_allowed_use[m] = 0;
					}
				}


			}



			// Read internal_from_node_id

			// Read length

			// Read capacity

			Link[k].FreeTravelTime = Link[k].length / free_speed * 60.0;

			parser_link.GetValueByFieldName("VDF_alpha", Link[k].VDF_Alpha);
			parser_link.GetValueByFieldName("VDF_beta", Link[k].VDF_Beta);
			parser_link.GetValueByFieldName("VDF_plf", Link[k].VDF_plf);


			for (int m = 1; m <= number_of_modes; m++)
			{
				char CSV_field_name[50];
				sprintf(CSV_field_name, "toll_%s", g_mode_type_vector[m].mode_type.c_str());
				parser_link.GetValueByFieldName(CSV_field_name, Link[k].mode_Toll[m], false, false);

				Link[k].mode_AdditionalCost[m] = Link[k].mode_Toll[m] / g_mode_type_vector[m].vot * 60.0;
			}



			if (capacity > 0)
				Link[k].BoverC = Link[k].VDF_Alpha / pow(capacity * lanes, Link[k].VDF_Beta);
			else
				Link[k].BoverC = 0;


			parser_link.GetValueByFieldName("VDF_cp", Link[k].Q_cp);
			parser_link.GetValueByFieldName("VDF_cd", Link[k].Q_cd);
			parser_link.GetValueByFieldName("VDF_n", Link[k].Q_n);
			parser_link.GetValueByFieldName("VDF_s", Link[k].Q_s);

			Link[k].free_speed = free_speed;
			Link[k].Cutoff_Speed = free_speed * 0.75;  // use 0.75 as default ratio, when free_speed = 70 mph, Cutoff_Speed = 52.8 mph in I-10 data set
			parser_link.GetValueByFieldName("geometry", Link[k].geometry, false);
			k++;
		}

		printf("total_base_link_volume = %f\n", total_base_link_volume);

		if (total_base_link_volume > 0)
			baselinkvolume_loaded_flag = 1;
		parser_link.CloseCSVFile();
	}



}

static void InitLinkPointers(char* LinksFileName)
{
	int k, Node, internal_from_node_id;
	// Node is the internal node id
	FirstLinkFrom = (int*)Alloc_1D(no_nodes, sizeof(int));
	LastLinkFrom = (int*)Alloc_1D(no_nodes, sizeof(int));

	FirstLinkFrom[1] = 1;
	Node = 1;

	for (k = 1; k <= number_of_links; k++)
	{
		internal_from_node_id = Link[k].internal_from_node_id;
		if (internal_from_node_id == Node)
			continue;

		else if (internal_from_node_id >= Node + 1)
		{
			LastLinkFrom[Node] = k - 1;
			Node = internal_from_node_id;
			FirstLinkFrom[Node] = k;
		}

		else if (internal_from_node_id < Node)
		{
			// ExitMessage("Sort error in link file '%s': a link from node %d was found after "
			//	"a link from node %d \n", LinksFileName, internal_from_node_id, Node);
		}
		else if (internal_from_node_id > Node + 1)
		{
			// InputWarning("link file '%s' has no links out from "
			//	"nodes %d through %d. \n", LinksFileName, Node + 1, internal_from_node_id - 1);
			LastLinkFrom[Node] = k - 1;
			for (Node++; Node < internal_from_node_id; Node++)
			{
				FirstLinkFrom[Node] = 0;
				LastLinkFrom[Node] = -1;
			}
			FirstLinkFrom[Node] = k; /* Node equals internal_from_node_id now. */
		}
	}

	if (Node == no_nodes)
	{
		LastLinkFrom[Node] = number_of_links; /* Now Node equals no_nodes in any case */
	}
	else
	{
		// InputWarning("link file '%s' has no links out from "
		//	"nodes %d through %d. \n", LinksFileName, Node + 1, no_nodes);
		LastLinkFrom[Node] = k - 1;
		for (Node++; Node <= no_nodes; Node++)
		{
			FirstLinkFrom[Node] = 0;
			LastLinkFrom[Node] = -1;
		}
	}

	if (g_tap_log_file == 1)
	{
		for (Node = 1; Node <= no_nodes; Node++)
		{ 
			fprintf(logfile, "node_id = %d, FirstLinkFrom = %d, LastLinkFrom = %d \n", g_map_node_seq_no_2_external_node_id[Node], FirstLinkFrom[Node], FirstLinkFrom[Node]);

		}
	}

}

void FindLinksTo(void)
{
	int Node, k;

	LinksTo = (sorted_list*)Alloc_1D(no_nodes, sizeof(sorted_list*));
	for (Node = 1; Node <= no_nodes; Node++)
		LinksTo[Node] = NULL;

	for (k = 1; k <= number_of_links; k++)
		AddToSortedList(k, &(LinksTo[Link[k].internal_to_node_id]));
}
void InitLinks()
{
	char LinksFileName[100] = "link.csv";
	char* InterFlowFileName;
	FILE* LinksFile;
	int k;

	Link = (struct link_record*)Alloc_1D(number_of_links, sizeof(struct link_record));
	ReadLinks();
	FindLinksTo();
	InitLinkPointers(LinksFileName);
	UpdateLinkAdditionalCost();
}



void StatusMessage(const char* group, const char* format, ...)
{
	double new_time;
}

int Read_ODtable(double*** ODtable, double*** DiffODtable, int no_zones)
{
	char ch, Coloumn[2], Semicoloumn[2]; /* Reserve room for the '\0' in the fscanf_s. */
	int Orig, Dest, NewOrig, NewDest;
	double Value;

	// current OD demand 

	for (int m = 1; m <= number_of_modes; m++)
	{
		FILE* file;
		fopen_s(&file, g_mode_type_vector[m].demand_file.c_str(), "r");
		printf("read demand file %s\n", g_mode_type_vector[m].demand_file.c_str());
		fprintf(summary_log_file, "read demand file %s\n", g_mode_type_vector[m].demand_file.c_str());

		if (file == NULL)
		{
			if (g_mode_type_vector[m].demand_file.length() > 0)
			{
			printf("Failed to open demand file %s\n", g_mode_type_vector[m].demand_file.c_str());
			}
			return 0;
		}

		int o_zone_id, d_zone_id;
		double volume;

		// Skip the header line
		char header[100];
		if (fgets(header, sizeof(header), file) == NULL)
		{
			printf("Failed to read header\n");
			return 0;
		}

		int line_count = 0;
		double total_volume = 0;
		// Read the data
		int result;
		while ((result = fscanf(file, "%d,%d,%lf", &o_zone_id, &d_zone_id, &volume)) != EOF)
		{
			if (o_zone_id > no_zones)
			{
				printf("o_zone_id: %d, d_zone_id: %d, volume: %.4lf\n", o_zone_id, d_zone_id,
					volume);
				printf("Error o_zone_id %d  > # of zones %d \n", o_zone_id, no_zones);
				break;
			}
			if (d_zone_id > no_zones)
			{
				printf("o_zone_id: %d, d_zone_id: %d, volume: %.4lf\n", o_zone_id, d_zone_id,
					volume);
				printf("Error d_zone_id %d  > # of zones %d \n", d_zone_id, no_zones);
				break;
			}


			if (result == 3)  // we have read all the 3 values correctly
			{
				if (line_count <= 3)
				{
					printf("o_zone_id: %d, d_zone_id: %d, volume: %.4lf\n", o_zone_id, d_zone_id,
						volume);

				}
				ODtable[m][o_zone_id][d_zone_id] = volume;
				DiffODtable[m][o_zone_id][d_zone_id] = volume; // by default in case there is no baseline demand being specified 
				total_volume += volume;
				line_count++;
			}
			else
			{
				printf("Error reading line %d\n", line_count);
				break;
			}
		}

		printf(" mode type = %s, total_volume = %f\n", g_mode_type_vector[m].mode_type.c_str(), total_volume);
		fprintf(summary_log_file, " mode type = %s, total_volume = %f\n", g_mode_type_vector[m].mode_type.c_str(), total_volume);

		fclose(file);

	}

	// baseline OD demand --. diff OD  = current OD - baseline ODdemand 

	if (g_base_demand_mode == 0 || baselinkvolume_loaded_flag == 0)
		return 0; // skip 

	for (int m = 1; m <= number_of_modes; m++)
	{
		FILE* file;
		std::string original_filename = g_mode_type_vector[m].demand_file;
		std::string modified_filename;

		// Check if the original filename ends with ".csv" and replace it with "_base.csv"
		size_t pos = original_filename.find(".csv");
		if (pos != std::string::npos) {
			modified_filename = original_filename.substr(0, pos) + "_base.csv";
		}
		else {
			// If the file does not have ".csv", append "_base"
			modified_filename = original_filename + "_base";
		}

		fopen_s(&file, modified_filename.c_str(), "r");
		printf("read demand file %s\n", modified_filename.c_str());

		if (file == NULL)
		{
			// by default, we can skip this requirement, but if we load baseline link volume we should have base OD demand for consistency 
			break;
		}

		int o_zone_id, d_zone_id;
		double volume;

		// Skip the header line
		char header[100];
		if (fgets(header, sizeof(header), file) == NULL)
		{
			printf("Failed to read header\n");
			return 0;
		}

		int line_count = 0;
		double total_volume_diff = 0;
		// Read the data
		int result;
		while ((result = fscanf(file, "%d,%d,%lf", &o_zone_id, &d_zone_id, &volume)) != EOF)
		{
			if (result == 3)  // we have read all the 3 values correctly
			{
				if (line_count <= 3)
				{
					printf("o_zone_id: %d, d_zone_id: %d, volume: %.4lf\n", o_zone_id, d_zone_id,
						volume);

				}
				DiffODtable[m][o_zone_id][d_zone_id] = ODtable[m][o_zone_id][d_zone_id] - volume;  // diff OD demand  = current demand - base demand 
				total_volume_diff += ODtable[m][o_zone_id][d_zone_id] - volume;
				line_count++;
			}
			else
			{
				printf("Error reading line %d\n", line_count);
				break;
			}
		}

		printf(" mode type = %s, total_volume_diff = %f\n", g_mode_type_vector[m].mode_type.c_str(), total_volume_diff);

		fclose(file);

	}
	return 1;
}

double Link_Travel_Time(int k, double* Volume)
{
	double IncomingDemand = Volume[k] / fmax(0.01, Link[k].lanes) / fmax(0.001, demand_period_ending_hours - demand_period_starting_hours) / fmax(0.0001, Link[k].VDF_plf);

	Link[k].Travel_time =
		Link[k].FreeTravelTime * (1.0 + Link[k].VDF_Alpha * (pow(IncomingDemand / fmax(0.1, Link[k].Link_Capacity), Link[k].VDF_Beta)));

	if (g_ODME_mode == 1 && Link[k].Obs_volume >= 0)
	{
		Link[k].Travel_time += 2 * g_ODME_link_volume_penalty * (Volume[k] - Link[k].Obs_volume);

	}

	if (g_ODME_mode == 1 && g_ODME_obs_VMT > 1 && Link[k].link_type >= 1 && g_System_VMT >= 1)  //p[hsical links
	{
		Link[k].Travel_time += 2 * g_ODME_VMT_penalty / number_of_links * (g_System_VMT - g_ODME_obs_VMT) * Link[k].length;



	}
	if (Link[k].Travel_time < 0)
		Link[k].Travel_time = 0;
	Link[k].BPR_TT = Link[k].Travel_time;

	return (Link[k].Travel_time);
}


double Link_QueueVDF(int k, double Volume, double& IncomingDemand, double& DOC, double& P, double& t0, double& t2, double& t3, double& vt2, double& Q_mu, double& Q_gamma, double& congestion_ref_speed,
	double& avg_queue_speed, double& avg_QVDF_period_speed, double& Severe_Congestion_P, double model_speed[300])
{

	IncomingDemand = Volume / fmax(0.01, Link[k].lanes) / fmax(0.001, demand_period_ending_hours - demand_period_starting_hours) / fmax(0.0001, Link[k].VDF_plf);
	DOC = IncomingDemand / fmax(0.1, Link[k].Lane_Capacity);

	double Travel_time =
		Link[k].FreeTravelTime * (1.0 + Link[k].VDF_Alpha * (pow(DOC, Link[k].VDF_Beta)));

	congestion_ref_speed = Link[k].Cutoff_Speed;
	if (DOC < 1)
		congestion_ref_speed = (1 - DOC) * Link[k].free_speed + DOC * Link[k].Cutoff_Speed;


	//step 3.2 calculate speed from VDF based on D/C ratio
	avg_queue_speed = congestion_ref_speed / (1.0 + Link[k].VDF_Alpha * pow(DOC, Link[k].VDF_Beta));


	P = Link[k].Q_cd * pow(DOC, Link[k].Q_n);  // applifed for both uncongested and congested conditions

	double H = demand_period_ending_hours - demand_period_starting_hours;

	if (P > H)
		avg_QVDF_period_speed = avg_queue_speed;
	else
		avg_QVDF_period_speed = P / H * avg_queue_speed + (1.0 - P / H) * (congestion_ref_speed + Link[k].free_speed) / 2.0;


	Link[k].QVDF_TT = Link[k].length / fmax(0.1, avg_QVDF_period_speed) * 60.0;






	double base = Link[k].Q_cp * pow(P, Link[k].Q_s) + 1.0;
	vt2 = Link[k].Cutoff_Speed / fmax(0.001, base);

	t2 = (demand_period_starting_hours + demand_period_ending_hours) / 2.0;
	t0 = t2 - 0.5 * P;
	t3 = t2 + 0.5 * P;

	Q_mu = std::min(Link[k].Lane_Capacity, IncomingDemand / std::max(0.01, P));

	//use  as the lower speed compared to 8/15 values for the congested states
	double RTT = Link[k].length / fmax(0.01, congestion_ref_speed);
	double wt2 = Link[k].length / vt2 - RTT; // in hour


	//step 5 compute gamma parameter is controlled by the maximum queue
	Q_gamma = wt2 * 64 * Q_mu / pow(P, 4);  // because q_tw = w*mu =1/4 * gamma (P/2)^4, => 1/vt2 * mu = 1/4 * gamma  * (P/2)^4


	double td_w = 0;
	//step scan the entire analysis period
	Severe_Congestion_P = 0;


	for (int t_in_min = demand_period_starting_hours * 60; t_in_min <= demand_period_ending_hours * 60; t_in_min += 5)  // 5 min interval
	{
		int t_interval = t_in_min / 5;
		double t = t_in_min / 60.0;  // t in hour
		double td_queue = 0;
		double td_speed = 0;
		model_speed[t_interval] = Link[k].free_speed;

		if (t0 <= t && t <= t3)  // within congestion duration P
		{
			//1/4*gamma*(t-t0)^2(t-t3)^2
			td_queue = 0.25 * Q_gamma * pow((t - t0), 2) * pow((t - t3), 2);
			td_w = td_queue / fmax(0.001, Q_mu);
			//L/[(w(t)+RTT_in_hour]
			td_speed = Link[k].length / (td_w + RTT);
		}
		else if (t < t0) //outside
		{
			td_queue = 0;
			double factor = (t - demand_period_starting_hours) / fmax(0.001, t0 - demand_period_starting_hours);
			td_speed = (1 - factor) * Link[k].free_speed + factor * fmax(congestion_ref_speed, avg_queue_speed);
		}
		else  // t> t3
		{
			td_queue = 0;
			double factor = (t - t3) / fmax(0.001, demand_period_ending_hours - t3);
			td_speed = (1 - factor) * fmax(congestion_ref_speed, avg_queue_speed) + (factor)*Link[k].free_speed;
		}

		// dtalog.output() << "td_queue t" << t << " =  " << td_queue << ", speed =" << td_speed << '\n';
		// g_DTA_log_file << "td_queue t" << t << " =  " << td_queue << ", speed =" << td_speed << '\n';



		if (t_in_min <= 410)
		{
			int idebug = 1;
		}
		double td_flow = 0; // default: get_volume_from_speed(td_speed, vf, k_critical, s3_m);
		model_speed[t_interval] = td_speed;

		if (td_speed < Link[k].free_speed * 0.5)
			Severe_Congestion_P += 5.0 / 60;  // 5 min interval
	}

	return P;
}

double Link_Travel_Time_Integral(int k, double* Volume)
{
	double IncomingDemand = Volume[k] / fmax(0.001, demand_period_ending_hours - demand_period_starting_hours) / fmax(0.0001, Link[k].VDF_plf);
	double integral = 0;
	if (Link[k].VDF_Beta >= 0.0)
		integral += IncomingDemand + (Volume[k] * Link[k].FreeTravelTime *
			(1.0 + (Link[k].BoverC / (Link[k].VDF_Beta + 1)) * pow(IncomingDemand, Link[k].VDF_Beta + 1)));

	return integral;
}

double Link_Travel_Time_Der(int k, double* Volume)
{
	if (Link[k].VDF_Beta == 0.0)
		return 0.0;
	else
		return (Link[k].FreeTravelTime * Link[k].BoverC * Link[k].VDF_Beta *
			pow(Volume[k], (Link[k].VDF_Beta - 1)));
}

double AdditionalCost(int k, int m)
{
	double AddCost = 0;

	AddCost = Link[k].mode_Toll[m] / g_mode_type_vector[m].vot * 60.0;

	return AddCost;
}

double Link_GenCost(int k, double* Volume)
{
	return (Link[k].mode_AdditionalCost[1] + Link_Travel_Time(k, Volume));
}

double LinkCost_Integral(int k, double* Volume)
{
	return (Link[k].mode_AdditionalCost[1] * Volume[k] + Link_Travel_Time_Integral(k, Volume));
}

double Link_GenCostDer(int k, double* Volume)
{
	return (Link_Travel_Time_Der(k, Volume));
}

/* External functions */

void ClearVolume(double* VolumeArray)
{
	int k;
	for (k = 1; k <= number_of_links; k++)
		VolumeArray[k] = 0.0;
}

void VolumeDifference(double* Volume1, double* Volume2, double* Difference)
{  // SubVolume, MainVolume
	int k;
	for (k = 1; k <= number_of_links; k++)
	{
		Difference[k] = Volume1[k] - Volume2[k];

		for (int m = 1; m <= number_of_modes; m++)
		{
			Link[k].mode_SDVolume[m] = Link[k].mode_SubVolume[m] - Link[k].mode_MainVolume[m];

			if (fabs(Difference[k] - Link[k].mode_SDVolume[m]) > 0.01)
			{
				printf("");
			}
		}
	}
}

void UpdateVolume(double* MainVolume, double* SDVolume, double Lambda)
{
	int k;
	for (k = 1; k <= number_of_links; k++)
	{
		MainVolume[k] += Lambda * SDVolume[k];
	}


	for (int k = 1; k <= number_of_links; k++)
	{
		for (int m = 1; m <= number_of_modes; m++)
		{
			Link[k].mode_MainVolume[m] += Lambda * Link[k].mode_SDVolume[m];
		}
	}
}

void UpdateLinkAdditionalCost(void)
{
	int k;

	for (k = 1; k <= number_of_links; k++)
		for (int m = 1; m <= number_of_modes; m++)
			AdditionalCost(k, m);
}

double UpdateLinkCost(double* MainVolume)
{
	int k;
	double system_wide_travel_time = 0;

	for (k = 1; k <= number_of_links; k++)
	{
		Link[k].Travel_time = Link_Travel_Time(k, MainVolume);

		Link[k].GenCost = Link_GenCost(k, MainVolume);
		system_wide_travel_time += (MainVolume[k] * Link[k].Travel_time);
	}

	return system_wide_travel_time;
}

void UpdateLinkCostDer(double* MainVolume)
{
	int k;

	for (k = 1; k <= number_of_links; k++)
	{
		Link[k].GenCostDer = Link_GenCostDer(k, MainVolume);
	}
}

void GetLinkTravelTimes(double* Volume, double* TravelTime)
{
	int k;

	for (k = 1; k <= number_of_links; k++)
	{
		TravelTime[k] = Link_Travel_Time(k, Volume);
	}
}

double TotalLinkCost(double* Volume)
{
	int k;
	double Sum = 0;

	for (k = 1; k <= number_of_links; k++)
		Sum += Link[k].GenCost * Volume[k];
	return Sum;
}

double OFscale = 1.0;

double OF_Links(double* MainVolume)
{
	int k;
	double Sum = 0;

	for (k = 1; k <= number_of_links; k++)
		Sum += LinkCost_Integral(k, MainVolume);

	return Sum / OFscale;
}

double OF_LinksDirectionalDerivative(double* MainVolume, double* SDVolume, double Lambda)
{
	//
	//purpose:
	//    this function calculates the directional derivative of the objective function with respect to the step size lambda.in optimization, the directional derivative indicates how the objective function changes as you move in a specific direction.
	//
	//        parameters :
	//        mainvolume : current flow volumes on each link.
	//        sdvolume : search direction volumes(difference between aon assignment and current flows).
	//        lambda : step size parameter indicating how far to move along the search direction.

		//Return Normalized Directional Derivative :
		//Normalize the sum by OFscale to get the directional derivative :
		//DirectionalDerivative
		//    =
		//    LinkCostSum
		//    OFscale
		//    DirectionalDerivative =
		//    OFscale
		//    LinkCostSum
		//    ?

		//    Mathematical Interpretation :
		//The directional derivative represents how the total system cost changes as you adjust the flow along the search direction by Lambda.In the context of traffic assignment, it helps determine whether increasing or decreasing Lambda will reduce the overall congestion and travel time.

	int k;
	double* Volume;
	double LinkCostSum = 0;

	Volume = (double*)Alloc_1D(number_of_links, sizeof(double));

	for (k = 1; k <= number_of_links; k++)
	{
		Volume[k] = MainVolume[k] + Lambda * SDVolume[k];
	}
	for (k = 1; k <= number_of_links; k++)
	{
		LinkCostSum += Link_GenCost(k, Volume) * SDVolume[k];

		//   LinkCostSum += Link_Travel_Time_Integral(k, Volume) * SDVolume[k];

	}

	free(Volume);
	return LinkCostSum / OFscale;
}

double Sum_ODtable(double*** ODtable, double* total_o_table, int no_zones)
{
	int Orig, Dest;
	double sum = 0.0;

	for (Orig = 1; Orig <= no_zones; Orig++)
		total_o_table[Orig] = 0;

	for (int m = 1; m <= number_of_modes; m++)
		for (Orig = 1; Orig <= no_zones; Orig++)
			for (Dest = 1; Dest <= no_zones; Dest++)
			{
				total_o_table[Orig] += ODtable[m][Orig][Dest];
				sum += ODtable[m][Orig][Dest];
			}


	return (sum);
}

int Read_ODflow(double* TotalODflow, int* number_of_modes, int* no_zones)
{
	FILE* ODflowFile;
	double RealTotal, InputTotal;

	MDODflow = (double***)Alloc_3D(*number_of_modes, *no_zones, *no_zones, sizeof(double));
	TotalOFlow = (double*)Alloc_1D(*no_zones, sizeof(double));

	MDDiffODflow = (double***)Alloc_3D(*number_of_modes, *no_zones, *no_zones, sizeof(double));

	MDRouteCost = (double***)Alloc_3D(*number_of_modes, *no_zones, *no_zones, sizeof(double));

	int with_basedemand = Read_ODtable(MDODflow, MDDiffODflow, *no_zones);

	RealTotal = (double)Sum_ODtable(MDODflow, TotalOFlow, *no_zones);


	zone_outbound_link_size = (int*)Alloc_1D(*no_zones, sizeof(int));

	for (int n = 1; n <= *no_zones; n++)
	{
		zone_outbound_link_size[n] = 0;
	}
	for (int k = 1; k <= number_of_links; k++)
	{
		int from_node_id = Link[k].external_from_node_id;
		if (from_node_id <= *no_zones)
		{
			zone_outbound_link_size[from_node_id] += 1;  // from_node_id is the zone_id; 
		}
	}

	

	// checking 
	int total_infeasible_outbound_zones = 0; 
	float total_infeasible_outbound_zone_demand = 0;
	float total_zone_demand = 0;
	for (int z = 1; z < *no_zones; z++)
	{
		if (zone_outbound_link_size[z] == 0 && TotalOFlow[z]>0.01)
		{
			printf("Error: There is no outbound link from zone %d with positive demand %f\n", z, TotalOFlow[z]);
			total_infeasible_outbound_zones++; 
			total_infeasible_outbound_zone_demand += TotalOFlow[z];
		}
		total_zone_demand += TotalOFlow[z];
	}
	printf("Error: %d zones have no outbound link with positive demand %f, = %f percentage of total demand\n", total_infeasible_outbound_zones,
		total_infeasible_outbound_zone_demand, total_infeasible_outbound_zone_demand*100/fmax(0.01, total_zone_demand));


	*TotalODflow = (double)RealTotal;

	return with_basedemand;

}

#include <cstdarg>
#include <cstdio>
#include <cstdlib>

void ExitMessage(const char* format, ...) {
    va_list ap;

    va_start(ap, format); // Initialize the va_list with the format parameter
    vfprintf(stderr, format, ap); // Print the formatted error message to stderr
    va_end(ap); // Clean up the va_list

    exit(EXIT_FAILURE); // Terminate the program
}


void CloseLinks(void)
{
	int Node;
	free(zone_outbound_link_size);
	free(Link);
	free(FirstLinkFrom);
	free(LastLinkFrom);
	for (Node = 1; Node <= no_nodes; Node++)
		FreeSortedList(LinksTo[Node]);
	free(LinksTo);
}


double LinksSDLineSearch(double* MainVolume, double* SDVolume) {
	int n;
	double lambdaleft = 0, lambdaright = 1, lambda = 0.5;
	double grad;

	// Initial check at lambda = 0
	grad = OF_LinksDirectionalDerivative(MainVolume, SDVolume, 0.0);
	StatusMessage("Line search step", "0.0");
	StatusMessage("Line search grad", "%lf", grad);
	if (grad >= 0) {
		LastLambda = 0.0;
		return 0.0;
	}

	// Check at lambda = 1
	grad = OF_LinksDirectionalDerivative(MainVolume, SDVolume, 1.0);
	StatusMessage("Line search step", "1.0");
	StatusMessage("Line search grad", "%lf", grad);
	if (grad <= 0) {
		LastLambda = 1.0;
		return 1.0;
	}

	// Bisection method for line search within [0, 1]
	for (n = 1; n <= MinLineSearchIterations; n++) {
		grad = OF_LinksDirectionalDerivative(MainVolume, SDVolume, lambda);

		if (grad <= 0.0) {
			lambdaleft = lambda;
		}
		else {
			lambdaright = lambda;
		}

		lambda = 0.5 * (lambdaleft + lambdaright);
	}

	// Additional iterations to ensure convergence, if necessary
	for (; lambdaleft == 0 && n <= MAX_NO_BISECTITERATION; n++) {
		grad = OF_LinksDirectionalDerivative(MainVolume, SDVolume, lambda);

		if (grad <= 0.0) {
			lambdaleft = lambda;
		}
		else {
			lambdaright = lambda;
		}

		lambda = 0.5 * (lambdaleft + lambdaright);
	}

	ActualIterations = n - 1;
	LastLambda = lambdaleft;
	return lambdaleft;
}

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>

using namespace std;

// Define constants and structures
const int T = 100; // Total simulation time
const int MAX_LINKS = 100; // Maximum number of links

struct Vehicle {
	int id;
	vector<int> path; // Sequence of links
	int currentLink;
	int TA; // Arrival time at the current link
	int TD; // Departure time from the current link
	int passengers; // Number of passengers in the vehicle
	int capacity; // Maximum passenger capacity
};

struct Link {
	int id;
	int FFTT; // Free-flow travel time
	double capacity; // Current available capacity
	queue<int> waitingQueue; // Vehicles waiting to enter the link
	vector<int> A; // Cumulative arrival counts
	vector<int> D; // Cumulative departure counts
};

// Function to simulate traffic
//void simulateTraffic(vector<Vehicle>& vehicles, vector<Link>& links, int totalTime) {
//	for (int t = 0; t < totalTime; ++t) { // Time loop
//		for (auto& link : links) { // Link loop
//			// Update link capacity based on arrivals and departures
//			if (t > 0) {
//				link.capacity -= (link.A[t - 1] - link.D[t - 1]);
//			}
//
//			for (auto& vehicle : vehicles) { // Vehicle loop
//				if (vehicle.currentLink == link.id && vehicle.TD == t) {
//					// Check if vehicle can move
//					if (link.capacity >= 1.0) {
//						// Move vehicle to the next link in its path
//						if (vehicle.currentLink + 1 < vehicle.path.size()) {
//							int nextLinkId = vehicle.path[vehicle.currentLink + 1];
//							auto& nextLink = links[nextLinkId];
//
//							vehicle.TA = t;
//							vehicle.TD = t + nextLink.FFTT;
//
//							// Update link capacities
//							link.capacity -= 1.0; // Assuming 1 PCE for simplicity
//							nextLink.capacity -= 1.0;
//
//							// Update cumulative flow counts
//							link.D[t]++;
//							nextLink.A[t]++;
//
//							// Update vehicle state
//							vehicle.currentLink++;
//						}
//					}
//					else {
//						// Vehicle waits if no capacity is available
//						vehicle.TD = t + 1;
//					}
//				}
//			}
//		}
//	}
//}

//int Simu_main() {
//	// Create sample vehicles
//	vector<Vehicle> vehicles = {
//		{0, {0, 1, 2}, 0, 0, 0, 4, 4}, // id, path, currentLink, TA, TD, passengers, capacity
//		{1, {1, 2, 3}, 0, 0, 0, 3, 4}
//	};
//
//	// Create sample links
//	vector<Link> links(MAX_LINKS);
//	for (int i = 0; i < MAX_LINKS; ++i) {
//		links[i] = { i, 5, 10.0, queue<int>(), vector<int>(T, 0), vector<int>(T, 0) };
//	}
//
//	// Run the simulation
//	simulateTraffic(vehicles, links, T);
//
//	// Print results
//	for (const auto& link : links) {
//		cout << "Link " << link.id << ": " << endl;
//		cout << "Arrival Counts: ";
//		for (int t = 0; t < T; ++t) cout << link.A[t] << " ";
//		cout << endl;
//
//		cout << "Departure Counts: ";
//		for (int t = 0; t < T; ++t) cout << link.D[t] << " ";
//		cout << endl;
//	}
//
//	return 0;
//}
