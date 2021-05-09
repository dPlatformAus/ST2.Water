/*
	A program to analyse bonding changes in simulated water data provided by Peter Pool
	Duncan Andrews (thanks also to Toby Hudson, Peter Harrowell and Peter Pool)
*/
#include <algorithm>
#include <vector>

#include <sys/types.h>

#include <time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <dirent.h>
#define NAMLEN(dirent) strlen((dirent)->d_name)

using namespace std;

#define DIM 3 // the number of dimensions
#define N 1728 // the total number of water molecules in a frame
#define NE 4000 // max number of events per tween
#define NC 150 // max number of clusters per tween
#define NCM 1000 // max number of molecules per cluster

#define NBONDS   8  // max number of bonds per molecule

char str_time_stamp[15];
char input_directory[15]="input/";
char output_directory[15]="output/";

double cutoff_radius_sq, cutoff_radius = 0.35; // the maximum radius within which oxygens could be considered bonded
double cutoff_hdist_factor = 0.1; // the fraction of the bond distance within which a hydrogen could be considered to be "in the bond"
int include_secondary = 1;
int timestamp_files = 0;
int read_velocities = 0;
int alternate_header = 0;
int skip_files = 0;

int errno=0;

#include "misc.h"
#include "water.h"
#include "frame.h"
#include "tween.h"

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void load_analysis_configuration()
{
	/*D: Tries to load main configuration settings.*/
   FILE *fp;

   if ( (fp=fopen("config/analysis_config.txt", "r")) != NULL) {
		fscanf(fp, "cutoff_radius=%lf\n", &cutoff_radius);
		fscanf(fp, "cutoff_hdist_factor=%lf\n", &cutoff_hdist_factor);
		fscanf(fp, "include_secondary=%d\n", &include_secondary);
		fclose(fp);
	} else {
		fprintf(stderr, "Error: failed to load config file\n");
		//system("PAUSE");
		exit(1);
   }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector <std::string> read_directory( const std::string& path = std::string() )
  {
  std::vector <std::string> result;
  dirent* de;
  DIR* dp;
  errno = 0;
  dp = opendir( path.empty() ? "." : path.c_str() );
  if (dp)
    {
    while (true)
      {
      errno = 0;
      de = readdir( dp );
      if (de == NULL) break;
      if ( !( (de->d_name[0] == '.' && NAMLEN(de) == 1) || (de->d_name[1] == '.' && NAMLEN(de) == 2) ) ) result.push_back( std::string( de->d_name ) );
      }
    closedir( dp );
    std::sort( result.begin(), result.end() );
    }
  return result;
  }
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void count_input_files(int &num_files)
{
	DIR *dp;
	struct dirent *ep;
	num_files=0;
	dp = opendir (input_directory);
	if (dp != NULL) {
		while ( (ep=readdir(dp)) != NULL){
			if ( !( (ep->d_name[0] == '.' && NAMLEN(ep) == 1) || (ep->d_name[1] == '.' && NAMLEN(ep) == 2) ) ) num_files++;
		}
		(void) closedir (dp);
	} else {
         perror ("Couldn't open the input directory\n");
		//system("PAUSE");
		exit(1);
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void load_input_files(char str_folder_files[][50])
{
	int i=0;
	DIR *dp;
	struct dirent *ep;
	dp = opendir (input_directory);
	if (dp != NULL) {
		while ( (ep=readdir(dp)) != NULL){
			if ( !( (ep->d_name[0] == '.' && NAMLEN(ep) == 1) || (ep->d_name[1] == '.' && NAMLEN(ep) == 2) ) ) {
				sprintf(str_folder_files[i], "%s", ep->d_name);
				i++;
			}
		}
		closedir (dp);
	} else {
         perror ("Couldn't open the input directory\n");
		//system("PAUSE");
		exit(1);
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void prune_input_files(int skip_files, int num_files, const std::vector<std::string>& str_folder_files, std::vector<std::string>& str_input_files)
{
	int i, j, skip_count;
	i=j=skip_count=0;
	str_input_files.push_back( str_folder_files[i] );
	//sprintf(str_input_files[j], "%s", str_folder_files[i]);
	for (i=1; i<num_files; i++){
		if (skip_count==skip_files){
			skip_count=-1;
			str_input_files.push_back( str_folder_files[i] );
			//sprintf(str_input_files[++j], "%s", str_folder_files[i]);
		}
		skip_count++;
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void compare_frames(Frame *prev_frame, Frame *current_frame, Tween *the_tween, int i)
{//compare frames ... record new bonds and lost bonds between each frame
	int j, k, l, m, n, exists;

	//fprintf(stderr, "\nComparing Frames--------:");
		the_tween->set_frames(i-1, i);
		for (j=0; j<N; j++){ //for each water molecule
			//find bonds in current_frame not in prev_frame (created)
			for (k=0; k<current_frame->molecules[j].num_bonds; k++){
				for (m=0; m<current_frame->molecules[j].bonds[k].max_rank_t; m++){
					exists=0;
					for (l=0; l<prev_frame->molecules[j].num_bonds; l++){
						if (current_frame->molecules[j].bonds[k].water_id == prev_frame->molecules[j].bonds[l].water_id){ //this water molecule is still bonded to the "other" water molecule this bond refers to
							for (n=0; n<prev_frame->molecules[j].bonds[l].max_rank_t; n++){
								if (current_frame->molecules[j].bonds[k].hydrogen_id[m] == prev_frame->molecules[j].bonds[l].hydrogen_id[n]){ //it is the same hydrogen in both versions of the bond
									exists=1;
								}
							}
						}
					}
					if (!exists){ // this bond did not exist in the previous frame, it has been created
						the_tween->add_event(1, j, current_frame->molecules[j].bonds[k].water_id, current_frame->molecules[j].bonds[k].hydrogen_id[m], m+1, current_frame->molecules[j].bonds[k].max_rank_t, current_frame->molecules[j].bonds[k].rank_o[m], current_frame->molecules[j].bonds[k].max_rank_o[m]);
					}
				}
			}
			//find bonds in prev_frame not in current_frame (broken)
			for (k=0; k<prev_frame->molecules[j].num_bonds; k++){
				for (m=0; m<prev_frame->molecules[j].bonds[k].max_rank_t; m++){
					exists=0;
					for (l=0; l<current_frame->molecules[j].num_bonds; l++){
						if (prev_frame->molecules[j].bonds[k].water_id == current_frame->molecules[j].bonds[l].water_id){ //this water molecule is still bonded to the "other" water molecule this bond refers to
							for (n=0; n<current_frame->molecules[j].bonds[l].max_rank_t; n++){
								if (prev_frame->molecules[j].bonds[k].hydrogen_id[m] == current_frame->molecules[j].bonds[l].hydrogen_id[n]){ //it is the same hydrogen in both versions of the bond
									exists=1;
								}
							}
						}
					}
					if (!exists){ // this does not exist in the next frame, it has been broken
						the_tween->add_event(0, j, prev_frame->molecules[j].bonds[k].water_id, prev_frame->molecules[j].bonds[k].hydrogen_id[m], m+1, prev_frame->molecules[j].bonds[k].max_rank_t, prev_frame->molecules[j].bonds[k].rank_o[m], prev_frame->molecules[j].bonds[k].max_rank_o[m]);
					}
				}
			}
		}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void identify_clusters(Frame *prev_frame, Frame *current_frame, Tween *the_tween)
{
	int j, k, l;

	//fprintf(stderr, "\nIdentifying Clusters------:");

		the_tween->get_clusters();
		//analyse clusters
		for (j=1; j<=the_tween->num_clusters; j++){
			for (k=0; k<N; k++){   //count molecules in cluster
				if (j==the_tween->cluster_ids[k]){ //molecule is assigned to this cluster
					the_tween->clusters[j].add_molecule(k, prev_frame->molecules[k].num_bonds, current_frame->molecules[k].num_bonds);
					//get num bonds broken and created																																			********************************************************************************************************** after this cluster still appears just no bonds broken or created
					for (l=0; l<the_tween->num_events; l++){																															///////here we have to be more carefull about secondaries to the same molecule *********************************************
						if (k==the_tween->events[l].water_id_1 && the_tween->events[l].hydrogen_id>0 && the_tween->events[l].rank_t==1){ //if the event is for a hydrogen on this molecule (don't count events twice) (******count secondaries? ... if so check not already counting a primary to the same molecule!!******)
																																																															// If event.rank_t >1 then there is another hydrogen participating in this bond ... so I think we don't count If event.rank_t >1
							if (the_tween->events[l].type==1){
								the_tween->clusters[j].bonds_created++;
							}else{
								the_tween->clusters[j].bonds_broken++;
							}
						}
					}
				}
			}
		}

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void output_statistics(FILE *the_file, Tween *the_tween, int i)
{
	int j;

	//fprintf(stderr, "\nOutput Statistics:");


        if (the_tween->num_clusters>NC) {
          fprintf(stderr, "error: numclusters NC exceeded on tween \n");
        }

		for (j=1; j<=the_tween->num_clusters; j++){
			fprintf(the_file, "%d, %d, %d, %d, %d, %d, %d", i-1, i, j, the_tween->clusters[j].num_molecules, the_tween->clusters[j].bonds_broken, the_tween->clusters[j].bonds_created, the_tween->clusters[j].valence_change);
			fprintf(the_file, ", %d, %d, %d, %d, %d, %d", the_tween->clusters[j].valence_gt4, the_tween->clusters[j].valence_lt4, the_tween->clusters[j].valence_change_gt4, the_tween->clusters[j].valence_change_lt4, the_tween->clusters[j].valence_remains_gt4, the_tween->clusters[j].valence_remains_lt4);
			fprintf(the_file, ", %d, %d, %d\n", the_tween->clusters[j].valence_defect_sum_start, the_tween->clusters[j].valence_defect_sum_end, the_tween->clusters[j].valence_defect_sum_full);

            if (the_tween->clusters[j].num_molecules>NCM) {
                fprintf(stderr, "error: num_molecules in a cluster exceeds NCM i=%d j=%d\n", i, j);
            }
		}

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void save_previous_defect_lifetime(int prev_valence, int consecutive_frames, int gt4[], int lt4[])
{
	if (consecutive_frames>0){
		//fprintf(stderr, "\n%d",consecutive_frames);
		if (prev_valence<4){
			lt4[consecutive_frames-1]++;
		}
		if (prev_valence>4){
			gt4[consecutive_frames-1]++;
		}
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void output_defect_lifetime_2(int num_frames)
{
	int i, j, frame_id, num_bonds;
	Output_File f_out, f_in;
	int *prev_valence=new int[N];//=4;
	int *consecutive_frames=new int[N];//=0;
	int *gt4=new int[num_frames];
	int *lt4=new int[num_frames];

	for (i=0; i<num_frames; i++) gt4[i]=lt4[i]=0;
	for (i=0; i<N; i++){ //for each water molecule
		prev_valence[i]=4;
		consecutive_frames[i]=0;
	}

	fprintf(stderr, "\nGather Defect Lifetimes: %d \n", num_frames);

	f_in.read("defect_lifetime_1.csv");
	f_out.open("defect_lifetimes_2.csv");

	//read data
	for (i=0; i<num_frames; i++) {
		fscanf(f_in.fp, "%d", &frame_id);
		if (frame_id != i) fprintf(stderr, "\nDefect Lifetimes Frame_ID miss match!");
		for (j=0; j<N; j++){ //for each water molecule
			fscanf(f_in.fp, ",%d", &num_bonds);
			if (num_bonds != 4){
				if (num_bonds < 4){
					if (prev_valence[j] < 4){ //still lt4 so just increment counter
						consecutive_frames[j]++;
					} else { // prev_valence >= 4 so set ...
						save_previous_defect_lifetime(prev_valence[j], consecutive_frames[j], gt4, lt4);
						consecutive_frames[j]=0;
						prev_valence[j]=num_bonds;
						consecutive_frames[j]++;
					}
				} else { // ----------------- num bonds > 4
					if (prev_valence[j] > 4){ //still gt4 so just increment counter
						consecutive_frames[j]++;
					} else { // prev_valence <= 4 so set ...
						save_previous_defect_lifetime(prev_valence[j], consecutive_frames[j], gt4, lt4);
						consecutive_frames[j]=0;
						prev_valence[j]=num_bonds;
						consecutive_frames[j]++;
					}
				}
			} else { // num bonds = 4
				save_previous_defect_lifetime(prev_valence[j], consecutive_frames[j], gt4, lt4);
				consecutive_frames[j]=0;
				prev_valence[j]=num_bonds;
			}
		}
	}
	for (j=0; j<N; j++){ //for each water molecule
		save_previous_defect_lifetime(prev_valence[j], consecutive_frames[j], gt4, lt4);
		consecutive_frames[j]=0;
	}


	fprintf(f_out.fp, "consecutive frames gt4, count\n");
	for (i=0; i<num_frames; i++) fprintf(f_out.fp, "%d, %d\n", i+1, gt4[i]);
	fprintf(f_out.fp, "consecutive frames lt4, count\n");
	for (i=0; i<num_frames; i++) fprintf(f_out.fp, "%d, %d\n", i+1, lt4[i]);

	f_out.close();
	f_in.close();

	delete [] prev_valence;
	delete [] consecutive_frames;
	delete [] gt4;
	delete [] lt4;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void output_clusters_per_frame(int num_frames)
{
	int i, j, count_clusters, end_frame, in_val;
	int num_clusters[NC+1];
	Output_File f_out, f_in;
	for (i=0; i<=NC; i++) num_clusters[i]=0;

	fprintf(stderr, "\nCount Clusters Per Frame: %d \n", num_frames);

	f_in.read("cluster_data_table.csv");
	fscanf(f_in.fp, "start_frame, end_frame, cluster_id, num_molecules, bonds_broken, bonds_created, valence_change, valence_gt4, valence_lt4, valence_change_gt4, valence_change_lt4, valence_remains_gt4, valence_remains_lt4, valence_defect_sum_start, valence_defect_sum_end, valence_defect_sum_full\n");

	for (i=1; i<num_frames; i++){
		count_clusters=0;
		do {
			for (j=0; j<16; j++){
				fscanf(f_in.fp, "%d", &in_val);
				fscanf(f_in.fp, ", ");
				if (j==1) end_frame = in_val;
			}
			count_clusters++;
		} while (end_frame == i && !feof(f_in.fp));
		num_clusters[count_clusters]++;
	}

	f_out.open("clusters_per_frame.csv");
	fprintf(f_out.fp, "number of clusters per frame, count\n");
	for (i=0; i<=NC; i++){
		if (num_clusters[i]>0) fprintf(f_out.fp, "%d, %d\n", i, num_clusters[i]);
	}


	f_in.close();
	f_out.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void output_valence_distribution_averages(int num_frames)
{
	int i, j;
	double in_val;
	double valence_count[NBONDS];
	Output_File f_out, f_in;
	for (i=0; i<NBONDS; i++) valence_count[i]=0;

	fprintf(stderr, "\nAverage Valence Distributions: %d \n", num_frames);

	f_in.read("valence_distribution.csv");
	fscanf(f_in.fp, ",0 Bonds, 1 Bonds, 2 Bonds, 3 Bonds, 4 Bonds, 5 Bonds, 6 Bonds, 7 Bonds\n"); //all over the place ... this should be ... for i = 0 to NBONDS !!!!
	f_out.open("average_valence_distribution.csv");
	fprintf(f_out.fp, ",0 Bonds, 1 Bonds, 2 Bonds, 3 Bonds, 4 Bonds, 5 Bonds, 6 Bonds, 7 Bonds\n");

	for (i=0; i<num_frames; i++){
		fscanf(f_in.fp, "percent");
		for (j=0; j<NBONDS; j++){
			fscanf(f_in.fp, ", %lf", &in_val);
			valence_count[j]+=in_val;
		}
		fscanf(f_in.fp, "\n");
	}
	fprintf(f_out.fp, "percent");
	for (i=0; i<NBONDS; i++) {
		fprintf(f_out.fp, ", %3.4f",  (valence_count[i]/num_frames) );
	}

	f_in.close();
	f_out.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	int i, start_time, finish_time, num_files, num_frames;

	start_time = time(NULL);

	if (timestamp_files) make_time_stamp_string();

	cutoff_radius_sq = cutoff_radius*cutoff_radius;

	load_analysis_configuration();
/* No longer reading input directory this way because you get files in the order they come ... the new way below is better and sorts them preventing errors
	count_input_files(num_files);
	char str_folder_files[num_files][50];
	load_input_files(str_folder_files);*/

    //new method of reading input files - includes a sort
	std::vector<std::string> str_folder_files;
	str_folder_files=read_directory(input_directory);
	if (errno) {
         fprintf(stderr, "Couldn't open the input directory, errno=%d\n", errno);
		//system("PAUSE");
		exit(1);
	}
	num_files = (int) str_folder_files.size();

	//prune skipped files
	num_frames=(int) ceil(num_files/(1.0+skip_files));
	std::vector<std::string> str_input_files;
	//char str_input_files[num_frames][50];  the old way
	prune_input_files(skip_files, num_files, str_folder_files, str_input_files);

	fprintf(stderr, "Number of input files: %d \n", num_files);
	fprintf(stderr, "Skip frames: %d \n", skip_files);
	fprintf(stderr, "Frames to process: %d \n", num_frames);
	fprintf(stderr, "Vectors to process: %d \n", (int) str_input_files.size());


	//	exit(1);

	Frame *prev_frame=new Frame;
	Frame *current_frame=new Frame;
	Tween *the_tween=new Tween;

	//first create the output files and write any header info
	//valence_distribution
	Output_File f_valence_dist;
	f_valence_dist.open("valence_distribution.csv");
	fprintf(f_valence_dist.fp, ",0 Bonds, 1 Bonds, 2 Bonds, 3 Bonds, 4 Bonds, 5 Bonds, 6 Bonds, 7 Bonds\n");

	//defect_lifetimes1
	Output_File f_defect_life;
	f_defect_life.open("defect_lifetime_1.csv");

	//cluster_data_table
	Output_File f_data_table;
	f_data_table.open("cluster_data_table.csv");
	fprintf(f_data_table.fp, "start_frame, end_frame, cluster_id, num_molecules, bonds_broken, bonds_created, valence_change");
	fprintf(f_data_table.fp, ", valence_gt4, valence_lt4, valence_change_gt4, valence_change_lt4, valence_remains_gt4, valence_remains_lt4");
	fprintf(f_data_table.fp, ", valence_defect_sum_start, valence_defect_sum_end, valence_defect_sum_full\n");



	//load first frame as prev
	//    fprintf(stderr, "before!!!");
	prev_frame->load_xyz_coordinates(num_frames, 0, str_input_files[0]);
	//    fprintf(stderr, "after!!!");
	prev_frame->analyse();
	prev_frame->write_valence_distribution(f_valence_dist.fp, 1);
	prev_frame->write_defect_lifetime_1(f_defect_life.fp);

	for (i=1; i<num_frames; i++){
	 //   fprintf(stderr, "before!!!");
		current_frame->load_xyz_coordinates(num_frames, i, str_input_files[i]);
	//    fprintf(stderr, "after!!!");
		current_frame->analyse();
		current_frame->write_valence_distribution(f_valence_dist.fp, 1);
		current_frame->write_defect_lifetime_1(f_defect_life.fp);

		compare_frames(prev_frame, current_frame, the_tween, i);
		identify_clusters(prev_frame, current_frame, the_tween);
		output_statistics(f_data_table.fp, the_tween, i);

		delete the_tween;
		the_tween=new Tween;
		if (the_tween == NULL){
			fprintf(stderr, "Error: No memory for the_tween %d\n", i);
			//system("PAUSE");
			exit(1);
		}

		*prev_frame = *current_frame;
		delete current_frame;
		current_frame=new Frame;
		if (current_frame == NULL){
			fprintf(stderr, "Error: No memory for current_frame %d\n", i);
			//system("PAUSE");
			exit(1);
		}
	}

	//close the output files
	f_data_table.close();
	f_defect_life.close();
	f_valence_dist.close();

	//post processing (read output files and analyse further)
	output_defect_lifetime_2(num_frames);
	output_clusters_per_frame(num_frames);
	output_valence_distribution_averages(num_frames);



	fprintf(stderr, "\n");

	delete prev_frame;
	delete current_frame;
	delete the_tween;

	finish_time = time(NULL);
	fprintf(stderr, "Execution time: %d\n", finish_time-start_time);

	//system("PAUSE");
	return 0;
}
