struct Event
{
	int type; //0=bond broken 1=bond created
	int water_id_1;
	int water_id_2;
	int hydrogen_id; /* the hydrogen atom participating in the bond
										values 1 and 2 are h1 and h2 of water_id_1,
										values -1 and -2 are h1 and h2 of water_id_2 */
	int rank_t;
	int max_rank_t;
	int rank_o;
	int max_rank_o;
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Cluster_molecule
{
	int water_id;
	int start_valence;
	int end_valence;
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Cluster
{
	public:
		int num_molecules;
		Cluster_molecule molecules[NCM]; //index is just 0,1,2,3 ... not water_id (that's stored as a property of the Cluster_molecule struct)
		int bonds_broken;
		int bonds_created;
		int valence_gt4;
		int valence_lt4;
		int valence_change;
		int valence_change_gt4;
		int valence_change_lt4;
		int valence_remains_gt4;
		int valence_remains_lt4;
		int valence_defect_sum_start;
		int valence_defect_sum_end;
		int valence_defect_sum_full; //if bonds_broken==bonds_created and valence_change then ( < 0 defect pair created, 0 defects moved, > 0 defect recombined ) number is always even

		Cluster();
		void calculate_valence_defect_sums();
		void add_molecule(int water_id, int start_bonds, int end_bonds);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Cluster::Cluster() {
	num_molecules = 0;
	bonds_broken = 0;
	bonds_created = 0;
	valence_change = 0;
	valence_gt4 = 0;
	valence_lt4 = 0;
	valence_change_gt4 = 0;
	valence_change_lt4 = 0;
	valence_remains_gt4 = 0;
	valence_remains_lt4 = 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Cluster::calculate_valence_defect_sums() {
	int i;
	valence_defect_sum_start=valence_defect_sum_end=0;
	for (i=0; i<num_molecules; i++){
		valence_defect_sum_start+=abs(molecules[i].start_valence-4);
		valence_defect_sum_end+=abs(molecules[i].end_valence-4);
	}
	valence_defect_sum_full=valence_defect_sum_start-valence_defect_sum_end;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Cluster::add_molecule(int water_id, int start_bonds, int end_bonds) {
	if (num_molecules<NCM){
		molecules[num_molecules].water_id = water_id;
		molecules[num_molecules].start_valence = start_bonds;
		molecules[num_molecules].end_valence = end_bonds;
		if (start_bonds>4 || end_bonds>4) valence_gt4=1; //these 2 only record that the valence is <> 4, doesn't check whether it changes
		if (start_bonds<4 || end_bonds<4) valence_lt4=1;
		if (start_bonds != end_bonds){
			valence_change=1;
			if (start_bonds>4 || end_bonds>4)  valence_change_gt4=1; //these 2 record that the valence is <> 4 and the valence change
			if (start_bonds<4 || end_bonds<4) valence_change_lt4=1;
		} else {
			if (start_bonds>4 || end_bonds>4) valence_remains_gt4=1; //these 2 record that the valence is <> 4 and the valence does not change
			if (start_bonds<4 || end_bonds<4) valence_remains_lt4=1;
		}
		num_molecules++;
		calculate_valence_defect_sums();
	} else {
		fprintf(stderr, "\nToo many molecules > %d in cluster", NCM);
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

class Tween
{
	public:
		int previous_frame_id;
		int next_frame_id;
		int num_events;
		int num_lost;
		int num_made;
		Event events[NE]; //some tweens have LOTS of events between frames (each event is stored twice) (would prefer to use a linked list for this)
		int cluster_ids[N]; //index is same as water_id, value is cluster_id, used to recursively identify clusters
		int num_clusters;
		Cluster clusters[NC+1]; //all relevent data about clusters, index is cluster_id value is a cluster object

		Tween();
		void clear_events();
		void set_frames(int pf, int nf);
		void add_event(int t, int w1, int w2, int h, int rt, int mrt, int ro, int mro);
		void write_tween(FILE *the_file);
		void write_events(int num_molecules, int water_id[], FILE *the_file);
		void output_data();
		void assign_cluster_id(int water_id);
		void get_clusters();
		void output_clusters();
		void write_clusters(FILE *the_file);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Tween::Tween() {
	int i;
	num_events = 0;
	num_lost = 0;
	num_made = 0;
	num_clusters = 0;
	for (i=0; i<N; i++) {
		cluster_ids[i]=0;
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Tween::clear_events() {
	num_events = 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Tween::set_frames(int pf, int nf) {
	previous_frame_id=pf;
	next_frame_id=nf;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Tween::add_event(int t, int w1, int w2, int h, int rt, int mrt, int ro, int mro) {
	if (num_events<NE){
		events[num_events].type = t;
		events[num_events].water_id_1 = w1;
		events[num_events].water_id_2 = w2;
		events[num_events].hydrogen_id = h;
		events[num_events].rank_t = rt;
		events[num_events].max_rank_t = mrt;
		events[num_events].rank_o = ro;
		events[num_events].max_rank_o = mro;
		num_events++;
		if (t){ //made
			num_made++;
		} else { //lost
			num_lost++;
		}
	} else {
		fprintf(stderr, "\nToo many events > %d between frames %d & %d", NE, previous_frame_id, next_frame_id);
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 void Tween::write_tween(FILE *the_file)
{
	int i;

	fprintf(the_file, "Bonds created, %d, Bonds broken, %d\n", num_made/2, num_lost/2); //each event is counted twise (once in the context of each molecule in the event) so divide by two to get the real number
	fprintf(the_file, "Type, Water ID 1, Water ID 2, Hydrogen ID, Hydrogen Rank, Max Hydrogen Rank, Other Rank, Max Other Rank\n");
	for (i=0; i<num_events; i++){ //for each event
		fprintf(the_file, "%d,%d,%d,%d,%d,%d,%d,%d\n", events[i].type, events[i].water_id_1, events[i].water_id_2, events[i].hydrogen_id, events[i].rank_t, events[i].max_rank_t, events[i].rank_o, events[i].max_rank_o);
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 void Tween::write_events(int num_molecules, int water_id[], FILE *the_file)
{
	int i, j;

	fprintf(the_file, "Type, Water ID 1, Water ID 2, Hydrogen ID, Hydrogen Rank, Max Hydrogen Rank, Other Rank, Max Other Rank\n");
	for (i=0; i<num_events; i++){ //for each event
		for (j=0; j<num_molecules; j++){
			if (water_id[j] == events[i].water_id_1){
				fprintf(the_file, "%d,%d,%d,%d,%d,%d,%d,%d\n", events[i].type, events[i].water_id_1, events[i].water_id_2, events[i].hydrogen_id, events[i].rank_t, events[i].max_rank_t, events[i].rank_o, events[i].max_rank_o);
			}
		}
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Tween::output_data()
{
	Output_File f_out;
	char file_suffix[30];
	sprintf(file_suffix, "frames_%d-%d_changes.csv", previous_frame_id, next_frame_id);
	f_out.open(file_suffix);
	write_tween(f_out.fp);
	f_out.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Tween::assign_cluster_id(int water_id)
{
	int i;
	cluster_ids[water_id]=num_clusters;
	for (i=0; i<num_events; i++){
		if (events[i].water_id_1 == water_id && cluster_ids[events[i].water_id_2]==0 && events[i].rank_t==1){ //had to add ************************************* && events[i].rank_t==1 ... if rank_t > 1 that means another hydrogen is better suited to this bond ... that one is counted ... this one shouldn't be! ... can get the nothing happening clusters
			assign_cluster_id(events[i].water_id_2);
		}
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Tween::get_clusters()
{
	int i;
	for (i=0; i<num_events; i++){						//had to add ************************************* && events[i].rank_t==1 ... if rank_t > 1 that means another hydrogen is better suited to this bond ... that one is counted ... this one shouldn't be! ... can get the nothing happening clusters
		if (cluster_ids[events[i].water_id_1]==0 && events[i].rank_t==1){ //the 1st water id in this event has not been assigned to a cluster -- could add select primary only, or include secondary here ***
			if (num_clusters<NC){
				num_clusters++;
				assign_cluster_id(events[i].water_id_1);
			} else{
				fprintf(stderr, "\nToo many clusters > %d", NC);
			}
		}
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Tween::output_clusters()
{
	Output_File f_out;
	char file_suffix[30];
	sprintf(file_suffix, "frames_%d-%d_clusters.csv", previous_frame_id, next_frame_id);
	f_out.open(file_suffix);
	write_clusters(f_out.fp);
	f_out.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 void Tween::write_clusters(FILE *the_file)
{
	int i;

	fprintf(the_file, "Water ID, Cluster ID");
	for (i=0; i<N; i++){
		fprintf(the_file, "%d,%d", i, cluster_ids[i]);
		fprintf(the_file, "\n");
	}
}
