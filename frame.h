class Frame
{
	public:
		int frame_id;
		Water molecules[N];
		char str_input_coordinates_file[50];
		int num_bond_counts[2][NBONDS]; //holds the number of molecules with a certain number of bonds [0=primary, 1=all]
		double box[DIM];
		double NFI ,NACC1, EGES, NACC2, TEMP, PRES, RHO, DELTA, RHOP, EGESP;

		Frame();
		void sort_bond_lengths_asc(BondLen sort_data[], int data_length);
		double MinImg(double x1, double x2, double period);
		double inter_particle_distance_sq(const Atom &a1, const Atom &a2);
		void get_bond_mid_point(const Atom &a1, const Atom &a2, Atom &midpoint);
		void load_xyz_coordinates(int num_frames, int id, const std::string& str_input);
		void write_molecule_linkage_header(FILE *the_file);
		void write_molecule_linkage(int mi, FILE *the_file);
		void write_linkages(int num_molecules, int water_id[], FILE *the_file);
		void write_valence_distribution(FILE *the_file, int percent);
		void write_bond_network(FILE *the_file);
		int molecules_bonded(const Water &w1, const Water &w2, BondLen h_id[], int &num_h_id);
		void get_bonds();
		void analyse_bonds();
		void count_bonds();
		void analyse();
		void output_data();
		void write_defect_lifetime_1(FILE *the_file);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Frame::Frame()
{
    int i;
	for (i=0; i<NBONDS; i++){
		num_bond_counts[0][i]=0;
		num_bond_counts[1][i]=0;
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Frame::sort_bond_lengths_asc(BondLen sort_data[], int data_length)
{//shell sort
	int swapped = 1, d, i;
	BondLen temp;
	d=data_length;
	while( swapped || (d>1)){
		swapped = 0;
		d = (d+1) / 2;
		for (i = 0; i < (data_length - d); i++){
			if (sort_data[i + d].bond_length < sort_data[i].bond_length){
				temp = sort_data[i + d];
				sort_data[i + d] = sort_data[i];
				sort_data[i] = temp;
				swapped = 1;
			}
		}
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Frame::MinImg(double x1, double x2, double period)
{
  /* Given the x positions of two atoms (or the y's or z's), and the size of the periodic box
  "period", this function calculates how far their nearest images are apart.  If the two particles
  are at opposite ends of the box, there is actually a pair of images that are quite close to one
  another.  */
  double dist, excess;
  dist = x2-x1;
  excess = period * (int) (dist/(period/2)); //excess is the number of box widths to remove to bring the distance back to the closest image (which must be within half a box width)
  return dist-excess;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Frame::inter_particle_distance_sq(const Atom &a1, const Atom &a2)
{
	/* Given two atoms, this calculates the square of the inter-particle distance.
	It uses the nearest periodic images of those atoms. */
	double distsq=0.0, dist;
	int d;
	for (d=0; d<DIM; d++) {
		dist = MinImg(a1.coord[d], a2.coord[d], box[d]);
		distsq += dist * dist;
	}
	return distsq;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Frame::get_bond_mid_point(const Atom &a1, const Atom &a2, Atom &midpoint)
{
	int d;
	for (d=0; d<DIM; d++) {
		midpoint.coord[d]=a1.coord[d]+MinImg(a1.coord[d], a2.coord[d], box[d])/2;
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Frame::load_xyz_coordinates(int num_frames, int id, const std::string& str_input )
{
	/* Tries to load the coordinates file specified in config/analysis_config.txt*/
	int i;
	FILE *fp;
	char file_name[56];
	frame_id=id;
	sprintf(str_input_coordinates_file, "%s", str_input.c_str());
	sprintf(file_name, "%s%s", input_directory, str_input.c_str());
	fprintf(stderr, "Frame %d of %d loading file: %s\n", frame_id+1, num_frames, file_name);

	if ( (fp=fopen(file_name, "r")) != NULL) {
		//load header parameters
		if (alternate_header){
			fscanf(fp, "%lf%lf%lf\n", &NFI ,&NACC1, &box[0]);
			fscanf(fp, "%lf%lf%lf%lf\n", &box[1], &box[2], &EGES, &NACC2);
			fscanf(fp, "%lf%lf%lf\n", &TEMP, &PRES, &RHO);
			fscanf(fp, "%lf%lf%lf\n", &DELTA, &RHOP, &EGESP);
		} else {
			fscanf(fp, "%lf%lf%lf%lf\n", &NFI ,&NACC1, &box[0], &box[1]);
			fscanf(fp, "%lf%lf%lf\n", &box[2], &EGES, &NACC2);
			fscanf(fp, "%lf%lf%lf\n", &TEMP, &PRES, &RHO);
			fscanf(fp, "%lf%lf%lf\n", &DELTA, &RHOP, &EGESP);
		}

		//load atoms

		//read in positions
		for (i=0; i<N; i++){ //for each water molecule
			molecules[i].id = i;
			fscanf(fp, "%lf%lf%lf\n", &molecules[i].h1.coord[0], &molecules[i].h1.coord[1], &molecules[i].h1.coord[2]);
			fscanf(fp, "%lf%lf%lf\n", &molecules[i].h2.coord[0], &molecules[i].h2.coord[1], &molecules[i].h2.coord[2]);
			fscanf(fp, "%lf%lf%lf\n", &molecules[i].o.coord[0], &molecules[i].o.coord[1], &molecules[i].o.coord[2]);
		}
		if (read_velocities){
			//read in velocities
			for (i=0; i<N; i++){ //for each water molecule
				fscanf(fp, "%lf%lf%lf\n", &molecules[i].h1.vel[0], &molecules[i].h1.vel[1], &molecules[i].h1.vel[2]);
				fscanf(fp, "%lf%lf%lf\n", &molecules[i].h2.vel[0], &molecules[i].h2.vel[1], &molecules[i].h2.vel[2]);
				fscanf(fp, "%lf%lf%lf\n", &molecules[i].o.vel[0], &molecules[i].o.vel[1], &molecules[i].o.vel[2]);
			}
		}
		fclose(fp);
   } else {
		fprintf(stderr, "Error: failed to load xyz configuration\n");
		//system("PAUSE");
		exit(1);
   }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 void Frame::write_molecule_linkage_header(FILE *the_file)
{
	fprintf(the_file, "Molecule ID, Num Bonds, Num Primary, Bond 1 Molecule ID, Bond 1 Num Eligible Hydrogen, Hydrogen ID, Hydrogen Rank, Hydrogen Max Rank ...\n");
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 void Frame::write_molecule_linkage(int mi, FILE *the_file)
{
	int j, k;
		fprintf(the_file, "%d,%d,%d", mi, molecules[mi].num_bonds, molecules[mi].num_primary_bonds);
		for (j=0; j<molecules[mi].num_bonds; j++) {
			fprintf(the_file, ",%d,%d", molecules[mi].bonds[j].water_id, molecules[mi].bonds[j].max_rank_t);
			for (k=0; k<molecules[mi].bonds[j].max_rank_t; k++){
				fprintf(the_file, ",%d,%d,%d", molecules[mi].bonds[j].hydrogen_id[k], molecules[mi].bonds[j].rank_o[k], molecules[mi].bonds[j].max_rank_o[k]);
			}
		}
		fprintf(the_file, "\n");
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 void Frame::write_linkages(int num_molecules, int water_id[], FILE *the_file)
{
	int i;
	write_molecule_linkage_header(the_file);
	for (i=0; i<num_molecules; i++){ //for each water molecule
		write_molecule_linkage(water_id[i], the_file);
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 void Frame::write_valence_distribution(FILE *the_file, int percent)
{
	int i;
	if (percent){
		fprintf(the_file, "percent");
		for (i=0; i<NBONDS; i++){
			fprintf(the_file, ", %3.4f",  ( (num_bond_counts[0][i] * 100.0)/ N) );
		}
	} else {
		fprintf(the_file, "count");
		for (i=0; i<NBONDS; i++){
			fprintf(the_file, ", %d", num_bond_counts[0][i]);
		}
	}
	fprintf(the_file, "\n");
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 void Frame::write_bond_network(FILE *the_file)
{
	int i;
	fprintf(the_file, ",0 Bonds, 1 Bonds, 2 Bonds, 3 Bonds, 4 Bonds, 5 Bonds, 6 Bonds, 7 Bonds\n");
	write_valence_distribution(the_file, 0);
	write_valence_distribution(the_file, 1);
	write_molecule_linkage_header(the_file);
	for (i=0; i<N; i++){ //for each water molecule
		write_molecule_linkage(i, the_file);
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int Frame::molecules_bonded(const Water &w1, const Water &w2, BondLen h_id[], int &num_h_id)
{
	/* check if two molecules are bonded. If bonded return 1. If not bonded return 0 */
	int i, result = 0;
	double rsq, cutoff_hdist;
	Atom midpoint; // the bond midpoint coordinates, Atom object is used to pass to inter_particle_distance_sq(const Atom &a1, const Atom &a2)
	double hdist[4]; //squared distances of the hydrogen atoms from the midpoint, 0=w1.h1, 1=w1.h2, 2=w2.h1, 3=w2.h2
	num_h_id=0;
	if (w1.id != w2.id) {
		rsq = inter_particle_distance_sq(w1.o, w2.o);
		if (rsq < cutoff_radius_sq) {
			get_bond_mid_point(w1.o, w2.o, midpoint);
			hdist[0]=inter_particle_distance_sq(w1.h1, midpoint);
			hdist[1]=inter_particle_distance_sq(w1.h2, midpoint);
			hdist[2]=inter_particle_distance_sq(w2.h1, midpoint);
			hdist[3]=inter_particle_distance_sq(w2.h2, midpoint);
			cutoff_hdist=cutoff_hdist_factor*rsq;	//calculate maximum distance a hydrogen can be from the mid point and still be considered "in the bond"
			for (i=0; i<4; i++){
				if (hdist[i] < cutoff_hdist){	//hydrogen is within the cutoff length
					result = 1;		//there is a bond
					h_id[num_h_id].id = i+1;
					if (h_id[num_h_id].id>2) h_id[num_h_id].id=2-h_id[num_h_id].id; //1=1, 2=2, 3=-1, 4=-2
					h_id[num_h_id].bond_length = hdist[i];
					num_h_id++;
				}
			}
			sort_bond_lengths_asc(h_id, num_h_id);
		}
	}
	return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Frame::get_bonds()
{
	int i, j, k, num_h_id;
	BondLen h_id[4];
	for (i=0; i<N; i++){ //for each water molecule
		for (j=i+1; j<N; j++){ //for each water molecule
			if (molecules_bonded(molecules[i], molecules[j], h_id, num_h_id)) {
				molecules[i].add_bond(j, h_id, num_h_id);
				for (k=0; k<4; k++){
					h_id[k].id = -1*h_id[k].id;		// molecules_bonded returns negative value for h if H belongs to j, so the opposite value is true for j's version of the bond
				}
				molecules[j].add_bond(i, h_id, num_h_id);
			}
		}
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Frame::count_bonds()
{
	int i, j, k, primary_count, has_primary;
	for (i=0; i<N; i++){ //for each water molecule
		primary_count=0;
		for (j=0; j<molecules[i].num_bonds; j++){
			has_primary=0;
			for (k=0; k<molecules[i].bonds[j].max_rank_t; k++){
				if (molecules[i].bonds[j].rank_o[k] == 1){
					has_primary=1;
				}
			}
			primary_count++; /// should this be if (has_primary)  ***!!!!!!
		}
        if (primary_count>=NBONDS) {
           fprintf(stderr, "error: primary count too high on molecule %d\n", i);
        }
		num_bond_counts[0][primary_count]++;

		if (molecules[i].num_bonds>=NBONDS) {
           fprintf(stderr, "error: bond count too high on molecule %d:  %d\n", i, molecules[i].num_bonds);
        }
		num_bond_counts[1][molecules[i].num_bonds]++;
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Frame::analyse_bonds()
{
	int i, j, k, l, bw, new_rank;
	int nh1, nh2, ch1, ch2;
	BondLen h1_bid[NBONDS], h2_bid[NBONDS];

	for (i=0; i<N; i++){ //for each water molecule, do local hydrogen in multiple bonds routine (make function)
		nh1=nh2=0;
		for (j=0; j<molecules[i].num_bonds; j++){ //count number of bonds local hydrogens are involved in
			for (k=0; k<molecules[i].bonds[j].max_rank_t; k++){
				if (molecules[i].bonds[j].hydrogen_id[k]==1) nh1++;
				if (molecules[i].bonds[j].hydrogen_id[k]==2) nh2++;
			}
		}
		if (nh1>1 || nh2>1){ //if h1 or h2 involved in more than one bond
			//record bond lengths so primary, secondary etc can be identified
			ch1=ch2=0;
			for (j=0; j<molecules[i].num_bonds; j++){
				bw=molecules[i].bonds[j].water_id; //the water_id stored in the bond (the other molecule)
				if (nh1>1){ //record bond ids and distances for h1
					for (k=0; k<molecules[i].bonds[j].max_rank_t; k++){
						if (molecules[i].bonds[j].hydrogen_id[k]==1){
							if (ch1>=NBONDS){
								fprintf(stderr, "ch1= %d!!\n", ch1);
								//system("PAUSE");
								exit(1);
							}
							h1_bid[ch1].id=j;
							h1_bid[ch1].bond_length=inter_particle_distance_sq(molecules[i].h1, molecules[bw].o);
							ch1++;
						}
					}
				}
				if (nh2>1){//record bond ids and distances for h2
					for (k=0; k<molecules[i].bonds[j].max_rank_t; k++){
						if (molecules[i].bonds[j].hydrogen_id[k]==2){
							if (ch2>=NBONDS){
								fprintf(stderr, "ch2= %d!!\n", ch1);
								//system("PAUSE");
								exit(1);
							}
							h2_bid[ch2].id=j;
							h2_bid[ch2].bond_length=inter_particle_distance_sq(molecules[i].h2, molecules[bw].o);
							ch2++;
						}
					}
				}
			}
		}
		//sort bond lengths and record ranks
		if (nh1>1){ //if hydrogen 1 has multiple bonds
			sort_bond_lengths_asc(h1_bid, nh1);
			for (j=0; j<nh1; j++){
				new_rank=j+1;
				for (k=0; k<molecules[i].bonds[h1_bid[j].id].max_rank_t; k++){
					if (molecules[i].bonds[h1_bid[j].id].hydrogen_id[k]==1){
						molecules[i].bonds[h1_bid[j].id].rank_o[k]=new_rank;
						molecules[i].bonds[h1_bid[j].id].max_rank_o[k]=nh1;
					}
				}
				//store in foreign molecule's bond too
				bw=molecules[i].bonds[h1_bid[j].id].water_id; //the water_id stored in the bond (the other molecule)
				for (k=0; k<molecules[bw].num_bonds; k++){ //for each bond on the other molecule
					if (molecules[bw].bonds[k].water_id==i){ //if the bond is to this molecule, record the rank data to the other molecule's bond record too
						for (l=0; l<molecules[bw].bonds[k].max_rank_t; l++){
							if (molecules[bw].bonds[k].hydrogen_id[l]==-1){
								molecules[bw].bonds[k].rank_o[l]=new_rank;
								molecules[bw].bonds[k].max_rank_o[l]=nh1;
							}
						}
					}
				}
			}
		}
		if (nh2>1){ //if hydrogen 2 has multiple bonds
			sort_bond_lengths_asc(h2_bid, nh2);
			for (j=0; j<nh2; j++){
				new_rank=j+1;
				for (k=0; k<molecules[i].bonds[h2_bid[j].id].max_rank_t; k++){
					if (molecules[i].bonds[h2_bid[j].id].hydrogen_id[k]==2){
						molecules[i].bonds[h2_bid[j].id].rank_o[k]=new_rank;
						molecules[i].bonds[h2_bid[j].id].max_rank_o[k]=nh2;
					}
				}
				//store in foreign molecule's bond too
				bw=molecules[i].bonds[h2_bid[j].id].water_id; //the water_id stored in the bond (the other molecule)
				for (k=0; k<molecules[bw].num_bonds; k++){ //for each bond on the other molecule
					if (molecules[bw].bonds[k].water_id==i){ //if the bond is to this molecule, record the rank data to the other molecule's bond record too
						for (l=0; l<molecules[bw].bonds[k].max_rank_t; l++){
							if (molecules[bw].bonds[k].hydrogen_id[l]==-2){
								molecules[bw].bonds[k].rank_o[l]=new_rank;
								molecules[bw].bonds[k].max_rank_o[l]=nh2;
							}
						}
					}
				}
			}
		}
		molecules[i].count_primary_bonds();
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Frame::analyse()
{
	get_bonds();
	analyse_bonds();
	count_bonds();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Frame::output_data()
{
	Output_File f_out;
	char file_suffix[30];
	sprintf(file_suffix, "frame_%d_bond_network.csv", frame_id);
	f_out.open(file_suffix);
	//write_bond_network(stderr);
	write_bond_network(f_out.fp);
	f_out.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Frame::write_defect_lifetime_1(FILE *the_file)
{
	int i;
	fprintf(the_file, "%d", frame_id);
	for (i=0; i<N; i++){
		if (molecules[i].num_bonds != 4) {
			fprintf(the_file, ", %d", molecules[i].num_bonds);
		} else {
			fprintf(the_file, ", 4");
		}
	}
	fprintf(the_file, "\n");
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


