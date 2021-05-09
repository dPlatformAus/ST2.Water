struct BondLen
{
	int id;
	double bond_length;
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Atom
{
		double coord[3]; // the position coordinates, 0=x, 1=y, 2=z
		double vel[3];	// the velocity in the x,y,z directions, 0=x, 1=y, 2=z
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Bond //data on the potential link between 2 water molecules
{
		int water_id;							// the id of the other water molecule participating in the bond
		int hydrogen_id[4];		/* the hydrogen atoms that could be participating in the bond, ordered by rank (closest to furthest from mid point, ie 0=closest, then 1 then 2 then 3)
														values 1 and 2 are h1 and h2 of this molecule,
														values -1 and -2 are h1 and h2 of the other molecule */
		//consider using array for hydrogen_id storing all hydrogens that could possibly be participating in this bond ... rank_o max_rank_o would only appliy to hydrogen_id[0], or store an array for those too and analyse all of them
		int max_rank_t;						//number of hids in hydrogen_id[4]
		//for each hydrogen in hydrogen_id[4]
		int rank_o[4];						//1=primary, 2=secondary, etc (secondary means another bond to a different water molecule uses the same hydrogen as this one and the other bond's water molecule is closer) (if hydrogen belongs to other molecule then stores rank of bond in context of other molecule, not this one)
		int max_rank_o[4];				/*	1=this is the only bond this hydrogen is involved in (if hydrogen belongs to other molecule then this stores max rank for the other molecule's hydrogen)
														2=there is also a secondary bond associated with this hydrogen ... etc*/
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Water
{
	public:
		int id;
		Atom h1;
		Atom h2;
		Atom o;
		int num_bonds;	 // the number of water molecules bonded to this one (primary and secondary)
		int num_primary_bonds;	 // the number of water molecules with primary bonds to this one
		Bond bonds[10]; // the ids, & etc, of the molecules and hydrogens involved in the bonds

		Water();
		void clear_bonds();
		void add_bond(int w_id, BondLen h_ids[], int num_h_id);
		void count_primary_bonds();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Water::Water() {
	num_bonds = 0;
	num_primary_bonds = 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Water::clear_bonds() {
	num_bonds = 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Water::add_bond(int w_id, BondLen h_ids[], int num_h_id) {
	int i;
	bonds[num_bonds].water_id = w_id;
	bonds[num_bonds].max_rank_t = num_h_id;
	for (i=0; i<num_h_id; i++){
		bonds[num_bonds].hydrogen_id[i] = h_ids[i].id;
		bonds[num_bonds].rank_o[i] = 1;			//initially all bonds are considered primary with respect to additional water molecules ... the bond analysis function (in frame) will alter these values later
		bonds[num_bonds].max_rank_o[i] = 1;
	}
	num_bonds++;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Water::count_primary_bonds() { // called at end of the bond analysis function (in frame)
	int i;
	for (i=0; i<num_bonds; i++){
		if (bonds[i].rank_o[0] == 1){ //if the highest ranked hydrogen participating in this bond (so far there have only ever been 1 ... so if somthing changes that may need to revisit) is in it as a primary bond
			num_primary_bonds++;
		}
	}
}


