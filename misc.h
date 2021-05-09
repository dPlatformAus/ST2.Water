//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void make_time_stamp_string()
{
	time_t current_time;
	time(&current_time);
	struct tm *ptm = localtime(&current_time);
	int the_year, the_month, the_day, the_hour, the_minute, the_second;
	char strmonth[3], strday[3], strhr[3], strmin[3], strsec[3];
	the_year = 1900 + ptm->tm_year;
	the_month = 1 + ptm->tm_mon;

	if (the_month>9){
	  sprintf(strmonth, "%d", the_month);
	}else{
	  sprintf(strmonth, "%s%d", "0", the_month);
	}

	the_day = ptm->tm_mday;
	if (the_day>9){
	  sprintf(strday, "%d", the_day);
	}else{
	  sprintf(strday, "%s%d", "0", the_day);
	}

	the_hour = ptm->tm_hour;
	if (the_hour>9){
	  sprintf(strhr, "%d", the_hour);
	}else{
	  sprintf(strhr, "%s%d", "0", the_hour);
	}

	the_minute = ptm->tm_min;
	if (the_minute>9){
	  sprintf(strmin, "%d", the_minute);
	}else{
	  sprintf(strmin, "%s%d", "0", the_minute);
	}

	the_second = ptm->tm_sec;
	if (the_second>9){
	  sprintf(strsec, "%d", the_second);
	}else{
	  sprintf(strsec, "%s%d", "0", the_second);
	}

	sprintf(str_time_stamp, "%d%s%s%s%s%s_", the_year, strmonth, strday, strhr, strmin, strsec);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Output_File
{
	public:
		FILE *fp;
		char file_name[50];
		int file_open;

		Output_File();
		~Output_File();
		void open(const char suffix[30]);
		void read(const char suffix[30]);
		void close();
		void warning_close();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Output_File::Output_File() {
	file_open=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Output_File::~Output_File() {
	warning_close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Output_File::open(const char suffix[30]) {
	warning_close();
	sprintf(file_name, "%s%s%s", output_directory, str_time_stamp, suffix);
	fp = fopen(file_name, "w");
	file_open=1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Output_File::read(const char suffix[30]) {
	warning_close();
	sprintf(file_name, "%s%s%s", output_directory, str_time_stamp, suffix);
	if ( (fp=fopen(file_name, "r")) != NULL) {
		file_open=1;
	} else {
		fprintf(stderr, "Error: failed to load file: %s\n", file_name);
		//system("PAUSE");
		exit(1);
   }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Output_File::close() {
	if (file_open){
		fclose(fp);
		file_open=0;
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Output_File::warning_close() {
	if (file_open) {
		close();
		fprintf(stderr, "\nwarning: file \"%s\" closed implicitly\n", file_name);
	}
}
